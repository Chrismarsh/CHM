#include "MS_wind.hpp"

MS_wind::MS_wind(config_file cfg)
        :module_base(parallel::domain)

{
    depends_from_met("U_R");
    depends_from_met("vw_dir");

    provides("U_R");
    provides("vw_dir");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

//Calculates the curvature required
void MS_wind::init(mesh domain)
{
    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
        d->interp_smoothing.init(interp_alg::tpspline,3,{ {"reuse_LU","true"}});

    }




}


void MS_wind::run(mesh domain)
{
    double PI = 3.14159;

    // omega_s needs to be scaled on [-0.5,0.5]
    double max_omega_s = -99999.0;

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        std::vector<boost::tuple<double, double, double> > u;
        std::vector<boost::tuple<double, double, double> > v;
        for (auto &s : global_param->get_stations(face->get_x(), face->get_y()))
        {
            if (is_nan(s->get("U_R")) || is_nan(s->get("vw_dir")))
                continue;

            double theta = s->get("vw_dir") * M_PI / 180.;

            auto f = domain->find_closest_face(s->x(),s->y());
            //figure out which lookup map we need
            int d = int(theta*180/M_PI/45.);
            if (d == 0) d = 8;
            double speedup = f->get_parameter("MS"+std::to_string(d));

            double W = s->get("U_R") / speedup;
            W = std::max(W, 0.1);


//            double phi = math::gis::bearing_to_polar(s->get("vw_dir") );

            double zonal_u = -W * sin(theta);
            double zonal_v = -W * cos(theta);

            u.push_back(boost::make_tuple(s->x(), s->y(), zonal_u));
            v.push_back(boost::make_tuple(s->x(), s->y(), zonal_v));
        }
//http://mst.nerc.ac.uk/wind_vect_convs.html

        auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
        double zonal_u = face->get_module_data<data>(ID)->interp(u, query);
        double zonal_v = face->get_module_data<data>(ID)->interp(v, query);

        double theta = 3.0 * M_PI * 0.5 - atan2(zonal_v, zonal_u);

        if (theta > 2.0 * M_PI)
            theta = theta - 2.0 * M_PI;

        //eqn 15
        double omega_s = face->slope() * cos(theta - face->aspect());

        if (fabs(omega_s) > max_omega_s)
            max_omega_s = fabs(omega_s);

        double W = sqrt(zonal_u * zonal_u + zonal_v * zonal_v);

        face->get_module_data<data>(ID)->corrected_theta = theta;
        face->get_module_data<data>(ID)->W = W;
    }
    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        double theta= face->get_module_data<data>(ID)->corrected_theta;
        double W= face->get_module_data<data>(ID)->W;

        //what liston calls 'wind slope'
        double omega_s = face->slope() * cos(theta - face->aspect());

        //scale between [-0.5,0.5]
        omega_s = omega_s / (max_omega_s * 2.0);

        //scale wind here via MS

        double aspect = face->aspect();
        double dirdiff = 0;

        double d2r = M_PI/180.0;


        //Ryan terain divergence
        // follow LIston's micromet implimentation to handle 0-360
        if (aspect > 270.0*d2r &&  theta < 90.0*d2r)
            dirdiff = aspect - theta - 360.0*d2r;
        else if (aspect < 90.0*d2r && theta > 270.0*d2r)
            dirdiff = aspect - theta + 360.0*d2r;
        else
            dirdiff = aspect - theta;

        if (std::abs(dirdiff) < 90.0*d2r)
        {
            theta = theta - 0.5 * std::min(omega_s, 45.0*d2r) * sin((2.0 * dirdiff));
            if (theta > 2.0*M_PI)
                theta = theta - 2*M_PI;
            else if (theta < 0.0)
                theta = theta + 2.0*M_PI;

        }

        //figure out which lookup map we need
        int d = int(theta*180/M_PI/45.);
        if (d == 0) d = 8;

        double speedup = face->get_parameter("MS"+std::to_string(d));
        W = W*speedup;

        W = std::max(W,0.1);
        W = std::min(W,30.0);
        face->set_face_data("U_R", W);
        face->set_face_data("vw_dir", theta * 180.0 / 3.14159);


        Vector_2 v = math::gis::bearing_to_cartesian(theta* 180.0 / 3.14159);
        Vector_3 v3(-v.x(),-v.y(), 0); //negate as direction it's blowing instead of where it is from!!

        face->set_face_vector("wind_direction",v3);

    }


#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

        auto face = domain->face(i);
        std::vector<boost::tuple<double, double, double> > u;
        for (size_t j = 0; j < 3; j++)
        {
            auto neigh = face->neighbor(j);
            if (neigh != nullptr)
                u.push_back(boost::make_tuple(neigh->get_x(), neigh->get_y(), neigh->face_data("U_R")));
        }

        auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());

        double new_u = face->get_module_data<data>(ID)->interp_smoothing(u, query);
        face->get_module_data<data>(ID)->temp_u = new_u;
    }
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        face->set_face_data("U_R",face->get_module_data<data>(ID)->temp_u );
    }
}

MS_wind::~MS_wind()
{

}