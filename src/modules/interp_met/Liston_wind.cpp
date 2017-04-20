#include "Liston_wind.hpp"

Liston_wind::Liston_wind(config_file cfg)
        :module_base(parallel::domain)

{
    depends_from_met("U_R");
    depends_from_met("vw_dir");

    provides("U_R");
    provides("vw_dir");
    provides("vw_dir_divergence");
    distance = cfg.get<double>("distance",300);//300.0;

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

//Calculates the curvature required
void Liston_wind::init(mesh domain)
{

    ys = cfg.get("ys",0.5);
    yc = cfg.get("yc",0.5);

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<lwinddata>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
        face->coloured = false;
    }

    if (domain->face(0)->has_parameter("Liston curvature"))
    {
        LOG_DEBUG << "Liston curvature available as parameter, using that.";
        return;
    }

    double curmax = -9999.0;

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
            auto face = domain->face(i);

            Point_3 me = face->center();
            Delaunay::Face_handle north;
            Delaunay::Face_handle south;
            Delaunay::Face_handle east;
            Delaunay::Face_handle west;

            Delaunay::Face_handle northeast;
            Delaunay::Face_handle northwest;
            Delaunay::Face_handle southeast;
            Delaunay::Face_handle southwest;

            north = domain->find_closest_face( math::gis::point_from_bearing(me,0,distance) ) ; // me.x(), me.y() + distance
            south = domain->find_closest_face( math::gis::point_from_bearing(me,180,distance) ); //me.x(), me.y() - distance
            west = domain->find_closest_face(  math::gis::point_from_bearing(me,270,distance) ); //me.x() - distance, me.y()
            east = domain->find_closest_face(  math::gis::point_from_bearing(me,90,distance)  ); //me.x() + distance, me.y()

            double z = face->get_z();
            double zw = west->get_z();
            double ze = east->get_z();
            double zs = south->get_z();
            double zn = north->get_z();

            double znw = 0.;
            double zne = 0.;
            double zse = 0.;
            double zsw = 0.;

            northeast = domain->find_closest_face(math::gis::point_from_bearing(me,45,distance)); //me.x() + distance, me.y() + distance
            zne = northeast->get_z();

            northwest = domain->find_closest_face(math::gis::point_from_bearing(me,315,distance)); //me.x() - distance, me.y() + distance
            znw = northwest->get_z();

            southeast = domain->find_closest_face(math::gis::point_from_bearing(me,135,distance)); //me.x() + distance, me.y() - distance
            zse = southeast->get_z();

            southwest = domain->find_closest_face(math::gis::point_from_bearing(me,225,distance)); //me.x() - distance, me.y() - distance
            zsw = southwest->get_z();

            double curve = .25 * ((z - .5 * (zw + ze)) / (2 * distance) + (z - .5 * (zs + zn)) / (2 * distance) +
                                  (z - .5 * (zsw + zne)) / (2 * sqrt(2 * distance)) +
                                  (z - .5 * (znw + zse)) / (2 * sqrt(2 * distance)));

            auto* c = face->get_module_data<lwinddata>(ID);
            c->curvature = curve;

            if (fabs(curve) > curmax)
            {
                curmax = fabs(curve);
            }


        }


    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto* c = face->get_module_data<lwinddata>(ID);

        double value = c->curvature / curmax / 2.0;//rescale to [-0.5,+0.5];

        //with very coarse meshes, with very few total triangles,
        // there are edge cases where curmax=0 and makes curvature NAN. Just set it to 0, no curvature, and don't do silly speedup/down
        c->curvature = std::isnan(value) ? 0 : value;

        face->set_parameter("Liston curvature",  c->curvature);
    }


    if ( cfg.get("serialize",true) )
    {
        LOG_DEBUG << "Serializing liston curvature";
        domain->serialize_parameter(cfg.get("serialize_output", "liston_curvature.mesh"),
                                    "Liston curvature");
    }

}


void Liston_wind::run(mesh domain)
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

            double W = s->get("U_R");
            W = std::max(W, 0.1);

            double theta = s->get("vw_dir") * 3.14159 / 180.;
            double phi = math::gis::bearing_to_polar(s->get("vw_dir") );

            double zonal_u = -W * sin(theta);
            double zonal_v = -W * cos(theta);
//            double zonal_u = -W * sin(phi); //negate as it needs to be the direction the wind is *going*
//            double zonal_v = -W * cos(phi);

            u.push_back(boost::make_tuple(s->x(), s->y(), zonal_u));
            v.push_back(boost::make_tuple(s->x(), s->y(), zonal_v));
        }
//http://mst.nerc.ac.uk/wind_vect_convs.html

        auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
        double zonal_u = face->get_module_data<lwinddata>(ID)->interp(u, query);
        double zonal_v = face->get_module_data<lwinddata>(ID)->interp(v, query);

        double theta = 3.0 * PI * 0.5 - atan2(zonal_v, zonal_u);
//        double theta = atan2(-zonal_v, -zonal_u);
        if (theta > 2 * PI)
            theta = theta - 2 * PI;

        //eqn 15
        double omega_s = face->slope() * cos(theta - face->aspect());

        if (fabs(omega_s) > max_omega_s)
            max_omega_s = fabs(omega_s);

        double W = sqrt(zonal_u * zonal_u + zonal_v * zonal_v);

        face->get_module_data<lwinddata>(ID)->corrected_theta = theta;
        face->get_module_data<lwinddata>(ID)->W = W;
    }
    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        double theta= face->get_module_data<lwinddata>(ID)->corrected_theta;
        double W= face->get_module_data<lwinddata>(ID)->W;

        //what liston calls 'wind slope'
        double omega_s = face->slope() * cos(theta - face->aspect());

        //scale between [-0.5,0.5]
        omega_s = omega_s / max_omega_s / 2.0;


        double omega_c = face->get_parameter("Liston curvature");

//        double ys = 0.5;
//        double yc = 0.5;
        double ys = 1;
        double yc = 1;

        double Ww = 1 + ys * omega_s + yc * omega_c;
        W = W * Ww;

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
            if (theta > 2*M_PI)
                theta = theta - 2*M_PI;
            else if (theta < 0.0)
                theta = theta + 2*M_PI;

        }

        W = std::max(W,0.1);
        face->set_face_data("U_R", W);
        face->set_face_data("vw_dir", theta * 180.0 / 3.14159);


        Vector_2 v = math::gis::bearing_to_cartesian(theta* 180.0 / 3.14159);
        Vector_3 v3(-v.x(),-v.y(), 0); //negate as direction it's blowing instead of where it is from!!

        face->set_face_data("vw_dir_divergence",dirdiff* 180.0 / 3.14159);
        face->set_face_vector("wind_direction",v3);

    }


}

Liston_wind::~Liston_wind()
{

}