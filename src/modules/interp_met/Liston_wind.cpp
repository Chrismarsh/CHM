#include "Liston_wind.hpp"

Liston_wind::Liston_wind(config_file cfg)
        :module_base(parallel::domain)

{
    depends_from_met("u");
    depends_from_met("vw_dir");

    provides("vw");
    provides("vw_dir");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

//Calculates the curvature required
void Liston_wind::init(mesh domain, boost::shared_ptr<global> global_param)
{

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<lwinddata>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->stations.size());
        face->coloured = false;
    }

    if (domain->face(0)->has_parameter("Liston curvature"))
    {
        LOG_DEBUG << "Liston curvature available as parameter, using that.";
        return ;
    }

    double curmax = -9999.0;

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
            auto face = domain->face(i);

            std::set<Point_3> myPoints;
            std::stack<Delaunay::Face_handle> neighbours;
            std::vector<Delaunay::Face_handle> neighbours_to_uncolor;
            int cur_depth = 0;
            int max_depth = 10;

            neighbours.push(face);

            while (!neighbours.empty() && (cur_depth <= max_depth))
            {
                Delaunay::Face_handle nface = neighbours.top();
                neighbours.pop();
                //this triangle
                //origin of fitting coord system = first input data point
                myPoints.insert(nface->center());
                for (int i = 0; i < 3; i++)
                {
                    myPoints.insert(nface->vertex(i)->point());
                }

                //add all of f's neighbours
                for (int k = 0; k < 3; k++)
                {
                    auto n = nface->neighbor(k);
//                    if (!domain->is_infinite(n) && !n->coloured)
                    if (n!=nullptr && !n->coloured)
                    {
                        n->coloured = true;
                        neighbours_to_uncolor.push_back(n);
                        neighbours.push(n);
                    }

                }
                cur_depth++;
            }

            for (auto it:neighbours_to_uncolor)
            {
                it->coloured = false;
            }
            neighbours_to_uncolor.clear();

            Point_3 me = face->center();
            Delaunay::Face_handle north;
            Delaunay::Face_handle south;
            Delaunay::Face_handle east;
            Delaunay::Face_handle west;

            Delaunay::Face_handle northeast;
            Delaunay::Face_handle northwest;
            Delaunay::Face_handle southeast;
            Delaunay::Face_handle southwest;


            double distance = 100;//300.0;


            north = domain->locate_face(me.x(), me.y() + distance);
            if (north == nullptr )// || domain->is_infinite(north))
                north = face;

            south = domain->locate_face(me.x(), me.y() - distance);
            if (south == nullptr)// || domain->is_infinite(south))
                south = face;

            west = domain->locate_face(me.x() - distance, me.y());
            if (west == nullptr)//  || domain->is_infinite(west))
                west = face;

            east = domain->locate_face(me.x() + distance, me.y());
            if (east == nullptr)//  || domain->is_infinite(east))
                east = face;

            double z = face->get_z();
            double zw = west->get_z();
            double ze = east->get_z();
            double zs = south->get_z();
            double zn = north->get_z();

            double znw = 0.;
            double zne = 0.;
            double zse = 0.;
            double zsw = 0.;

            northeast = domain->locate_face(me.x() + distance, me.y() + distance);
            if (northeast == nullptr )//|| domain->is_infinite(northeast))
                zne = (ze + zn) / 2.;
            else
                zne = northeast->get_z();

            northwest = domain->locate_face(me.x() - distance, me.y() + distance);
            if (northwest == nullptr)// || domain->is_infinite(south))
                znw = (zw + zn) / 2.;
            else
                znw = northwest->get_z();

            southeast = domain->locate_face(me.x() + distance, me.y() - distance);
            if (southeast == nullptr)//|| domain->is_infinite(west))
                zse = (ze + zs) / 2.;
            else
                zse = southeast->get_z();

            southwest = domain->locate_face(me.x() - distance, me.y() - distance);
            if (southwest == nullptr)// || domain->is_infinite(east))
                zsw = (zw + zs) / 2.;
            else
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
        c->curvature = value;
        face->set_parameter("Liston curvature", value);
    }


    if ( cfg.get("serialize",true) )
    {
        LOG_DEBUG << "Serializing liston curvature";
        domain->serialize_parameter(cfg.get("serialize_output", "liston_curvature.mesh"),
                                    "Liston curvature");
    }

}


void Liston_wind::run(mesh domain, boost::shared_ptr<global> global_param)
{
    double PI = 3.14159;

    std::vector< boost::tuple<double, double, double> > u;
    std::vector< boost::tuple<double, double, double> > v;
    for (auto& s : global_param->stations)
    {
        if( is_nan(s->get("u")) || is_nan(s->get("vw_dir")))
            continue;

        double W = s->get("u");

        double theta = s->get("vw_dir") * 3.14159/180.;
        double zonal_u = -W * sin(theta);
        double zonal_v = -W * cos(theta);
        u.push_back(boost::make_tuple(s->x(), s->y(), zonal_u ) );
        v.push_back(boost::make_tuple(s->x(), s->y(), zonal_v ) );
    }

    // omega_s needs to be scaled on [-0.5,0.5]
    double max_omega_s = -99999.0;


    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto elem = domain->face(i);
        auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
        double zonal_u = elem->get_module_data<lwinddata>(ID)->interp(u, query);
        double zonal_v = elem->get_module_data<lwinddata>(ID)->interp(v, query);

        double theta = 3.0 * PI * 0.5 - atan2(zonal_v , zonal_u);

        if (theta > 2*PI)
            theta = theta - 2*PI;

        //eqn 15
        double omega_s = elem->slope() * cos(theta - elem->aspect());

        if( fabs(omega_s) > max_omega_s)
            max_omega_s = fabs(omega_s);
    }

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto elem = domain->face(i);
        auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
        double zonal_u = elem->get_module_data<lwinddata>(ID)->interp(u, query);
        double zonal_v = elem->get_module_data<lwinddata>(ID)->interp(v, query);

        double W = sqrt(zonal_u * zonal_u + zonal_v * zonal_v);
        double corrected_theta = 3.0 * PI * 0.5 - atan2(zonal_v , zonal_u);

        if (corrected_theta > 2*PI)
            corrected_theta = corrected_theta - 2*PI;


        double omega_s = elem->slope() * cos(corrected_theta - elem->aspect());

        omega_s = omega_s / max_omega_s / 2.0;

        double omega_c = elem->get_parameter("Liston curvature");//elem->get_module_data<lwinddata>(ID)->curvature;

        double ys = 0.5;
        double yc = 0.5;

        double Ww = 1 + ys * omega_s + yc * omega_c;

        W = W * Ww;

        double theta_d = -0.5 * omega_s * sin(2 * (elem->aspect() - corrected_theta));
        corrected_theta = theta_d + corrected_theta;

        elem->set_face_data("vw", W);
        elem->set_face_data("vw_dir", corrected_theta * 180.0 / 3.14159);
    }


}

Liston_wind::~Liston_wind()
{

}