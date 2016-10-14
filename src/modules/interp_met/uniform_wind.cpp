#include "uniform_wind.hpp"

uniform_wind::uniform_wind(config_file cfg)
        :module_base(parallel::domain)

{
    depends_from_met("U_R");
    depends_from_met("vw_dir");

    provides("U_R");
    provides("vw_dir");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

//Calculates the curvature required
void uniform_wind::init(mesh domain)
{

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<lwinddata>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
        face->coloured = false;
    }

}


void uniform_wind::run(mesh domain)
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
            double zonal_u = -W * sin(theta);
            double zonal_v = -W * cos(theta);
            u.push_back(boost::make_tuple(s->x(), s->y(), zonal_u));
            v.push_back(boost::make_tuple(s->x(), s->y(), zonal_v));
        }

        // Interp over stations
        auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
        double zonal_u = face->get_module_data<lwinddata>(ID)->interp(u, query);
        double zonal_v = face->get_module_data<lwinddata>(ID)->interp(v, query);

	// Convert back to direction and magnitude
        double theta = 3.0 * PI * 0.5 - atan2(zonal_v, zonal_u);

        if (theta > 2 * PI)
            theta = theta - 2 * PI;

        double W = sqrt(zonal_u * zonal_u + zonal_v * zonal_v);
        double corrected_theta = 3.0 * PI * 0.5 - atan2(zonal_v, zonal_u);

        if (corrected_theta > 2 * PI)
            corrected_theta = corrected_theta - 2 * PI;

        face->get_module_data<lwinddata>(ID)->corrected_theta = corrected_theta;
        face->get_module_data<lwinddata>(ID)->W = W;
    }
    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        double corrected_theta= face->get_module_data<lwinddata>(ID)->corrected_theta;
        double W= face->get_module_data<lwinddata>(ID)->W;

        W = std::max(W,0.1);
        face->set_face_data("U_R", W);
        face->set_face_data("vw_dir", corrected_theta * 180.0 / 3.14159);
    }


}

uniform_wind::~uniform_wind()
{

}
