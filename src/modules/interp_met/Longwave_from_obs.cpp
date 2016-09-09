#include "Longwave_from_obs.hpp"

Longwave_from_obs::Longwave_from_obs(config_file cfg)
        :module_base(parallel::data)

{
    provides("ilwr");
    //provides("lw_lapse_rate");

    depends_from_met("Qli");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

Longwave_from_obs::~Longwave_from_obs()
{

}
void Longwave_from_obs::init(mesh domain, boost::shared_ptr<global> global_param)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
    }
}
void Longwave_from_obs::run(mesh_elem& face, boost::shared_ptr<global> global_param)
{

    // Use constant annual LW lapse rate based on emperical study in Alps
    double lapse_rate = 2.8/100; // 2.8 W/m^2 / 100 meters (Marty et al. 2002)

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : global_param->get_stations( face->get_x(), face->get_y()))
    {
        if( is_nan(s->get("Qli")))
            continue;
        double v = s->get("Qli") - lapse_rate * (0.0 - s->z());
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), v ) );
    }


    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double value = face->get_module_data<data>(ID)->interp(lowered_values, query);

    //raise value back up to the face's elevation from sea level
    value =  value + lapse_rate * (0.0 - face->get_z());

    face->set_face_data("ilwr",value);

    //face->set_face_data("lw_lapse_rate",lapse_rate);

}
