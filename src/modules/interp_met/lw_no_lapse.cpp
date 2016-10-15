#include "lw_no_lapse.hpp"

lw_no_lapse::lw_no_lapse(config_file cfg)
        :module_base(parallel::data)

{
    provides("ilwr");
    //provides("lw_lapse_rate");

    depends_from_met("Qli");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

lw_no_lapse::~lw_no_lapse()
{

}
void lw_no_lapse::init(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<lw_no_lapse::data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
    }}
void lw_no_lapse::run(mesh_elem& face)
{

    // Use no lapse rate
    double lapse_rate = 0.0;

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

}
