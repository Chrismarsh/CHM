#include "Liston_monthly_llra_ta.hpp"

Liston_monthly_llra_ta::Liston_monthly_llra_ta(config_file cfg)
        :module_base(parallel::data)

{
    provides("t");
    provides("t_lapse_rate");

    depends_from_met("t");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

Liston_monthly_llra_ta::~Liston_monthly_llra_ta()
{

}
void Liston_monthly_llra_ta::init(mesh domain, boost::shared_ptr<global> global_param)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
    }
}
void Liston_monthly_llra_ta::run(mesh_elem& face, boost::shared_ptr<global> global_param)
{

    double lapse_rate = -9999;

    switch(global_param->month())
    {
        case 1:
            lapse_rate=0.0044;
            break;
        case 2:
            lapse_rate=0.0059;
            break;
        case 3:
            lapse_rate=0.0071;
            break;
        case 4:
            lapse_rate=0.0078;
            break;
        case 5:
            lapse_rate=0.0081;
            break;
        case 6:
            lapse_rate=0.0082;
            break;
        case 7:
            lapse_rate=0.0081;
            break;
        case 8:
            lapse_rate=0.0081;
            break;
        case 9:
            lapse_rate=0.0077;
            break;
        case 10:
            lapse_rate=0.0068;
            break;
        case 11:
            lapse_rate=0.0055;
            break;
        case 12:
            lapse_rate=0.0047;
            break;

    }


    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : global_param->get_stations( face->get_x(), face->get_y()))
    {
        if( is_nan(s->get("t")))
            continue;
        double v = s->get("t") - lapse_rate * (0.0 - s->z());
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), v ) );
    }


    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double value = face->get_module_data<data>(ID)->interp(lowered_values, query);

    //raise value back up to the face's elevation from sea level
    value =  value + lapse_rate * (0.0 - face->get_z());

    face->set_face_data("t",value);

    face->set_face_data("t_lapse_rate",lapse_rate);

}