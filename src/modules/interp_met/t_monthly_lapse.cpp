#include "t_monthly_lapse.hpp"

t_monthly_lapse::t_monthly_lapse(config_file cfg)
        :module_base(parallel::data)

{
    provides("t");
    provides("t_lapse_rate");

    depends_from_met("t");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

t_monthly_lapse::~t_monthly_lapse()
{

}
void t_monthly_lapse::init(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
    }

    MLR[0]=cfg.get("MLR_1",0.0049);
    MLR[1]=cfg.get("MLR_2",0.0049);
    MLR[2]=cfg.get("MLR_3",0.0060);
    MLR[3]=cfg.get("MLR_4",0.0060);
    MLR[4]=cfg.get("MLR_5",0.0060);
    MLR[5]=cfg.get("MLR_6",0.0053);
    MLR[6]=cfg.get("MLR_7",0.0053);
    MLR[7]=cfg.get("MLR_8",0.0053);
    MLR[8]=cfg.get("MLR_9",0.0046);
    MLR[9]=cfg.get("MLR_10",0.0046);
    MLR[10]=cfg.get("MLR_11",0.0049);
    MLR[11]=cfg.get("MLR_12",0.0049);

}
void t_monthly_lapse::run(mesh_elem& face)
{

    double lapse_rate = MLR[global_param->month()-1];

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