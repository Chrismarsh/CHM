#include "Liston_monthly_llra_rh.hpp"


Liston_monthly_llra_rh::Liston_monthly_llra_rh()
        :module_base(parallel::data)

{
    provides("rh");
    provides("Td_lapse_rate");

    depends("t");
    depends_from_met("rh");


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

Liston_monthly_llra_rh::~Liston_monthly_llra_rh()
{

}
void Liston_monthly_llra_rh::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{
    size_t ID = elem->_debug_ID;
    // 1/m
    double lambda_table[] = {0.00041,
                           0.00042,
                           0.00040,
                           0.00039,
                           0.00038,
                           0.00036,
                           0.00033,
                           0.00033,
                           0.00036,
                           0.00037,
                           0.00040,
                           0.00040};
    double lambda = lambda_table[ global_param->month() - 1 ];

    //taken from mio
    const double Aw = 611.21, Bw = 17.502, Cw = 240.97; //parameters for water
    const double Ai = 611.15, Bi = 22.452, Ci = 272.55; //parameters for ice

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : global_param->stations)
    {

        double t = s->get("t")+273.15;
        double rh = s->get("rh")/100.;

        double dewPointTemp = mio::Atmosphere::RhtoDewPoint(rh,t,false); // K
        double f = t < 273.15 ?  Ci / Bi : Cw/Bw; // T<0 -> use w.r.t ice
        double dewPointLapseRate = lambda * f;

        double newTd = dewPointTemp - dewPointLapseRate * (s->z() - 0.0);

        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), newTd ) );
    }

    interp_base* interp=nullptr;
    std::string interp_method = "spline";
    if(interp_method == "spline")
        interp = new thin_plate_spline();

    auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
    double value = (*interp)(lowered_values, query);//K

    //raise value back up to the face's elevation from sea level
    //use water for the moment as per Liston, Elder
    double t = elem->face_data("t")+273.15;
    double f = t < 273.15 ?  Ci / Bi : Cw/Bw; // T<0 -> use w.r.t ice
    double dewPointLapseRate = lambda * f;

    double Td = value - dewPointLapseRate*(0.0 - elem->get_z());
    double rh = mio::Atmosphere::DewPointtoRh(Td,t,false);

    elem->set_face_data("rh",rh*100.0);
    elem->set_face_data("Td_lapse_rate",dewPointLapseRate);

    delete interp;
}