#include "Kunkel_monthlyTd_rh.hpp"


Kunkel_monthlyTd_rh::Kunkel_monthlyTd_rh()
        :module_base(parallel::data)

{
    provides("rh");
    provides("Td_lapse_rate");

    depends("t");
    depends_from_met("rh");


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

Kunkel_monthlyTd_rh::~Kunkel_monthlyTd_rh()
{

}
void Kunkel_monthlyTd_rh::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{
//    size_t ID = elem->_debug_ID;
    // 1/km
    double lapse_rates[] = {
            0.41,
            0.42,
            0.40,
            0.39,
            0.38,
            0.36,
            0.33,
            0.33,
            0.36,
            0.37,
            0.4,
            0.4
    } ;

    double lapse = lapse_rates[ global_param->month() - 1 ] / 1000.; // -> 1/m

    //taken from mio
    const double  Bw = 17.502, Cw = 240.97; //parameters for water
    const double Bi = 22.452, Ci = 272.55; //parameters for ice

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : global_param->stations)
    {

        double t = s->get("t")+273.15;
        double rh = s->get("rh")/100.;

        double Tdz0 = mio::Atmosphere::RhtoDewPoint(rh,t,false) - 273.15; // K

        double C = t < 273.15 ? Ci : Cw;
        double B = t < 273.15 ? Bi : Bw;

        double z = 0.;
        double z0 = elem->get_z();
//        double am = lapse;
        double Td_z = -lapse*C*(z-z0) / B + Tdz0;
//        double Td_z = (-am*(z-z0)*(C+Tdz0)/B+Tdz0)/(1+am*(z-z0)*(C+Tdz0)/(B*C));
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), Td_z ) );
    }

    interp_base* interp=nullptr;
    std::string interp_method = "spline";
    if(interp_method == "spline")
        interp = new thin_plate_spline();

    auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
    double Tdz0 = (*interp)(lowered_values, query);//C

    //raise value back up to the face's elevation from sea level
    double t = elem->face_data("t") + 273.15;
    double C = t < 273.15 ? Ci : Cw;
    double B = t < 273.15 ? Bi : Bw;

    double z0 = 0.;
    double z = elem->get_z();
//    double am = lapse;
    double Td_z = -lapse*C*(z-z0) / B + Tdz0;
//    double Td_z = (-am*(z-z0)*(C+Tdz0)/B+Tdz0)/(1+am*(z-z0)*(C+Tdz0)/(B*C));

    double rh = mio::Atmosphere::DewPointtoRh(Td_z+273.15,t,false);

    elem->set_face_data("rh",rh*100.0);
    elem->set_face_data("Td_lapse_rate",lapse);

    delete interp;
}