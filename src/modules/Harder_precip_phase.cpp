#include "Harder_precip_phase.hpp"
Harder_precip_phase::Harder_precip_phase(std::string ID)
{
    _depends->push_back("t");
    _depends->push_back("rh");
    _depends->push_back("ea");

    _provides->push_back("p_snow");
    _provides->push_back("p_rain");
    _provides->push_back("Ti");
    _provides->push_back("frac_precip_rain");
    _provides->push_back("frac_precip_snow");

    this->ID = ID;
    _parallel_type = parallel::data;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}
Harder_precip_phase::~Harder_precip_phase()
{

}
void Harder_precip_phase::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{

    double Ta = elem->face_data("t")+273.15; //K
    double T =  elem->face_data("t");
    double RH = elem->face_data("rh");
    //double ea = elem->face_data("ea")*0.001; // Pa -> kPa
    double ea = RH/100 * 0.611*exp( (17.3*T) / (237.3+T));

    // (A.6)
    double D = 2.06 * pow(10,-5) * pow(Ta/273.15,1.75);

    // (A.9)
    double lambda_t = 0.000063 * Ta + 0.00673;

    // (A.10) (A.11)
    double L;
    if(T < 0.0)
    {
        L = 1000.0 * (2834.1 - 0.29 *T - 0.004*T*T);
    }
    else
    {
        L = 1000.0 * (2501.0 - (2.361 * T));
    }

    /*
     * The *1000 and /1000 are important unit conversions. Doesn't quite match the harder paper, but Phil assures me it is correct.
     */
    double mw = 0.01801528 * 1000.0; //[kgmol-1]
    double R = 8.31441 /1000.0; // [J mol-1 K-1]
    double rho = (mw * ea) / (R*Ta);


    auto fx = [=](double Ti)
    {
        return boost::math::make_tuple(
                T+D*L*(rho/(1000.0)-.611*mw*exp(17.3*Ti/(237.3+Ti))/(R*(Ti+273.15)*(1000.0)))/lambda_t-Ti,
                D*L*(-0.6110000000e-3*mw*(17.3/(237.3+Ti)-17.3*Ti/pow(237.3+Ti,2))*exp(17.3*Ti/(237.3+Ti))/(R*(Ti+273.15))+0.6110000000e-3*mw*exp(17.3*Ti/(237.3+Ti))/(R*pow(Ti+273.15,2)))/lambda_t-1);
    };

    double guess = T;
    double min = -50;
    double max = 50;
    double digits = 6;

    double Ti = boost::math::tools::newton_raphson_iterate(fx, guess, min, max, digits);

    double b=2.630006;
    double c=0.09336;
    double frTi = 1 / (1+b*pow(c,Ti));

    elem->set_face_data("Ti",Ti);
    elem->set_face_data("frac_precip_rain",frTi);
    elem->set_face_data("frac_precip_snow",1-frTi);

}