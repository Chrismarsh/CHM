#include "Harder_precip_phase.hpp"
Harder_precip_phase::Harder_precip_phase(config_file cfg)
        :module_base(parallel::data)
{
    depends("t");
    depends("rh");
    depends("p");

    provides("Ti");
    provides("frac_precip_rain");
    provides("frac_precip_snow");
    provides("p_snow");
    provides("p_rain");

    //     default values:
    //     b=2.630006;
    //     c=0.09336;
    b = cfg.get("const.b",2.630006);
    c = cfg.get("const.c",0.09336);

    LOG_DEBUG << "Successfully instantiated module " << this->ID;


}
Harder_precip_phase::~Harder_precip_phase()
{

}
void Harder_precip_phase::run(mesh_elem& face, boost::shared_ptr<global> global_param)
{

    double Ta = face->face_data("t")+273.15; //K
    double T =  face->face_data("t");
    double RH = face->face_data("rh");
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

    double frTi = 1 / (1+b*pow(c,Ti));

    frTi = std::trunc(100.0*frTi) / 100.0; //truncate to 2 decimal positions

    face->set_face_data("Ti",Ti);
    face->set_face_data("frac_precip_rain",frTi);
    face->set_face_data("frac_precip_snow",1-frTi);

    double p = face->face_data("p");

    face->set_face_data("p_rain", p * frTi);
    face->set_face_data("p_snow", p * (1-frTi));

}