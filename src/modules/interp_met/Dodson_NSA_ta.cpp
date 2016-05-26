#include "Dodson_NSA_ta.h"

Dodson_NSA_ta::Dodson_NSA_ta(config_file cfg)
:module_base(parallel::data)
{

    depends_from_met("t");
    provides("t");
    provides("t_lapse_rate");
}
Dodson_NSA_ta::~Dodson_NSA_ta()
{

}
void Dodson_NSA_ta::init(mesh domain, boost::shared_ptr<global> global_param)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->stations.size());
    }
}
void Dodson_NSA_ta::run(mesh_elem &elem, boost::shared_ptr<global> global_param)
{

    double Po = 100000.0; //sea level pressure (Pa)
    double Cp = 1005.0; //specific heat of dry air J/[KgK]
    double R = 8.3143; //gas constant J/(molK)
    double m = 0.02897; //molecular weight of dry air kg/mol
    double Tb = 300.0; //assumed sea level temperature 300K
    double g = 9.810616; //gravity m/s/s
    double lapse = 0.0065; //K/m

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;


    for (auto& s : global_param->stations)
    {
        if( is_nan(s->get("t")))
            continue;

        double elev = s->z(); //station elevation
        //pressure at station's elevation
        double Pz = Po * pow(Tb/(Tb+lapse*elev),(m*g)/(lapse*R));
        double ta = s->get("t") + 273.15; //to K

        //calculate virtual temp (eqn 1)
        double ratio = (Po/Pz);
        double exp = R/(m*Cp);
        double theta = ta * pow(ratio,exp);

        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), theta ) );
    }

    auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());

    //interpolated virtual temp, now go back to station
    double theta = elem->get_module_data<data>(ID)->interp(lowered_values, query);
    double elev = elem->get_z();
    double Pz = Po * pow(Tb/(Tb+ (-lapse)*elev),(m*g)/((-lapse)*R));
    double ratio = (Po/Pz);
    double exp = R/(m*Cp);

    double Ta = ( theta/pow(ratio,exp) );
    Ta -= 273.15;

    elem->set_face_data("t",Ta);
    elem->set_face_data("t_lapse_rate",lapse);
}