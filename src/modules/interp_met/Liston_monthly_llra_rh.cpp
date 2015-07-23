#include "Liston_monthly_llra_rh.hpp"


Liston_monthly_llra_rh::Liston_monthly_llra_rh( std::string ID)
{
    _provides->push_back("rh");
    _provides->push_back("ea");
    _provides->push_back("es");
    _provides->push_back("Td_lapse_rate");
    _depends->push_back("t");

    _depends_from_met->push_back("rh");

    this->ID = ID;
    _parallel_type = parallel::data;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

Liston_monthly_llra_rh::~Liston_monthly_llra_rh()
{

}
void Liston_monthly_llra_rh::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{


    double lambda = 0.0;

    switch(global_param->month())
    {
        case 1:
            lambda = 0.00041;
            break;
        case 2:
            lambda = 0.00042;
            break;
        case 3:
            lambda = 0.00040;
            break;
        case 4:
            lambda = 0.00039;
            break;
        case 5:
            lambda = 0.00038;
            break;
        case 6:
            lambda = 0.00036;
            break;
        case 7:
            lambda = 0.00033;
            break;
        case 8:
            lambda = 0.00033;
            break;
        case 9:
            lambda = 0.00036;
            break;
        case 10:
            lambda = 0.00037;
            break;
        case 11:
            lambda = 0.00040;
            break;
        case 12:
            lambda = 0.00040;
            break;

    }

    if(lambda == 0.0)
        BOOST_THROW_EXCEPTION(interp_error() << errstr_info("rh lapse rate == 0"));

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : global_param->stations)
    {
        //use water for the moment as per Liston, Elder
        double a = 611.21;
        double b = 17.502;
        double c = 240.97;
        double temp = s->get(global_param->get_variable("Tair"));
        double rh = s->get(global_param->get_variable("RH"));

        //because boom otherwise
        if (rh <= 0.0)
        {
            rh = 1.0;
        }
        //solve RH ~= 100* e/es for e
        double es = a * exp((b * temp) / (c + temp));
        double e = rh / 100.0 * es;

        double dewPointTemp = (c * log(e / a)) / (b - log(e / a));

        double dewPointLapseRate = lambda * c / b;

        double newTd = dewPointTemp - dewPointLapseRate * (0.0 - s->z());

        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), newTd ) );
    }

    interp_base* interp=nullptr;
    std::string interp_method = "spline";
    if(interp_method == "spline")
        interp = new thin_plate_spline();

    auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
    double value = (*interp)(lowered_values, query);

    //raise value back up to the face's elevation from sea level
    //use water for the moment as per Liston, Elder
    double a = 611.21;
    double b = 17.502;
    double c = 240.97;
    double dewPointLapseRate = lambda * c / b;
    double Td = value - (-dewPointLapseRate)*(0 - elem->get_z());
    double e = a * exp((b * Td) / (c + Td));
    double temp = elem->face_data("t");
    double es = a * exp((b * temp) / (c + temp));

    //RH value replaces Tdew value
    double rh = 100.0 * e / es;
    if (rh > 100.0)
        rh = 100.0;

    elem->set_face_data(global_param->get_variable("RH"),rh);
    elem->set_face_data("es",es);
    elem->set_face_data("ea",e);
    elem->set_face_data("Td_lapse_rate",dewPointLapseRate);

    delete interp;
}