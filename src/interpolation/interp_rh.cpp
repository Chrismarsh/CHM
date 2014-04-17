#include "interp_rh.hpp"

interp_rh::interp_rh()
{
    
}
interp_rh::~interp_rh()
{

}

void interp_rh::operator()(std::string method, mesh_elem& m, station_list& stations)
{
    boost::shared_ptr<interp_visitor> visitor;
    if (method == "LLRA_rh_var")
    {
        visitor = boost::make_shared<LLRA_rh_var>();
    } 
    else
    {
        BOOST_THROW_EXCEPTION(interpolation_error()
                << errstr_info("Unknown visitor requested: " + method));
    }
    
    interp_2d interp;
    double rh = interp("idw",m,stations,visitor);
    m.add_face_data("Rh",rh); //TODO: fix hardcoded variable string

}


//********************************************
//  Algorithms
///*******************************************


double LLRA_rh_var::get_lambda_rate(int month)
{
    double lambda;
    switch (month)
    {
        case 1:
            lambda = 0.00041F;
            break;
        case 2:
            lambda = 0.00042F;
            break;
        case 3:
            lambda = 0.00040F;
            break;
        case 4:
            lambda = 0.00039F;
            break;
        case 5:
            lambda = 0.00038F;
            break;
        case 6:
            lambda = 0.00036F;
            break;
        case 7:
            lambda = 0.00033F;
            break;
        case 8:
            lambda = 0.00033F;
            break;
        case 9:
            lambda = 0.00036F;
            break;
        case 10:
            lambda = 0.00037F;
            break;
        case 11:
            lambda = 0.00040F;
            break;
        case 12:
            lambda = 0.00040F;
            break;

    }

    return lambda;
}

double LLRA_rh_var::lower(mesh_elem& m, boost::shared_ptr<station> s)
{
    //use water for the moment as per Liston, Elder
    double a = 611.21;
    double b = 17.502;
    double c = 240.97;
    double temp = s->now().get(TAIR);
    double rh = s->now().get(RH);
    
    //because boom otherwise
    if (rh <= 0.0)
    {
        rh = 1.0;
    }
    //solve RH ~= 100* e/es for e
    double es = a * exp((b * temp) / (c + temp));
    double e = rh / 100 * es;

    double dewPointTemp = (c * log(e / a)) / (b - log(e / a));

    _month = s->now().month();
    double lambda = get_lambda_rate(_month);

    double dewPointLapseRate = lambda * c / b;

    double newTd = dewPointTemp - dewPointLapseRate * (0.0 - m.get_z());
    return newTd;

}

double LLRA_rh_var::raise(double value, mesh_elem& m)
{
    
    //use water for the moment as per Liston, Elder
    double a = 611.21;
    double b = 17.502;
    double c = 240.97;
    double lambda = get_lambda_rate(_month);
    double dewPointLapseRate = lambda * c / b;
    double Td = value - (-dewPointLapseRate)*(0 - m.get_z());
    double e = a * exp((b * Td) / (c + Td));
    double temp = m.get_face_data(TAIR);
    double es = a * exp((b * temp) / (c + temp));

    //RH value replaces Tdew value
    double rh = 100 * e / es;
    if (rh > 100)
        rh = 100;
    
    return rh;
}