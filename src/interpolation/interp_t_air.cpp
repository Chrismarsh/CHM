#include "interp_t_air.hpp"

interp_t_air::interp_t_air()
{
    
}
interp_t_air::~interp_t_air()
{

}

void interp_t_air::operator()(std::string method, mesh_elem& m, station_list& stations, boost::shared_ptr<global> global_param)
{
    boost::shared_ptr<interp_visitor> visitor;
    if (method == "LLRA_const")
    {
        visitor = boost::make_shared<LLRA_const>();
    } 
    else if (method == "LLRA_var")
    {
        visitor = boost::make_shared<LLRA_var>();
    }
    else
    {
        BOOST_THROW_EXCEPTION(interpolation_error()
                << errstr_info("Unknown visitor requested: " + method));
    }
    
    interp_2d interp;
    double tair = interp("spline",m,stations,visitor,global_param);
    m.add_face_data(TAIR,tair);

}


//********************************************
//  Algorithms
///*******************************************

double LLRA_const::lower(mesh_elem& m, boost::shared_ptr<station>  s, boost::shared_ptr<global> global_param)
{

    double lapse_rate = 0.0065;
    double v = s->now().get(TAIR) - lapse_rate * (0.0 - s->get_elevation());
   
   return v;
    
}

double LLRA_const::raise(double value, mesh_elem& m, boost::shared_ptr<global> global_param)
{
    double lapse_rate = 0.0065;

    double v =  value + lapse_rate * (0.0 - m.get_z());
   return v;

}

double LLRA_var::get_lapse_rate(int month)
{
    double lapse_rate;
    switch(month)
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
    
    return lapse_rate;
}
double LLRA_var::lower(mesh_elem& m, boost::shared_ptr<station>  s, boost::shared_ptr<global> global_param)
{

//    LOG_DEBUG << "Month: " << _month;
    double lapse_rate = get_lapse_rate(global_param->month());
    double v = s->now().get(TAIR) - lapse_rate * (0.0 - s->get_elevation());
    
    return v;

}

double LLRA_var::raise(double value, mesh_elem& m, boost::shared_ptr<global> global_param)
{
    double lapse_rate = get_lapse_rate(global_param->month());

    double v =  value + lapse_rate * (0.0 - m.get_z());
    return v;

}


