
#include "interp_2d.hpp"

double interp_2d::operator()(std::string method, mesh_elem& elem, station_list&  stations,  boost::shared_ptr<interp_visitor> visitor )
{
    boost::shared_ptr<interp_alg_base> interp;
            
    if(method == "idw")
    {
        interp =  boost::make_shared<inv_dist>();
    }
//    else if (method == "spline")
//    {
//        interp =  boost::make_shared<spline>();        
//    }
    else 
    {
        BOOST_THROW_EXCEPTION(interp_unknown_type() << errstr_info("Unknown type '" + method + "'"));
    }
    
    return (*interp)(stations,elem,visitor);
    
}


interp_2d::interp_2d()
{
    
}
interp_2d::~interp_2d()
{
    
}