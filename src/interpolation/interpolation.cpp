
#include "interpolation.hpp"

double interpolation::operator()(std::string method, station_list const&  stations, mesh_elem& elem, std::string variable)
{
    boost::shared_ptr<interpolation_base> interp;
            
    if(method == "idw")
    {
        interp =  boost::make_shared<inv_dist>();
    }
    else if (method == "spline")
    {
        interp =  boost::make_shared<spline>();        
    }
    else 
    {
        BOOST_THROW_EXCEPTION(interp_unknown_type() << errstr_info("Unknow type '" + method + "'"));
    }
    
    return (*interp)(stations,elem,variable);
    
}


interpolation::interpolation()
{
    
}
interpolation::~interpolation()
{
    
}