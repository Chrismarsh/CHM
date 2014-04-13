#include "inv_dist.hpp"

inv_dist::inv_dist()
{
    
}

inv_dist::~inv_dist()
{
    
}

double inv_dist::operator()(station_list&  stations, mesh_elem& elem, std::string variable, boost::shared_ptr<interp_visitor> visitor)
{
    
    double numerator = 0.0;
    double denominator = 0.0;

    double z0 = 0;

    if (stations.size() == 1)
    {
        BOOST_THROW_EXCEPTION( interpolation_error()
                                << errstr_info("IDW requires >=2 stations"));
    }
    for(auto& itr : stations)
    {
        
        double z = visitor->lower(elem,variable,itr);
        
//        itr->now()->get<double>(variable);

        double sx = itr->get_x();
        double sy = itr->get_y();
        
        double ex = elem.get_x();
        double ey = elem.get_y();
        
        double xdiff = (sx  - ex);
        double ydiff = (sy  - ey);
        double di = pow(sqrt(
                             pow(xdiff,2.0) + pow(ydiff,2.0)
                            )
                        ,2.0
                    );
        if(di == 0)
        {
                numerator = z;
                denominator = 1.0;
                //break;
        }
        else
        {
                numerator += z/ di;
                denominator += 1 / di; 

        }
        
        
    }
   z0 = (numerator/denominator);
   z0 = visitor->raise(z0,elem,variable);
    
    return z0;

}
