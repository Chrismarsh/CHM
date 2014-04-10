#include "inv_dist.hpp"

inv_dist::inv_dist()
{
    
}

inv_dist::~inv_dist()
{
    
}

double inv_dist::operator()(station_list const&  stations, mesh_elem& elem, std::string variable)
{
    
    double numerator = 0.0;
    double denominator = 0.0;

    double z0 = 0;


    for(auto& itr : stations)
    {
        
        double z = itr->now().get<double>(variable);
        double xdiff = ( itr->get_x() - elem.get_x());
        double ydiff = ( itr->get_y() - elem.get_y());
        double di = pow(sqrt(
                             pow(xdiff,2.0) + pow(ydiff,2.0)
                            )
                        ,2.0
                    );
        if(di <= 0)
        {
                numerator = z;
                denominator = 1.0;
                break;
        }
        else
        {
                numerator += z/ di;
                denominator += 1 / di; 

        }
         z0 += (numerator/denominator);
    }

   
    z0 /= stations.size();
    
    return z0;

}
