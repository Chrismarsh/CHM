#include "inv_dist.hpp"

inv_dist::inv_dist()
{
    
}

inv_dist::~inv_dist()
{
    
}

double inv_dist::operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_points)
{
    
    double numerator = 0.0;
    double denominator = 0.0;

    double z0 = 0;
    if (sample_points.size() == 1)
    {
        BOOST_THROW_EXCEPTION( interpolation_error()
                                << errstr_info("IDW requires >=2 stations"));
    }
    for(size_t i=0;i<sample_points.size();i++)
    {
        double z = sample_points.at(i).get<2>();

        double sx = sample_points.at(i).get<0>();
        double sy = sample_points.at(i).get<1>();
        
        double ex = query_points.get<0>();
        double ey = query_points.get<1>();
        
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
        }
        else
        {
                numerator += z/ di;
                denominator += 1 / di; 
        }
    }

   z0 = (numerator/denominator);

   return z0;

}
