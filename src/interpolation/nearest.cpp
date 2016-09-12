#include "nearest.hpp"

nearest::nearest()
{
    
}

nearest::~nearest()
{
    
}

double nearest::operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point)
{

    if (sample_points.size() > 1)
    {
        BOOST_THROW_EXCEPTION( interpolation_error()
                                << errstr_info("nearest requires exactly 1 station"));
    }

    return sample_points.at(0).get<2>();

}
