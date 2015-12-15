#include "interpolation.h"

interpolation::interpolation(interp_alg ia, size_t size)
{
    init(ia,size);
}

void interpolation::init(interp_alg ia, size_t size)
{
    switch(ia)
    {
        case interp_alg::tpspline:
            if(size)
                base = boost::movelib::make_unique<thin_plate_spline>(size);
            else
                base = boost::movelib::make_unique<thin_plate_spline>();
        break;
    };
}

interpolation::interpolation()
{
    base = nullptr;
}
interpolation::~interpolation()
{

}

double interpolation::operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point)
{
    base->operator()(sample_points,query_point);
}
