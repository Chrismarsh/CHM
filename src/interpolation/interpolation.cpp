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
//            if(size)
//                base = boost::movelib::make_unique<thin_plate_spline>(size);
//            else
//                base = boost::movelib::make_unique<thin_plate_spline>();
            if(size)
                base = boost::make_shared<thin_plate_spline>(size);
            else
                base = boost::make_shared<thin_plate_spline>();
        break;
    };

    this->size = size;
    this->ia = ia;
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
//    std::vector< boost::tuple<double,double,double> > isfinite_sample_points;
//    for(auto& itr: sample_points)
//    {
//        if (itr.get<2>() != -9999.)
//        {
//            isfinite_sample_points.push_back(itr);
//        }
//    }
    if (sample_points.size() == 0)
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info("Interpolation sample point length = 0."));
    }

    if(sample_points.size() != this->size)
    {

        this->init(ia,sample_points.size());
    }

    base->operator()(sample_points,query_point);
}
