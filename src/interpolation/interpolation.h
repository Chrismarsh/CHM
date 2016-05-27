#pragma once

#include "interp_base.hpp"
#include "inv_dist.hpp"
#include "TPSpline.hpp"

#include <vector>
#include <boost/tuple/tuple.hpp>
//#include <boost/move/unique_ptr.hpp>
//#include <boost/move/make_unique.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
enum interp_alg
{
    tpspline,
    idw
};

class interpolation
{
public:
    /*
     * Some of the interpolators need to allocate memory, and if called multiple times, it is faster to preinit and set size
     * as required where size is the number of items used to interpolate from. E.g., # stations.
     */
    interpolation(interp_alg ia, size_t size=0);
    interpolation();
    ~interpolation();

    void init(interp_alg ia, size_t size=0);

    double operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point);

private:
//    boost::movelib::unique_ptr<interp_base> base;
    boost::shared_ptr<interp_base> base;
    size_t size;
    interp_alg ia;
};
