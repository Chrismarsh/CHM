#pragma once
#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"
#include <gsl/gsl_fit.h>
#include <meteoio/MeteoIO.h>
class rh_from_obs : public module_base
{
public:
    rh_from_obs();
    ~rh_from_obs();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);
};
