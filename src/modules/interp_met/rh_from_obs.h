#pragma once
#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include <gsl/gsl_fit.h>
#include <meteoio/MeteoIO.h>
class rh_from_obs : public module_base
{
public:
    rh_from_obs(config_file cfg);
    ~rh_from_obs();
    virtual void run(mesh_elem& face);
    virtual void init(mesh domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};
