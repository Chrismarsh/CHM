#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"
#include <meteoio/MeteoIO.h>

#include "sno.h"
#include "snomacros.h"
class snodata : public face_info
{
public:
    sno data;
    double sum_runoff;
    double sum_melt;
    int dead;

};
class snobal : public module_base
{
public:
    snobal(config_file cfg);

    ~snobal();

    double drift_density; // if we have blowing snow, this is the density of those particles

    virtual void run(mesh_elem &face);
    virtual void init(mesh domain);

};



