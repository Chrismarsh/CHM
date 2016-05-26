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

    virtual void run(mesh_elem &elem, boost::shared_ptr <global> global_param);

    virtual void run(mesh domain, boost::shared_ptr <global> global_param);

    virtual void init(mesh domain, boost::shared_ptr<global> global);

};



