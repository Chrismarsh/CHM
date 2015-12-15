#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"
#include <meteoio/MeteoIO.h>
//#include "_snobal.h"
//#include "do_data_tstep.h"
//#include "init_snow.h"
#include "sno.h"
#include "snomacros.h"
class snodata : public face_info
{
public:
    sno data;
    int dead;

};
class snobal : public module_base
{
public:
    snobal();

    ~snobal();

    virtual void run(mesh_elem &elem, boost::shared_ptr <global> global_param);

    virtual void run(mesh domain, boost::shared_ptr <global> global_param);

    virtual void init(mesh domain, boost::shared_ptr<global> global);

};



