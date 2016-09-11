#pragma once



#include "module_base.hpp"
#include <meteoio/MeteoIO.h>

/**
 * Dodson, R. and Marks, D.: Daily air temperature interpolated at high spatial resolution over a large mountainous region, Clim. Res., 8(Myers 1994), 1–20, doi:10.3354/cr008001, 1997.
 */
class Dodson_NSA_ta : public module_base
{
public:
    Dodson_NSA_ta(config_file cfg);
    ~Dodson_NSA_ta();
    void run(mesh_elem &face, boost::shared_ptr<global> global_param);
    virtual void init(mesh domain, boost::shared_ptr<global> global_param);
    struct data : public face_info
    {
        interpolation interp;
    };
};


