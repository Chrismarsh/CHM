#pragma once

#include <boost/shared_ptr.hpp>

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

#include <string>
class snow_slide : public module_base
{
public:
    snow_slide(config_file cfg);

    ~snow_slide();

    virtual void run(mesh domain);

    virtual void init(mesh domain);


    struct data : public face_info
    {
        double maxDepth; // m
        double snowdepthavg_copy; // m
        double swe_copy; // m (Note: swe units outside of snowslide are still mm)
        double delta_avalanche_snowdepth; // m^3
        double delta_avalanche_mass; // m^3
    };

};



