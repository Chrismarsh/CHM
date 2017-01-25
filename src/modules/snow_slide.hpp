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
        double maxDepth;
        double snowdepthavg_copy;
        double swe_copy;
    };

};



