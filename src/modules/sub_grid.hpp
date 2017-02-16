#pragma once

#include <boost/shared_ptr.hpp>
#include "module_base.hpp"
#include "logger.hpp"
#include <string>
#include <math.h>

class sub_grid : public module_base {
public:
    sub_grid(config_file cfg);
    
    ~sub_grid();

    virtual void run(mesh_elem &face);

};