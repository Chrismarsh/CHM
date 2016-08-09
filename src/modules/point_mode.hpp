#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"

/*
 * When using the model in point mode. Doesn't do any interp, but rather uses a specific input station as input.
 */
class point_mode : public module_base
{
public:
    point_mode(config_file cfg);

    ~point_mode();

    virtual void run(mesh_elem &face, boost::shared_ptr <global> global_param);


    bool t ;
    bool rh ;
    bool vw ;
    bool p ;
    bool ilwr ;
    bool iswr ;
    bool vw_dir ;
};



