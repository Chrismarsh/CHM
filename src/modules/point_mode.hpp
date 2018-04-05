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

    virtual void run(mesh_elem &face);


    bool t ;
    bool rh ;
    bool U_R ;
    bool p ;
    bool ilwr ;
    bool iswr ;
    bool vw_dir ;
    bool iswr_diffuse ;
    bool iswr_direct ;
    bool U_2m_above_srf ;
    bool T_g;
};



