#pragma once
#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"


/**
 * Proff of concept for deformable mesh
 */
class deform_mesh : public module_base
{
public:
    deform_mesh(config_file cfg);

    ~deform_mesh();

    void run(mesh domain);

};