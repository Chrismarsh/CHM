#pragma once
#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"



/**
 * Proof of concept for crop_rotation
 */
class crop_rotation : public module_base
{
public:
    crop_rotation(config_file cfg);

    ~crop_rotation();

    void run(mesh_elem &face);

};