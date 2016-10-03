#pragma once

#include "module_base.hpp"
#include <ogr_spatialref.h>


/**
 * \addtogroup modules
 * @{
 * \class solar
 * \brief Calculates solar position. Deals with UTM/geographic meshes.
 * This could have be it's own function, however was put into a module so-as to be able to cache the results if it is a UTM grid
 *
 * Depends:
 *
 */
class solar : public module_base
{

public:

    //if we have a UTM mesh, cache the calculated lat and long
    struct data : public face_info
    {
        double lat;
        double lng;
    };

    solar(config_file cfg);
    ~solar();
    void run(mesh_elem &face);
    void init(mesh domain);
};
