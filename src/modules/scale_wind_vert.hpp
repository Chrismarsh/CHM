#pragma once

#include <boost/shared_ptr.hpp>
#include <constants/Atmosphere.h>
#include "module_base.hpp"
#include <boost/shared_ptr.hpp>
#include "logger.hpp"
#include <string>
#include <math.h>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/ilu.hpp>
/**
 * \addtogroup modules
 * @{
 * \class scale_wind_vert
 * \brief Scales wind speed from reference height to defined heights
 *
 * Depends:
 * - U_R [m/s]
 * - Z_R [m] - Height of wind speed measurement/model layer
 *
 */

class scale_wind_vert : public module_base {
public:
    scale_wind_vert(config_file cfg);
    
    ~scale_wind_vert();
    virtual void init(mesh domain);

    //this module can swap between a domain parallel and a data parallel state
    ///domain parallel allows for blending through vegetation to avoid sharp gradietns that can complicate blowing snow, &c.
    virtual void run(mesh domain);
    virtual void run(mesh_elem &face);

    //scales the windspeed at a single triangle
    void point_scale(mesh_elem &face);

    bool ignore_canopy;
    //virtual void init(mesh domain);
    struct d: public face_info
    {
        double temp_u;
        interpolation interp;
    };
};