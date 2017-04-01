#pragma once

#include <boost/shared_ptr.hpp>
#include <constants/Atmosphere.h>
#include "module_base.hpp"
#include <boost/shared_ptr.hpp>
#include "logger.hpp"
#include <string>
#include <math.h>

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

    double get_U_2m_above_srf(mesh_elem face);

    ~scale_wind_vert();

    virtual void run(mesh_elem &face);

    //virtual void init(mesh domain);
};