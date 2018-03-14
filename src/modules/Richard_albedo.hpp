#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include <meteoio/MeteoIO.h>

/**
 * Eqn 4, 5
 * Essery, R., and P. Etchevers (2004), Parameter sensitivity in simulations of snowmelt, J. Geophys. Res., 109(D20111), 1–15, doi:10.1029/2004JD005036.
 */
class Richard_albedo : public module_base
{
public:
    struct data : public face_info
    {
        double albedo;

    };

    Richard_albedo(config_file cfg);
    ~Richard_albedo();
    void run(mesh_elem& face);
    void init(mesh domain);
    void checkpoint(mesh domain,  netcdf& chkpt);
    void load_checkpoint(mesh domain,  netcdf& chkpt);

    double amin;
    double amax;
    double a1;
    double a2;
    double albedo_snow;
    double albedo_bare;
    double min_swe_refresh;
};


