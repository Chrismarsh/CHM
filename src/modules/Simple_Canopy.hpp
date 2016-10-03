#pragma once

#include <boost/shared_ptr.hpp>
#include "logger.hpp"
#include "module_base.hpp"
#include <meteoio/MeteoIO.h>
#include <constants/Atmosphere.h>
#include <constants/PhysConst.h>
#include <constants/Vegetation.h>
#include <string>
#include "face.hpp"

/*
 * @Brief Solves energy and mass equations for canopy states using parameterizations based on LAI and canopy closure.
 *
 */

class Simple_Canopy : public module_base
{
public:
    Simple_Canopy(config_file cfg);

    ~Simple_Canopy();

    virtual void run(mesh_elem &elem);

    virtual void init(mesh domain);

    double delta(double ta);

    double lambda(double ta);

    double gamma(double air_pressure, double ta);

    double Qs(double air_pressure, double ta);

    struct data : public face_info {

        // parameters
        double LAI;
        double CanopyHeight;
        int canopyType;

        // model states
        double rain_load;
        double Snow_load;

        // Output diagnostic variables
        double cum_net_snow;
        double cum_net_rain;
        double cum_Subl_Cpy;
        double cum_intcp_evap;
        double cum_SUnload_H2O;
    };



};



