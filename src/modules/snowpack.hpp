#pragma once

#include <boost/shared_ptr.hpp>

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"



#include <meteoio/MeteoIO.h>
#include <snowpack/libsnowpack.h>

#include <string>
class Lehning_snowpack : public module_base
{
public:
    Lehning_snowpack(config_file cfg);

    ~Lehning_snowpack();

    virtual void run(mesh_elem &elem, boost::shared_ptr <global> global_param);

    virtual void init(mesh domain, boost::shared_ptr <global> global_param);


    struct data : public face_info
    {
        //main snowpack model
        boost::shared_ptr<Snowpack> sp;

        /*
         * This is the PRIMARY data structure of the SNOWPACK program \n
         * It is used extensively not only during the finite element solution but also to control
         */
        boost::shared_ptr<SnowStation> Xdata;

        boost::shared_ptr<SnowpackConfig> Spackconfig;
        boost::shared_ptr<Meteo> meteo;
        boost::shared_ptr<Stability> stability;
        mio::Config config;
        double cum_precip;
    };

    double sn_dt; // calculation step length

};



