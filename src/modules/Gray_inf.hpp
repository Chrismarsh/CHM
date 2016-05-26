#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"
#include <cmath>



class Gray_inf : public module_base
{
public:
    Gray_inf(config_file cfg);

    ~Gray_inf();

    void run(mesh_elem &elem, boost::shared_ptr <global> global_param);
    void init(mesh domain, boost::shared_ptr<global> global);

    class data : public face_info
    {
    public:
        double storage;  //mm
        double max_storage;
        double porosity; // [-]
        double soil_depth; //mm

        double opportunity_time; //h
        double last_ts_potential_inf; // last time steps potential infiltration
        double total_inf;
        double total_excess;

    };
};



