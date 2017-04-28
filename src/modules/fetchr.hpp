#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"


class fetchr : public module_base
{
public:
    fetchr(config_file cfg);

    ~fetchr();

    virtual void run(mesh_elem& face);

//number of steps along the search vector to check for a higher point
    int steps;
    //max distance to search
    double max_distance;

    //size of the step to take
    double size_of_step;

    bool incl_veg;

    //Obstacle heigh increment (m/m)
    //0.06 m/m corresponds to prarie shelter belts
    double I;

};



