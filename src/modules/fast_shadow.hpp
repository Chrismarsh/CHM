#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"

// Calculates horizon shadows using an adaptation of Dozier and Frew 1990 to an unstructured mesh
//Dozier, J., & Frew, J. (1990).
// Rapid calculation of terrain parameters for radiation modeling from digital elevation data.
// IEEE Transactions on Geoscience and Remote, 28(5), 963–969.
class fast_shadow : public module_base
{
public:
    fast_shadow(config_file cfg);

    ~fast_shadow();

    virtual void run(mesh_elem& face);

//number of steps along the search vector to check for a higher point
    int steps;
    //max distance to search
    double max_distance;

    //size of the step to take
    double size_of_step;

};



