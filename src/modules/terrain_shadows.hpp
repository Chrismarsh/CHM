#pragma once

#include "../logger.h"
#include "triangulation.h"
#include "module_base.hpp"

#include <cstdlib>
#include <string>
#include <utility> //for pair
#include <vector>
#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

struct module_shadow_face_info :  face_info
{
    module_shadow_face_info()
    {
        shadow = 0;
    }
    std::vector< std::pair<size_t,Point_3> > v;
    int shadow;
};

struct vertex_flag : vertex_info
{
    vertex_flag()
    {
        visited = false;
    }
    bool visited;
    Point_3 prj_vertex;
    Point_3 org_vertex;
};
class terrain_shadow : public module_base
{
    public:
        terrain_shadow(std::string ID);
        ~terrain_shadow();
        virtual void run(mesh domain, boost::shared_ptr<global> global_param);
        

};

