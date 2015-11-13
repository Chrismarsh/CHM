#pragma once



#include "module_base.hpp"
#include <meteoio/MeteoIO.h>
#include "TPSpline.hpp"

class Dodson_NSA_ta : public module_base
{
public:
    Dodson_NSA_ta();
    ~Dodson_NSA_ta();
    void run(mesh_elem &elem, boost::shared_ptr<global> global_param);
};


