#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"
#include <meteoio/MeteoIO.h>


class Iqbal_iswr : public module_base
{
public:
    Iqbal_iswr(config_file cfg);
    ~Iqbal_iswr();
    void run(mesh_elem& elem, boost::shared_ptr<global> global_param);
};
