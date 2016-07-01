//#pragma once

#include <boost/shared_ptr.hpp>

#include "logger.hpp"
//#include "triangulation.hpp"
#include "module_base.hpp"
//#include "TPSpline.hpp"

#include <string>
class Simple_Canopy : public module_base
{
public:
    Simple_Canopy(config_file cfg);

    ~Simple_Canopy();

    virtual void run(mesh_elem &elem, boost::shared_ptr <global> global_param);

    virtual void init(mesh domain, boost::shared_ptr <global> global_param);


};



