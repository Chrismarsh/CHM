#pragma once


#include "../logger.h"
#include "module_base.hpp"

class degree_day : public module_base
{
    public:
        degree_day(std::string ID);
        ~degree_day();
        virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);

};