#pragma once

#include "soil_base.hpp"
#include "soil_DTO.hpp"
class soil_ET : public soil_base
{
public:
    soil_ET(soil_ET_DTO& _DTO);
    ~soil_ET();

    void run() override;

private:

    soil_ET_DTO& DTO;
    
    double set_ET_layer(double et, double& percent, int& type);
};

