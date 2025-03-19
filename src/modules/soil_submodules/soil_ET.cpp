#include "soil_ET.hpp"
#include <iostream>
soil_ET::soil_ET(soil_ET_DTO& _DTO) : DTO(_DTO)
{

};

soil_ET::~soil_ET()
{

};

void soil_ET::run()
{  
    
    DTO.actual_ET = 0.0;

    if (DTO.swe > 0.0)
       return; 

    std::cout << "in soil ET" << std::endl;
    double available_to_evap = DTO.potential_ET;
    if (DTO.depression_storage + DTO.soil_storage > 0.0)
    {
        available_to_evap *= DTO.depression_storage / 
                (DTO.depression_storage + DTO.soil_storage);
    }
    else
        available_to_evap = 0.0;

    if (DTO.depression_storage > 0.0 && available_to_evap > 0.0)
    {
        if (DTO.depression_storage > available_to_evap)
            DTO.depression_storage = std::max(0.0,DTO.depression_storage - available_to_evap);
        else
        {
            available_to_evap = DTO.depression_storage;
            DTO.depression_storage = 0.0;
        }
        DTO.actual_ET = available_to_evap;
    }
    else
        available_to_evap = 0.0;
    
    std::cout << "in soil ET" << std::endl;
    available_to_evap = DTO.potential_ET - available_to_evap;

    if (available_to_evap > 0.0 && DTO.soil_storage > 0.0 && DTO.ground_cover_type > 0)
    {
        double percent_available_lower;
        double percent_available_rechr;
        double ET_lower;
        double ET_rechr;
        double soil_lower_storage = DTO.soil_storage - DTO.soil_rechr_storage;
        double soil_lower_max = DTO.soil_storage_max - DTO.soil_rechr_max;

        if ( soil_lower_max > 0.0 ) // soil_lower > 0.0
            percent_available_lower = soil_lower_storage / soil_lower_max;
        else
            percent_available_lower = 0.0;

        percent_available_rechr = DTO.soil_rechr_storage / DTO.soil_rechr_max;

        ET_rechr = available_to_evap;
        
        ET_rechr = set_ET_layer(ET_rechr,percent_available_rechr,DTO.soil_type_rechr);

        if (ET_rechr > available_to_evap) 
        {
            ET_lower = 0.0;
            ET_rechr = available_to_evap;
        }
        else
            ET_lower = available_to_evap - ET_rechr;

        if (ET_lower > 0.0)
        {
            ET_lower = set_ET_layer(ET_lower, percent_available_lower, DTO.soil_type_lower);
        }

        double ET = 0.0;
        std::cout << "in soil ET" << std::endl; 
        switch (DTO.ground_cover_type)
        {
        case 0: // bare soil, no evap
            break; 
        case 1: // recharge layer only, crops 
            if (ET_rechr > DTO.soil_rechr_storage)
            {
                ET = DTO.soil_rechr_storage;
                DTO.soil_rechr_storage = 0.0;
            }
            else
            {
                DTO.soil_rechr_storage -= ET_rechr;
                ET = ET_rechr;
            }
            DTO.soil_storage -= ET_rechr;
            break;
        case 2: // all soil moisture, grasses & shrubs
            if (ET_rechr + ET_lower >= DTO.soil_storage)
            {
                ET = DTO.soil_storage;
                DTO.soil_storage = 0.0;
                DTO.soil_rechr_storage = 0.0;
            }
            else
            {
                ET = ET_rechr + ET_lower;
                DTO.soil_storage -= ET;

                DTO.soil_rechr_storage = std::max(DTO.soil_rechr_storage - ET_rechr, 0.0);
            }
            break;
        }

        DTO.actual_ET += ET;
        std::cout << "in soil ET" << std::endl;
        if (DTO.is_lake(DTO))
        {
            std::cout << "out" <<  std::endl;
            DTO.actual_ET = DTO.potential_ET;
        }
    };
};

double soil_ET::set_ET_layer(double et, double& percent, int& type)
{
    switch (type)
    {
    case 1: // sandy soil
        if (percent < 0.25)
            et = 0.5 * percent * et;
        break;
    case 2: // loam soil
        if (percent < 0.5)
            et = percent * et;
        break;
    case 3: // clay soil
        if (percent <= 0.33)
            et = 0.5 * percent * et;
        else if (percent < 0.67)
            et = percent * et;
        break;
    default: //organic soil
        // do nothing, use default above
        break;
    }

    return et;
};


