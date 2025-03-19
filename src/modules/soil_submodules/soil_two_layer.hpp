#pragma once

#include "soil_base.hpp"
#include "soil_two_layer.hpp"
#include "soil_DTO.hpp"
#include "I_K_estimate.hpp"
// Temporary forward declaration
// TODO add include after writing this
// class k_estimate;

class soil_two_layer : public soil_base
{
public: // methods
    soil_two_layer(two_layer_DTO& _DTO,I_K_estimate& estimate) : DTO(_DTO), k_estimator(estimate) {};
    ~soil_two_layer() {};

    void run() override;

private: //methods
    
    two_layer_DTO& DTO; 
    I_K_estimate& k_estimator; 
    
    void initialize_single_step_vars(); 
    void set_K_values();
    void set_layer_thaw_fraction(); // set the new value in a layer, ensure values are between 0 and 1.
    void set_condensation();
    void organize_soil_layers();
    void manage_detention();
    void manage_depression();
    void manage_groundwater();
    void manage_subsurface_runoff();


    void _push_excess_down(double& layer_storage, double& layer_max, double& layer_down); 
    double transfer_min(double& val1, double& val2);
};
