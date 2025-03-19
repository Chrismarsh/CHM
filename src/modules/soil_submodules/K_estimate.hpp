#pragma once

#include "I_K_estimate.hpp"
#include "soil_DTO.hpp"
#include <iostream>
#include <memory>
// using an interface is pointless for now, since both are exposed, but this is more flexible in the future
// and would allow for dependecy injection from the parent module.
class I_Darcy_Vels 
{
public:
    virtual ~I_Darcy_Vels() = default;
    double lateral_rechr = 0.0;

    double lateral_lower = 0.0;

    double vertical_depression = 0.0;

    double vertical_lower = 0.0;

    double lateral_ground_water = 0.0;

    double lateral_detention = 0.0;

    virtual void calculate() = 0;

    I_Darcy_Vels(two_layer_DTO& _DTO) : DTO(_DTO) {};

    void init_vels();
protected:
    two_layer_DTO& DTO;
};

class Darcy_Vels : public I_Darcy_Vels
{
public:
    explicit Darcy_Vels(two_layer_DTO& _DTO) : I_Darcy_Vels(_DTO), exponent(3.0 + 2.0/DTO.pore_size_dist), exponent_organic(3.0 + 2.0/DTO.pore_size_dist_organic) {};
    ~Darcy_Vels() override {};
    virtual void calculate() override;

private:
    // water density * gravity acceleration / dynamic viscosity of water
    // only used exactly in this manner
    const double factor = 1000.0 * 9.8 * 0.001787;
    const double exponent;
    const double exponent_organic; 
    void set_snow();
    void set_clear();
    double get_lateral_lower();
    double get_reused(); // this is an expression reused many times, stored in a function
    double get_lateral_ground_water();
    double get_detention_snow();
    double get_detention_organic();

};
    
class K_estimate : public I_K_estimate
{
public:
    explicit K_estimate(two_layer_DTO& _DTO); 
    ~K_estimate() override {};

    void run(void) override;

private:

    two_layer_DTO& DTO;
    void check_soil_zeros();
    void set_K_values(I_Darcy_Vels& Vels);
    std::unique_ptr<I_Darcy_Vels> Vels;

     
};
