#pragma once
#include "soil_DTO.hpp"
#include <algorithm>
#include <iostream>
class soil_class_base
{
public:
	soil_class_base(two_layer_DTO& DTO) : d(DTO) {};

    ~soil_class_base() {};
protected:
	two_layer_DTO& d;

    void _push_excess_down(double& layer_storage, double& layer_max,
            double& layer_down);

};

class initializer : public soil_class_base
{
public:
    initializer(two_layer_DTO& DTO) : soil_class_base(DTO) {};
    ~initializer() {};

    void zero_single_step_vars();
    void layer_thaw_fraction();

};

class condensator : public soil_class_base
{
public:
	condensator(two_layer_DTO& DTO) : soil_class_base(DTO) {};
	~condensator() {};

	void set();

};

class infiltrator : public soil_class_base
{
public:
	infiltrator(two_layer_DTO& DTO) : soil_class_base(DTO) {};
	~infiltrator() {};

	void distribute();
	
};

class detention_layer : public soil_class_base
{
public:
    detention_layer(two_layer_DTO& DTO) : soil_class_base(DTO) {};
    ~detention_layer() {};
    
    void manage();
    
};

class depression_layer : public soil_class_base
{
public:
    depression_layer(two_layer_DTO& DTO) : soil_class_base(DTO) {};
    ~depression_layer() {};
    void manage();

};

class groundwater_layer : public soil_class_base
{
public:
    groundwater_layer(two_layer_DTO& DTO) : soil_class_base(DTO) {};
    ~groundwater_layer() {};
    void manage();

};

class subsurface_runoff : public soil_class_base
{
public:
    subsurface_runoff(two_layer_DTO& DTO) : soil_class_base(DTO) {};
    ~subsurface_runoff() {};

    void manage();

};
