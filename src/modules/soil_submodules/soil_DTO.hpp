#pragma once

struct shared_DTO
{
    // shared between more than one submodule
    double soil_storage_max = 0.0;
    double swe = 0.0; //in - snobal TODO should this include input from PBSM3D?
    double soil_storage = 0.0;
    double soil_rechr_storage = 0.0;
    double soil_rechr_max = 0.0;
    double depression_storage = 0.0; 
    double actual_ET = 0.0; 
};

struct soil_ET_DTO : virtual shared_DTO
{

    double actual_soil_ET = 0.0; //out
    int ground_cover_type = 0;
    int soil_type_rechr = 0;
    int soil_type_lower = 0;

    virtual bool is_lake(soil_ET_DTO& DTO) = 0;
    
};

struct two_layer_DTO : virtual shared_DTO
{
    double thaw_front_depth = 0.0; //in - XG
    double freeze_front_depth = 0.0; //in - XG
    double freeze_thaw_first_front = 0.0; //in - XG
    double infil = 0.0; //in - infil_all
    double runoff = 0.0; //in - infil_all 
    double routing_residual = 0.0; //in - NetRoute

    double condensation = 0.0; //out (and used)
    double soil_excess_to_runoff = 0.0; //out (and used)
	double runoff_to_depression = 0.0; //out
    double soil_excess_to_gw = 0.0; //out (and used)
    double ground_water_out = 0.0; //out (and used)
    double soil_to_ssr = 0.0; //out (and used)
    double rechr_to_ssr = 0.0;
    double excess = 0.0;
    
    bool allow_runoff_from_infiltration;
			// All K_ stuff are in
    double K_soil_to_gw = 0.0; // All these K's could be swapped to held from previous and could be global/system wide
    double K_rechr_to_ssr = 0.0;
    double K_lower_to_ssr = 0.0;
    double K_detention_to_runoff = 0.0;
    double K_depression_to_ssr = 0.0;
    double K_depression_to_gw = 0.0;
    double K_ground_water_out = 0.0;
    double Ksaturated_rechr = 0.0;
    double Ksaturated_lower = 0.0;
    double Ksaturated_ground_water = 0.0;
    double Ksaturated_snow = 0.0;
    double Ksaturated_organic = 0.0;

    // held from previous 
    double thaw_fraction_rechr = 0.0;
    double thaw_fraction_lower = 0.0;
    double porosity = 0.0;
    bool excess_to_ssr = true; 
    double detention_max = 0.0;
    double detention_storage = 0.0;
    double detention_snow_max = 0.0;
    double detention_organic_max = 0.0;
    double detention_snow_init = 0.0;
    double detention_organic_init = 0.0;
    double depression_max = 0.0;
    double depression_to_gw = 0.0;
    double ground_water_storage = 0.0;
    double ground_water_max = 0.0;
    double snow_covered_threshold = 0.0;

    // K_estimate
    double pore_size_dist = 0.0;
    double pore_size_dist_organic = 0.0;
    double soil_index = 0.0;
    double local_slope = 0.0;
    double snow_grain_diameter = 0.0; // divided by 1000 in K_estimate, likely converting from mm to m, this var is mm
    double snow_density = 0.0;

    virtual int get_dt(two_layer_DTO& DTO) = 0;
    //virtual bool is_day_start(two_layer_DTO& DTO) = 0;
    //virtual void call_model_error(two_layer_DTO& DTO) = 0;
};

struct main_DTO : two_layer_DTO, soil_ET_DTO
{
};
