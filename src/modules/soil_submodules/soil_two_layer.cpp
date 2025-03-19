#include "soil_two_layer.hpp"

void soil_two_layer::run() 
{

    DTO.routing_residual = 0.0; // TODO zero'd because routing doesn't exist yet.
                                // also would need to remove or use face_area

    k_estimator.run();

    initialize_single_step_vars(); 

    //set_K_values();

    set_layer_thaw_fraction();

    set_condensation();

    organize_soil_layers();

    manage_detention();

    manage_depression();

    manage_groundwater();

    manage_subsurface_runoff();
    
};

void soil_two_layer::initialize_single_step_vars()
{
    DTO.condensation = 0.0;
    DTO.soil_excess_to_runoff = 0.0;
    DTO.soil_excess_to_gw = 0.0;
    DTO.ground_water_out = 0.0;
    DTO.soil_to_ssr = 0.0;
    DTO.rechr_to_ssr = 0.0;
    DTO.excess = 0.0;
    
};

void soil_two_layer::set_K_values()
{
    //K_estimate.calculate();
};

void soil_two_layer::set_layer_thaw_fraction()
{
    DTO.thaw_fraction_rechr = 0.0;
    DTO.thaw_fraction_lower = 0.0;

    double rechr_depth = 0.0;
    double soil_depth = 0.0;
    
        // TODO this might not be right, but the above is what CRHM does: 
        // depth = storage / porosity, because porosity = storage / depth;
    if (DTO.porosity > 0.0)  
    {    
    // TODO, porosity in rechr is the same as lower    
        rechr_depth = DTO.soil_rechr_max / DTO.porosity / 1000.0;
        soil_depth = DTO.soil_storage_max / DTO.porosity / 1000.0;
    }
    

    if (DTO.allow_runoff_from_infiltration || DTO.freeze_thaw_first_front == 0.0)
    {
        DTO.thaw_fraction_rechr = 1.0;
        DTO.thaw_fraction_lower = 1.0;
    }
    else
    {  
        DTO.thaw_fraction_rechr = 0.0;
        DTO.thaw_fraction_lower = 0.0;
    }

    if (DTO.soil_storage_max > 0.0 && DTO.allow_runoff_from_infiltration  &&
             DTO.freeze_thaw_first_front > 0.0) 
    {

        // TODO Verify this calculation
        if (DTO.thaw_front_depth < rechr_depth)
        {
            DTO.thaw_fraction_rechr = DTO.thaw_front_depth / rechr_depth;
        }
        else if ( DTO.thaw_front_depth < soil_depth)
        {
            DTO.thaw_fraction_rechr = 1.0;
            DTO.thaw_fraction_lower = (DTO.thaw_front_depth - rechr_depth) 
                / (soil_depth - rechr_depth);
        }
        else 
        {
            DTO.thaw_fraction_rechr = 1.0;
            DTO.thaw_fraction_lower = 1.0;        
        }
           
    }
    
};

void soil_two_layer::set_condensation()
{
    if (DTO.actual_ET < 0.0 && DTO.swe == 0.0)
    {
        DTO.condensation = -1.0 * DTO.actual_ET;
        DTO.actual_ET = 0.0;
    }
};

void soil_two_layer::organize_soil_layers()
{
    if (DTO.soil_storage_max > 0.0)
    {
        double soil_lower_storage = DTO.soil_storage - DTO.soil_rechr_storage;

        double potential = DTO.infil + DTO.condensation;

        double possible = DTO.thaw_fraction_rechr * (DTO.soil_rechr_max - DTO.soil_rechr_storage);
        if (possible > potential || !DTO.allow_runoff_from_infiltration)
            possible = potential;
        else
            DTO.soil_excess_to_runoff = potential - possible;

        DTO.soil_rechr_storage += possible;

        if (DTO.soil_rechr_storage > DTO.soil_rechr_max)
            _push_excess_down(DTO.soil_rechr_storage,DTO.soil_rechr_max,soil_lower_storage);

        DTO.soil_storage = soil_lower_storage + DTO.soil_rechr_storage;

        if (DTO.soil_storage > DTO.soil_storage_max)
            _push_excess_down(DTO.soil_storage,DTO.soil_storage_max,DTO.soil_excess_to_gw);
        
        if (DTO.swe == 0.0) // if there is no snowcover
        {
            DTO.rechr_to_ssr = DTO.soil_rechr_storage / DTO.soil_rechr_max * DTO.K_rechr_to_ssr * DTO.thaw_fraction_rechr;
            DTO.rechr_to_ssr = std::min(DTO.rechr_to_ssr,DTO.soil_rechr_storage * DTO.thaw_fraction_rechr);

            DTO.soil_rechr_storage = std::max(0.0, DTO.soil_rechr_storage - DTO.rechr_to_ssr);

            DTO.soil_storage -= DTO.rechr_to_ssr;
            DTO.soil_to_ssr = DTO.rechr_to_ssr;

        }

        if (DTO.soil_excess_to_gw > DTO.K_soil_to_gw * DTO.thaw_fraction_lower)
        {
            double excess_to_gw_max = DTO.K_soil_to_gw * DTO.thaw_fraction_lower;
            _push_excess_down(DTO.soil_excess_to_gw,excess_to_gw_max,DTO.excess);
        }

        // Line 607 of SoilX crhmcommetns branch, comment says upper layer but code says lower-layer, ask logan about this.
        if (DTO.excess_to_ssr && DTO.excess > 0.0)
        {
            double excess_to_ssr_max = DTO.excess * (1.0 - DTO.thaw_fraction_lower);
            _push_excess_down(DTO.excess,excess_to_ssr_max,DTO.soil_to_ssr);
        }

        
    }    
    else
    {
        DTO.excess = DTO.infil + DTO.condensation;
    }

        

};

void soil_two_layer::manage_detention()
{
    double face_area = 1.0; // Later will be pulled from face, but since routine_residual is always zero, ignoring it here.
    DTO.soil_excess_to_runoff += DTO.runoff + DTO.excess + DTO.routing_residual / face_area; // routing_residual comes from the crhm varaible redirected_residual which has units of mm*km^2/int (not sure why), so face_area is there for now for consistency.
    
    if (DTO.soil_excess_to_runoff > 0.0)
    {
        if (DTO.swe == 0.0)
            DTO.detention_max = DTO.detention_snow_max;
        else
            DTO.detention_max = DTO.detention_organic_max;

        double detention_space = DTO.detention_max - DTO.detention_storage;
        
        if (detention_space > 0.0)
        {
            if (DTO.soil_excess_to_runoff > detention_space)
            {
                DTO.soil_excess_to_runoff = std::max(0.0,DTO.soil_excess_to_runoff - detention_space); 
                DTO.detention_storage += detention_space;
            }
            else
            {
                DTO.detention_storage += DTO.soil_excess_to_runoff;
                DTO.soil_excess_to_runoff = 0.0;
            }
        }
    }

    if (DTO.detention_storage > 0.0 && DTO.K_detention_to_runoff > 0.0)
    {
        double transfer = std::min(DTO.detention_storage,DTO.K_detention_to_runoff);
        DTO.soil_excess_to_runoff += transfer;
        DTO.detention_storage -= transfer;

        if (DTO.detention_storage < 0.0001) // from CRHM, for safety and to drop any Floating-point errors
            DTO.detention_storage = 0.0;
    }
};

void soil_two_layer::manage_depression()
{
    if (DTO.soil_excess_to_runoff > 0.0 && DTO.depression_max > 0.0)
    {
        double exponent = -1.0 * std::min(12.0,DTO.soil_excess_to_runoff / DTO.depression_max);

        double depression_space = (DTO.depression_max - DTO.depression_storage) * (1 - exp(exponent));        

        if (DTO.soil_storage_max == 0.0)
            depression_space = DTO.depression_max - DTO.depression_storage;

        if (depression_space > 0.0)
        {
            if (DTO.soil_excess_to_runoff > depression_space)
            {
                DTO.soil_excess_to_runoff = std::max(0.0,DTO.soil_excess_to_runoff - depression_space);
                DTO.depression_storage += depression_space;
                // TODO add total tracker
            }
            else
            {
                DTO.depression_storage += DTO.soil_excess_to_runoff;
                // TODO add total tracker
                DTO.soil_excess_to_runoff = 0.0;
            }

        }
    }

    if (DTO.depression_storage > 0.0 && DTO.K_depression_to_gw > 0.0)
    {
        double value = transfer_min(DTO.depression_storage,DTO.K_depression_to_gw);
        DTO.depression_to_gw += value;
        DTO.depression_storage -= value;
        if (DTO.depression_storage < 0.0) // floatig point error safety?
            DTO.depression_storage = 0.0;
    }


};



void soil_two_layer::manage_groundwater()
{
    DTO.soil_excess_to_gw += DTO.depression_to_gw;
    DTO.depression_to_gw = 0.0;
    DTO.ground_water_storage += DTO.soil_excess_to_gw;

    if (DTO.ground_water_storage > DTO.ground_water_max)
    {
        _push_excess_down(DTO.ground_water_storage,DTO.ground_water_max,DTO.ground_water_out);
    }

    if (DTO.ground_water_max > 0.0) // divide by zero safety
    {
        double spilled = DTO.ground_water_storage / DTO.ground_water_max * DTO.K_ground_water_out;
        DTO.ground_water_storage -= spilled;
        DTO.ground_water_out += spilled;
    }

    // TODO daily and simulation totals, like in CRHM?

};

void soil_two_layer::manage_subsurface_runoff()
{
    if (DTO.depression_storage > 0.0 && DTO.K_depression_to_ssr > 0.0)
    {
        double value = transfer_min(DTO.depression_storage,DTO.K_depression_to_ssr);
        DTO.soil_to_ssr += value;
        DTO.depression_storage -= value;
        if (DTO.depression_storage < 0.0)
            DTO.depression_storage = 0.0;
    }

    if (DTO.K_lower_to_ssr > 0.0)
    {
        double available = DTO.soil_storage - DTO.soil_rechr_storage;
        double lower_unfrozen = DTO.K_lower_to_ssr * DTO.thaw_fraction_lower;
        double value = transfer_min(lower_unfrozen,available);
        DTO.soil_storage -= value;
        DTO.soil_to_ssr += value;
    }

};


void soil_two_layer::_push_excess_down(double& layer_storage, double& layer_max, double& layer_down)
{
    double excess = layer_storage - layer_max;

    layer_storage = layer_max;

    layer_down += excess;
};
    
double soil_two_layer::transfer_min(double& val1, double& val2)
{
    return std::min(val1,val2);
};




