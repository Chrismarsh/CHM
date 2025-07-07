#include "soil_classes.hpp"
void soil_class_base::_push_excess_down(double& layer_storage, double& layer_max,
            double& layer_down)
{
    double excess = layer_storage - layer_max;

    layer_storage = layer_max;

    layer_down += excess;
};

void initializer::zero_single_step_vars()
{
    d.condensation = 0.0;
    d.soil_excess_to_runoff = 0.0;
    d.soil_excess_to_gw = 0.0;
    d.ground_water_out = 0.0;
    d.soil_to_ssr = 0.0;
    d.rechr_to_ssr = 0.0;
    d.excess = 0.0;
	d.runoff_to_depression = 0.0;
    d.routing_residual = 0.0; //For now, no routing so zeroed here, TODO
};

void initializer::layer_thaw_fraction()
{
    //rerun once daily
    if (!d.get_new_day())
        return;


    d.thaw_fraction_rechr = 0.0;
    d.thaw_fraction_lower = 0.0;
	
	// reference so intention is more clear
	// If we are allowing runoff from infiltration we are also allowing the soil to freeze
	bool& allow_soil_to_freeze = d.allow_runoff_from_infiltration;
    double rechr_depth = 0.0;
    double soil_depth = 0.0;
    
        // TODO this might not be right, but the above is what CRHM does 
        // depth = storage / porosity, because porosity = storage / depth;
    if (d.soil_storage_max ==0.0)
    {
        d.thaw_fraction_lower = 0.0;
        d.thaw_fraction_rechr = 0.0;
        return;
    }

    if (d.porosity > 0.0)  
    {    
    // TODO, porosity in rechr is the same as lower    
        rechr_depth = d.soil_rechr_max / d.porosity / 1000.0;
        soil_depth = d.soil_storage_max / d.porosity / 1000.0;
    }
    else
    {
        d.thaw_fraction_rechr = 0.0;
        d.thaw_fraction_lower = 0.0;
        return;
    };

    if (!allow_soil_to_freeze || d.freeze_thaw_first_front == 0.0)
    {
        d.thaw_fraction_rechr = 1.0;
        d.thaw_fraction_lower = 1.0;
    }
    else
    {  
        d.thaw_fraction_rechr = 0.0;
        d.thaw_fraction_lower = 0.0;
    }

    if (d.soil_storage_max > 0.0 && allow_soil_to_freeze  &&
             d.freeze_thaw_first_front > 0.0) 
    {

        // TODO Verify this calculation
        if (d.thaw_front_depth < rechr_depth)
        {   
            d.thaw_fraction_rechr = d.thaw_front_depth / rechr_depth;
        }
        else if ( d.thaw_front_depth < soil_depth)
        {
            d.thaw_fraction_rechr = 1.0;
            d.thaw_fraction_lower = (d.thaw_front_depth - rechr_depth) 
                / (soil_depth - rechr_depth);
        }
        else 
        {
            d.thaw_fraction_rechr = 1.0;
            d.thaw_fraction_lower = 1.0;        
        }
           
    }
};

void condensator::set()
{
    if (d.actual_ET < 0.0 && d.swe == 0.0)
    {
        d.condensation = -1.0 * d.actual_ET;
        d.actual_ET = 0.0;
    }
};

void infiltrator::distribute()
{
	if (d.soil_storage_max > 0.0)
	{
		double soil_lower_storage = d.soil_storage - d.soil_rechr_storage;

		double potential = d.infil + d.condensation;

		double possible = d.thaw_fraction_rechr * (d.soil_rechr_max - d.soil_rechr_storage);
		if (possible > potential || !d.allow_runoff_from_infiltration)
			possible = potential;
		else
			d.soil_excess_to_runoff = potential - possible;
		d.soil_rechr_storage += possible;

		if (d.soil_rechr_storage > d.soil_rechr_max)
			_push_excess_down(d.soil_rechr_storage,d.soil_rechr_max,soil_lower_storage);

		d.soil_storage = soil_lower_storage + d.soil_rechr_storage;

		if (d.soil_storage > d.soil_storage_max)
			_push_excess_down(d.soil_storage,d.soil_storage_max,d.soil_excess_to_gw);
		
		if (d.swe == 0.0) // if there is no snowcover
		{
			d.rechr_to_ssr = d.soil_rechr_storage / d.soil_rechr_max * d.K_rechr_to_ssr * d.thaw_fraction_rechr;
			d.rechr_to_ssr = std::min(d.rechr_to_ssr,d.soil_rechr_storage * d.thaw_fraction_rechr);

			d.soil_rechr_storage = std::max(0.0, d.soil_rechr_storage - d.rechr_to_ssr);

			d.soil_storage -= d.rechr_to_ssr;
			d.soil_to_ssr = d.rechr_to_ssr;

		}

		if (d.soil_excess_to_gw > d.K_soil_to_gw * d.thaw_fraction_lower)
		{
			double excess_to_gw_max = d.K_soil_to_gw * d.thaw_fraction_lower;
			_push_excess_down(d.soil_excess_to_gw,excess_to_gw_max,d.excess);
		}

		// Line 607 of SoilX crhmcommetns branch, comment says upper layer but code says lower-layer, ask logan about this.
		if (d.excess_to_ssr && d.excess > 0.0)
		{
			double excess_to_ssr_max = d.excess * (1.0 - d.thaw_fraction_lower);
			_push_excess_down(d.excess,excess_to_ssr_max,d.soil_to_ssr);
		}

		
	}    
	else
	{
		d.excess = d.infil + d.condensation;
	}
};

void detention_layer::manage()
{
    double face_area = 1.0; // Later will be pulled from face, but since routine_residual is always zero, ignoring it here.
    d.soil_excess_to_runoff += d.runoff + d.excess + d.routing_residual / face_area; // routing_residual comes from the crhm varaible redirected_residual which has units of mm*km^2/int (not sure why), so face_area is there for now for consistency.
    
    if (d.soil_excess_to_runoff > 0.0)
    {
        if (d.swe == 0.0)
            d.detention_max = d.detention_organic_max;
        else
            d.detention_max = d.detention_snow_max;

        double detention_space = d.detention_max - d.detention_storage;
        
        if (detention_space > 0.0)
        {
            if (d.soil_excess_to_runoff > detention_space)
            {
                d.soil_excess_to_runoff = std::max(0.0,d.soil_excess_to_runoff - detention_space); 
                d.detention_storage += detention_space;
            }
            else
            {
                d.detention_storage += d.soil_excess_to_runoff;
                d.soil_excess_to_runoff = 0.0;
            }
        }
    }

    if (d.detention_storage > 0.0 && d.K_detention_to_runoff > 0.0)
    {
        double transfer = std::min(d.detention_storage,d.K_detention_to_runoff);
        d.soil_excess_to_runoff += transfer;
        d.detention_storage -= transfer;

        if (d.detention_storage < 0.0001) // from CRHM, for safety and to drop any Floating-point errors
            d.detention_storage = 0.0;
    }
};

void depression_layer::manage()
{
    if (d.soil_excess_to_runoff > 0.0 && d.depression_max > 0.0)
        {
            double exponent = -1.0 * std::min(12.0,d.soil_excess_to_runoff / d.depression_max);

            double depression_space = (d.depression_max - d.depression_storage) * (1 - exp(exponent));        

            if (d.soil_storage_max == 0.0)
                depression_space = d.depression_max - d.depression_storage;

            if (depression_space > 0.0)
            {
                if (d.soil_excess_to_runoff > depression_space)
                {
                    d.soil_excess_to_runoff = std::max(0.0,d.soil_excess_to_runoff - depression_space);
                    d.depression_storage += depression_space;
					d.runoff_to_depression += depression_space; 
					// TODO add total tracker
                }
                else
                {
                    d.depression_storage += d.soil_excess_to_runoff;
                    // TODO add total tracker
                    d.soil_excess_to_runoff = 0.0;
                }

            }
        }

        if (d.depression_storage > 0.0 && d.K_depression_to_gw > 0.0)
        {
            double amount_to_move = std::min(d.depression_storage,d.K_depression_to_gw);
            d.depression_to_gw += amount_to_move;
            d.depression_storage -= amount_to_move;
            if (d.depression_storage < 0.0) // floatig point error safety?
                d.depression_storage = 0.0;
        }
};

void groundwater_layer::manage()
{
    d.soil_excess_to_gw += d.depression_to_gw;
    d.depression_to_gw = 0.0;
    d.ground_water_storage += d.soil_excess_to_gw;
                                                                                                   
    if (d.ground_water_storage > d.ground_water_max)
    {
        _push_excess_down(d.ground_water_storage,d.ground_water_max,d.ground_water_out);
    }
                                                                                                   
    if (d.ground_water_max > 0.0) // divide by zero safety
    {
        double spilled = d.ground_water_storage / d.ground_water_max * d.K_ground_water_out;
        d.ground_water_storage -= spilled;
        d.ground_water_out += spilled;
    }
};

void subsurface_runoff::manage()
{
    if (d.depression_storage > 0.0 && d.K_depression_to_ssr > 0.0)
    {
        double amount_to_move = std::min(d.depression_storage,d.K_depression_to_ssr);
        d.soil_to_ssr += amount_to_move;
        d.depression_storage -= amount_to_move;
        if (d.depression_storage < 0.0)
            d.depression_storage = 0.0;
    }

    if (d.K_lower_to_ssr > 0.0)
    {
        double available = d.soil_storage - d.soil_rechr_storage;
        double lower_unfrozen = d.K_lower_to_ssr * d.thaw_fraction_lower;
        double amount_to_move = std::min(lower_unfrozen,available);
        d.soil_storage -= amount_to_move;
        d.soil_to_ssr += amount_to_move;
    }   
};
