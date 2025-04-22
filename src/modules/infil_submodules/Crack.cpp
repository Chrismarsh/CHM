#include "Crack.hpp"

Crack::Crack(double& _major, double& _min_swe_to_freeze, unsigned int& _infDays, bool& _AllowPriorInf, double& _lenstemp,double& _steps_per_day, Crack::info& _d) : major(_major), min_swe_to_freeze(_min_swe_to_freeze), infDays(_infDays), AllowPriorInf(_AllowPriorInf), lenstemp(_lenstemp), steps_per_day(_steps_per_day), d(_d)
{

};

void Crack::init_inputs(double _snowmelt,double _rainfall, double _swe, double _soil_storage_at_freeze, double _airtemp, bool _newday)
{
    snowmelt = _snowmelt;
    rainfall = _rainfall;
    swe = _swe;
    soil_storage_at_freeze = _soil_storage_at_freeze;
    airtemp = _airtemp;
    is_newday = _newday;
};

void Crack::run()
{

    d.daily_rain_total += rainfall;
    if (is_newday) // Gray's infiltration, 1985
    {
        

        if (d.daily_melt_total > 0.0)
        {
            
            if (soil_storage_at_freeze == 0) // Unlimited
            {
                inf += d.daily_melt_total;
                d.major_melt_count = 1; 
            }
            else if (soil_storage_at_freeze > 0 && soil_storage_at_freeze < 100) // Limited
            {
                
                increment_major_count();
                
                Check_for_ice_lens();

                if (is_first_major())
                {
                    SPDLOG_DEBUG("First Major");
                    Calc_Index();
                    Calc_Actual_Inf();
                    //increment_major_count();
                }
                else if (is_limited_phase())
                {
                    SPDLOG_DEBUG("Limited Phase");
                    Calc_Actual_Inf();
                    
                    //increment_major_count();
                }
                else if (is_prior_first_major() && AllowPriorInf)
                {
                    SPDLOG_DEBUG("Prior");
                    inf = d.daily_melt_total;
                }

            }
            else if (soil_storage_at_freeze == 100) // Restricted
            {
                inf = 0.;
                d.major_melt_count = 1;
            }

           
            runoff = d.daily_melt_total - inf;
            // melt_runoff and snowinf only track melt related quantities
            // total runoff and infiltrated amounts from ANY source are stored in
            // inf and runoff
            
            // This is a weird function, if there is any snowinf, then the rain on the snow also infiltrates
            // this is ported directly from CRHM module crack
            if (inf > 0.0)
            {
                inf += d.daily_rain_total;
            }
            else
            {
                runoff += d.daily_rain_total;
            }
            melt_runoff = runoff;
            snow_inf = inf;


        }
        d.yesterday_melt = d.daily_melt_total;
        d.current_inf = inf;
        d.current_snow_inf = snow_inf;
        d.current_runoff = runoff;
        d.current_melt_runoff = melt_runoff;
        
        rain_on_snow = d.daily_rain_total;
        d.daily_melt_total = 0.0;
        d.daily_rain_total = 0.0;

    }
    else
    {
        runoff = d.current_runoff;
        melt_runoff = d.current_melt_runoff; 
        inf = d.current_inf;
        snow_inf = d.current_snow_inf;
    }
    
    if (!is_CRHM_compare_test)    
        d.daily_melt_total += snowmelt;
    update_t_max();
};


void Crack::Calc_Index() {
    d.index = 5 * (1 - soil_storage_at_freeze/100.0) * std::pow(swe,0.584);
    // d.major_major_per_melt is obtained by dividing d.index by the 
    // total number of time steps to get to d.index
    // This only works if 86400 / dt is a fraction which turns infDays into an integer
    // Example: if dt is 345600 (4 days in seconds) and infDays is 6 days. 
    // the denominator is 1.5, which then requires 2 major melts for it to stop, not 1.5
    // this is actually OK behaviour, but difficult to understand.
     
    d.max_major_per_melt = d.index / infDays;
    d.index = std::min(d.index / swe,1.0);
    d.init_SWE = swe;
};

void Crack::Calc_Actual_Inf() {
    inf = d.daily_melt_total * d.index;
    if (inf > d.max_major_per_melt && is_major_melt()) {
        inf = d.max_major_per_melt;
    }
};


void Crack::update_t_max() 
{
    if (is_newday)
        d.tmax = airtemp;
    else
        d.tmax = std::max(d.tmax,airtemp);
};

void Crack::Check_for_ice_lens()
{
    if (d.major_melt_count > 0 && d.tmax < lenstemp)
    {
        SPDLOG_DEBUG("Ice lens found");
        d.major_melt_count = infDays + 4;
    }
    d.tmax = 0.0;
};

bool Crack::is_first_major()
{
    // TODO This isn't really `is_first_major()`, only the part to the left of OR. Shoudl split this into two functions
    return (is_major_melt() && swe >= d.init_SWE && 
            (is_prior_first_major() || is_limited_phase()));//( (d.major_melt_count == 0) & (is_major_melt()) ) || ( (swe >= d.init_SWE) & (is_limited_phase()));
};

bool Crack::is_major_melt()
{
    return d.daily_melt_total > major;
};

void Crack::increment_major_count()
{
    if (is_major_melt())
    {
        d.major_melt_count++;
    }
}; 


bool Crack::is_limited_phase()
{
    return d.major_melt_count > 0 && d.major_melt_count <= infDays;
};

bool Crack::is_prior_first_major()
{
    return d.major_melt_count == 0;
};
