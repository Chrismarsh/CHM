//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#include "Infil_All.hpp"
REGISTER_MODULE_CPP(Infil_All);

Infil_All::Infil_All(config_file cfg) : module_base("Infil_All", parallel::data, cfg)
{

    depends("swe");
    depends("snowmelt_int");
    depends("rainfall_int"); // NEW
    depends("soil_storage_at_freeze"); // NEW, depends on Volumetric model, equivalent to fallstat in crhm
    depends("soil_storage");
    depends("t");

    provides("inf");
    provides("total_inf");
    provides("snowinf"); // NEW
    provides("total_snowinf"); // NEW
    provides("total_excess");
    provides("total_meltexcess"); // NEW
    provides("runoff");
    provides("melt_runoff"); // NEW
    provides("total_rain_on_snow"); // NEW
    provides("rain_on_snow"); // NEW
    provides("frozen");
    provides("major_melt_count");
}

Infil_All::~Infil_All()
{

}

void Infil_All::init(mesh& domain)
{
    //store all of snobals global variables from this timestep to be used as ICs for the next timestep
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<Infil_All::data>(ID);

        // Triangle Totals
        d.total_inf = 0.;
        d.total_snowinf = 0.; // NEW
        d.total_excess = 0.;
        d.total_meltexcess = 0.; // NEW
        d.total_rain_on_snow = 0.; // NEW

        // Triangle variables
        d.frozen = false; // NEW, Maybe initial condition, not always necessary because of SWE check to freeze the ground
	    d.major_melt_count = 0; // NEW, For Gray frozen soil routine, counts number of major melts
        d.index = 0;
        d.max_major_per_melt = 0.;
        d.init_SWE = 0.;
        d.daily_melt_total = 0.;
        d.soil_storage = face->soil_attribute<double>("soil_storage");
        d.current_day_is_major = false;
        d.last_day = 0;
        d.tmax = 0.0;

        // Model Parameters
        infDays = cfg.get("max_inf_days",6);
        min_swe_to_freeze = cfg.get("min_swe_to_freeze",25);
        major = cfg.get("major",5); 
        AllowPriorInf = cfg.get("AllowPriorInf",true);
        thaw_type = cfg.get("thaw_type",0); // Default is Ayers

        SoilDataObj = std::make_unique<Soil::soils_na>();
        if (thaw_type == AYERS)
        {    
            d.texture = face->soil_attribute<std::string>("soil_texture","soils");
            d.ground_cover = face->soil_attribute<std::string>("soil_groundcover","soils");
        }
        else if (thaw_type == GREENAMPT)
        {
            d.soil_type = face->soil_attribute<std::string>("soil_type","soils");
            d.ksaturated = SoilDataObj->saturated_conductivity(d.soil_type);
        }
        lenstemp = cfg.get("temperature_ice_lens",-10.0);


        d.soil_storage_max = face->parameter("soil_storage_max"_s);


   }
}
void Infil_All::run(mesh_elem &face)
{
    // TODO if its water it should probably take all rain as "infil", there will be no snowmelt... sorta, snow melting and leaking into the water under the ice in spring??
    if(is_water(face))
    {
        set_all_nan_on_skip(face);
        return;
    }
    
    auto& d = face->get_module_data<Infil_All::data>(ID);

    auto id = face->cell_local_id;

    // CRHM does total infil for snow and total infil separately, wonder if I should do this
    double runoff = 0.;
    double melt_runoff = 0.;
    double inf = 0.;
    double snowinf = 0.;
    double rain_on_snow = 0.;

    double snowmelt = (*face)["snowmelt_int"_s];
    double rainfall = (*face)["rainfall_int"_s]; // NEW
    double swe = (*face)["swe"_s]; 
    double soil_storage_at_freeze = (*face)["soil_storage_at_freeze"_s];
    double airtemp = (*face)["t"_s];

    if (swe > min_swe_to_freeze && !d.frozen)
    {
        d.frozen = true; // Initiate frozen soil at 25 mm depth (as in CRHM)

        d.index = 0.;
        d.max_major_per_melt = 0.;
        d.init_SWE = 0.;
    }
    else if (swe <= 0.0 && d.major_melt_count > 0)
    {
        d.frozen = false;
        d.major_melt_count = 0;
    }

    if (d.frozen) // Gray's infiltration, 1985
    {
        if (rainfall > 0.0)
        {
            rain_on_snow = rainfall;
        }

        if (snowmelt > 0.0)
        {
            
            if (soil_storage_at_freeze == 0) // Unlimited
            {
                inf += snowmelt;
                d.major_melt_count = 1; 
            }
            else if (soil_storage_at_freeze > 0 && soil_storage_at_freeze < 100) // Limited
            {
                
                Check_for_ice_lens(d,airtemp);

                daily_melt_increment(d,snowmelt);
                increment_major_count(d);
                if (is_first_major(d,snowmelt,swe))
                {
                    SPDLOG_DEBUG("First Major");
                    Calc_Index(d,swe,soil_storage_at_freeze);
                    snowinf = Calc_Actual_Inf(d,snowmelt);
   
                    increment_major_count(d);
                }
                else if (is_limited_phase(d))
                {
                    SPDLOG_DEBUG("Limited Phase");
                    snowinf = Calc_Actual_Inf(d,snowmelt);
                    
                    increment_major_count(d);
                }
                else if (is_prior_first_major(d))
                {
                    SPDLOG_DEBUG("Prior");
                    snowinf = snowmelt;
                }

            }
            else if (soil_storage_at_freeze == 100) // Restricted
            {
                snowinf = 0.;
                d.major_melt_count = 1;
            }

           
            melt_runoff = snowmelt - snowinf;

            // melt_runoff and snowinf only track melt related quantities
            // total runoff and infiltrated amounts from ANY source are stored in
            // inf and runoff
            
            // This is a weird function, if there is any snowinf, then the rain on the snow also infiltrates
            // this is ported directly from CRHM module crack
            if (snowinf > 0.0)
            {
                inf += rain_on_snow;
            }
            else
            {
                runoff += rain_on_snow;
            }
            runoff += melt_runoff;
            inf += snowinf;


            Increment_Totals(d,runoff,melt_runoff,inf,snowinf,rain_on_snow);
            if (is_new_day(d))
            {
                d.last_day = global_param->day();
            }      

        }
    }
    else if (thaw_type == AYERS) // if not frozen, do Ayers
    {
        if (rainfall > 0.0)
        {
            // TODO set maxinfil at the beginning
            double maxinfil = SoilDataObj->ayers_texture(d.texture,d.ground_cover); 
            if (maxinfil > rainfall)
            {
                inf = rainfall;
            }
            else
            {
                inf = maxinfil;
                runoff = rainfall - maxinfil;
            }
        }
        
        melt_to_infil(inf,snowinf,snowmelt);
        
        // Increment totals
        
        Increment_Totals(d,runoff,melt_runoff,inf,snowinf,rain_on_snow);
    }
    else if (thaw_type == GREENAMPT) // if not frozen, do GreenAmpt
    {
        d.GA_temp = std::make_unique<data::tempvars>();

        if(rainfall > 0.0) {
            d.GA_temp->intensity = convert_to_rate_hourly(rainfall);

            if(d.soil_type == "pavement"){ // TODO Not a real option, handle this
                                           // ,this is a string handle pavement separately
                runoff = rainfall;
            }
            else if(is_space_in_dry_soil(d.soil_storage,d.soil_storage_max,rainfall)){
                inf =  rainfall;
            }
            else {
                
                //double F0 = d.soil_storage;
                Initialize_GA_Variables(d);

                // TODO what about ponding
                if (d.GA_temp->intensity > d.GA_temp->initial_rate) { // ponding is ongoing

                    d.GA_temp->final_storage = d.GA_temp->initial_storage + rainfall;
                    
                    double ponding_time = global_param->dt() / 3600.0;
                    
                    find_final_storage(d,d.GA_temp,d.GA_temp->initial_storage,ponding_time);    
                    
                    d.GA_temp->pond = rainfall - (d.GA_temp->final_storage - d.GA_temp->initial_storage); 

                    // TODO CRHM version calls find_final_storage again here, but appears unchanged. 
                    // Note that during tests that this could cause a difference.
                }
                else {

                    d.GA_temp->final_storage = d.GA_temp->initial_storage + rainfall;
                    d.GA_temp->final_rate = calc_GA_infiltration_rate(d,d.GA_temp,d.GA_temp->final_storage); //TODO calcf1 not a function anymore

                    if (d.GA_temp->intensity > d.GA_temp->final_rate) { // ponding starts midway through the time step
                        initialize_ponding_vars(d,d.GA_temp);

                        double ponding_time = global_param->dt() / 3600.0 - d.GA_temp->time_to_ponding;
                        
                        find_final_storage(d,d.GA_temp,d.GA_temp->storage_at_ponding,ponding_time);
                        
                        d.GA_temp->pond = rainfall - (d.GA_temp->final_storage - d.GA_temp->initial_storage); 
                    }

                }


                inf = d.GA_temp->final_storage - d.soil_storage; 
                if(d.GA_temp->pond > 0.0){
                    runoff = d.GA_temp->pond; 
                }
            }


            // Increment totals
            Increment_Totals(d,runoff,melt_runoff,inf,snowinf,rain_on_snow);
            d.soil_storage += d.GA_temp->final_storage;  

            melt_to_infil(inf,snowinf,snowmelt);

            d.GA_temp.reset();

        } // if(net_rain[hh] + net_snow[hh] > 0.0) greenampt routine
    }  



    // set variables to face
    (*face)["total_excess"_s]=d.total_excess;
    (*face)["total_meltexcess"_s]=d.total_meltexcess;
    (*face)["total_inf"_s]=d.total_inf;
    (*face)["total_snowinf"_s]=d.total_snowinf;
    (*face)["total_rain_on_snow"_s]=d.total_rain_on_snow;

    (*face)["runoff"_s]=runoff;
    (*face)["inf"_s]=inf;
    (*face)["rain_on_snow"_s]=rain_on_snow;
    (*face)["snowinf"_s]=snowinf;
    (*face)["melt_runoff"_s]=melt_runoff;
    (*face)["frozen"_s]=static_cast<int>(d.frozen);
    (*face)["major_melt_count"_s]=d.major_melt_count;
}

//General Functions
void Infil_All::Increment_Totals(Infil_All::data &d, double &runoff, double &melt_runoff, double &inf, double &snowinf, double &rain_on_snow) {
    d.total_inf += inf;
    d.total_excess += runoff;
    d.total_snowinf += snowinf;
    d.total_meltexcess += melt_runoff;
    d.total_rain_on_snow += rain_on_snow;
}      

void Infil_All::melt_to_infil(double& inf,double& snowinf,double& snowmelt)
{
    // This function determines what to do with the excess melt after swe = 0 when the crack model shuts off
    // Logan suggested setting it as runoff
    if (snowmelt > 0.0)
    {
        inf += snowmelt;             
        snowinf += snowmelt;
    }
};

// Crack Functions
void Infil_All::Calc_Index(Infil_All::data &d, double &swe, double &theta) {
    d.index = 5 * (1 - theta/100.0) * std::pow(swe,0.584);
    // d.major_major_per_melt is obtained by dividing d.index by the 
    // total number of time steps to get to d.index
    // This only works if 86400 / dt is a fraction which turns infDays into an integer
    // Example: if dt is 345600 (4 days in seconds) and infDays is 6 days. 
    // the denominator is 1.5, which then requires 2 major melts for it to stop, not 1.5
    // this is actually OK behaviour, but difficult to understand.
     
    d.max_major_per_melt = d.index / (infDays * 86400.0 / global_param->dt() );
    d.index = d.index / swe;
    d.init_SWE = swe;
}

double Infil_All::Calc_Actual_Inf(Infil_All::data &d, double &melt) {
    double inf = melt * d.index;
    if (inf > d.max_major_per_melt) {
        inf = d.max_major_per_melt;
    }
    return inf;
}


void Infil_All::Check_for_ice_lens(Infil_All::data &d, double &t) 
{
    d.tmax = std::max(d.tmax,t);

    if (is_new_day(d))
    {
        if (d.major_melt_count > 0 && d.tmax < lenstemp)
        {
            SPDLOG_DEBUG("Ice lens found"); 
            d.major_melt_count = infDays + 4;
        }
        d.tmax = 0.0;
    }
}

bool Infil_All::is_first_major(Infil_All::data& d, double& snowmelt, double& swe)
{
    return ( (d.major_melt_count == 0) & (is_major_melt(d)) ) || ( (swe >= d.init_SWE) & (is_limited_phase(d)));
};

bool Infil_All::is_major_melt(Infil_All::data& d)
{
    return d.daily_melt_total > major;
};

void Infil_All::increment_major_count(Infil_All::data& d)
{
    if (is_major_melt(d) && !d.current_day_is_major)
    {
        d.major_melt_count++;
        d.current_day_is_major = true;
    }
}; 


bool Infil_All::is_limited_phase(Infil_All::data& d)
{
    return d.major_melt_count > 0 && d.major_melt_count <= infDays;
};

bool Infil_All::is_prior_first_major(Infil_All::data& d)
{
    return d.major_melt_count == 0 and AllowPriorInf;
};

bool Infil_All::is_new_day(Infil_All::data& d)
{
    int current_day = global_param->day();
    
    if (current_day == d.last_day)
        return false;
    else 
    {
        return true;
    }
};

void Infil_All::daily_melt_increment(Infil_All::data& d, double& snowmelt)
{

    if (!is_new_day(d))
        d.daily_melt_total += snowmelt;
    else
    {
        d.daily_melt_total = snowmelt; 
        d.current_day_is_major = false;
    }
};
// Ayers

// Green-Ampt Functions

double Infil_All::convert_to_rate_hourly(double &rainfall) {
    return rainfall / (global_param->dt() / 3600.0);
}

bool Infil_All::is_space_in_dry_soil(double &moist, double &max, double &rainfall) {
    return moist == 0.0 && max >= rainfall;
}

void Infil_All::Initialize_GA_Variables(Infil_All::data &d) {
    // This function requires d.soil_storage so the full object d must be passed.
    // For simple reading, defined GA pointer to be consistent with other functions that use GA 
    // rather than GA_temp
    //
    // TODO Make this a constructor for the tempvars struct
    std::unique_ptr<Infil_All::data::tempvars> &GA = d.GA_temp;
    
    GA->soil_storage_deficit = (1.0 - d.soil_storage/d.soil_storage_max); // TODO GA in Dingman is porosity - pore space filed
                                                                        // Here: 1.0 means we've filled all the pores
                                                                        // 0.4 - 0.2 = 0.2 (porosity)
                                                                        // 1.0 - 0.5/1.0 = 0.5 (current)
                                                                        // Is this a problem?
    GA->initial_rate = calc_GA_infiltration_rate(d,GA,d.soil_storage);
    GA->initial_storage = d.soil_storage;
    GA->final_storage = GA->initial_storage;
    GA->final_rate = GA->initial_rate;
    GA->capillary_suction = SoilDataObj->capillary_suction(d.soil_type)
        * GA->soil_storage_deficit;
}

void Infil_All::initialize_ponding_vars(Infil_All::data& d,std::unique_ptr<Infil_All::data::tempvars> &GA) {
    GA->storage_at_ponding = d.ksaturated * GA->capillary_suction / (GA->intensity - d.ksaturated); 
    GA->time_to_ponding = (GA->storage_at_ponding - GA->initial_storage)/GA->intensity;
}

void Infil_All::find_final_storage(Infil_All::data& d,std::unique_ptr<Infil_All::data::tempvars> &GA, \
        double &initial_storage, double &dt) {
    
    double LastF1;

    do {
    
        LastF1 = GA->final_storage;
    
        GA->final_storage = initial_storage + d.ksaturated*dt + GA->capillary_suction * \
                            log((GA->final_storage + GA->capillary_suction) \
                            / (initial_storage + GA->capillary_suction));
    
    } while(fabs(LastF1 - GA->final_storage) > 0.001);

}

double Infil_All::calc_GA_infiltration_rate(Infil_All::data& d,std::unique_ptr<Infil_All::data::tempvars> &GA, double &F){

    return d.ksaturated*(GA->capillary_suction/F + 1.0);

}

//End Green-Ampt functions
