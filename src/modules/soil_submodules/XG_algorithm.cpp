#include "XG_algorithm.hpp"
XG_algorithm::XG_algorithm(const double& t, const double& _moist, 
            const double& _rechr, state& _S, const params& _P) : 
        surface_temp(t), soil_moist(_moist), soil_rechr(_rechr), 
		S(_S), P(_P), depths(this->P.getdepths())
{
	depths.back() = 100.0;	
};



void XG_algorithm::run()
{
    if (S.last_step_new_day)
    {   
        S.reset_degree_day_counter();
    }

    
    S.accumulate_degree_days(surface_temp);
    
        
    if (S.is_newday)
    {
        S.determine_freeze_thaw_idle();

        S.set_thermal_conductivities(P,soil_moist,soil_rechr)
            .set_freezethaw_ratios(P);
        
        if(S.freezing()) // handle freezing
        {
            if(S.net_negative_degree_days())
            {
                freeze(); // XG-Algorithm - Freezing

                // check for thaw front
                if(S.thaw_front_overtaken())
                {
                    if(S.exist_excess_fronts())
                    {
                        S.merge_freeze_fronts(this);
                    }
                    else
                    {
                        S.Zdt = 0.0;
                        S.Bth = 0.0;
                        S.Zd_front[1] = 0.0;
                    }
                    
                }
            }
            S.Zd_front[0] = -S.Zdf;
        } // freezing handled
        else if(S.thawing()) // Surface thawing lower ground frozen
        {
            S.Bth += S.B;  // Accumulate thawing degree-days

            if(S.Bth <= 0.0)
            {
                S.Bfr = S.Bth;
                S.Zdt = 0.0;
                S.Bth = 0.0;
            }
            else
            {
                thaw(); // XG-Algorithm - Thawing
                
                // check for freeze front
                if(S.freeze_front_overtaken())
                {
                    if(S.exist_excess_fronts())
                        S.merge_thaw_fronts(this);
                    else
                    {
                        S.Zdf = 0.0;
                        S.Bfr = 0.0;
                        S.Zdt = 0.0;
                        S.Bth = 0.0;
                        S.TrigState = 0;
                        S.Zd_front[1] = 0.0;
                    }

                    
                }
                
            } // if check freeze front
            
            S.Zd_front[0] = S.Zdt;
        } // thawing handled
        
        S.last_step_new_day = true;
    }
    thaw_front_depth = S.Zdt;
    freeze_front_depth = S.Zdf;
    first_front_depth = S.Zd_front[0];
};

void XG_algorithm::freeze(void)
{
    size_t layer = 1;
    double Za;
    double L = 335000;

    S.Zdf = 0.0;
    
    double ftc;
    if (P.k_update == 2)
        ftc = Interpolated_ftc(S.Zdf, layer);
    else
        ftc = S.ftc[layer-1];

    Za = stefan_equation(S.Bfr,ftc,layer);

    while (Za > depths[layer-1] && layer < P.N_Soil_layers)
    {
        S.Zdf += depths[layer-1];

        Za = (Za - depths[layer-1])/S.pf[layer];
        ++layer;

        if(P.k_update > 0 && P.freeze_kw_ki_update && layer > S.Fz_low)
        {
            S.ftc[S.Fz_low-1] = S.get_ttc(S.Fz_low-1,P);
            S.ftc_contents[S.Fz_low-1] = 1;
            S.pf[S.Fz_low] = std::sqrt(S.ftc[S.Fz_low-1]*S.layer_h2o[S.Fz_low]/(S.ftc[S.Fz_low]*S.layer_h2o[S.Fz_low-1]));
            S.Fz_low = layer;
            S.tc_composite[layer-1] = S.ftc[layer-1];
        }
    }
    S.Zdf += Za;

    S.Zdf = std::min(S.Zdf,P.Zpf_init);
    
};

void XG_algorithm::thaw(void)
{
    size_t layer = 1;
    double Za;
    double L = 335000;

    S.Zdt = 0.0;

    double ttc;
    if (P.k_update == 2)
        ttc = Interpolated_ftc(S.Zdf, layer);
    else
        ttc = S.ttc[layer-1];

    Za = stefan_equation(S.Bth,ttc,layer);

    while (Za > depths[layer-1] && layer < P.N_Soil_layers)
    {
        S.Zdt += depths[layer-1];

        Za = (Za - depths[layer-1])/S.pt[layer];
        ++layer;

        if(P.k_update > 0 && P.thaw_ki_kw_update && layer > S.Th_low)
        {
            S.ttc[S.Th_low-1] = S.get_ftc(S.Th_low-1,P);
            S.ttc_contents[S.Th_low-1] = 1;
            S.pt[S.Th_low] = std::sqrt(S.ttc[S.Th_low-1]*S.layer_h2o[S.Th_low]/(S.ttc[S.Th_low]*S.layer_h2o[S.Th_low-1]));
            S.Th_low = layer;
            S.tc_composite[layer-1] = S.ttc[layer-1];
        }
    }
    S.Zdt += Za;

    S.Zdt = std::min(S.Zdt,P.Zpf_init);

};

double XG_algorithm::stefan_equation(double& surface_index, double& thermal_conductivity,size_t& layer)
{
    const double L = 335000;
    // TODO Don't know where the 211 and 8640011 comes from. 
    // Examining the original equation from (Xie and Gough, 2013)
    // A factor of 86400 is needed
    assert(S.layer_h2o[layer-1] != 0 && "About to trip divide by zero in stefan equation!");
    return std::sqrt(2.0*86400.0*thermal_conductivity * surface_index / 
            (S.layer_h2o[layer-1]*L));

};

double XG_algorithm::Interpolated_ttc(double Za, size_t layer)
{ 
  if(!P.thaw_ki_kw_update)
    return (S.ttc[layer-1]);

  double split = (Za - S.Zdt)/depths[layer-1];
  if(split >= 1.0)
    split = 1.0;

  double combination = S.ttc[layer-1] - split*(S.ttc[layer-1] - S.ftc[layer-1]); // thawed(18k) to frozen (4k)

  S.tc_composite2[layer-1] = combination;

  return (combination);
};

double XG_algorithm::Interpolated_ftc(double Za, size_t layer) { //

  if(!P.freeze_kw_ki_update)
    return (S.ftc[layer-1]);

  double split = (Za - S.Zdf)/depths[layer-1];
  if(split >= 1.0)
    split = 1.0;

  double combination = S.ftc[layer-1] + split*(S.ttc[layer-1] - S.ftc[layer-1]); // frozen (4k) to thawed(18k)

  S.tc_composite2[layer-1] = combination;

  return (combination);
};

#define tolerance 0.000001
void XG_algorithm::find_thaw_D(double dt) { // XG-Algorithm - Thawing - used by init
// solve for Bth from Zdt using Bisection method
    
    if(dt == 0)
        return;

    auto solution = [&dt](double Zdt)
    { return Zdt - dt; };

    double low = 0.0;
    double high = 50000.0;

    do {
        double mid = (high + low)/2;
        S.Bth = mid;
        thaw();
        if (solution(S.Zdt) > 0)
            high = mid;
        else
            low = mid;

    } while (high - low > tolerance);
    
  //TODO Throw CHM exception here, indicates that Zdt is too large
};

void XG_algorithm::find_freeze_D(double df) { // XG-Algorithm - Thawing - used by init
// solve for Bfr from Zdt using Bisection method
    if(df == 0)
        return;
    
    auto solution = [&df](double Zdf)
    { return Zdf - df; };

    double low = 0.0;
    double high = 50000.0;

    do {
        double mid = (high + low)/2;
        S.Bfr = mid;
        freeze();
        if (solution(S.Zdf) > 0)
            high = mid;
        else
            low = mid;

    } while (high - low > tolerance);
    //TODO Throw CHM exception here, indicates Zdf is too large
};

void XG_algorithm::init_freezethaw_degreedays(const double& Zdf_init, const double& Zdt_init, const double& Zpf_init)
{
	if (Zdf_init > 0.0 && S.Bfr == 0.0)
		find_freeze_D(Zdf_init);

	if (Zdt_init > 0.0 && Zdt_init < Zpf_init && S.Bth == 0.0)
		find_thaw_D(Zdt_init);
};

void XG_algorithm::state::push_front(double D) {

 // if(nfront >= front_size-3){ // space to allocate plus Zdf/Zdt(2 slots) plus top of stack indicator
 //   string S = string("'") + Name + " (XG)' too many fronts in hru = " + to_string(hh+1).c_str();
 //   CRHMException TExcept(S.c_str(), TExcept::TERMINATE);
 //   LogError(TExcept);
 //   throw TExcept;
 // }
 // TODO, above if statement is needed for some edge cases, but I don't know what front_size even is. I have this on my XG note TODO


  for(size_t ii = nfront+1; 2 <= ii ; --ii) // move contents up
    Zd_front[ii+1] = Zd_front[ii];

  ++nfront;
  Zd_front[2] = D; // add new entry
};

double XG_algorithm::state::pop_front(void) {

    double D = std::fabs(Zd_front[2]); // always positive

    for(size_t ii = 2; ii < nfront+1; ++ii) // move contents down
        Zd_front[ii] = Zd_front[ii+1];

    Zd_front[nfront+1] = 0.0; // clear memory

    --nfront;

    return D;
};

double XG_algorithm::state::last_front(void){

  if(!nfront)
    return 0.0;
  else
    return (Zd_front[2]);
};


double XG_algorithm::state::get_ftc(size_t layer,const params& P){ // unfrozen(thawed) soil to be frozen
  if(P.calc_conductivity){
    return (P.soil_solid_km_kw[layer] - P.soil_solid_km[layer])*std::pow(this->layer_h2o[layer]/(1000.0*P.por[layer]),2) + P.soil_solid_km[layer];
  }
  else
    return (1.0 - P.por[layer])*P.soil_solid_km[layer] + this->layer_h2o[layer]/1000.0*kw + (P.por[layer] - this->layer_h2o[layer]/1000.0)*ka;
};

double XG_algorithm::state::get_ttc(size_t layer,const params& P){ // frozen soil to be unfrozen(thawed)
  if(P.calc_conductivity){
    return  P.soil_solid_km[layer]*std::pow(P.soil_solid_km_ki[layer]/P.soil_solid_km[layer], this->layer_h2o[layer]/(1000.0*P.por[layer]));
  }
  else
    return (1.0 - P.por[layer])*P.soil_solid_km[layer] + this->layer_h2o[layer]/1000.0*ki + (P.por[layer] - this->layer_h2o[layer]/1000.0)*ka;
};

void XG_algorithm::state::determine_freeze_thaw_idle()
{
    if (TrigAcc > P.Trigthrhld)
        TrigAcc = P.Trigthrhld;
    else if (TrigAcc < -P.Trigthrhld)
        TrigAcc = -P.Trigthrhld;
    
    // Start Thaw starting from Idle
    if (TrigAcc >= P.Trigthrhld / 2.0 && TrigState == 0 && (Zdf > 0.0 || nfront > 0))
    {
        TrigState = 1;
        Zd_front[1] = -Zdf;
        t_trend = 0.0;
    }            

    // Start Freeze starting from Idle
    if (TrigAcc <= -P.Trigthrhld / 2.0 && TrigState == 0)
    {
        TrigState = -1;
        Zd_front[1] = Zdt;
        t_trend = 0.0;
    }

    // Start Idle starting from Freeze
    if (TrigState == -1 && TrigAcc >= P.Trigthrhld/2.0 && t_trend > 0.0)
    {
        TrigState = 0;

        if (Zdt > 0.0 && Zdf > 0.0)
        {
            if (Zdt > Zdf)
            {
                push_front(Zdt);
                Zdt = 0.0;
                Bth = 0.0;
                Zd_front[0] = 0.0;
                Zd_front[1] = -Zdf;
            }
        }
    }

    // Start Idle starting from Thaw
    if (TrigState == 1 && TrigAcc <= -P.Trigthrhld / 2.0 && t_trend < 0.0)
    {
        TrigState = 0;

        if (Zdf > 0.0 && Zdt > 0.0)
        {
            if (Zdf > Zdt)
            {
                push_front(-Zdf);
                Zdf = 0.0;
                Bfr = 0.0;
                Zd_front[0] = 0.0;
                Zd_front[1] = Zdt;
            }
        }
    }
};

bool XG_algorithm::state::freezing()
{
    return TrigState < 0;
};

bool XG_algorithm::state::thawing()
{
    return TrigState > 0;
};

bool XG_algorithm::state::net_negative_degree_days()
{
    Bfr -= B;
    return Bfr > 0.0;
};

void XG_algorithm::state::reset_degree_day_counter()
{
    last_step_new_day = false;
    B = 0.0;
};

void XG_algorithm::state::merge_freeze_fronts(XG_algorithm* XG)
{
    double Last = last_front();
    
    if(Last < 0.0) // frozen front
    {
        Zdf = pop_front();
        XG->find_freeze_D(Zdf);
        double Last = last_front();
        
        if(Last > 0.0) // thaw front
        {
            Zdt = pop_front();
            XG->find_thaw_D(Zdt);
            Zd_front[1] = Zdt;
        }
        else if(Last < 0.0) // never two frozen fronts together
        {
           //CRHM has thos throw an error, never two frozen fronts 
        }
        else // no thaw front
        {
            Zdt = 0.0;
            Bth = 0.0;
            Zd_front[1] = 0.0;
        }
    }
    else if(Last < 0.0) // never two freeze fronts together
    {
        // CRHM had another throw here, but I dont think its possible to reach this ever.
    }
    else // no thaw layer
    {
        Zdt = 0.0;
        Bth = 0.0;
        Zd_front[1] = 0.0;
    }

};

void XG_algorithm::state::merge_thaw_fronts(XG_algorithm* XG)
{
    double Last = last_front();
    
    if(Last > 0.0) // thaw front
    {
        Zdt = pop_front();
        XG->find_thaw_D(Zdt);
        double Last = last_front();
        
        if(Last < 0.0) // frozen front
        {
            Zdf = pop_front();
            XG->find_freeze_D(Zdf);
            Zd_front[1] = -Zdf;
        }
        else if(Last > 0.0) // never two thaw fronts together
        {
           //CRHM has thos throw an error, never two frozen fronts 
        }
        else // no frozen front
        {
            Zdf = 0.0;
            Bfr = 0.0;
            Zd_front[1] = 0.0;
        }
    }
    else if(Last < 0.0) // never two freeze fronts together
    {
        // CRHM had another throw here, but I dont think its possible to reach this ever.
    }
    else // no thaw layer
    {
        Zdf = 0.0;
        Bfr = 0.0;
        Zd_front[1] = 0.0;
    }
};

bool XG_algorithm::state::exist_excess_fronts()
{
    return nfront > 0;
};

bool XG_algorithm::state::freeze_front_overtaken()
{
    return Zdf > 0.0 && Zdt >= Zdf;
};

bool XG_algorithm::state::thaw_front_overtaken()
{
    return Zdt > 0.0 && Zdf >= Zdt;
};

XG_algorithm::state& XG_algorithm::state::set_layer_moisture_maximums(const XG_algorithm::params& P)
{
    //XG_algorithm xg(0.0,0.0,0.0,*this,P);
    //xg.printdebug();
    double sum = 0.0;
    std::vector<double> depths = P.getdepths();
    for (double val : depths)
    {
        sum += val;
    }
    if (sum < this->Zdf || sum < this->Zdt)
    {
        //TODO Add A CHM exception to say that the total soil depth is less than initial Zdt,Zdf
    } 
    double rechrmax = P.soil_rechr_max;
    double soilmax = P.soil_moist_max;
    //xg.printdebug();
    for (size_t layer = 0; layer < P.N_Soil_layers; ++layer)
    {
        
        //xg.printdebug();
        this->XG_max[layer] = 0;
        //xg.printdebug();
        this->XG_max[layer] = P.por[layer] * depths[layer] * 1000.0;
        //xg.printdebug();
        
        this->theta[layer] = P.theta_default[layer];
        //xg.printdebug();
        if (rechrmax > 0.0)
        {
            //xg.printdebug();
            if (rechrmax > this->XG_max[layer])
            {
                //xg.printdebug();
                this->XG_rechr_d += depths[layer];
                this->rechr_fract[layer] = 1.0;
                rechrmax -= this->XG_max[layer];
            }
            else
            {
                //xg.printdebug();
                const double amount = rechrmax / this->XG_max[layer];
                this->rechr_fract[layer] = rechrmax / this->XG_max[layer];

                this->XG_rechr_d += depths[layer] * amount;
                const double amount_remaining = 1.0 - amount;
                //xg.printdebug();
                if (soilmax >= this->XG_max[layer]*amount_remaining)
                {
                    this->moist_fract[layer] = amount_remaining;
                    soilmax -= this->XG_max[layer] * amount_remaining;
                    this->XG_moist_d += depths[layer];
                }
                else
                {
                    this->moist_fract[layer] = (soilmax -  rechrmax) / this->XG_max[layer];
                    const double used = this->rechr_fract[layer] + this->moist_fract[layer];
                    this->default_fract[layer] = 1.0 - used;
                    this->XG_moist_d += this->XG_rechr_d + depths[layer] * used;
                    soilmax = 0.0;
                }
                rechrmax = 0.0;
            }
        }
        else if (soilmax > 0.0)
        {
            //xg.printdebug();
            if (soilmax >= this->XG_max[layer]) {
                this->XG_moist_d += depths[layer];
                this->moist_fract[layer] = 1.0;
                soilmax -= this->XG_max[layer];
            }
            else
            {
                const double amount = soilmax / this->XG_max[layer];
                this->XG_moist_d += depths[layer] * amount;
                this->moist_fract[layer] = amount;
                this->default_fract[layer] = 1.0 - amount;
                soilmax = 0.0;
            }
        }
        else
        {
            this->default_fract[layer] = 1.0;
        }
    }
    
    //xg.printdebug();
    if (rechrmax != 0.0 || soilmax != 0.0)
    {
        // put CHM exception here
    }
    return *this;
};

XG_algorithm::state& XG_algorithm::state::set_thermal_conductivities(const XG_algorithm::params& P,const double& soil_moist, const double& soil_rechr)
{
    Th_low = 1;
    Fz_low = 1;
    check_XG_moist = 0.0;
    XG_moist.resize(P.N_Soil_layers,0.0);

    // Process each soil layer
    for (size_t layer = 0; layer < P.N_Soil_layers; ++layer) 
    {
        // Calculate moisture content
        if (P.soil_moist_max > 0.0) 
        {  // handle soil_moist_max = 0.0 (slough case)
            XG_moist[layer] = 
                rechr_fract[layer] * XG_max[layer] * (soil_rechr / P.soil_rechr_max) +
                moist_fract[layer] * XG_max[layer] * 
                (soil_moist - soil_rechr) / (P.soil_moist_max - P.soil_rechr_max);
        }
        else
        {
            XG_moist[layer] = 0.0;
        }

        // Update moisture checks and adjustments
        check_XG_moist += XG_moist[layer];
        XG_moist[layer] += 
            default_fract[layer] * XG_max[layer] * P.theta_default[layer];


        // Calculate and validate theta
        theta[layer] = XG_moist[layer] / XG_max[layer];
        theta[layer] = (theta[layer] < P.theta_min) ? P.theta_min : theta[layer];
        //if (theta[layer] <= P.theta_min) 
        //{
        //    theta[layer] = P.theta_min;  // enforce minimum value
        //}

        // Convert to water content (kg/m³)
        layer_h2o[layer] = theta[layer] * P.por[layer] * 1000.0;
        // Update thermal conductivities
        if (P.k_update) 
        {  // dynamic update mode
            ftc[layer] = (ftc_contents[layer] == 1) 
                ? get_ttc(layer,P) 
                : get_ftc(layer,P);
            
            ttc[layer] = (ttc_contents[layer] == 1) 
                ? get_ftc(layer,P) 
                : get_ttc(layer,P);

        } 
        else 
        {
            ftc[layer] = get_ftc(layer,P);
            ttc[layer] = get_ttc(layer,P);
            ftc_contents[layer] = 0;
            ttc_contents[layer] = 0;
        }
    }
    return *this;
};

XG_algorithm::state& XG_algorithm::state::set_freezethaw_ratios(const XG_algorithm::params& P)
{
    for (size_t layer = 1; layer < P.N_Soil_layers; ++layer)
    {
        pf[layer] = std::sqrt(
                (ftc[layer-1] / layer_h2o[layer-1]) /
                (ftc[layer] / layer_h2o[layer])
                );

        pt[layer] = std::sqrt(
                (ttc[layer-1] / layer_h2o[layer-1]) /
                (ttc[layer] / layer_h2o[layer])
                );
    };

    return *this;
};

void XG_algorithm::state::accumulate_degree_days(const double& surface_temp)
{
    B += surface_temp / P.time_step_per_day;
    TrigAcc += B;

    t_trend -= t_trend / 192;
    t_trend += B/192;
};
//void XG_algorithm::init_freezethaw_front_depths(const double& Zdf_init, const double& Zdt_init, const double& Zpf_init)
//{
//	if (Zdf_init > 0.0 && Bfr == 0.0)
//		find_freeze_D(Zdf_init);
//
//	if (Zdt_init > 0.0 && Zdt_init < Zpf_init && Bth == 0.0)
//		find_thaw_D(Zdt_init);
//};

//void XG_algorithm::run()
//{
//
//    if (DTO.is_day_start(DTO))
//        B = 0.0;
//
//    B += DTO.t_surface;
//
//    if (STO.is_day_start(DTO))
//    {
//        // idle to freeze/thaw or freeze/thaw to idle
//        
//        // Keep T_total within +/- T_threshold
//        // TODO might not be needed, might be diagnositic like B
//        if (T_total > T_threshold)
//            T_total = T_threshold;
//        else if (T_total < -T_threshold)
//            T_total = -T_threshold;
//
//        if (std::abs(T_total) >= T_threshold/2.0 && freeze_thaw_state == 0)
//        {
//            if (freeze_front_depth < 0.0 || num_fronts)
//                freeze_thaw_state = 1;
//            else
//                freeze_thaw_state = -1;
//
//            temp_trend = 0.0;
//        }
//
//        if (freeze_thaw_state && std::abs(T_total) >= T_threshold/2.0 && std::abs(temp_trend) > 0.0)
//        {
//            freeze_thaw_state = 0;
//
//            if (freeze_front_depth > 0.0 && thaw_front_depth > 0.0)
//            {
//                if ( is_freezing() )
//                    switch_to_idle(thaw_front_depth);
//                else
//                    switch_to_idle(-freeze_front_depth);
//            }
//        }
//
//        // Calculate thermal conductivities
//
//        // Issues 
//
//        thaw_low_layer = 1;
//        freeze_low_layer = 1; 
//         
//        double thaw_k_T_rechr = get_k_T(DTO.soil_rechr_storage);
//        double freeze_k_T_rechr = ;
//        double thaw_k_T_lower = ;
//        double freeze_k_T_lower = ;
//         
//        for (int ii = 0; ii < N_layers; ++ii)
//        {
//            if (DTO.soil_storage_max > 0.0)
//            {
//                 
//            }
//
//        }    
//    }
//
//
//
//
//
//
//            
//};
//
//bool XG_algorithm::is_freezing()
//{
//    if (freeze_thaw_state == -1 && T_total > 0.0 && temp_trend > 0.0)
//        return true;
//    else if (freeze_thaw_state == 1 && T_total < 0.0 && temp_trend < 0.0)
//        return false;
//}
//
//void XG_algorithm::switch_to_idle(double front_depth)
//{
//    push_front(front_depth);// TODO write this function
//    // TODO in crhm both parts modify Zd_front_array (something I've decided not to include in XG) come back if necessary
//    if (front_depth > 0.0 && thaw_front_depth > freeze_front_depth)
//    {
//        thaw_front_depth = 0.0;
//        thaw_degree_days = 0.0;
//    }
//    else (front_depth < 0.0 && freeze_front_depth > thaw_front_depth)
//    {
//        freeze_front_depth = 0.0;
//        freeze_degree_days = 0.0;
//    }
//
//};
//
//
//
//
