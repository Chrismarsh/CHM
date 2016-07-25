#include "Simple_Canopy.hpp"


Simple_Canopy::Simple_Canopy(config_file cfg)
        : module_base(parallel::data)
{
    depends("p_rain");
    depends("p_snow");
    depends("iswr");
    depends("rh");
    depends("t");
    depends("vw"); // Wind speed at reference height Z_U_R [m/s]

    depends("vw_dir");
    depends("ilwr");
    depends("snowdepthavg");
    depends("snow_albedo"); //
    //depends("air_pressure"); // TODO: add to met interp (hard coded at 915 mbar for now)

    provides("Snow_load");
    provides("rain_load");
    provides("ta_subcanopy");
    provides("rh_subcanopy");
    provides("vw_subcanopy");
    //provides("wdir_subcanopy");
    provides("p_subcanopy");
    provides("p_rain_subcanopy");
    provides("p_snow_subcanopy");
    provides("frac_precip_rain_subcanopy");
    provides("frac_precip_snow_subcanopy");
    provides("iswr_subcanopy");
    provides("ilwr_subcanopy");
    //provides("diff_subcanopy");
    //provides("ir_h_subcanopy");

}

Simple_Canopy::~Simple_Canopy()
{

}

void Simple_Canopy::run(mesh_elem &elem, boost::shared_ptr <global> global_param)
{

    auto data = elem->get_module_data<Simple_Canopy::data>(ID);

    // Get meteorological data for current face
    double ta           = elem->face_data("t");
    double rh           = elem->face_data("rh");
    double vw           = elem->face_data("vw");
    double wdir         = elem->face_data("vw_dir");
    double iswr         = elem->face_data("iswr"); // SW in above canopy
    //double diff   = elem->face_data("iswr_diffuse"); // not used currently
    //double ir_h   = elem->face_data("iswr_direct"); // not used currently
    double ilwr         = elem->face_data("ilwr"); // LW in above canopy
    double p_rain       = elem->face_data("p_rain"); // rain (mm/timestep) above canopy
    double p_snow       = elem->face_data("p_snow"); // snow (mm/timestep) above canopy
    double snowdepthavg = elem->face_data("snowdepthavg");
    double Albedo       = elem->face_data("snow_albedo"); // Broad band snow albedo TODO: should be GROUND albedo (soil or snow), fractional as well
    double air_pressure = 915; //elem->face_data("air_pressure"); //"Average surface pressure", "(kPa)" TODO: Get from face_data

    // Checks on boundary conditions

    // Albedo
    if (Albedo == -9999) // If it is not defined TODO: current hack, should be initialized in mesher
        Albedo=0.1; // Assume no snow

    // Snow depth
    if (snowdepthavg == -9999) // If it is not defined TODO: current hack, should be initialized in mesher
        snowdepthavg = 0;

    // TODO:  Canopy parameters should come from mesh
    double Ht       = data->CanopyHeight; //", NHRU, "[0.1, 0.25, 1.0]", "0.001", "100.0", "forest/vegetation height", "(m)", &Ht);
    double LAI      = data->LAI; //", NHRU, "[2.2]", "0.1", "20.0", "leaf-area-index", "()", &LAI);

    /* Unsure what simple_canopy wind needs, read papers then code below. For now just give it wind at Z_U_R
    // Scale Wind speed
    if(Atmosphere::Z_U_R >= Ht)
    { // Forcing wind speed height greater than top of canopy height
        // Scale wind speed down to canopy top with logarithmic profile
        double vw_canopy = log_scale_wind(vw, Atmosphere::Z_U_R, Ht); // Wind speed at canopy top
        double Z_vw = Ht; // Update current height of vw

        // Scale through canopy using exponential profile to middle of canopy (one layer canopy) TODO: Scale to middle of canopy (need truck space)
        //vw_canopy = exp_scale_wind(vw);

    } else
    { // Forcing wind speed height less than than top of canopy height

        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Forcing height lower than canopy height, currently can't handel this situation."));

    }
     */

    // End Checks

    // Options for canopy representation TODO: get from cfg



    //ClassCRHMCanopy::run(void)

    // declared observations

    double Ts; //", NHRU, "snow surface temperature IN CANOPY", "(°C)", &Ts);

    //double Qnsn; //", NHRU, "net all-wave at snow surface", "(W/m^2)", &Qnsn);

    double Qnsn_Var; //", NHRU, "net all-wave at snow surface", "(W/m^2*int)", &Qnsn_Var);

    double Qsisn; //", NHRU, "incident short-wave at surface", "(W/m^2)", &Qsisn);

    double Qlisn; //", NHRU, "incident long-wave at surface", "(W/m^2)", &Qlisn);

    double Qlosn; //", NHRU, "reflected long-wave at surface", "(W/m^2)", &Qlosn);

    // declared variables

    double k; //"", NHRU, "extinction coefficient", "()", &k);

    double Tauc; //"", NHRU, "short-wave transmissivity", "(W/m^2)", &Tauc);

    double ra; //"", NHRU, "", "(s/m)", &ra);

    double drip_Cpy; //"", NHRU, "canopy drip", "(mm/int)", &drip_Cpy);

    double direct_rain; //", NHRU, "direct rainfall through canopy", "(mm/int)", &direct_rain);

    double net_rain; //"", NHRU, " direct_rain + drip", "(mm/int)", &net_rain);

    double cum_net_rain; //"", NHRU, " direct_rain + drip", "(mm)", &cum_net_rain);

    double Subl_Cpy; //"", NHRU, "canopy snow sublimation", "(mm/int)", &Subl_Cpy);

    //double cum_Subl_Cpy; //"", NHRU, "canopy snow sublimation", "(mm)", &cum_Subl_Cpy);

    double Pevap; //"", NHRU, "used when ground is snow covered to calculate canopy evaporation (Priestley-Taylor)", "(mm)", &Pevap);

    //double rain_load; //"", NHRU, "canopy rain load", "(mm)", &data->rain_load);

    //double Snow_load; //"", NHRU, "canopy snow load (timetep start)", "(mm)", &Snow_load);

    double direct_snow; //\", NHRU, "snow 'direct' Thru", "(mm/int)", &direct_snow);

    double SUnload; //"", NHRU, "unloaded canopy snow", "(mm)", &SUnload);

    double SUnload_H2O; //"", NHRU, "unloaded canopy snow as water", "(mm)", &SUnload_H2O);

   // double cum_SUnload_H2O; //"", NHRU, "Cummulative unloaded canopy snow as water", "(mm)", &cum_SUnload_H2O);

    double net_snow; //"", NHRU, "hru_snow minus interception", "(mm/int)", &net_snow);

    //double cum_net_snow; //"", NHRU, "Cummulative Canopy unload ", "(mm)", &cum_net_snow);

    double net_p; //"", NHRU, "total precipitation after interception", "(mm/int)", &net_p);

    double u_FHt; //"", NHRU, "wind speed at forest top (z = FHt)", "(m/s)", &u_FHt);

    double Cc; //"", NHRU, "Canopy coverage", "()", &Cc); UNITS???

    double intcp_evap; //"", NHRU, "HRU Evaporation from interception", "(mm/int)", &intcp_evap);

    //double cum_intcp_evap; //", NHRU, "HRU Evaporation from interception", "(mm)", &cum_intcp_evap);

    double Kstar_H;

    // Parameters used

    // Default values used for now
    double Alpha_c          = Vegetation::alb_c; // "canopy albedo" 0.05-0.2
    double B_canopy         = 0.038; //TODO: What is this? Where does it come from?", NHRU, "[0.038]", "0.0", "0.2", "canopy enhancement parameter. Suggestions are Colorado - 0.23 and Alberta - 0.038", "()", &B_canopy);
    double Zref             = 2; //", "0.01", "100.0", "temperature measurement height", "(m)", &Zref); TODO: Take from config
    double Zwind            = Atmosphere::Z_U_R; //", "0.01", "100.0", "wind measurement height", "(m)", &Zwind); // Set as defined above
    double Z0snow           = Snow::Z0_SNOW; //", "0.0001", "0.01", "snow roughness length", "(m)", &Z0snow);
    double Sbar             = 6.6; //", "0.0", "100.0", "maximum canopy snow interception load", "(kg/m^2)", &Sbar);
    double Zvent            = 0.75; //", "0.0", "1.0", "ventilation wind speed height (z/Ht)", "()", &Zvent);
    double unload_t         = 1.0; //", "-10.0", "20.0", "if ice-bulb temp >= t : canopy snow is unloaded as snow", "(°C)", &unload_t);
    double unload_t_water   = 4.0; //", "-10.0", "20.0", "if ice-bulb temp >= t: canopy snow is unloaded as water", "(°C)", &unload_t_water);

    double SolAng           = global_param->solar_el() * mio::Cst::to_rad; // degrees to radians

    // Initialize
    net_rain = 0.0;
    direct_rain = 0.0;
    drip_Cpy = 0.0;
    intcp_evap = 0.0;
    direct_snow = 0.0;
    SUnload = 0.0;
    SUnload_H2O = 0.0;
    Subl_Cpy = 0.0;

    // Canopy temperature is first approximated by the air temperature.
    double T1 = ta + mio::Cst::t_water_freezing_pt; // Canopy temperature (C to K)

    // Get Exposure (how much is canopy above snowdepth)
    double Exposure = Ht - snowdepthavg; // (m)
    if(Exposure < 0.0)
        Exposure = 0.0;


    double LAI_ = LAI*Exposure/Ht; // Rescaling LAI???

    // terrain view factor (equivalent to 1-Vf), where Vf is the sky view factory. // Where does this equation come from?
    double Vf = 0.45 - 0.29*log(LAI);

    // Rescaleing Vf ???
    double Vf_ = Vf + (1.0 - Vf)*sin((Ht - Exposure)/Ht*M_PI_2); // Where does equation come from?

    // What is this doing??? // TODO: add slope in
    if(SolAng > 0.001) { // radians
        k    = 1.081*SolAng*cos(SolAng)/sin(SolAng); // "extinction coefficient"
        Tauc = exp(-k*LAI_); // "short-wave transmissivity", "(W/m^2)"
    }
    else {
        k    = 0.0; // "extinction coefficient"
        Tauc = 0.0; // "short-wave transmissivity", "(W/m^2)"
    }


    Kstar_H = iswr*(1.0 - Alpha_c - Tauc*(1.0 - Albedo)); //  what is Kstar_H???

    // Incident long-wave at surface, "(W/m^2)"
    Qlisn = ilwr*Vf_ + (1.0 - Vf_)*Vegetation::emiss_c*PhysConst::sbc*pow(T1, 4.0) + B_canopy*Kstar_H;

    // Incident short-wave at surface, "(W/m^2)"
    Qsisn = iswr*Tauc;


    double rho = air_pressure*1000/(PhysConst::Rgas*T1); // density of Air (pressure kPa to Pa = *1000)

    double U1 = vw; // Wind speed (m/s) at height Z_vw [m] (top of canopy)

    // Aerodynamic resistance of canopy
    ra = (log(Zref/Z0snow)*log(Zwind/Z0snow))/pow(PhysConst::kappa,2)/U1; // (s/m)

    double deltaX = 0.622*PhysConst::Ls*Qs(air_pressure, T1)/(PhysConst::Rgas*(pow(T1,2))); // Must be (kg K-1)

    double q = (rh/100)*Qs(air_pressure, T1); // specific humidity (kg/kg)

    // snow surface temperature of snow in canopy TODO: Check units
    Ts = T1 + (Snow::emiss*(ilwr - PhysConst::sbc*pow(T1, 4.0)) + PhysConst::Ls*(q - Qs(air_pressure, T1))*rho/ra)/
                  (4.0*Snow::emiss*PhysConst::sbc*pow(T1, 3.0) + (PhysConst::Cp + PhysConst::Ls*deltaX)*rho/ra);

    Ts -= mio::Cst::t_water_freezing_pt; // K to C

    // Check if Ts is above freezing or there is no snow TODO: don't understand logic here, why does it matter if there is snow on ground? negative snow depth?
    if(Ts > 0.0 || snowdepthavg <= 0.0)
        Ts = 0.0;

    // Emitted long-wave from snowpack in canopy (downwards and upwards?)
    Qlosn = Snow::emiss*PhysConst::sbc*pow(Ts + mio::Cst::t_water_freezing_pt, 4.0);

    // Net radiation for ground snowpack surface
    //Qnsn = Qlisn - Qlosn + Qsisn*(1.0 - Albedo); // Not used here

//==============================================================================
// coupled forest snow interception and sublimation routine:
// after Hedstom & Pomeroy 1998?/ Parviainen & Pomeroy 2000:
// calculate maximum canopy snow load (L*):

    if(data->Snow_load > 0.0 || p_snow > 0.0) { // handle snow
        double RhoS = 67.92 + 51.25 * exp(ta / 2.59);
        double LStar = Sbar * (0.27 + 46.0 / RhoS) * LAI;

        if (data->Snow_load > LStar) { // after increase in temperature
            direct_snow = data->Snow_load - LStar;
            data->Snow_load = LStar;
        }

        // calculate intercepted snowload

        // Calculate wind speed at canopy top // TODO: use U1 instead of vw here?
        if (Ht - 2.0 / 3.0 * Zwind > 1.0) // Find source of equations
            u_FHt = vw * log((Ht - 2.0 / 3.0 * Zwind) / 0.123 * Zwind) /
                    log((Zwind - 2.0 / 3.0 * Zwind) / 0.123 * Zwind);
        else
            u_FHt = 0.0;

        double I1 = 0.0; // new Interecption ? [mm?]

        // calculate horizontal canopy-coverage (Cc):

        Cc = 0.29 * log(LAI) + 0.55; //TODO: Where do hard coded param comes from?
        if (Cc <= 0.0) {
            Cc = 0.0;
        } else if (Cc > 1.0) {
            Cc = 1.0;
        }

        if(p_snow > 0.0 && fabs(p_snow/LStar) < 50.0){ //TODO: Where do hard coded param comes from?
            if (u_FHt <= 1.0)  // if wind speed at canopy top > 1 m/s //TODO: Where do hard coded param comes from?
                I1 = (LStar-data->Snow_load)*(1.0-exp(-Cc*p_snow/LStar));
            else
                I1 = (LStar-data->Snow_load)*(1.0-exp(-p_snow/LStar));

            if(I1 <= 0)
                I1 = 0;

            data->Snow_load += I1;

            // calculate canopy snow throughfall before unloading:

            direct_snow += (p_snow - I1);
        }

        // calculate snow ventilation windspeed:

        const double gamma = 1.15; // What is this param?
        double xi2 = 1-Zvent;
        double windExt2 = (gamma * LAI * xi2);

        double uVent = u_FHt * exp(-1 * windExt2);


// calculate sublimation of intercepted snow from ideal intercepted ice sphere (500 microns diameter):

        double Alpha, A1, B1, C1, J, D, Lamb, Mpm, Nu, Nr, SStar, Sigma2;

        double Es = 611.15 * exp(22.452*ta/(ta + 273.0));  // {sat pressure}

        double SvDens = Es*PhysConst::M/(PhysConst::R*(ta + 273.0)); // {sat density}

        Lamb = 6.3e-4*(ta+273.0) + 0.0673;  // thermal conductivity of atmosphere
        Nr = 2.0 * Snow::Radius * uVent / Atmosphere::KinVisc;  // Reynolds number
        Nu = 1.79 + 0.606 * sqrt(Nr); // Nusselt number
        SStar = M_PI * pow(Snow::Radius,2) * (1.0 - Snow::AlbedoIce) * iswr;  // SW to snow particle !!!! changed
        A1 = Lamb * (ta + 273) * Nu;
        B1 = PhysConst::Ls * PhysConst::M /(PhysConst::R * (ta + 273.0))- 1.0;
        J = B1/A1;
        Sigma2 = rh/100 - 1;
        D = 2.06e-5* pow((ta+273.0)/273.0, -1.75); // diffusivity of water vapour
        C1 = 1.0/(D*SvDens*Nu);

        Alpha = 5.0;
        Mpm = 4.0/3.0 * M_PI * PhysConst::DICE * pow(Snow::Radius,3) *(1.0 + 3.0/Alpha + 2.0/pow(Alpha,2));

// sublimation rate of single 'ideal' ice sphere:

        double Vs = (2.0* M_PI* Snow::Radius*Sigma2 - SStar* J)/(PhysConst::Ls* J + C1)/Mpm;

// snow exposure coefficient (Ce):

        double Ce;
        if ((data->Snow_load/LStar) <= 0.0)
            Ce = 0.07;
        else
            Ce = Vegetation::ks* pow((data->Snow_load/LStar), -Vegetation::Fract);

// calculate 'potential' canopy sublimation:

        double Vi = Vs*Ce;

// calculate 'ice-bulb' temperature of intercepted snow:

        double IceBulbT = ta - (Vi* PhysConst::Ls/1e6/PhysConst::Ci);

// determine whether canopy snow is unloaded:

        if(IceBulbT >= unload_t){
            if(IceBulbT >= unload_t_water){
                drip_Cpy = data->Snow_load;
                SUnload_H2O = data->Snow_load;
            }
            else {
                SUnload = data->Snow_load*(IceBulbT - unload_t)/(unload_t_water - unload_t);
                drip_Cpy = data->Snow_load - SUnload;
                SUnload_H2O = drip_Cpy;
            }

            data->Snow_load = 0.0;
            data->cum_SUnload_H2O += SUnload_H2O;
        }

// limit sublimation to canopy snow available and take sublimated snow away from canopy snow at timestep start

        //Subl_Cpy = -data->Snow_load*Vi*Hs*Global::Interval*24*3600/Hs; // make W/m2 (original in CRHM)
        Subl_Cpy = -data->Snow_load*Vi*PhysConst::Ls*global_param->dt()/PhysConst::Ls; // make W/m2 TODO: check Interval is same as dt() (in seconds
        // TODO: Hs/HS = 1 !!!
        
        if(Subl_Cpy > data->Snow_load){
            Subl_Cpy = data->Snow_load;
            data->Snow_load = 0.0;
        }
        else{
            data->Snow_load -= Subl_Cpy;
            if(data->Snow_load < 0.0)
                data->Snow_load = 0.0;
        }

// calculate total sub-canopy snow:

        net_snow = direct_snow + SUnload;

    } // handle snow



    double smax = Cc*LAI*0.2;

//  Forest rain interception and evaporation model
// 'sparse' Rutter interception model (i.e. Valente 1997):

// calculate direct throughfall:

    if(p_rain > 0.0){

        direct_rain = p_rain * (1-Cc);

        // calculate rain accumulation on canopy before evap loss:

        if (data->rain_load + p_rain*Cc > smax){
            drip_Cpy += (data->rain_load + p_rain*Cc - smax);
            data->rain_load = smax;
        }
        else
            data->rain_load += p_rain*Cc;

    }

// calculate 'actual evap' of water from canopy and canopy storage after evaporation::
    // TODO: I don't understand why two different potential evaporation calcs are used for snow in canopy or no snow in canopy?
    // For now just use PT method for all cases because I don't want to include the Classevap CRHM class


    // If Liquid water exists in canopy
    if(data->rain_load > 0.0){
        /*if(data->Snow_load == 0){ // use Granger when no snowcover IN CANOPY // Changed to use Snow_load to check if canopy has snow
            if(data->rain_load >= hru_evap*Cc){ // (evaporation in mm)
                intcp_evap = hru_evap*Cc;  //
                data->rain_load -= hru_evap*Cc;
            }
            else{
                intcp_evap = data->rain_load;
                data->rain_load = 0.0;
            }
        }
        else{ */// use Priestley-Taylor when snowcover IN CANOPY
            //double Q = iswr*86400/Global::Freq/1e6/lambda(ta); // convert w/m2 to mm/m^2/int (original CRHM)
            double temp_Global_Freq = global_param->dt()/86400; // time steps per day (following CRHM convention)
            double Q = iswr*86400/temp_Global_Freq/1e6/lambda(ta); // convert w/m2 to mm/m^2/int TODO: Units don't make sense here (missing density of water??)

            if(iswr > 0.0)
                Pevap = 1.26*delta(ta)*Q/(delta(ta) + gamma(air_pressure, ta));
            else
                Pevap = 0.0;

            if(data->rain_load >= Pevap*Cc){  // (evaporation in mm)
                intcp_evap = Pevap*Cc;  // check
                data->rain_load -= Pevap*Cc;
            }
            else{
                intcp_evap = data->rain_load; // check
                data->rain_load = 0.0;
            }
        /*}*/
    } // if data->rain_load > 0.0

    // cumulative amounts....

    net_rain = direct_rain + drip_Cpy;
    net_p = net_rain + net_snow;

    data->cum_net_rain += net_rain;
    data->cum_net_snow += net_snow;
    data->cum_intcp_evap += intcp_evap;
    data->cum_Subl_Cpy += Subl_Cpy;

    // Output computed canopy states and fluxes downward to snowpack and upward to atmosphere
    elem->set_face_data("Snow_load",data->Snow_load);
    elem->set_face_data("rain_load",data->rain_load);
    //elem->set_face_data("cum_Subl_Cpy",data->cum_Subl_Cpy);
    elem->set_face_data("ta_subcanopy",ta);
    elem->set_face_data("rh_subcanopy",rh);
    elem->set_face_data("vw_subcanopy",vw);
    elem->set_face_data("iswr_subcanopy",Qsisn); // (W/m^2)
    elem->set_face_data("ilwr_subcanopy",Qlisn); // (W/m^2)
    elem->set_face_data("p_rain_subcanopy",net_rain); // (mm/int)
    elem->set_face_data("p_snow_subcanopy",net_snow); // (mm/int)
    elem->set_face_data("p_subcanopy",net_p); // Total precip (mm/int)
    elem->set_face_data("frac_precip_rain_subcanopy",net_rain/net_p); // Fraction rain (-)
    elem->set_face_data("frac_precip_snow_subcanopy",net_snow/net_p); // Fraction snow (-)

}

void Simple_Canopy::init(mesh domain, boost::shared_ptr <global> global_param)
{
    // TODO: Parallel call needed here?
    //#pragma omp parallel for


    // For each face
    for(size_t i=0;i<domain->size_faces();i++)
    {
        // Get current face
        auto face = domain->face(i);

        auto d = face->make_module_data<Simple_Canopy::data>(ID);

        // Check if Canopy exists at this face/triangle
        double LC; // landcover type

        //TODO: Why is LC set to zero? why is has_parameter("landcover") false?

        LC = face->get_parameter("landcover");
        bool LCexists = face->has_parameter("landcover");
        if(face->has_parameter("landcover")) {
            LC = face->get_parameter("landcover");
            std::string LC_string = std::to_string(LC);
            d->HasCanopy    = global_param->parameters.get<bool>("landcover." + std::to_string(LC) + ".canopy");
            if(d->HasCanopy) {
                d->CanopyHeight = global_param->parameters.get<double>(
                        "landcover." + std::to_string(LC) + ".CanopyHeight");
                d->LAI = global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".LAI");


                // ClassCRHMCanopy
                d->rain_load = 0.0;
                d->Snow_load = 0.0;
                d->cum_net_snow = 0.0; // "Cumulative Canopy unload ", "(mm)"
                d->cum_net_rain = 0.0; // " direct_rain + drip", "(mm)"
                d->cum_Subl_Cpy = 0.0; //  "canopy snow sublimation", "(mm)"
                d->cum_intcp_evap = 0.0; // "HRU Evaporation from interception", "(mm)"
                d->cum_SUnload_H2O = 0.0; // "Cumulative unloaded canopy snow as water", "(mm)"


            }
        } else {
            // TODO: throw error here
        }


        //d->LAI          = cfg.get<double>("canopy.CanopyLeafAreaIndex");
        //d->CanopyHeight = cfg.get<double>("canopy.CanopyHeight");


        // Initialize canopy state variables
        // For each canopy layer
        // Snow in branches average bulk temperature




    }
}

double Simple_Canopy::delta(double ta) // Slope of sat vap p vs t, kPa/°C
{
    if (ta > 0.0)
        return(2504.0*exp( 17.27 * ta/(ta+237.3)) / pow(ta+237.3,2));
    else
        return(3549.0*exp( 21.88 * ta/(ta+265.5)) / pow(ta+265.5,2));
}

double Simple_Canopy::lambda(double ta) // Latent heat of vaporization (mJ/(kg °C))
{
    return( 2.501 - 0.002361 * ta );
}

double Simple_Canopy::gamma(double air_pressure, double ta) // Psychrometric constant (kPa/°C)
{
    return( 0.00163 * air_pressure / lambda(ta)); // lambda (mJ/(kg °C))
}

double Simple_Canopy::Qs(double air_pressure, double T1) {
    /* INPUT
    air_pressure    - unadjusted air pressure in Pa
    T1              - Air temperature in K
    // OUTPUT
    Qs              - Saturated mixing ratio (kg/kg)
     */
    T1 = T1 - mio::Cst::t_water_freezing_pt; // K to C
    double es = 611.213*exp(22.4422*T1/(272.186+T1)); // Pa
    return(0.622 * ( es / (air_pressure - es) )); // kg/kg TODO: Check Qs calc
}



