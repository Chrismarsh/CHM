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

#include "snowpack.hpp"
REGISTER_MODULE_CPP(Lehning_snowpack);

Lehning_snowpack::Lehning_snowpack(config_file cfg)
        : module_base("Lehning_snowpack", parallel::data, cfg)
{
    depends("iswr");
    depends("ilwr");
    depends("rh");
    depends("t");
    depends("U_2m_above_srf");
    depends("p");
    depends("frac_precip_rain");

    depends("snow_albedo");

    // Optional subcanopy variables if a canopy module is included (used if exist)
    optional("ta_subcanopy");
    optional("rh_subcanopy");
    optional("p_subcanopy");
    optional("frac_precip_rain_subcanopy");
    optional("iswr_subcanopy");
    optional("ilwr_subcanopy");
    optional("drift_mass");

    optional("T_g");

    provides("dQ");
    provides("swe");
    provides("T_s");
    provides("T_s_0");
    provides("n_nodes");
    provides("n_elem");
    provides("snowdepthavg");

//    if(!has_optional("snow_albedo"))
//        provides("snow_albedo");

    provides("H");
    provides("E");
    provides("G");
    provides("ilwr_out");
    provides("iswr_out");
    provides("R_n");
    provides("runoff");
    provides("mass_snowpack_removed");
    provides("sum_runoff");
    provides("sum_subl");
    provides("sublimation");
    provides("evap");

    provides("MS_SWE");
    provides("MS_WATER");
    provides("MS_TOTALMASS");
    provides("MS_SOIL_RUNOFF");
    
    provides("snow_sphericity");
    provides("snow_grain_size");
    provides("frac_ice_content");
    provides("Sliq");
    provides("Tsnow");

}

Lehning_snowpack::~Lehning_snowpack()
{

}

void Lehning_snowpack::run(mesh_elem &face)
{
    if(is_water(face))
    {
        set_all_nan_on_skip(face);
        return;
    }
    auto& data = face->get_module_data<Lehning_snowpack::data>(ID);

    /**
     * Builds this timestep's meteo data
     */
    CurrentMeteo Mdata(data.config);
    Mdata.date   =  mio::Date( global_param->year(), global_param->month(), global_param->day(),global_param->hour(),global_param->min(),-6 );
    // Optional inputs if there is a canopy or not
    if(has_optional("ta_subcanopy")) {
        Mdata.ta     =  (*face)["ta_subcanopy"_s]+mio::Cst::t_water_freezing_pt;
    } else {
        Mdata.ta     =  (*face)["t"_s]+mio::Cst::t_water_freezing_pt;
    }

    if(has_optional("rh_subcanopy")) {
        Mdata.rh     =  (*face)["rh_subcanopy"_s]/100.;
    } else {
        Mdata.rh     =  (*face)["rh"_s]/100.;
    }


    Mdata.vw     =  (*face)["U_2m_above_srf"_s];

    Mdata.vw_max = mio::IOUtils::nodata;// TODO: fix md(MeteoData::VW_MAX);

    Mdata.vw_drift = Mdata.vw ;//mio::IOUtils::nodata;
    Mdata.dw_drift = Mdata.dw;//0;


    if(has_optional("iswr_subcanopy")) {
        Mdata.iswr     =  std::max(0.0,(*face)["iswr_subcanopy"_s]);
    } else {
        Mdata.iswr     =  std::max(0.0,(*face)["iswr"_s]);
    }

    // If  Snowpack, "SW_MODE" : "BOTH"  then rswr and iswr needs to be definined.
    // If an external albedo is not used, a parametrized one is used. but rswr and iswr both must be defined.
    // If "SW_MODE" : "INCOMING", is used, then mAlbedo needs to be undefined to trigger the appropriate rswr and albedo calculations.
    // rswr must still be set, but we use the garbage that is in Xdata instead.
//    if(has_optional("snow_albedo"))
//    {
        //measured albedo in snowpack will be fed from an albedo model
        //in the config 'both' will enable this
        Mdata.mAlbedo   =  (*face)["snow_albedo"_s];
        Mdata.rswr      =  std::max(0.0,(*face)["snow_albedo"_s] * Mdata.iswr);

//    }
//    else
//    {
//        Mdata.rswr = std::max(0.0,data.Xdata->Albedo * Mdata.iswr);
//        Mdata.mAlbedo = Constants::undefined; //this will trigger calculating a paramaterized albedo
//    }


    if(has_optional("ilwr_subcanopy")) {
        Mdata.ilwr     =  (*face)["ilwr_subcanopy"_s];
    } else {
        Mdata.ilwr     =  (*face)["ilwr"_s];
    }

    Mdata.ea =  mio::Atmosphere::blkBody_Emissivity(Mdata.ilwr, Mdata.ta); //follows apline3d
//    Mdata.ea  = SnLaws::AirEmissivity(Mdata.ilwr, Mdata.ta, "default"); //atmospheric emissivity!
    // Note this is used later throughout the code line 673 in Snowpack.cc, dealing with the emissivity of the snowpack! (might be a bug).
    // Left the ea calculation here for now - NIC

    //double thresh_rain = 2 + mio::Cst::t_water_freezing_pt;

    // Define fraction rain and total precipitaiton (soild and liquid)
    if(has_optional("p_subcanopy")) {
        Mdata.psum_ph = (*face)["frac_precip_rain_subcanopy"_s]; //  0 = snow, 1 = rain
        Mdata.psum = (*face)["p_subcanopy"_s];
    } else {
        Mdata.psum_ph = (*face)["frac_precip_rain"_s]; //  0 = snow, 1 = rain
        Mdata.psum = (*face)["p"_s];
    }

    Mdata.rho_hn = 100;

    //setup a single ground temp measurement
    mio::MeteoData soil_meas;
//    soil_meas.addParameter("HTS1");
//    soil_meas.addParameter("TS1");
//    soil_meas("HTS1") = -0.1;
//    soil_meas("TS1") = 269;
//    Mdata.setMeasTempParameters(soil_meas);
//    Mdata.ts.push_back(soil_meas("TS1"));
//    Mdata.zv_ts.push_back(soil_meas("HTS1"));

    Mdata.tss = data.Xdata->Ndata[data.Xdata->getNumberOfElements()].T;  //we use previous timestep value//mio::IOUtils::nodata; //Constants::undefined;;//
    //setting this to tss is inline with Alpine3d if there is no soil node. However, it might make more sense to use a const ground temp?
    if(has_optional("T_g"))
        Mdata.ts0 = (*face)["T_g"_s] + 273.15;
    else
        Mdata.ts0 =  const_T_g + 273.15; //Mdata.tss <- is in line with Alpine3D, but this makes no sense. So use a const temp.

    Mdata.hs = mio::IOUtils::nodata; //follows alpine3d
    Mdata.elev      = (*face)["solar_el"_s]*mio::Cst::to_rad;

    data.cum_precip  += Mdata.psum; //running sum of the precip. snowpack removes the rain component for us.
    data.meteo->compMeteo(Mdata,*(data.Xdata),false,false); // no canopy model

    double mass_erode = 0;

    if(has_optional("drift_mass"))
    {

        double mass = (*face)["drift_mass"_s];
        mass = is_nan(mass) ? 0 : mass;

        if(mass > 0)
        {

            Mdata.psum += mass;
            data.cum_precip  += Mdata.psum;
            Mdata.psum_ph = 0;
        }
        else
        {
            mass_erode = mass; // snowpack expects the mass erode to be negative
        }
    }

    // To collect surface exchange data for output
    SurfaceFluxes surface_fluxes;
//    surface_fluxes.reset(false);
//    surface_fluxes.drift = 0.;
//    surface_fluxes.mass[SurfaceFluxes::MS_WIND] = 0.;


    // Boundary condition (fluxes)
    BoundCond Bdata;

    try
    {
        data.sp->runSnowpackModel(Mdata, *(data.Xdata), data.cum_precip, Bdata,surface_fluxes,mass_erode);
        surface_fluxes.collectSurfaceFluxes(Bdata, *(data.Xdata), Mdata);
    }catch(...)
    {
        if (data.Xdata->swe > 3)
        {

            auto details = "(" + std::to_string(face->center().x()) +
                           "," + std::to_string(face->center().y()) +
                           "," + std::to_string(face->center().z())
                           + ") ID = " + std::to_string(face->cell_local_id);
            CHM_THROW_EXCEPTION(module_error, "Snowpack died. Triangle center = " + details);
        }
    }


    (*face)["sublimation"_s]=surface_fluxes.mass[SurfaceFluxes::MS_SUBLIMATION];

    if(data.Xdata->swe > 0)
    {

        double bulk_T_s=0;
        for(size_t i = 0; i < data.Xdata->getNumberOfElements(); ++i)
        {
            bulk_T_s += data.Xdata->Ndata[i].T;
        }

        bulk_T_s /= data.Xdata->getNumberOfElements();

        (*face)["T_s"_s]=bulk_T_s;
        (*face)["T_s_0"_s]=Mdata.tss;
        (*face)["n_nodes"_s]=data.Xdata->getNumberOfNodes();
        (*face)["n_elem"_s]=data.Xdata->getNumberOfElements();
        (*face)["H"_s]=surface_fluxes.qs;
        (*face)["E"_s]=surface_fluxes.ql;


        (*face)["G"_s]=surface_fluxes.qg0 == -999 ? 0 : surface_fluxes.qg0; //qg0 is the correct ground heatflux to match snowpack output. qg is just uninit

        (*face)["ilwr_out"_s]=-surface_fluxes.lw_out; //these are actually positive!
        (*face)["iswr_out"_s]=-surface_fluxes.sw_out;
        (*face)["R_n"_s]= surface_fluxes.lw_net + (surface_fluxes.sw_in-surface_fluxes.sw_out );
        (*face)["dQ"_s]=surface_fluxes.dIntEnergy;
//        if(!has_optional("snow_albedo"))
//        {
            (*face)["snow_albedo"_s]=data.Xdata->Albedo;  //even if we have a measured albedo, Xdata will reflect this. //surface_fluxes.pAlbedo);
//        }

        // New parameters
        if(data.Xdata->getNumberOfElements() > 0)
        {
            // Top element is at index nElems-1 (elements are ordered from ground (0) to surface)
            const auto& top_element = data.Xdata->Edata[data.Xdata->getNumberOfElements() - 1];
            (*face)["snow_sphericity"_s] = top_element.sp;
            (*face)["snow_grain_size"_s] = top_element.rg;
            (*face)["frac_ice_content"_s] = top_element.theta[ICE]; // ICE is defined as 1 in DataClasses.h
            (*face)["Sliq"_s] = top_element.theta[WATER]; // WATER is defined as 2 in DataClasses.h
        }
        
        if(data.Xdata->getNumberOfNodes() > 0)
        {
            // Top node is at index nNodes-1 (nodes are ordered from ground (0) to surface)
            const auto& top_node = data.Xdata->Ndata[data.Xdata->getNumberOfNodes() - 1];
            (*face)["Tsnow"_s] = top_node.T - mio::Cst::t_water_freezing_pt; // Convert from K to °C
        }

    } else{
       set_all_nan_on_skip(face);

    }

    //always write out 0 swe regardless of amount of swee
    (*face)["swe"_s]=data.Xdata->swe;
    (*face)["mass_snowpack_removed"_s]=data.Xdata->ErosionMass;
    (*face)["snowdepthavg"_s]=data.Xdata->cH - data.Xdata->Ground; // cH includes soil depth if SNP_SOIL == 1, hence subtracting Ground height
    (*face)["runoff"_s]=surface_fluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF];
    (*face)["evap"_s]=surface_fluxes.mass[SurfaceFluxes::MS_EVAPORATION];
//    (*face)["sum_runoff"_s]=  (*face["sum_runoff"_s] + surface_fluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF]);
    data.sum_subl +=surface_fluxes.mass[SurfaceFluxes::MS_SUBLIMATION];
    (*face)["sum_subl"_s]=   data.sum_subl ;


    (*face)["MS_SWE"_s]=surface_fluxes.mass[SurfaceFluxes::MS_SWE];
    (*face)["MS_WATER"_s]=surface_fluxes.mass[SurfaceFluxes::MS_WATER];
    (*face)["MS_TOTALMASS"_s]=surface_fluxes.mass[SurfaceFluxes::MS_TOTALMASS];
    (*face)["MS_SOIL_RUNOFF"_s]=surface_fluxes.mass[SurfaceFluxes::MS_SOIL_RUNOFF];
}

void Lehning_snowpack::init(mesh& domain)
{
    const_T_g = cfg.get("const_T_g",-4.0);

    for(size_t i=0;i<domain->size_local_faces();i++)
    {
        auto face = domain->face(i);

        auto& d = face->make_module_data<Lehning_snowpack::data>(ID);


        //setup critical keys.
        //overwrite the user if a dangerous key is set
        d.config.addKey("METEO_STEP_LENGTH", "Snowpack", std::to_string( 3600.0 / global_param->dt())); // Hz. Number of met per hour
        d.config.addKey("MEAS_TSS", "Snowpack", "false");

        //specified as minutes, snowpack will convert to s for us. CHM dt is in s
        d.config.addKey("CALCULATION_STEP_LENGTH","Snowpack", std::to_string(global_param->dt()  / 60 ) );
        //default values for
        //	"Snowpack": { }

        d.config.addKey("MEAS_TSS","Snowpack","false");
        d.config.addKey("ENFORCE_MEASURED_SNOW_HEIGHTS","Snowpack","false");
        d.config.addKey("SW_MODE","Snowpack","BOTH");
        d.config.addKey("HEIGHT_OF_WIND_VALUE","Snowpack","2");
        d.config.addKey("HEIGHT_OF_METEO_VALUES","Snowpack","2");
        d.config.addKey("ATMOSPHERIC_STABILITY","Snowpack","MO_MICHLMAYR");
        d.config.addKey("ROUGHNESS_LENGTH","Snowpack","0.001");
        d.config.addKey("CHANGE_BC","Snowpack","false");
        d.config.addKey("THRESH_CHANGE_BC","Snowpack","-1.0");
        d.config.addKey("SNP_SOIL","Snowpack","false");
        d.config.addKey("SOIL_FLUX","Snowpack","false");
        d.config.addKey("GEO_HEAT","Snowpack","0.06");
        d.config.addKey("CANOPY","Snowpack","false");

        //default values for
        //	"SnowpackAdvanced": { }
        d.config.addKey("MAX_NUMBER_MEAS_TEMPERATURES","SnowpackAdvanced","1");
        d.config.addKey("ALPINE3D","SnowpackAdvanced","true"); //must be true for any blowing snow module
        d.config.addKey("SNOW_EROSION","SnowpackAdvanced","false");
        d.config.addKey("MEAS_INCOMING_LONGWAVE","SnowpackAdvanced","true");
        d.config.addKey("THRESH_RAIN","SnowpackAdvanced","2");
        d.config.addKey("THRESH_RAIN_RANGE","SnowpackAdvanced","2");
        d.config.addKey("WATERTRANSPORTMODEL_SNOW","SnowpackAdvanced","BUCKET");
        d.config.addKey("VARIANT","SnowpackAdvanced","DEFAULT");
        d.config.addKey("ADJUST_HEIGHT_OF_WIND_VALUE","SnowpackAdvanced","false"); // we always provide a 2m wind, even if there is snowcover
        d.config.addKey("HN_DENSITY","SnowpackAdvanced","MEASURED"); //We can then set it in at run time. Do it this way so we can have temporally variable if we want.

        d.config.addKey("ADVECTIVE_HEAT","SnowpackAdvanced","TRUE"); // This enables heat transport caused by moving water. Error if true/false not specified.
        
        d.config.addKey("COMBINE_ELEMENTS","SnowpackAdvanced","true"); //Defines whether joining elements will be considered at all
        //Activates algorithm to reduce the number of elements deeper in the snowpack AND to split elements again when they come back to the surface
        //Only works when COMBINE_ELEMENTS == TRUE.
        d.config.addKey("REDUCE_N_ELEMENTS","SnowpackAdvanced","true");

        double grooming_week_start = cfg.get("GROOMING_WEEK_START", 40); //First week of grooming
        double grooming_week_end = cfg.get("GROOMING_WEEK_END", 17); //Last week of grooming
        double grooming_hour = cfg.get("GROOMING_HOUR", 21); //Hour at which grooming is performed
        double grooming_depth_start = cfg.get("GROOMING_DEPTH_START", 0.4); // How much snow on the ground to start grooming
        double grooming_depth_impact = cfg.get("GROOMING_DEPTH_IMPACT", 0.4); //maximum depth of snow impacted by grooming
        // Check for grooming parameter in domain - if true, enable SNOW_GROOMING
        bool is_grooming = false;
        if(face->has_parameter("grooming")){
            is_grooming = face->parameter("grooming"_s);
        } else if(face->has_parameter("Resort")){ // Check if Resort parameter exists and use it to look up grooming in parameter mapping (similar to landcover for SimpleCanopy, look SnowCast config)
            int resort_type = face->parameter("Resort"_s);
            try {
                is_grooming = global_param->parameters.get<bool>("Resort." + std::to_string(resort_type) + ".grooming");
            } catch(const boost::property_tree::ptree_bad_path& e) {
                SPDLOG_ERROR("No grooming parameter defined for Resort type: {}", resort_type);
            }
        }
        
        if(is_grooming)
        {
            d.config.addKey("SNOW_GROOMING","TechSnow","true");
            d.config.addKey("GROOMING_WEEK_START","TechSnow",std::to_string(grooming_week_start));
            d.config.addKey("GROOMING_WEEK_END","TechSnow",std::to_string(grooming_week_end));
            d.config.addKey("GROOMING_HOUR","TechSnow",std::to_string(grooming_hour));
            d.config.addKey("GROOMING_DEPTH_START","TechSnow",std::to_string(grooming_depth_start));
            d.config.addKey("GROOMING_DEPTH_IMPACT","TechSnow",std::to_string(grooming_depth_impact));
        }



        // because we use our own config, we need to do the conversion
        //format is same key-val pairs that snowpack expects, case sensitive
        /**
         * [Snowpack]
         * [SnowpackAdvanced]
         */
        for(auto itr : cfg)
        {
            for(auto jtr : itr.second)
            {
                d.config.addKey(jtr.first.data(),itr.first.data(),jtr.second.data());
            }
        }
        

        d.Spackconfig = std::make_shared<SnowpackConfig>(d.config);

        d.cum_precip=0.;

        //addSpecial keys goes here to deal with Antarctica, canopy, and detect grass

        SN_SNOWSOIL_DATA SSdata;
        SSdata.SoilAlb = cfg.get<double>("sno.SoilAlbedo",0.09);
        SSdata.Albedo = SSdata.SoilAlb; // following snowpacks' no snow default.
        SSdata.BareSoil_z0 = cfg.get<double>("sno.BareSoil_z0",0.2);
        if (SSdata.BareSoil_z0 == 0.)
        {
            SPDLOG_WARN("[snowpack] BareSoil_z0 == 0, set to 0.2");
            SSdata.BareSoil_z0 = 0.2;
        }

        SSdata.WindScalingFactor= cfg.get<double>("sno.WindScalingFactor",1);
        SSdata.TimeCountDeltaHS = cfg.get<double>("sno.TimeCountDeltaHS",0.0);


        SSdata.meta.stationName = cfg.get<std::string>("sno.station_name","chm");
        SSdata.meta.position.setAltitude(face->get_z());

        SSdata.meta.position.setXY(face->get_x(),face->get_y(),face->get_z());
        SSdata.meta.setSlope(mio::IOUtils::nodata,mio::IOUtils::nodata);
//        SSdata.meta.setSlope(face->slope() * ,face->aspect());
//        SSdata.meta.setSlope(0,0);

        SSdata.HS_last = 0.; //cfg.get<double>("sno.HS_Last");

        //meta data in *sno files that we don't use
//        cfg.get<std::string>("sno.station_id");

//        cfg.get<double>("sno.latitude");
//        cfg.get<double>("sno.longitude");
//        cfg.get<double>("sno.altitude");
//        cfg.get<double>("sno.nodata");
//        cfg.get<double>("sno.tz");
//        cfg.get<std::string>("sno.source");
//        cfg.get<std::string>("sno.ProfileDate");



        //assumes no starting layers
        SSdata.nN = 1;
        SSdata.Height = 0.;

        SSdata.nLayers = 0;// cfg.get("sno.nSoilLayerData",0);
//        SSdata.nLayers += cfg.get("sno.nSnowLayerData",0);
//        SSdata.Ldata



        SSdata.Canopy_Height = cfg.get<double>("sno.CanopyHeight",0);
        SSdata.Canopy_LAI = cfg.get<double>("sno.CanopyLeafAreaIndex",0);
        SSdata.Canopy_Direct_Throughfall = cfg.get<double>("sno.CanopyDirectThroughfall",1);

        SSdata.ErosionLevel = cfg.get<double>("sno.ErosionLevel",0);

        d.Xdata = std::make_shared<SnowStation>(false,false);
        d.Xdata->initialize(SSdata,0);
//        d.Xdata->cos_sl = 1;
//        d.Xdata->windward = false;
//        d.Xdata->rho_hn = 0;
//        d.Xdata->hn = 0;
//        d.Xdata->mH = 0;

        d.sp = std::make_shared<Snowpack>(*(d.Spackconfig));
        d.meteo = std::make_shared<Meteo>( (d.config));
        d.stability = std::make_shared<Stability> ( (d.config), false);

        d.sum_subl = 0;


    }
}

// Maximum number of snow layers to support in checkpoint
// Based on SmetIO.cc writesnowcover/readsnowcover variables
static const size_t MAX_LAYERS = 200;

void Lehning_snowpack::checkpoint(mesh& domain, netcdf& chkpt)
{
    auto& nc = chkpt.get_ncfile();
    size_t nFaces = domain->size_local_faces();

    // First pass: find the actual maximum number of layers across all triangles
    size_t maxLayers = 0;
    for (size_t i = 0; i < nFaces; i++)
    {
        auto face = domain->face(i);
        auto& d = face->get_module_data<Lehning_snowpack::data>(ID);
        auto& Xdata = *(d.Xdata);
        
        size_t nElems = Xdata.getNumberOfElements();
        size_t nNodes = Xdata.getNumberOfNodes();
        maxLayers = std::max(maxLayers, nElems);
        maxLayers = std::max(maxLayers, nNodes);
    }
    
    maxLayers = std::max(size_t(1), std::min(maxLayers, MAX_LAYERS));

    // Create dimensions
    netCDF::NcDim triDim;
    netCDF::NcDim layerDim;
    try {
        triDim = nc.addDim("snowpack_tri_id", nFaces);
    } catch(netCDF::exceptions::NcNameInUse&) {
        triDim = nc.getDim("snowpack_tri_id");
    }
    try {
        layerDim = nc.addDim("snowpack_layer", maxLayers);
    } catch(netCDF::exceptions::NcNameInUse&) {
        layerDim = nc.getDim("snowpack_layer");
    }
    
    std::vector<netCDF::NcDim> dims2D = {triDim, layerDim};
    std::vector<netCDF::NcDim> dims1D = {triDim};
    
    // Helper lambda to create 1D variable
    auto create1D = [&](const std::string& name) {
        try {
            nc.addVar(name, netCDF::ncDouble, dims1D);
        } catch(netCDF::exceptions::NcNameInUse&) {}
    };
    
    // Helper lambda to create 2D variable
    auto create2D = [&](const std::string& name) {
        try {
            nc.addVar(name, netCDF::ncDouble, dims2D);
        } catch(netCDF::exceptions::NcNameInUse&) {}
    };
    
    // Create scalar variables (1D - one per triangle) - SmetIO style
    create1D("snowpack:nElems");
    create1D("snowpack:nNodes");
    create1D("snowpack:Albedo");
    create1D("snowpack:pAlbedo");
    create1D("snowpack:swe");
    create1D("snowpack:cH");
    create1D("snowpack:mH");
    create1D("snowpack:Ground");
    create1D("snowpack:mass_sum");
    create1D("snowpack:lwc_sum");
    create1D("snowpack:ColdContent");
    create1D("snowpack:dIntEnergy");
    create1D("snowpack:ErosionMass");
    create1D("snowpack:cum_precip");
    create1D("snowpack:sum_subl");
    create1D("snowpack:hn");
    create1D("snowpack:rho_hn");
    
    // Create element data variables (2D - per triangle per layer) - Exact SmetIO fields
    create2D("snowpack:E_L");                 // Layer_Thick [m]
    create2D("snowpack:E_Te");                // Layer temperature [K]
    create2D("snowpack:E_theta_ICE");         // Vol_Frac_I - ice content [0-1]
    create2D("snowpack:E_theta_i_reservoir");  // Vol_Frac_IR - ice reservoir [0-1]
    create2D("snowpack:E_theta_i_reservoir_cumul"); // Vol_Frac_CIR - cumulative ice reservoir [0-1]
    create2D("snowpack:E_theta_WATER");       // Vol_Frac_W - liquid water [0-1]
    create2D("snowpack:E_theta_WATER_PREF");  // Vol_Frac_WP - preferential flow water [0-1]
    create2D("snowpack:E_theta_AIR");         // Vol_Frac_V - voids/air [0-1]
    create2D("snowpack:E_theta_SOIL");        // Vol_Frac_S - soil content [0-1]
    create2D("snowpack:E_soil_rho");          // Rho_S - soil density [kg/m3]
    create2D("snowpack:E_soil_k");            // Conduc_S - soil conductivity [W/(mK)]
    create2D("snowpack:E_soil_c");            // HeatCapac_S - soil heat capacity [J/(kgK)]
    create2D("snowpack:E_rg");                // rg - grain radius [mm]
    create2D("snowpack:E_rb");                // rb - bond radius [mm]
    create2D("snowpack:E_dd");                // dd - dendricity [0-1]
    create2D("snowpack:E_sp");                // sp - sphericity [0-1]
    create2D("snowpack:E_mk");                // mk - grain marker
    create2D("snowpack:E_CDot");              // CDot - stress rate [Pa/s]
    create2D("snowpack:E_metamo");            // metamo - metamorphism state
    create2D("snowpack:E_dsm");               // dsm - dry snow metamorphism (NIED)
    create2D("snowpack:E_salinity");          // Sal - salinity [PSU]
    create2D("snowpack:E_h");                 // h - capillary pressure head [m]
    create2D("snowpack:E_depositionDate");    // depositionDate - layer deposition date as Julian date
    
    // Create node data variables (2D - per triangle per node) - Exact SmetIO fields
    create2D("snowpack:N_T");                 // Temperature [K]
    create2D("snowpack:N_hoar");              // mass_hoar - surface hoar mass
    
    double fillValue = mio::IOUtils::nodata;
    
    // Allocate buffers using the actual maxLayers, not MAX_LAYERS
    std::vector<double> buf_nElems(nFaces);
    std::vector<double> buf_nNodes(nFaces);
    std::vector<double> buf_Albedo(nFaces);
    std::vector<double> buf_pAlbedo(nFaces);
    std::vector<double> buf_swe(nFaces);
    std::vector<double> buf_cH(nFaces);
    std::vector<double> buf_mH(nFaces);
    std::vector<double> buf_Ground(nFaces);
    std::vector<double> buf_mass_sum(nFaces);
    std::vector<double> buf_lwc_sum(nFaces);
    std::vector<double> buf_ColdContent(nFaces);
    std::vector<double> buf_dIntEnergy(nFaces);
    std::vector<double> buf_ErosionMass(nFaces);
    std::vector<double> buf_cum_precip(nFaces);
    std::vector<double> buf_sum_subl(nFaces);
    std::vector<double> buf_hn(nFaces);
    std::vector<double> buf_rho_hn(nFaces);
    
    size_t total2D = nFaces * maxLayers;
    std::vector<double> buf_E_L(total2D, fillValue);
    std::vector<double> buf_E_Te(total2D, fillValue);
    std::vector<double> buf_E_theta_ICE(total2D, fillValue);
    std::vector<double> buf_E_theta_i_reservoir(total2D, fillValue);
    std::vector<double> buf_E_theta_i_reservoir_cumul(total2D, fillValue);
    std::vector<double> buf_E_theta_WATER(total2D, fillValue);
    std::vector<double> buf_E_theta_WATER_PREF(total2D, fillValue);
    std::vector<double> buf_E_theta_AIR(total2D, fillValue);
    std::vector<double> buf_E_theta_SOIL(total2D, fillValue);
    std::vector<double> buf_E_soil_rho(total2D, fillValue);
    std::vector<double> buf_E_soil_k(total2D, fillValue);
    std::vector<double> buf_E_soil_c(total2D, fillValue);
    std::vector<double> buf_E_rg(total2D, fillValue);
    std::vector<double> buf_E_rb(total2D, fillValue);
    std::vector<double> buf_E_dd(total2D, fillValue);
    std::vector<double> buf_E_sp(total2D, fillValue);
    std::vector<double> buf_E_mk(total2D, fillValue);
    std::vector<double> buf_E_CDot(total2D, fillValue);
    std::vector<double> buf_E_metamo(total2D, fillValue);
    std::vector<double> buf_E_dsm(total2D, fillValue);
    std::vector<double> buf_E_salinity(total2D, fillValue);
    std::vector<double> buf_E_h(total2D, fillValue);
    std::vector<double> buf_E_depositionDate(total2D, fillValue);
    
    std::vector<double> buf_N_T(total2D, fillValue);
    std::vector<double> buf_N_hoar(total2D, fillValue);
    
    // Fill buffers from mesh data
    for (size_t i = 0; i < nFaces; i++)
    {
        auto face = domain->face(i);
        auto& d = face->get_module_data<Lehning_snowpack::data>(ID);
        auto& Xdata = *(d.Xdata);
        
        size_t nElems = Xdata.getNumberOfElements();
        size_t nNodes = Xdata.getNumberOfNodes();
        
        buf_nElems[i] = static_cast<double>(nElems);
        buf_nNodes[i] = static_cast<double>(nNodes);
        buf_Albedo[i] = Xdata.Albedo;
        buf_pAlbedo[i] = Xdata.pAlbedo;
        buf_swe[i] = Xdata.swe;
        buf_cH[i] = Xdata.cH;
        buf_mH[i] = Xdata.mH;
        buf_Ground[i] = Xdata.Ground;
        buf_mass_sum[i] = Xdata.mass_sum;
        buf_lwc_sum[i] = Xdata.lwc_sum;
        buf_ColdContent[i] = Xdata.ColdContent;
        buf_dIntEnergy[i] = Xdata.dIntEnergy;
        buf_ErosionMass[i] = Xdata.ErosionMass;
        buf_cum_precip[i] = d.cum_precip;
        buf_sum_subl[i] = d.sum_subl;
        buf_hn[i] = Xdata.hn;
        buf_rho_hn[i] = Xdata.rho_hn;
        
        // Fill element data - Exact SmetIO fields
        for (size_t e = 0; e < nElems && e < maxLayers; e++)
        {
            size_t idx = i * maxLayers + e;
            auto& elem = Xdata.Edata[e];
            buf_E_L[idx] = elem.L;
            buf_E_Te[idx] = Xdata.Ndata[e+1].T; // Note: SmetIO uses node T for layer temperature
            buf_E_theta_ICE[idx] = elem.theta[ICE];
            buf_E_theta_i_reservoir[idx] = elem.theta_i_reservoir;
            buf_E_theta_i_reservoir_cumul[idx] = elem.theta_i_reservoir_cumul;
            buf_E_theta_WATER[idx] = elem.theta[WATER];
            buf_E_theta_WATER_PREF[idx] = elem.theta[WATER_PREF];
            buf_E_theta_AIR[idx] = elem.theta[AIR];
            buf_E_theta_SOIL[idx] = elem.theta[SOIL];
            buf_E_soil_rho[idx] = elem.soil[SOIL_RHO];
            buf_E_soil_k[idx] = elem.soil[SOIL_K];
            buf_E_soil_c[idx] = elem.soil[SOIL_C];
            buf_E_rg[idx] = elem.rg;
            buf_E_rb[idx] = elem.rb;
            buf_E_dd[idx] = elem.dd;
            buf_E_sp[idx] = elem.sp;
            buf_E_mk[idx] = static_cast<double>(elem.mk);
            buf_E_CDot[idx] = elem.CDot;
            buf_E_metamo[idx] = elem.metamo;
            buf_E_dsm[idx] = elem.dsm;
            buf_E_salinity[idx] = elem.salinity;
            buf_E_h[idx] = elem.h;
            // Store depositionDate as Julian date; use nodata if undefined
            buf_E_depositionDate[idx] = elem.depositionDate.isUndef() ? mio::IOUtils::nodata : elem.depositionDate.getJulian();
        }
        
        // Fill node data - Exact SmetIO fields
        for (size_t n = 0; n < nNodes && n < maxLayers; n++)
        {
            size_t idx = i * maxLayers + n;
            auto& node = Xdata.Ndata[n];
            buf_N_T[idx] = node.T;
            buf_N_hoar[idx] = node.hoar;
        }
    }
    
    // Batch write all 1D variables
    nc.getVar("snowpack:nElems").putVar(buf_nElems.data());
    nc.getVar("snowpack:nNodes").putVar(buf_nNodes.data());
    nc.getVar("snowpack:Albedo").putVar(buf_Albedo.data());
    nc.getVar("snowpack:pAlbedo").putVar(buf_pAlbedo.data());
    nc.getVar("snowpack:swe").putVar(buf_swe.data());
    nc.getVar("snowpack:cH").putVar(buf_cH.data());
    nc.getVar("snowpack:mH").putVar(buf_mH.data());
    nc.getVar("snowpack:Ground").putVar(buf_Ground.data());
    nc.getVar("snowpack:mass_sum").putVar(buf_mass_sum.data());
    nc.getVar("snowpack:lwc_sum").putVar(buf_lwc_sum.data());
    nc.getVar("snowpack:ColdContent").putVar(buf_ColdContent.data());
    nc.getVar("snowpack:dIntEnergy").putVar(buf_dIntEnergy.data());
    nc.getVar("snowpack:ErosionMass").putVar(buf_ErosionMass.data());
    nc.getVar("snowpack:cum_precip").putVar(buf_cum_precip.data());
    nc.getVar("snowpack:sum_subl").putVar(buf_sum_subl.data());
    nc.getVar("snowpack:hn").putVar(buf_hn.data());
    nc.getVar("snowpack:rho_hn").putVar(buf_rho_hn.data());
    
    // Batch write all 2D variables - Exact SmetIO fields
    nc.getVar("snowpack:E_L").putVar(buf_E_L.data());
    nc.getVar("snowpack:E_Te").putVar(buf_E_Te.data());
    nc.getVar("snowpack:E_theta_ICE").putVar(buf_E_theta_ICE.data());
    nc.getVar("snowpack:E_theta_i_reservoir").putVar(buf_E_theta_i_reservoir.data());
    nc.getVar("snowpack:E_theta_i_reservoir_cumul").putVar(buf_E_theta_i_reservoir_cumul.data());
    nc.getVar("snowpack:E_theta_WATER").putVar(buf_E_theta_WATER.data());
    nc.getVar("snowpack:E_theta_WATER_PREF").putVar(buf_E_theta_WATER_PREF.data());
    nc.getVar("snowpack:E_theta_AIR").putVar(buf_E_theta_AIR.data());
    nc.getVar("snowpack:E_theta_SOIL").putVar(buf_E_theta_SOIL.data());
    nc.getVar("snowpack:E_soil_rho").putVar(buf_E_soil_rho.data());
    nc.getVar("snowpack:E_soil_k").putVar(buf_E_soil_k.data());
    nc.getVar("snowpack:E_soil_c").putVar(buf_E_soil_c.data());
    nc.getVar("snowpack:E_rg").putVar(buf_E_rg.data());
    nc.getVar("snowpack:E_rb").putVar(buf_E_rb.data());
    nc.getVar("snowpack:E_dd").putVar(buf_E_dd.data());
    nc.getVar("snowpack:E_sp").putVar(buf_E_sp.data());
    nc.getVar("snowpack:E_mk").putVar(buf_E_mk.data());
    nc.getVar("snowpack:E_CDot").putVar(buf_E_CDot.data());
    nc.getVar("snowpack:E_metamo").putVar(buf_E_metamo.data());
    nc.getVar("snowpack:E_dsm").putVar(buf_E_dsm.data());
    nc.getVar("snowpack:E_salinity").putVar(buf_E_salinity.data());
    nc.getVar("snowpack:E_h").putVar(buf_E_h.data());
    nc.getVar("snowpack:E_depositionDate").putVar(buf_E_depositionDate.data());
    
    nc.getVar("snowpack:N_T").putVar(buf_N_T.data());
    nc.getVar("snowpack:N_hoar").putVar(buf_N_hoar.data());
}

void Lehning_snowpack::load_checkpoint(mesh& domain, netcdf& chkpt)
{
    auto& nc = chkpt.get_ncfile();
    size_t nFaces = domain->size_local_faces();
    
    // Get the actual layer dimension from the checkpoint file
    size_t maxLayers = nc.getDim("snowpack_layer").getSize();
    size_t total2D = nFaces * maxLayers;
    
    // Allocate buffers for batch reading
    std::vector<double> buf_nElems(nFaces);
    std::vector<double> buf_nNodes(nFaces);
    std::vector<double> buf_Albedo(nFaces);
    std::vector<double> buf_pAlbedo(nFaces);
    std::vector<double> buf_swe(nFaces);
    std::vector<double> buf_cH(nFaces);
    std::vector<double> buf_mH(nFaces);
    std::vector<double> buf_Ground(nFaces);
    std::vector<double> buf_mass_sum(nFaces);
    std::vector<double> buf_lwc_sum(nFaces);
    std::vector<double> buf_ColdContent(nFaces);
    std::vector<double> buf_dIntEnergy(nFaces);
    std::vector<double> buf_ErosionMass(nFaces);
    std::vector<double> buf_cum_precip(nFaces);
    std::vector<double> buf_sum_subl(nFaces);
    std::vector<double> buf_hn(nFaces);
    std::vector<double> buf_rho_hn(nFaces);
    
    std::vector<double> buf_E_L(total2D);
    std::vector<double> buf_E_Te(total2D);
    std::vector<double> buf_E_theta_ICE(total2D);
    std::vector<double> buf_E_theta_i_reservoir(total2D, 0.0);
    std::vector<double> buf_E_theta_i_reservoir_cumul(total2D, 0.0);
    std::vector<double> buf_E_theta_WATER(total2D);
    std::vector<double> buf_E_theta_WATER_PREF(total2D, 0.0);
    std::vector<double> buf_E_theta_AIR(total2D);
    std::vector<double> buf_E_theta_SOIL(total2D);
    std::vector<double> buf_E_soil_rho(total2D);
    std::vector<double> buf_E_soil_k(total2D);
    std::vector<double> buf_E_soil_c(total2D);
    std::vector<double> buf_E_rg(total2D);
    std::vector<double> buf_E_rb(total2D);
    std::vector<double> buf_E_dd(total2D);
    std::vector<double> buf_E_sp(total2D);
    std::vector<double> buf_E_mk(total2D);
    std::vector<double> buf_E_CDot(total2D, 0.0);
    std::vector<double> buf_E_metamo(total2D, 0.0);
    std::vector<double> buf_E_dsm(total2D, 0.0);
    std::vector<double> buf_E_salinity(total2D, 0.0);
    std::vector<double> buf_E_h(total2D, 0.0);
    std::vector<double> buf_E_depositionDate(total2D, mio::IOUtils::nodata);
    
    std::vector<double> buf_N_T(total2D);
    std::vector<double> buf_N_hoar(total2D, 0.0);
    
    // Batch read all 1D variables
    nc.getVar("snowpack:nElems").getVar(buf_nElems.data());
    nc.getVar("snowpack:nNodes").getVar(buf_nNodes.data());
    nc.getVar("snowpack:Albedo").getVar(buf_Albedo.data());
    nc.getVar("snowpack:pAlbedo").getVar(buf_pAlbedo.data());
    nc.getVar("snowpack:swe").getVar(buf_swe.data());
    nc.getVar("snowpack:cH").getVar(buf_cH.data());
    nc.getVar("snowpack:mH").getVar(buf_mH.data());
    nc.getVar("snowpack:Ground").getVar(buf_Ground.data());
    nc.getVar("snowpack:mass_sum").getVar(buf_mass_sum.data());
    nc.getVar("snowpack:lwc_sum").getVar(buf_lwc_sum.data());
    nc.getVar("snowpack:ColdContent").getVar(buf_ColdContent.data());
    nc.getVar("snowpack:dIntEnergy").getVar(buf_dIntEnergy.data());
    nc.getVar("snowpack:ErosionMass").getVar(buf_ErosionMass.data());
    nc.getVar("snowpack:cum_precip").getVar(buf_cum_precip.data());
    nc.getVar("snowpack:sum_subl").getVar(buf_sum_subl.data());
    nc.getVar("snowpack:hn").getVar(buf_hn.data());
    nc.getVar("snowpack:rho_hn").getVar(buf_rho_hn.data());
    
    // Batch read all 2D variables - Exact SmetIO fields
    nc.getVar("snowpack:E_L").getVar(buf_E_L.data());
    nc.getVar("snowpack:E_Te").getVar(buf_E_Te.data());
    nc.getVar("snowpack:E_theta_ICE").getVar(buf_E_theta_ICE.data());
    try { nc.getVar("snowpack:E_theta_i_reservoir").getVar(buf_E_theta_i_reservoir.data()); } catch(...) {}
    try { nc.getVar("snowpack:E_theta_i_reservoir_cumul").getVar(buf_E_theta_i_reservoir_cumul.data()); } catch(...) {}
    nc.getVar("snowpack:E_theta_WATER").getVar(buf_E_theta_WATER.data());
    try { nc.getVar("snowpack:E_theta_WATER_PREF").getVar(buf_E_theta_WATER_PREF.data()); } catch(...) {}
    nc.getVar("snowpack:E_theta_AIR").getVar(buf_E_theta_AIR.data());
    nc.getVar("snowpack:E_theta_SOIL").getVar(buf_E_theta_SOIL.data());
    nc.getVar("snowpack:E_soil_rho").getVar(buf_E_soil_rho.data());
    nc.getVar("snowpack:E_soil_k").getVar(buf_E_soil_k.data());
    nc.getVar("snowpack:E_soil_c").getVar(buf_E_soil_c.data());
    nc.getVar("snowpack:E_rg").getVar(buf_E_rg.data());
    nc.getVar("snowpack:E_rb").getVar(buf_E_rb.data());
    nc.getVar("snowpack:E_dd").getVar(buf_E_dd.data());
    nc.getVar("snowpack:E_sp").getVar(buf_E_sp.data());
    nc.getVar("snowpack:E_mk").getVar(buf_E_mk.data());
    try { nc.getVar("snowpack:E_CDot").getVar(buf_E_CDot.data()); } catch(...) {}
    try { nc.getVar("snowpack:E_metamo").getVar(buf_E_metamo.data()); } catch(...) {}
    try { nc.getVar("snowpack:E_dsm").getVar(buf_E_dsm.data()); } catch(...) {}
    try { nc.getVar("snowpack:E_salinity").getVar(buf_E_salinity.data()); } catch(...) {}
    try { nc.getVar("snowpack:E_h").getVar(buf_E_h.data()); } catch(...) {}
    try { nc.getVar("snowpack:E_depositionDate").getVar(buf_E_depositionDate.data()); } catch(...) {}
    
    nc.getVar("snowpack:N_T").getVar(buf_N_T.data());
    try { nc.getVar("snowpack:N_hoar").getVar(buf_N_hoar.data()); } catch(...) {}
    
    // Distribute data to mesh elements
    for (size_t i = 0; i < nFaces; i++)
    {
        auto face = domain->face(i);
        auto& d = face->get_module_data<Lehning_snowpack::data>(ID);
        auto& Xdata = *(d.Xdata);
        
        // Read scalar values from buffers
        size_t nElems = static_cast<size_t>(buf_nElems[i]);
        size_t nNodes = static_cast<size_t>(buf_nNodes[i]);
        
        // First, re-initialize the SnowStation with basic configuration (same as init())
        // This ensures all internal structures are properly set up
        SN_SNOWSOIL_DATA SSdata;
        SSdata.SoilAlb = cfg.get<double>("sno.SoilAlbedo", 0.09);
        SSdata.Albedo = SSdata.SoilAlb;
        SSdata.BareSoil_z0 = cfg.get<double>("sno.BareSoil_z0", 0.2);
        if (SSdata.BareSoil_z0 == 0.) {
            SSdata.BareSoil_z0 = 0.2;
        }
        SSdata.WindScalingFactor = cfg.get<double>("sno.WindScalingFactor", 1);
        SSdata.TimeCountDeltaHS = cfg.get<double>("sno.TimeCountDeltaHS", 0.0);
        SSdata.meta.stationName = cfg.get<std::string>("sno.station_name", "chm");
        SSdata.meta.position.setAltitude(face->get_z());
        SSdata.meta.position.setXY(face->get_x(), face->get_y(), face->get_z());
        SSdata.meta.setSlope(mio::IOUtils::nodata, mio::IOUtils::nodata);
        SSdata.HS_last = 0.;
        SSdata.nN = 1;
        SSdata.Height = 0.;
        SSdata.nLayers = 0;
        SSdata.Canopy_Height = cfg.get<double>("sno.CanopyHeight", 0);
        SSdata.Canopy_LAI = cfg.get<double>("sno.CanopyLeafAreaIndex", 0);
        SSdata.Canopy_Direct_Throughfall = cfg.get<double>("sno.CanopyDirectThroughfall", 1);
        SSdata.ErosionLevel = cfg.get<double>("sno.ErosionLevel", 0);
        
        Xdata.initialize(SSdata, 0);
        
        // Now overwrite with checkpoint data
        Xdata.Albedo = buf_Albedo[i];
        Xdata.pAlbedo = buf_pAlbedo[i];
        Xdata.swe = buf_swe[i];
        Xdata.cH = buf_cH[i];
        Xdata.mH = buf_mH[i];
        Xdata.Ground = buf_Ground[i];
        Xdata.mass_sum = buf_mass_sum[i];
        Xdata.lwc_sum = buf_lwc_sum[i];
        Xdata.ColdContent = buf_ColdContent[i];
        Xdata.dIntEnergy = buf_dIntEnergy[i];
        Xdata.ErosionMass = buf_ErosionMass[i];
        d.cum_precip = buf_cum_precip[i];
        d.sum_subl = buf_sum_subl[i];
        Xdata.hn = buf_hn[i];
        Xdata.rho_hn = buf_rho_hn[i];
        
        // Resize the vectors to accommodate the loaded data
        Xdata.resize(nElems);
        
        // Read element data from buffers - Exact SmetIO style
        // Use nElems (per-triangle actual count), not maxLayers, to avoid reading nodata
        for (size_t e = 0; e < nElems; e++)
        {
            size_t idx = i * maxLayers + e;
            auto& elem = Xdata.Edata[e];
            
            // Initialize all fields to reasonable defaults - matching SmetIO
            elem.L = buf_E_L[idx];
            elem.L0 = elem.L;
            elem.Te = buf_E_Te[idx];
            elem.gradT = 0.0;
            elem.meltfreeze_tk = mio::Cst::t_water_freezing_pt;
            elem.theta[ICE] = buf_E_theta_ICE[idx];
            elem.theta_i_reservoir = buf_E_theta_i_reservoir[idx];
            elem.theta_i_reservoir_cumul = buf_E_theta_i_reservoir_cumul[idx];
            elem.theta[WATER] = buf_E_theta_WATER[idx];
            elem.theta[WATER_PREF] = buf_E_theta_WATER_PREF[idx];
            elem.theta[AIR] = buf_E_theta_AIR[idx];
            elem.theta[SOIL] = buf_E_theta_SOIL[idx];
            elem.h = buf_E_h[idx];
            elem.soil[SOIL_RHO] = buf_E_soil_rho[idx];
            elem.soil[SOIL_K] = buf_E_soil_k[idx];
            elem.soil[SOIL_C] = buf_E_soil_c[idx];
            elem.Rho = (elem.theta[ICE] * 917.0) + (elem.theta[WATER] * 1000.0) + (elem.theta[SOIL] * elem.soil[SOIL_RHO]);
            elem.M = elem.Rho * elem.L;
            elem.sw_abs = 0.0;
            elem.rg = buf_E_rg[idx];
            elem.rb = buf_E_rb[idx];
            elem.dd = buf_E_dd[idx];
            elem.sp = buf_E_sp[idx];
            elem.ogs = elem.rg; // Default to grain radius
            elem.N3 = 4.0; // Typical coordination number
            elem.mk = static_cast<unsigned short>(buf_E_mk[idx] + 0.5);
            elem.type = 0;
            elem.metamo = buf_E_metamo[idx];
            elem.salinity = buf_E_salinity[idx];
            // Restore depositionDate from Julian date; keep as undefined (default) if nodata
            if (buf_E_depositionDate[idx] != mio::IOUtils::nodata) {
                elem.depositionDate.setDate(buf_E_depositionDate[idx], 0.0);
            }
            elem.dth_w = 0.0;
            // res_wat_cont will be computed by snowResidualWaterContent() below
            elem.Qmf = 0.0;
            elem.QIntmf = 0.0;
            elem.dEps = 0.0;
            elem.Eps = 0.0;
            elem.Eps_e = 0.0;
            elem.Eps_v = 0.0;
            elem.Eps_Dot = 0.0;
            elem.Eps_vDot = 0.0;
            elem.E = 0.0;
            elem.S = 0.0;
            elem.C = 0.0;
            elem.CDot = buf_E_CDot[idx];
            elem.ps2rb = 0.0;
            elem.s_strength = 0.0;
            elem.hard = 0.0;
            elem.S_dr = 0.0;
            elem.crit_cut_length = 0.0;
            elem.lwc_source = 0.0;
            elem.PrefFlowArea = 0.0;
            elem.theta_w_transfer = 0.0;
            elem.SlopeParFlux = 0.0;
            elem.Qph_up = 0.0;
            elem.Qph_down = 0.0;
            elem.dsm = buf_E_dsm[idx];
            elem.rime = 0.0;
            elem.rhov = 0.0;
            elem.Qmm = 0.0;
            elem.vapTrans_fluxDiff = 0.0;
            elem.vapTrans_snowDenChangeRate = 0.0;
            elem.vapTrans_cumulativeDenChange = 0.0;
            elem.vapTrans_underSaturationDegree = 0.0;
        }
        
        // These are normally set in SnowStation::initialize() but need to be explicitly
        // computed here after loading checkpoint data
        Xdata.SoilNode = 0;  // Will be computed below
        for (size_t e = 0; e < nElems; e++) {
            auto& elem = Xdata.Edata[e];
            
            // Compute residual water content based on ice content (required for boundary conditions)
            elem.snowResidualWaterContent();
            
            // Compute heat capacity (required for thermal matrix)
            elem.heatCapacity();
            
            // Update density from volumetric contents (ensures consistency)
            elem.updDensity();
            
            // Compute mass from density and length
            elem.M = elem.Rho * elem.L;
            
            // Count soil nodes for proper boundary condition handling
            if (elem.theta[SOIL] > 0.0) {
                Xdata.SoilNode++;
            }
        }
        
        // Read node data from buffers - Match SmetIO
        for (size_t n = 0; n < nNodes; n++)
        {
            size_t idx = i * maxLayers + n;
            auto& node = Xdata.Ndata[n];
            node.T = buf_N_T[idx];
            node.hoar = buf_N_hoar[idx];
            // Ensure all other fields are initialized (NodeData constructor already does this, but just in case)
            node.z = 0.0;
            node.u = 0.0;
            node.f = 0.0;
            node.udot = 0.0;
            node.S_n = 0.0;
            node.S_s = 0.0;
            node.ssi = 6.0; // Max stability
            node.dsm = 0.0;
            node.S_dsm = 0.0;
            node.Sigdsm = 0.0;
            node.rime = 0.0;
            node.water_flux = 0.0;
            node.rhov = 0.0;
        }
        
        // Ensure ground node has a valid temperature (not nodata)
        // When there's no snow, the ground node temperature is critical for stability
        if (Xdata.Ndata[0].T <= 0.0 || Xdata.Ndata[0].T > 400.0) {
            SPDLOG_DEBUG("Face {}: Invalid ground node temperature {}, resetting to freezing point",
                i, Xdata.Ndata[0].T);
            Xdata.Ndata[0].T = mio::Cst::t_water_freezing_pt;
        }
        
        // Ensure swe is consistent with nElems
        // If there are no elements, swe should be 0
        if (nElems == 0 && Xdata.swe > 0.0) {
            SPDLOG_DEBUG("Resetting swe from {} to 0 for face {} (no elements)", Xdata.swe, i);
            Xdata.swe = 0.0;
        }
        
        // Recompute node positions (z) from element thicknesses (L)
        // This is necessary for proper thermal calculations
        // Ndata[0] is at the bottom (ground), Ndata[nNodes-1] is at the surface
        Xdata.Ndata[0].z = 0.0;  // Ground level
        double computed_cH = Xdata.Ground;  // Start from ground level
        for (size_t e = 0; e < nElems; e++) {
            // Node e+1 is above element e
            Xdata.Ndata[e + 1].z = Xdata.Ndata[e].z + Xdata.Edata[e].L;
            computed_cH += Xdata.Edata[e].L;
        }
        // Ensure cH is consistent with the sum of element thicknesses plus ground
        // Use computed value if there's a mismatch (indicates checkpoint data inconsistency)
        if (nElems > 0 && std::abs(Xdata.cH - computed_cH) > 0.001) {
            SPDLOG_DEBUG("Adjusting cH from {} to {} for face {} (checkpoint inconsistency)",
                        Xdata.cH, computed_cH, i);
            Xdata.cH = computed_cH;
        }
        
        // Synchronize element temperatures with node temperatures for consistency
        // This ensures Te is the average of adjacent node temperatures
        for (size_t e = 0; e < nElems; e++) {
            Xdata.Edata[e].Te = (Xdata.Ndata[e].T + Xdata.Ndata[e + 1].T) / 2.;
        }
        
        // Set output variables from loaded state
        (*face)["swe"_s] = Xdata.swe;
        (*face)["snowdepthavg"_s] = Xdata.cH - Xdata.Ground;
        
        if (Xdata.swe > 0)
        {
            double bulk_T_s = 0;
            for (size_t e = 0; e < nElems; ++e)
            {
                bulk_T_s += Xdata.Edata[e].Te;
            }
            if (nElems > 0)
                bulk_T_s /= nElems;
            
            (*face)["T_s"_s] = bulk_T_s - mio::Cst::t_water_freezing_pt;
            (*face)["n_nodes"_s] = nNodes;
            (*face)["n_elem"_s] = nElems;
        }
        else
        {
            (*face)["T_s"_s] = mio::IOUtils::nodata;
            (*face)["n_nodes"_s] = 0;
            (*face)["n_elem"_s] = 0;
        }
        
        (*face)["mass_snowpack_removed"_s] = Xdata.ErosionMass;
        (*face)["snow_albedo"_s] = Xdata.Albedo;
        (*face)["sum_subl"_s] = d.sum_subl;
    }
    
    SPDLOG_INFO("Snowpack checkpoint loaded successfully for {} faces", nFaces);
}
