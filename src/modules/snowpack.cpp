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
    data.meteo->compMeteo(Mdata,*(data.Xdata),false); // no canopy model

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
            BOOST_THROW_EXCEPTION(module_error() << errstr_info("Snowpack died. Triangle center = " + details));
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

    for(size_t i=0;i<domain->size_faces();i++)
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

        d.config.addKey("COMBINE_ELEMENTS","SnowpackAdvanced","true"); //Defines whether joining elements will be considered at all
        //Activates algorithm to reduce the number of elements deeper in the snowpack AND to split elements again when they come back to the surface
        //Only works when COMBINE_ELEMENTS == TRUE.
        d.config.addKey("REDUCE_N_ELEMENTS","SnowpackAdvanced","true");


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
        

        d.Spackconfig = boost::make_shared<SnowpackConfig>(d.config);

        d.cum_precip=0.;

        //addSpecial keys goes here to deal with Antarctica, canopy, and detect grass

        SN_SNOWSOIL_DATA SSdata;
        SSdata.SoilAlb = cfg.get<double>("sno.SoilAlbedo",0.09);
        SSdata.Albedo = SSdata.SoilAlb; // following snowpacks' no snow default.
        SSdata.BareSoil_z0 = cfg.get<double>("sno.BareSoil_z0",0.2);
        if (SSdata.BareSoil_z0 == 0.)
        {
            LOG_WARNING << "[snowpack] BareSoil_z0 == 0, set to 0.2";
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

        d.Xdata = boost::make_shared<SnowStation>(false,false);
        d.Xdata->initialize(SSdata,0);
//        d.Xdata->cos_sl = 1;
//        d.Xdata->windward = false;
//        d.Xdata->rho_hn = 0;
//        d.Xdata->hn = 0;
//        d.Xdata->mH = 0;

        d.sp = boost::make_shared<Snowpack>(*(d.Spackconfig));
        d.meteo = boost::make_shared<Meteo>( (d.config));
        d.stability = boost::make_shared<Stability> ( (d.config), false);

        d.sum_subl = 0;


    }
}
