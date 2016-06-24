#include "snowpack.hpp"

Lehning_snowpack::Lehning_snowpack(config_file cfg)
        : module_base(parallel::data)
{
    depends("p_rain");
    depends("p_snow");
    depends("iswr");
    depends("rh");
    depends("t");
    depends("vw");
    depends("vw_dir");
    depends("p");
    depends("ilwr");

    optional("snow_albedo");

    provides("dQ");
    provides("swe");
    provides("T_s");
    provides("T_s_0");
    provides("n_nodes");
    provides("n_elem");
    provides("spack_albedo");
    provides("H");
    provides("E");
    provides("G");
    provides("ilwr_out");
    provides("iswr_out");
    provides("R_n");


}

Lehning_snowpack::~Lehning_snowpack()
{

}

void Lehning_snowpack::run(mesh_elem &elem, boost::shared_ptr <global> global_param)
{

    auto data = elem->get_module_data<Lehning_snowpack::data>(ID);

    /**
     * Builds this timestep's meteo data
     */
    CurrentMeteo Mdata(data->config);
    Mdata.date   =  mio::Date( global_param->year(), global_param->month(), global_param->day(),global_param->hour(),global_param->min(),-6 );
    Mdata.ta     =  elem->face_data("t")+273.15;
    Mdata.rh     =  elem->face_data("rh")/100.;
    Mdata.vw     =  elem->face_data("vw");
    Mdata.dw     =  elem->face_data("vw_dir");
    Mdata.vw_max = mio::IOUtils::nodata;//elem->face_data("vw");// TODO: fix md(MeteoData::VW_MAX);

    Mdata.vw_drift =  Mdata.vw ;//mio::IOUtils::nodata;
    Mdata.dw_drift = Mdata.dw;//0;



    Mdata.iswr      = elem->face_data("iswr");

    // If  Snowpack, "SW_MODE" : "BOTH"  then rswr and iswr needs to be definined.
    // If an external albedo is not used, a parametrized one is used. but rswr and iswr both must be defined.
    // If "SW_MODE" : "INCOMING", is used, then mAlbedo needs to be undefined to trigger the appropriate rswr and albedo calculations.
    // rswr must still be set, but we use the garbage that is in Xdata instead.
    if(has_optional("snow_albedo"))
    {
        //measured albedo in snowpack will be fed from an albedo model
        //in the config 'both' will enable this
        Mdata.mAlbedo   =  elem->face_data("snow_albedo");
        Mdata.rswr      =  elem->face_data("snow_albedo") * elem->face_data("iswr");

    }
    else
    {
        Mdata.rswr = data->Xdata->Albedo * Mdata.iswr;
        Mdata.mAlbedo = Constants::undefined; //this will trigger calculating a paramaterized albedo
    }



    double lw_in = elem->face_data("ilwr");
//    Mdata.ea  = SnLaws::AirEmissivity(lw_in, Mdata.ta, "default"); //atmospheric emissivity!
    Mdata.ea = mio::Atmosphere::blkBody_Emissivity(lw_in, Mdata.ta);
    double thresh_rain = 2 + 273.15;

//    Mdata.psum_ph = (Mdata.ta>= thresh_rain) ? 1. : 0.;
    Mdata.psum_ph = elem->face_data("frac_precip_rain"); //  0 = snow, 1 = rain
    Mdata.psum = elem->face_data("p");



    //setup a single ground temp measurement
//    mio::MeteoData soil_meas;
//    soil_meas.addParameter("HTS1");
//    soil_meas.addParameter("TS1");
//    soil_meas("HTS1") = -0.1;
//    soil_meas("TS1") = 269;
//    Mdata.setMeasTempParameters(soil_meas);
//    Mdata.ts.push_back(soil_meas("TS1"));
//    Mdata.zv_ts.push_back(soil_meas("HTS1"));

    Mdata.tss=  data->Xdata->Ndata[data->Xdata->getNumberOfElements()].T;  //we use previous timestep value//mio::IOUtils::nodata; //Constants::undefined;;//
    Mdata.ts0 = Mdata.tss;//mio::IOUtils::nodata;//273.15-4.;


//    Mdata.hs = mio::IOUtils::nodata;

    Mdata.diff      = elem->face_data("iswr_diffuse");
    Mdata.dir_h     = elem->face_data("iswr_direct");
    Mdata.elev      = global_param->solar_el()*mio::Cst::to_rad;

//    Mdata.tss_a12h = Constants::undefined;
//    Mdata.tss_a24h = Constants::undefined;
//    Mdata.hs_a3h = Constants::undefined;
//    Mdata.hs_rate = Constants::undefined;

    data->cum_precip  += Mdata.psum; //running sum of the precip. snowpack removes the rain component for us.
    data->meteo->compMeteo(Mdata,*(data->Xdata),false); // no canopy model

    // To collect surface exchange data for output
    SurfaceFluxes surface_fluxes;
    surface_fluxes.reset(false);
    surface_fluxes.drift = 0.;
    surface_fluxes.mass[SurfaceFluxes::MS_WIND] = 0.;


    // Boundary condition (fluxes)
    BoundCond Bdata;

    try
    {
        data->sp->runSnowpackModel(Mdata, *(data->Xdata), data->cum_precip, Bdata,surface_fluxes);
        surface_fluxes.collectSurfaceFluxes(Bdata, *(data->Xdata), Mdata);
    }catch(std::exception& e)
    {
        LOG_DEBUG << e.what();
        auto details = "("+std::to_string(elem->center().x()) + "," + std::to_string(elem->center().y())+","+std::to_string(elem->center().z())+")";
        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Snowpack died. Triangle center = "+details));
    }



    elem->set_face_data("swe",data->Xdata->swe);
    elem->set_face_data("T_s",Mdata.tss-273.15);
    elem->set_face_data("T_s_0",Mdata.tss-273.15);
    elem->set_face_data("n_nodes",data->Xdata->getNumberOfNodes());
    elem->set_face_data("n_elem",data->Xdata->getNumberOfElements());

    elem->set_face_data("H",surface_fluxes.qs);
    elem->set_face_data("E",surface_fluxes.ql);


    elem->set_face_data("G",surface_fluxes.qg0 == -999 ? 0 : surface_fluxes.qg0); //qg0 is the correct ground heatflux to match snowpack output. qg is just uninit

    elem->set_face_data("ilwr_out",surface_fluxes.lw_out);
    elem->set_face_data("iswr_out",surface_fluxes.sw_out);
    elem->set_face_data("dQ",surface_fluxes.dIntEnergy);
//    if(!has_optional("snow_albedo"))
//    {
        elem->set_face_data("spack_albedo",data->Xdata->Albedo);  //even if we have a measured albedo, Xdata will reflect this. //surface_fluxes.pAlbedo);
//    }


}

void Lehning_snowpack::init(mesh domain, boost::shared_ptr <global> global_param)
{
    for(size_t i=0;i<domain->size_faces();i++)
    {
        auto face = domain->face(i);

        auto d = face->make_module_data<Lehning_snowpack::data>(ID);

        // because we use our own config, we need to do the conversion
        //format is same key-val pairs that snowpack expects, case sensitive
        /**
         * [snowpack]
         * [snowpackadvanced]
         */
        for(auto itr : cfg)
        {
            for(auto jtr : itr.second)
            {
                d->config.addKey(jtr.first.data(),itr.first.data(),jtr.second.data());
            }
        }

        //setup critical keys.
        //overwrite the user if a dangerous key is set
//        d->config.addKey("METEO_STEP_LENGTH", "Snowpack", std::to_string(1)); // Hz. I think this should be just 1 for us.
        d->config.addKey("METEO_STEP_LENGTH", "Snowpack", std::to_string(4)); // Hz. I think this should be just 1 for us.
        d->config.addKey("MEAS_TSS", "Snowpack", "false");

        d->config.addKey("CALCULATION_STEP_LENGTH","Snowpack", std::to_string(global_param->dt() /60 ) ); //specified as  minutes

        d->Spackconfig = boost::make_shared<SnowpackConfig>(d->config);


        d->cum_precip=0.;

        //addSpecial keys goes here to deal with Antarctica, canopy, and detect grass

        SN_SNOWSOIL_DATA SSdata;
        SSdata.SoilAlb = cfg.get<double>("sno.SoilAlbedo");
        SSdata.Albedo = SSdata.SoilAlb; // following snowpacks' no snow default.
        SSdata.BareSoil_z0 = cfg.get<double>("sno.BareSoil_z0");
        if (SSdata.BareSoil_z0 == 0.)
        {
            LOG_WARNING << "[snowpack] BareSoil_z0 == 0, set to 0.2";
            SSdata.BareSoil_z0 = 0.2;
        }

        SSdata.WindScalingFactor= cfg.get<double>("sno.WindScalingFactor");
        SSdata.TimeCountDeltaHS = cfg.get<double>("sno.TimeCountDeltaHS");


        SSdata.meta.stationName = cfg.get<std::string>("sno.station_name");
        SSdata.meta.position.setAltitude(face->get_z());

        SSdata.meta.position.setXY(face->get_x(),face->get_y(),face->get_z());
        SSdata.meta.setSlope(face->slope() ,face->aspect());
        SSdata.meta.setSlope(face->slope() * 180./3.14159,face->aspect()* 180./3.14159);


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



        SSdata.Canopy_Height = cfg.get<double>("sno.CanopyHeight");
        SSdata.Canopy_LAI = cfg.get<double>("sno.CanopyLeafAreaIndex");
        SSdata.Canopy_Direct_Throughfall = cfg.get<double>("sno.CanopyDirectThroughfall");

        SSdata.ErosionLevel = cfg.get<double>("sno.ErosionLevel");

        d->Xdata = boost::make_shared<SnowStation>(false,false);
        d->Xdata->initialize(SSdata,0);
//        d->Xdata->windward = false;
//        d->Xdata->rho_hn = 0;
//        d->Xdata->hn = 0;
//        d->Xdata->mH = 0;

        d->sp = boost::make_shared<Snowpack>(*(d->Spackconfig));
        d->meteo = boost::make_shared<Meteo>( (d->config));
        d->stability = boost::make_shared<Stability> ( (d->config), false);

    }
}