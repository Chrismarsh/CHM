#include "Simple_Canopy.hpp"

Simple_Canopy::Simple_Canopy(config_file cfg)
        : module_base(parallel::data)
{
    depends("p_rain");
    depends("p_snow");
    depends("iswr");
    //depends("iswr_diffuse");
    //depends("iswr_direct");
    depends("rh");
    depends("t");
    depends("vw"); // TODO: Should this be scaled up to reference height first? (in wind speed module)
    depends("vw_dir");
    depends("p");
    depends("ilwr");

    provides("ta_subcanopy");
    provides("rh_subcanopy");
    provides("vw_subcanopy");
    provides("wdir_subcanopy");
    provides("p_rain_subcanopy");
    provides("p_snow_subcanopy");
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

    //auto data = elem->get_module_data<Simple_Canopy::data>(ID);

    /**
     * Builds this timestep's meterological data above canopy
     */
    //Mdata.date   = mio::Date( global_param->year(), global_param->month(), global_param->day(),global_param->hour(),global_param->min(),-6 );
    double ta     = elem->face_data("t");
    double rh     = elem->face_data("rh");
    double vw     = elem->face_data("vw");
    double wdir   = elem->face_data("vw_dir");
    double iswr  = elem->face_data("iswr"); // SW in above canopy
    //double diff   = elem->face_data("iswr_diffuse");
    //double ir_h   = elem->face_data("iswr_direct");
    double ilwr  = elem->face_data("ilwr"); // LW in above canopy
    double p_rain = elem->face_data("p_rain"); // rain (mm/timestep) above canopy
    double p_snow = elem->face_data("p_snow"); // snow (mm/timestep) above canopy

    // Do canopy stuff here

    // TEST missing canopy (simply pass on data for now)
    double ta_subcanopy     = ta;
    double rh_subcanopy     = rh;
    double vw_subcanopy     = vw;
    double wdir_subcanopy   = wdir;
    double iswr_subcanopy  = iswr;
    //double diff_subcanopy   = diff;
    //double ir_h_subcanopy   = ir_h;
    double ilwr_subcanopy  = ilwr;
    double p_rain_subcanopy = p_rain;
    double p_snow_subcanopy = p_snow;

    //cout << "we made it to Simple_Canopy";

    // Output computed canopy states and fluxes downward to snowpack and upward to atmosphere

    elem->set_face_data("ta_subcanopy",ta_subcanopy);
    elem->set_face_data("rh_subcanopy",rh_subcanopy);
    elem->set_face_data("vw_subcanopy",vw_subcanopy);
    elem->set_face_data("wdir_subcanopy",wdir_subcanopy);
    elem->set_face_data("iswr_subcanopy",iswr_subcanopy);
    //elem->set_face_data("diff_subcanopy",diff_subcanopy);
    //elem->set_face_data("ir_h_subcanopy",ir_h_subcanopy);
    elem->set_face_data("ilwr_subcanopy",ilwr_subcanopy);
    elem->set_face_data("p_rain_subcanopy",p_rain_subcanopy);
    elem->set_face_data("p_snow_subcanopy",p_snow_subcanopy);


}

void Simple_Canopy::init(mesh domain, boost::shared_ptr <global> global_param)
{
    // TODO: Parallel call needed here?
    // For each face
    for(size_t i=0;i<domain->size_faces();i++)
    {
        // Get current face
        auto face = domain->face(i);

        // What does this do?
        //auto d = face->make_module_data<Simple_Canopy::data>(ID);

        // Get canopy parameters for this face
        //LAI = cfg.get<double>("canopy.LAI");

        // Initialize canopy state variables
        // For each canopy layer
        // Snow in branches average bulk temperature

        // Snow mass (SWE) in branches




    }
}