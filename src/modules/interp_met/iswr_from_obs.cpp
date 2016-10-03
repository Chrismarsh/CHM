
#include "iswr_from_obs.h"


iswr_from_obs::iswr_from_obs(config_file cfg)
        : module_base(parallel::data)
{
    depends_from_met("Qsi");
    depends("solar_el");

    provides("iswr_direct");
    provides("iswr_diffuse");
    provides("iswr");
    provides("atm_trans");
}
iswr_from_obs::~iswr_from_obs()
{

}
void iswr_from_obs::init(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
    }
}
void iswr_from_obs::run(mesh_elem &face)
{
    //interpolate all the measured qsi

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : global_param->get_stations( face->get_x(), face->get_y()))
    {
        if( is_nan(s->get("Qsi")))
            continue;
        double v = s->get("Qsi");
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), v ) );
    }


    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double iswr_measured =face->get_module_data<data>(ID)->interp(lowered_values, query);

    double splitting_coef = 0.0;

    double elevation = face->face_data("solar_el");
    double iswr_modeled = 1375.0; //top of atmosphere
    double elevation_threshold = 3.0;

    if( elevation < elevation_threshold ) {
        //when the Sun is low above the horizon, Mt is getting abnormaly too large pretending
        // this is a clear sky day when almost all the radiation should be diffuse
        // no matter how the sky is
        splitting_coef = 1.0;
    } else {
        // clear sky index (ratio global measured to top of atmosphere radiation)
        const double Mt = iswr_measured / iswr_modeled; // should be <=1.2, aka clearness index Kt
        const double clear_sky = 0.147;

        // diffuse fraction: hourly ratio of diffuse to global radiation incident on a horizontal surface
        // splitting according to a combination of Reindl et al.(1990)'s models (Mt-model and Mt&Psolar->elev-model):
        if( Mt >= 0.78 ) { // Mt in [0.78;1] -> clear day
            splitting_coef = clear_sky;
        } else {
            if( Mt <= 0.3 ) { // Mt in [0;0.3] -> overcast
                splitting_coef = 1.02 - 0.248*Mt;
                if(splitting_coef>1.) splitting_coef=1.;
            } else {           // Mt in ]0.3;0.78[ -> cloudy
                splitting_coef = 1.4 - 1.749*Mt + 0.177*sin(elevation*mio::Cst::to_rad);
                if(splitting_coef>0.97) splitting_coef = 0.97;
                if(splitting_coef<clear_sky) splitting_coef = clear_sky;
            }
        }
    }


    double obs_dir = iswr_measured * splitting_coef;
    double obs_diff = iswr_measured * (1-splitting_coef);

    face->set_face_data("iswr_direct",obs_dir);
    face->set_face_data("iswr_diffuse",obs_diff);
    face->set_face_data("iswr",iswr_measured);
    face->set_face_data("atm_trans",iswr_measured/1375.);
}

