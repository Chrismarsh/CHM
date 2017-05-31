
#include "iswr_from_obs.h"


iswr_from_obs::iswr_from_obs(config_file cfg)
        : module_base(parallel::data)
{
    depends_from_met("Qsi");
    depends("solar_el");

    provides("iswr_direct_no_slope");
    provides("iswr_diffuse_no_slope");
    provides("iswr_observed");

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
    double iswr_observed =face->get_module_data<data>(ID)->interp(lowered_values, query);


    // This is what is used in SUMMA
    //! compute the fraction of direct radiation using the parameterization of Nijssen and Lettenmaier (1999)
    double Frad_direct = 0.7000;
    double directScale = 0.0900;
    double cosZenith = cos( M_PI/2.0 -  (face->face_data("solar_el")*mio::Cst::to_rad)) ; //zenith is from 90 vert -> 0 horz
    double scalarFractionDirect = 0;

    if(cosZenith > 0. )
        scalarFractionDirect = Frad_direct*cosZenith/(cosZenith + directScale);

    double split_dir = iswr_observed * scalarFractionDirect;
    double split_diff = iswr_observed - split_dir;

    face->set_face_data("iswr_direct_no_slope",split_dir);
    face->set_face_data("iswr_diffuse_no_slope",split_diff);
    face->set_face_data("iswr_observed",iswr_observed);

    face->set_face_data("atm_trans",split_dir/1375.);
}

