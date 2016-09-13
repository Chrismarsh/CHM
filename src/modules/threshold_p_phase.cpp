#include "threshold_p_phase.hpp"

threshold_p_phase::threshold_p_phase(config_file cfg)
    : module_base(parallel::data)
{
    depends("t");
    depends("p");

    t_thresh = cfg.get("threshold_temperature", 2.0) ;

    provides("frac_precip_rain");
    provides("frac_precip_snow");

    provides("p_snow");
    provides("p_rain");

}

threshold_p_phase::~threshold_p_phase()
{

}

void threshold_p_phase::run(mesh_elem &face)
{
    double p = face->face_data("p");

    if( face->face_data("t") >= t_thresh)
    {
        face->set_face_data("frac_precip_rain", 1);
        face->set_face_data("frac_precip_snow", 0);

        face->set_face_data("p_rain", p);
        face->set_face_data("p_snow", 0 );
    } else
    {
        face->set_face_data("frac_precip_rain", 0);
        face->set_face_data("frac_precip_snow", 1);

        face->set_face_data("p_rain", 0);
        face->set_face_data("p_snow", p );
    }



}