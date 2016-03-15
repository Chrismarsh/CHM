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

void threshold_p_phase::run(mesh_elem &elem, boost::shared_ptr <global> global_param)
{
    double p = elem->face_data("p");

    if( elem->face_data("t") >= t_thresh)
    {
        elem->set_face_data("frac_precip_rain", 1);
        elem->set_face_data("frac_precip_snow", 0);

        elem->set_face_data("p_rain", p);
        elem->set_face_data("p_snow", 0 );
    } else
    {
        elem->set_face_data("frac_precip_rain", 0);
        elem->set_face_data("frac_precip_snow", 1);

        elem->set_face_data("p_rain", 0);
        elem->set_face_data("p_snow", p );
    }



}