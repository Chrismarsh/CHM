#include "point_mode.hpp"

point_mode::point_mode(config_file cfg)
        : module_base(parallel::data)
{

     t = cfg.get<bool>("provide.t",true);
     rh = cfg.get<bool>("provide.rh",true);
     vw = cfg.get<bool>("provide.u",true);
     p = cfg.get<bool>("provide.p",true);
     ilwr = cfg.get<bool>("provide.ilwr",true);

    if(t)
    {
        depends_from_met("t");
        provides("t");
    }

    if(rh)
    {
        depends_from_met("rh");
        provides("rh");
    }

    if(vw)
    {
        depends_from_met("u");
        provides("vw");
        provides("vw_dir");
    }

    if(p)
    {
        depends_from_met("p");
        provides("p");
    }

    if(ilwr)
    {
        depends_from_met("Qli");
        provides("ilwr");
    }

}

point_mode::~point_mode()
{

}

void point_mode::run(mesh_elem &elem, boost::shared_ptr <global> global_param)
{
    // at this point, if the user has provided more than 1 station, they've been stripped out.
    // we can safetly take the 1st (and only) station.

    if(t)
    {
        double st = global_param->stations.at(0)->get("t");
        elem->set_face_data("t",st);
    }

    if(rh)
    {
        double srh = global_param->stations.at(0)->get("rh");
        elem->set_face_data("rh", srh);
    }
    if(vw)
    {
        double su = global_param->stations.at(0)->get("u");

        //make sure we don't have zero windpseeds
        su = std::max(su,0.5);
        elem->set_face_data("vw",su);
        elem->set_face_data("vw_dir",270.); //TODO: real wind direction
    }

    if(p)
    {
        double sp = global_param->stations.at(0)->get("p");
        elem->set_face_data("p", sp);
    }
    if(ilwr)
    {
        double silwr = global_param->stations.at(0)->get("Qli");
        elem->set_face_data("ilwr", silwr);

    }
}

