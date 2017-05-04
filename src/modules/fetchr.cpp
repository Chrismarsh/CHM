
#include "fetchr.hpp"

fetchr::fetchr(config_file cfg)
        : module_base(parallel::data)
{
    depends("vw_dir");

    provides("fetch");

    //number of steps along the search vector to check for a higher point
    steps = cfg.get("steps",10);
    //max distance to search
    max_distance = cfg.get("max_distance",1000.0);

    //size of the step to take
    size_of_step = max_distance / steps;

    I = cfg.get("I",0.06);

    incl_veg = cfg.get("incl_veg",true);
}

fetchr::~fetchr()
{


}

void fetchr::run(mesh_elem& face)
{

    face->set_face_data("fetch", max_distance);

    //direction it is from, need upwind fetch
    double wind_dir = face->face_data("vw_dir") ;

    //if we are using vegetation and the current face is covered in veg, set the fetch to 0
    if(incl_veg && face->has_parameter("landcover"))
    {
        int me_LC = face->get_parameter("landcover");

        double me_Z_CanTop = global_param->parameters.get<double>("landcover." + std::to_string(me_LC) + ".CanopyHeight");
        if(me_Z_CanTop > 1) // 1m might be too high?
        {
            face->set_face_data("fetch", 0);
            return;
        }

    }

    // search along wind_dir azimuth in j step increments
    for (int j = 1; j <= steps; ++j)
    {
        double distance = j * size_of_step;

        auto f = face->find_closest_face(wind_dir, distance);

        double Z_CanTop = 0;
        if (incl_veg && f->has_parameter("landcover"))
        {
            int LC = f->get_parameter("landcover");
            Z_CanTop = global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".CanopyHeight");
        }

        //include canopy height if available
        double Z_test =  f->center().z()+Z_CanTop;

        //equation 1, pg 771, Lapen and Martz 1993
        double Z_core = face->center().z() + distance*I;

        //reset the fetch if there is elevation, but also if we run into vegetation >1m
        if(Z_test >= Z_core /*|| Z_CanTop > 1*/)
        {
            face->set_face_data("fetch", distance);
            break;
        }
    }




}

