
#include "Burridge_iswr.h"

Burridge_iswr::Burridge_iswr()
        :module_base(parallel::data)
{

    depends("cloud_frac");

    provides("iswr");
    provides("iswr_diffuse");
    provides("iswr_direct");
    provides("atm_trans");
}

Burridge_iswr::~Burridge_iswr()
{

}

void Burridge_iswr::run(mesh_elem &elem, boost::shared_ptr<global> global_param)
{
    double cosZ = cos( (90.0-global_param->solar_el()) *mio::Cst::to_rad);

//    double aspect_south0 = elem->aspect() * mio::Cst::to_deg;
//    if (aspect_south0 >= 180.0)
//        aspect_south0 -=  180.0;
//    else
//        aspect_south0 += 180.0;
//    aspect_south0 *= mio::Cst::to_rad;
//
//    double slope = elem->slope();
//    double sun_az = global_param->solar_az(); //* mio::Cst::to_rad;
//    if (sun_az >= 180.0)
//        sun_az -=  180.0;
//    else
//        sun_az += 180.0;
//    sun_az *= mio::Cst::to_rad;
//
//
//    double sinZ = sqrt(1.0 - cosZ*cosZ);
//    double cosi = cos(slope) * cosZ +
//              sin(slope) * sinZ *
//              cos(sun_az - aspect_south0);
//
//
//    if (cosi < 0.0)
//        cosi = 0.0;
//    if(cosZ <= 0.0)
//        cosZ=0.0;

    double S = 1375.0;
    double cf = elem->face_data("cloud_frac");

    double dir = S  * (0.6+0.2*cosZ)*(1.0-cf);
    double diff = S * (0.3+0.1*cosZ)*(cf);

 //   dir = dir * cosi;
    diff = diff*cosZ;


    if (diff <0)
        diff = 0.0;
    if(dir <0)
        dir = 0.0;

    elem->set_face_data("iswr_diffuse",diff);
    elem->set_face_data("iswr_direct",dir);
    elem->set_face_data("iswr",dir+diff);
    elem->set_face_data("atm_trans", (dir+diff) / 1375.);
}
