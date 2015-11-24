#include <meteoio/meteolaws/Meteoconst.h>
#include "slope_iswr.hpp"

slope_iswr::slope_iswr()
        :module_base(parallel::data)
{
    depends("iswr");
    depends("iswr_diffuse");
    depends("iswr_direct");

    provides("iswr");
    provides("iswr_direct");
    provides("solar_angle");

    optional("shadowed");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}
void slope_iswr::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{
    

    double A = global_param->solar_az() * mio::Cst::to_rad;
    double E = global_param->solar_el() * mio::Cst::to_rad;
    
    //radiation data
    //solar vector
    //xyz cartesian
    arma::vec S;
    S << cos(E) * sin(A) << arma::endr
      << cos(E) * cos(A) << arma::endr
      << sin(E) << arma::endr;

    
    Vector_3 n = elem->normal();
        arma::vec N;
    
    N << n[0] << arma::endr
      << n[1] << arma::endr
      << n[2] << arma::endr;
    
    double angle = acos(arma::dot(S,N));
    angle = cos(angle);

    if(angle < 0.0 || E < 0.0523598776) //3deg -> rad
        angle = 0.0;

    elem->set_face_data("solar_angle",angle);

    //if we have remote shadowing
    if(has_optional("shadow"))
    {
        if(elem->face_data("shadow") == 1)
        {
            angle = 0;
        }
    }

    double swr =  angle * elem->face_data("iswr_direct");
    double diff = elem->face_data("iswr_diffuse");

    elem->set_face_data("iswr_direct",swr );
    elem->set_face_data("iswr", swr + diff );

}

slope_iswr::~slope_iswr()
{
    
    
}



