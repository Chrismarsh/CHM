
#include "solar.hpp"

Solar::Solar( std::string ID)
{
    LOG_DEBUG << "Successfully instatiated module " << ID;
}
void Solar::run(mesh_elem& elem)
{
    double A = 120* M_PI/180.0;
    double E = 30* M_PI/180.0;
    
    //radiation data
    //solar vector
    //xyz cartesian
    arma::vec S;
    S << cos(E) * sin(A) << arma::endr
        << cos(E) * cos(A) << arma::endr
        << sin(E) << arma::endr;
                                          
    double angle = acos(arma::dot(S,elem.get_facenormal()));
    angle = cos(angle);

    if(angle < 0.0)
        angle = 0.0;
    
    
    elem.add_face_data("solar_S_angle",angle);
}

Solar::~Solar()
{
    
    
}