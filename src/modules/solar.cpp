
#include "solar.hpp"

Solar::Solar( std::string ID)
{
    _provides->push_back("solar_S_angle");
    this->ID = ID;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}
void Solar::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{
    

    double A = global_param->solar_az() * M_PI/180.0;
    double E = global_param->solar_el() * M_PI/180.0;
    
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



