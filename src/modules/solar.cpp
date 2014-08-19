
#include "solar.hpp"

Solar::Solar( std::string ID)
{
    _provides->push_back("solar_S_angle");
    _provides->push_back("solar_short");
    
    _depends->push_back("z_prime");
    
    
    this->ID = ID;
    _parallel_type = parallel::data;
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

    
    Vector_3 n = elem->normal();
        arma::vec N;
    
    N << n[0] << arma::endr
      << n[1] << arma::endr
      << n[2] << arma::endr;
    
    double angle = acos(arma::dot(S,N));
    angle = cos(angle);

    if(angle < 0.0)
        angle = 0.0;
    
    
    elem->set_face_data("solar_S_angle",angle);
    
    double shadow = elem->face_data("shadowed");
    
    elem->set_face_data("solar_short", shadow == 1 ? 0 : angle ); //shadow == 1
    
    
}

Solar::~Solar()
{
    
    
}



