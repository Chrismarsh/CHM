
#include "solar.hpp"

Solar::Solar( std::string ID)
{
    _provides->push_back("solar_angle");
    _provides->push_back("Qsi");
    _provides->push_back("S0");
    
    _depends->push_back("shadowed");
    _depends->push_back("atm_trans");
    
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

    if(angle < 0.0 || E < 0.0523598776) //3deg -> rad
        angle = 0.0;
    
    
    elem->set_face_data("solar_angle",angle);
    
    double shadow = elem->face_data("shadowed");
    if(shadow == 1)
        angle = 0;
    double tau = elem->face_data("atm_trans");
    elem->set_face_data("Qsi", tau*angle*1375.0 );
    elem->set_face_data("S0",angle*1375.0 );

    
    
}

Solar::~Solar()
{
    
    
}



