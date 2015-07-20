#include "Sicart_ilwr.hpp"

Sicart_ilwr::Sicart_ilwr( std::string ID)
{

    _depends->push_back("t");
    _depends->push_back("ea");
    _depends->push_back("rh");
    _depends->push_back("atm_trans");

    _provides->push_back("Lin");


    this->ID = ID;
    _parallel_type = parallel::data;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

void Sicart_ilwr::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{

    double T = elem->face_data("t")+273.15; //C->K
    double tau = elem->face_data("atm_trans");
    double RH = elem->face_data("rh")/100.0;//requires fractional
    double e = elem->face_data("ea")*0.01; // pa->mb
    double sigma = 5.67*pow(10.0,-8.0); //boltzman

    double Lin = 1.24*pow(e/T,1.0/7.0)*(1.0+0.44*RH-0.18*tau)*sigma*pow(T,4.0);

    elem->set_face_data("Lin", Lin);
}

Sicart_ilwr::~Sicart_ilwr()
{

}