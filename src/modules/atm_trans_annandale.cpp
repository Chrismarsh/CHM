#include "atm_trans_annandale.hpp"

atm_trans_annandale::atm_trans_annandale( std::string ID)
{

    _depends->push_back("t");

    _provides->push_back("atm_trans");

    this->ID = ID;
    _parallel_type = parallel::data;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

void atm_trans_annandale::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{
    double z = elem->get_z();
    double k_rs = 0.16; //inland
    timeseries& ts = *(elem->get_underlying_timeseries());
    timeseries::iterator now = elem->now();

    double tmax = daily::max(ts,now,"t");
    double tmin = daily::min(ts,now,"t");
    double tau = k_rs *(1 + 2.7*pow(10,-5.0)*z)*pow(tmax-tmin,0.5);

    elem->set_face_data("atm_trans", tau);


}

atm_trans_annandale::~atm_trans_annandale()
{

}