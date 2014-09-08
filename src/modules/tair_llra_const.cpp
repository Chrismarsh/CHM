#include "tair_llra_const.hpp"

tair_llra_const::tair_llra_const( std::string ID)
{
    _provides->push_back("t");

    _provides->push_back("tair_llra_const");

    this->ID = ID;
    _parallel_type = parallel::data;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

tair_llra_const::~tair_llra_const()
{

}
void tair_llra_const::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{

    double lapse_rate = 0.0065;


    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : global_param->stations)
    {
        double v = s->get(global_param->get_variable("Tair")) - lapse_rate * (0.0 - s->z());
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), v ) );
    }

    interp_base* interp=nullptr;
    std::string interp_method = "spline";
    if(interp_method == "spline")
        interp = new thin_plate_spline();

    auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
    double value = (*interp)(lowered_values, query);

    //raise value back up to the face's elevation from sea level
    value =  value + lapse_rate * (0.0 - elem->get_z());

    elem->set_face_data(global_param->get_variable("Tair"),value);
    elem->set_face_data("tair_llra_const",value);

    delete interp;
}