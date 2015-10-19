#include "const_llra_ta.hpp"

const_llra_ta::const_llra_ta()
        :module_base(parallel::data)

{
    provides("t");
    provides("const_llra_ta");

    depends_from_met("t");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

const_llra_ta::~const_llra_ta()
{

}

void const_llra_ta::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
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
    elem->set_face_data("const_llra_ta",value);

    delete interp;
}


