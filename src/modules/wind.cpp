#include "wind.hpp"

wind::wind( std::string ID)
{

   _depends_from_met->push_back("u");

   _provides->push_back("u");


    this->ID = ID;
    _parallel_type = parallel::data;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

void wind::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{

    std::vector< boost::tuple<double, double, double> > values;
    for (auto& s : global_param->stations)
    {
        double u = s->get("u");
        values.push_back( boost::make_tuple(s->x(), s->y(), u ) );
    }

    interp_base* interp=nullptr;
    std::string interp_method = "spline";
    if(interp_method == "spline")
        interp = new thin_plate_spline();

    auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
    double value = (*interp)(values, query);

    elem->set_face_data("u", value);
}

wind::~wind()
{

}