#include "kunkel_rh.hpp"

kunkel_rh::kunkel_rh()
        : module_base(parallel::data)
{
    provides("rh");
    depends_from_met("rh");


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

kunkel_rh::~kunkel_rh()
{

}

void kunkel_rh::run(mesh_elem &elem, boost::shared_ptr <global> global_param)
{
    // 1/km
    double lapse_rates[] =
            {-0.09,
             0.0,
             0.09,
             0.11,
             0.11,
             0.12,
             0.14,
             0.15,
             0.11,
             0.07,
             -0.02,
             -0.07
            };

    double lapse = lapse_rates[global_param->month() - 1] / 1000.0; // -> 1/m
    std::vector<boost::tuple<double, double, double> > lowered_values;
    for (auto &s : global_param->stations)
    {

        double rh = s->get("rh");

        double rh_z = rh * exp(lapse * (0.0 - s->z()));

        lowered_values.push_back(boost::make_tuple(s->x(), s->y(), rh_z));
    }

    interp_base *interp = nullptr;
    std::string interp_method = "spline";
    if (interp_method == "spline")
        interp = new thin_plate_spline();

    auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
    double value = (*interp)(lowered_values, query);//C

    double rh = value * exp(lapse * (elem->get_z() - 0.0));

    rh = std::min(rh, 100.0);
    rh = std::max(10.0, rh);
    elem->set_face_data("rh", rh);

}

