#include "kunkel_rh.hpp"

kunkel_rh::kunkel_rh(config_file cfg)
        : module_base(parallel::data)
{
    provides("rh");
    depends_from_met("rh");


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

kunkel_rh::~kunkel_rh()
{

}
void kunkel_rh::init(mesh domain, boost::shared_ptr<global> global_param)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations_in_radius( face->get_x(), face->get_y(), global_param->station_search_radius).size());
    }
}
void kunkel_rh::run(mesh_elem &face, boost::shared_ptr <global> global_param)
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
    for (auto &s : global_param->get_stations_in_radius( face->get_x(), face->get_y(), global_param->station_search_radius))
    {
        if( is_nan(s->get("rh")))
            continue;
        double rh = s->get("rh");

        double rh_z = rh * exp(lapse * (0.0 - s->z()));

        lowered_values.push_back(boost::make_tuple(s->x(), s->y(), rh_z));
    }


    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double value = face->get_module_data<data>(ID)->interp(lowered_values, query);//C

    double rh = value * exp(lapse * (face->get_z() - 0.0));

    rh = std::min(rh, 100.0);
    rh = std::max(10.0, rh);
    face->set_face_data("rh", rh);

}

