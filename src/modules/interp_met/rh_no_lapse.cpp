#include "rh_no_lapse.hpp"

rh_no_lapse::rh_no_lapse(config_file cfg)
        : module_base(parallel::data)
{
    provides("rh");
    depends_from_met("rh");


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

rh_no_lapse::~rh_no_lapse()
{

}
void rh_no_lapse::init(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<rh_no_lapse::data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
    }
}
void rh_no_lapse::run(mesh_elem &face)
{

    std::vector<boost::tuple<double, double, double> > lowered_values;
    for (auto &s : global_param->get_stations( face->get_x(), face->get_y()))
    {
        if( is_nan(s->get("rh")))
            continue;
        double rh = s->get("rh");

        lowered_values.push_back(boost::make_tuple(s->x(), s->y(), rh));
    }

    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double rh = face->get_module_data<data>(ID)->interp(lowered_values, query);

    rh = std::min(rh, 100.0);
    rh = std::max(10.0, rh);
    face->set_face_data("rh", rh);

}

