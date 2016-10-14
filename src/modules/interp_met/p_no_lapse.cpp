#include "p_no_lapse.hpp"

p_no_lapse::p_no_lapse(config_file cfg)
        :module_base(parallel::data)
{

    depends_from_met("p");

    provides("p");


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}
void p_no_lapse::init(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<p_no_lapse::data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
    }
}
void p_no_lapse::run(mesh_elem& face)
{

    std::vector< boost::tuple<double, double, double> > ppt;
    std::vector< boost::tuple<double, double, double> > staion_z;
    for (auto& s : global_param->get_stations( face->get_x(), face->get_y()))
    {
        if( is_nan(s->get("p")))
            continue;
        double u = s->get("p");
        ppt.push_back( boost::make_tuple(s->x(), s->y(), u ) );
        staion_z.push_back( boost::make_tuple(s->x(), s->y(), s->z() ) );
    }

    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double p0 = face->get_module_data<data>(ID)->interp(ppt, query);

    face->set_face_data("p", std::max(0.0,p0));


}

p_no_lapse::~p_no_lapse()
{

}
