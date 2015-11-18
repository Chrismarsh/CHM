#include "Thornton_p.hpp"

Thornton_p::Thornton_p()
        :module_base(parallel::data)
{

    depends_from_met("p");

    provides("p");


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

void Thornton_p::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{
    //km^-1
    double monthly_factors[] = {0.35, 0.35, 0.35, 0.30, 0.25, 0.20, 0.20, 0.20, 0.20, 0.25, 0.30, 0.35};
    for(auto& mf: monthly_factors)
    {
        mf /= 1000.0; //to m^-1
    }
    std::vector< boost::tuple<double, double, double> > ppt;
    std::vector< boost::tuple<double, double, double> > staion_z;
    for (auto& s : global_param->stations)
    {
        double u = s->get("p");
        ppt.push_back( boost::make_tuple(s->x(), s->y(), u ) );
        staion_z.push_back( boost::make_tuple(s->x(), s->y(), s->z() ) );
    }

    thin_plate_spline interp;

    auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
    double p0 = interp(ppt, query);
    double z0 = interp(staion_z,query);
    double z = elem->get_z();

    int month = global_param->month()-1;
    double lapse = monthly_factors[month];

    double P = p0*( (1+lapse*(z-z0))/(1-lapse*(z-z0)));
    P = std::max(0.0,P);

    elem->set_face_data("p", P);


}

Thornton_p::~Thornton_p()
{

}