//
// Created by chris on 17/11/15.
//

#include "Thornton_var_p.h"

Thornton_var_p::Thornton_var_p(config_file cfg)
: module_base(parallel::data)
{
    depends_from_met("p");

    provides("p");
    provides("p_lapse");

}
Thornton_var_p::~Thornton_var_p()
{

}
template <typename T> int signum(T val) {
    return (T(0) < val) - (val < T(0));
}
void Thornton_var_p::init(mesh domain, boost::shared_ptr<global> global_param)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->stations.size());
    }
}
void Thornton_var_p::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{

    //generate lapse rates
    std::vector<double> sp;
    std::vector<double> sz;


    static boost::posix_time::ptime last_update;
    static double lapse=-999.0;

    //we have to run this at least once, use the first bool flag to do so
    //otherwise just check that the last update was at a different time step, if so calculate the lapse rate.
    //otherwise, just used the stored lapse rate
    if(last_update != global_param->posix_time() )
    {
        for (auto& s : global_param->stations)
        {
            double u = s->get("p");
            sp.push_back( u  );
            sz.push_back( s->z());
        }

        size_t n = sp.size();

        std::vector< boost::tuple<double,double> > combinations;
        gsl_combination* c = gsl_combination_calloc(n,2);
        do{
            auto tp = boost::make_tuple( c->data[0] , c->data[1]);
            combinations.push_back(tp);
        }while (gsl_combination_next (c) == GSL_SUCCESS);

        std::vector<double> normalize_precip;
        std::vector<double> z_diff;

        for(auto itr : combinations)
        {
            int i0 =itr.get<0>();
            int i1 =itr.get<1>();

            double p1 = sp[i0];
            double p2 = sp[i1];

            double np = (p1-p2)/(p1+p2);
            double z = sz[i0] - sz[i1];

            normalize_precip.push_back(np);
            z_diff.push_back(z);
        }
        // least squares linear fit to these points ( p v. z)

        double c0, c1, cov00, cov01, cov11, chisq;

        gsl_fit_linear (&z_diff[0], 1,&normalize_precip[0], 1, n,
                        &c0, &c1, &cov00, &cov01, &cov11,
                        &chisq);
        lapse = c1;//use the slope (y=mx+b c1 == m)
        last_update = global_param->posix_time();
    }

    elem->set_face_data("p_lapse",lapse);

    //now do the full interpolation
    std::vector< boost::tuple<double, double, double> > ppt;
    std::vector< boost::tuple<double, double, double> > station_z;
    for (auto& s : global_param->stations)
    {
        double u = s->get("p");
        ppt.push_back( boost::make_tuple(s->x(), s->y(), u ) );
        station_z.push_back( boost::make_tuple(s->x(), s->y(), s->z() ) );
    }



    auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
    double p0 = elem->get_module_data<data>(ID)->interp(ppt, query);
    double z0 = elem->get_module_data<data>(ID)->interp(station_z,query);
    double z = elem->get_z();

    double f = lapse*(z-z0);

    //as per Thorton 1997, if f > fmax, the error goes up, so constrain it here.
    double fmax = 0.95;
    if(fabs(f) >= fmax)
    {
        f = signum(f) * 0.95;
    }

    double P = p0*( (1+f)/(1-f));
    P = std::max(0.0,P);

    elem->set_face_data("p", P);


}