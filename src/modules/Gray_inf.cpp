#include "Gray_inf.hpp"

Gray_inf::Gray_inf(config_file cfg)
        : module_base(parallel::data)
{

    depends("swe");
    depends("snowmelt_int");

    provides("inf");
    provides("total_inf");
    provides("total_excess");
    provides("runoff");
    provides("soil_storage");
    provides("potential_inf");
    provides("opportunity_time");
    provides("available_storage");

}

Gray_inf::~Gray_inf()
{

}

void Gray_inf::init(mesh domain)
{
    //store all of snobals global variables from this timestep to be used as ICs for the next timestep
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto* d = face->make_module_data<Gray_inf::data>(ID);

        d->soil_depth = 400;
        d->porosity = .4;
        d->max_storage =  d->soil_depth * d->porosity;
//        d->storage =  d->max_storage  - (1 - d->max_storage * face->get_parameter("sm")/100.);
        d->storage =  d->max_storage * face->get_parameter("sm")/100.;

        d->last_ts_potential_inf = 0;
        d->opportunity_time=0.;
        d->total_inf = 0.;
        d->total_excess = 0.;

    }
}
void Gray_inf::run(mesh_elem &face)
{
    if(is_water(face))
    {
        set_all_nan_on_skip(face);
        return;
    }


    auto* d = face->get_module_data<Gray_inf::data>(ID);

    auto id = face->cell_id;
    double C = 2.;
    double S0 = 1;
    double SI = face->get_initial_condition("sm")/100.;

    double TI = 272.;

    double runoff = 0.;
    double inf = 0.;

    double snowmelt = face->face_data("snowmelt_int");

    double potential_inf = d->last_ts_potential_inf;
    double avail_storage = (d->max_storage - d->storage);

    if(snowmelt > 0)
    {
        d->opportunity_time += global_param->dt() / 3600.;
        double t0 = d->opportunity_time;
        potential_inf = C * pow(S0,2.92) * pow((1. - SI),1.64) * pow((273.15 - TI) / 273.15, -0.45) * pow(t0,0.44);


        //cap the total infiltration to be no more than our available storage
        if (potential_inf > avail_storage )
        {
            potential_inf = avail_storage;
        }

        d->last_ts_potential_inf = potential_inf;

        if( d->total_inf + snowmelt > potential_inf)
        {
            runoff = (d->total_inf + snowmelt) - potential_inf;
        }


        //infiltrate everything else
        inf    = snowmelt - runoff;


        d->storage += inf;

        if (d->storage > d->max_storage)
        {
            d->storage = d->max_storage;
        }


        d->total_inf += inf;
        d->total_excess += runoff;

    }


    if( !is_nan(face->get_initial_condition("sm")))
    {
        face->set_face_data("total_excess",d->total_excess);
        face->set_face_data("total_inf",d->total_inf);

        face->set_face_data("runoff",runoff);
        face->set_face_data("inf",inf);
        face->set_face_data("potential_inf",potential_inf);
        face->set_face_data("soil_storage", d->storage);
        face->set_face_data("opportunity_time",d->opportunity_time);
        face->set_face_data("available_storage",avail_storage);
    }
    else
    {
        set_all_nan_on_skip(face);
    }

}
