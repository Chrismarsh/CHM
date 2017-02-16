#include "Richard_albedo.hpp"

Richard_albedo::Richard_albedo(config_file cfg)
: module_base(parallel::data)
{

    depends("swe");
    depends("T_s"); // snow temp
    depends("p_snow");

    provides("snow_albedo");
    provides("melting_albedo");

}
Richard_albedo::~Richard_albedo()
{

}



void Richard_albedo::run(mesh_elem &face)
{
    if(is_water(face))
    {
        set_all_nan_on_skip(face);
        return;
    }
    double albedo = face->get_module_data<Richard_albedo::data>(ID)->albedo;

    double swe = face->face_data("swe");

    double dt = global_param->dt();
    if(global_param->first_time_step)
    {
        if(swe > 1.0)
        {
            albedo = albedo_snow;
        }
        else
        {
            albedo = albedo_bare;
        }
    }

    if ( swe > 0.)
    {
        if (face->face_data("T_s_0") >= 273.)  //melting snow, T_s_0???
        {
            albedo = (albedo - amin) * exp(-dt/a2) + amin;
            face->set_face_data("melting_albedo",1);
        }
        else
        {
            albedo = albedo - dt/a1; // cold snow decay
            face->set_face_data("melting_albedo",0);
        }

        double psnow = face->face_data("p_snow");
        albedo = albedo + (amax - albedo) * (psnow )/min_swe_refresh; //* dt

        albedo = std::max(albedo,amin);
        albedo = std::min(albedo,amax);
    }
    else
    {
        face->set_face_data("melting_albedo",0);
        albedo = albedo_bare;
    }


    face->set_face_data("snow_albedo",albedo);
    face->get_module_data<Richard_albedo::data>(ID)->albedo = albedo;
}

void Richard_albedo::init(mesh domain)
{

    //these
    amin = cfg.get("albedo_min",0.5);
    amax = cfg.get("albedo_max",0.84);
    a1 = cfg.get("a1",1.08e7); //cold snow decay const
    a2 = cfg.get("a2",7.2e5);  //melting snow decay const
    min_swe_refresh = cfg.get("min_swe_refresh",1.); //min swe to refresh albedo
    albedo_snow = cfg.get("init_albedo_snow",0.85); // intial snow albedo
    albedo_bare = cfg.get("init_albedo_bare",0.17); //initial bare ground albedo

#pragma omp parallel for
    for (auto i=0; i < domain->size_faces(); ++i)
    {
        auto face = domain->face(i);
        auto* d = face->make_module_data<data>(ID);
        d->albedo = albedo_bare;
    }
}
