#include "soil_module.hpp"

REGISTER_MODULE_CPP(soil_module);

soil_module::soil_module(config_file cfg) : module_base("soil_module", parallel::data, cfg)
{
    depends("swe");
    depends("thaw_front_depth"); 
    depends("freeze_front_depth");
    depends("first_front_depth");
    depends("ET");
    depends("inf");
    depends("runoff");
//    depends("routing_residual");

    provides("condensation");
    provides("actual_soil_ET");
    provides("soil_excess_to_runoff");
    provides("soil_excess_to_gw");
	provides("runoff_to_depression");
    provides("ground_water_out");
    provides("soil_to_ssr");
    provides("rechr_to_ssr");
    provides("soil_storage");
    provides("soil_rechr_storage");
    provides("depression_storage");
    provides("ground_water_storage");
    provides("detention_storage");
    provides("K_rechr_to_ssr");
    provides("K_lower_to_ssr");
    provides("K_detention_to_runoff");
    provides("K_depression_to_ssr");
    provides("K_depression_to_gw");
    provides("K_ground_water_out");
    provides("K_soil_to_gw");
};

soil_module::~soil_module()
{

};

void soil_module::init(mesh& domain)
{

    SoilDataObj = std::make_unique<Soil::soils_na>();

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<soil_module::data>(ID);
        // I do some evil things here to allow for the submodules to access module_base functions like is_water
        // A pointer to face is put in d, likewise a pointer to this instance of this class is also added, see the overridden functions
        // get_dt and is_lake below.
        // Changed is_lake to a variable from a function, so it stores the result of is_water rather than requiring a copy of face.
		//d.my_face = &face;
        set_local_module(d);
        set_soil_params(face,d);
   
        // dependency injection of K_estimation
        d.K_estimator = std::make_unique<K_estimate>(d);
        d.soil_layers = std::make_unique<soil_two_layer>(d,*d.K_estimator);
        // ET coud (should) have been dependecy injection. So TODO
        d.ET = std::make_unique<soil_ET>(d);  
        //if (d.soil_layers) Removed if statments for now because ET and soil modules are coupled. 
                //if (d.ET)
        set_ET_params(face,d);
    
        initial_soil_conditions(face,d);
    
        set_soil_outputs(face,d);

    }
};

void soil_module::run(mesh_elem& face)
{

    auto& d = face->make_module_data<soil_module::data>(ID);
    
    get_soil_inputs(face,d);


    // order of operations here is hard coded, but it wouldn't be physically wrong to impose ET before soil
    // This is fine because this soil module is run in this order.

    d.soil_layers->run();
    
    if (d.swe == 0.0)
        d.ET->run();
    else
        d.actual_soil_ET = 0.0;

    set_soil_outputs(face,d);


};

void soil_module::get_soil_inputs(mesh_elem& face,soil_module::data& d)
bool soil_module::is_new_day()
{
    // TODO This has hard coded elements, Chris suggested something different here: https://godbolt.org/z/3c51T1avT
    int td = global_param->posix_time().time_of_day().total_seconds();
    int time_to_midnight = 86400 - td;
    if (td >= 0 && td < global_param->dt()) //(time_to_midnight >= global_param->dt())
    {
        return true;
    }
    else
        return false;
};

{
    d.swe = (*face)["swe"_s];
    d.thaw_front_depth = (*face)["thaw_front_depth"_s];
    d.freeze_front_depth = (*face)["freeze_front_depth"_s];
    d.freeze_thaw_first_front = (*face)["first_front_depth"_s];
    d.actual_ET = (*face)["ET"_s];
    d.infil = (*face)["inf"_s];
    d.runoff = (*face)["runoff"_s];
    d.routing_residual = 0.0; //(*face)["routine_residual"_s];
	d.is_lake = is_water(face);
};

void soil_module::set_soil_outputs(mesh_elem& face,soil_module::data& d)
{
    (*face)["condensation"_s] = d.condensation;
    (*face)["actual_soil_ET"_s] = d.actual_soil_ET; 
    (*face)["soil_excess_to_runoff"_s] = d.soil_excess_to_runoff; 
    (*face)["soil_excess_to_gw"_s] = d.soil_excess_to_gw; 
    (*face)["runoff_to_depression"_s] = d.runoff_to_depression;
	(*face)["ground_water_out"_s] = d.ground_water_out; 
    (*face)["soil_to_ssr"_s] = d.soil_to_ssr;
    (*face)["rechr_to_ssr"_s] = d.rechr_to_ssr;
    (*face)["soil_storage"_s] = d.soil_storage;
    (*face)["soil_rechr_storage"_s] = d.soil_rechr_storage;
    (*face)["depression_storage"_s] = d.depression_storage;
    (*face)["ground_water_storage"_s] = d.ground_water_storage;
    (*face)["detention_storage"_s] = d.detention_storage;
    (*face)["K_rechr_to_ssr"_s] = d.K_rechr_to_ssr;
    (*face)["K_lower_to_ssr"_s] = d.K_lower_to_ssr;
    (*face)["K_detention_to_runoff"_s] = d.K_detention_to_runoff;
    (*face)["K_depression_to_ssr"_s] = d.K_depression_to_ssr;
    (*face)["K_depression_to_gw"_s] = d.K_depression_to_gw;
    (*face)["K_ground_water_out"_s] = d.K_ground_water_out;
    (*face)["K_soil_to_gw"_s] = d.K_soil_to_gw;
};

void soil_module::set_soil_params(mesh_elem& face, soil_module::data& d)
{
    // TODO actually connect to stuff
    if (face->has_soil())
    {
        d.soil_storage_max = face->soil_attribute<double>("soil_storage_max"_s);
        d.soil_rechr_max = face->soil_attribute<double>("soil_rechr_max"_s);
        d.excess_to_ssr = face->soil_attribute<bool>("excess_to_ssr"_s); 
        d.detention_snow_max = face->soil_attribute<double>("detention_snow_max"_s);
        d.detention_organic_max = face->soil_attribute<double>("detention_organic_max"_s);
        d.depression_max = face->soil_attribute<double>("depression_max"_s);
        d.ground_water_max = face->soil_attribute<double>("ground_water_max"_s);
        d.local_slope = face->soil_attribute<double>("local_slope"_s)*3.14159265/180;
        
        d.pore_size_dist = face->soil_attribute<double>("PSD_K_estimator");
        d.pore_size_dist_organic = face->soil_attribute<double>("PSD_K_organic");
        //const std::string soil_type = face->soil_attribute<std::string>("soil_type"_s,"soils");
        d.porosity = 0.5; //SoilDataObj->porosity(soil_type);
                          //TODO CRHM uses a fixed value for porosity, rather than searching in the SoilDataObj 
        d.soil_index = face->soil_attribute<double>("soil_index");
        d.snow_grain_diameter = face->soil_attribute<double>("snow_grain_diameter");

        d.Ksaturated_rechr = face->soil_attribute<double>("Ksaturated_rechr");
        d.Ksaturated_lower = face->soil_attribute<double>("Ksaturated_lower");
        d.Ksaturated_ground_water = face->soil_attribute<double>("Ksaturated_ground_water");
        d.Ksaturated_organic = face->soil_attribute<double>("Ksaturated_organic"); 
        // Note: Ksaturated_snow is computed in K_estimate   
        
        d.allow_runoff_from_infiltration = face->soil_attribute<bool>("allow_runoff_from_infiltration"_s);    
    }
    else
    {
        d.soil_storage_max = 0.0;
        d.soil_rechr_max = 0.0;
        d.excess_to_ssr = 0.0;
        d.detention_snow_max = 0.0;
        d.detention_organic_max = 0.0;
        d.depression_max = 0.0;
        d.ground_water_max = 0.0;
        d.local_slope = 0.0;

        d.pore_size_dist = 0.0;
        d.pore_size_dist_organic = 0.0;
        d.soil_index = 0.0;
        d.snow_grain_diameter =0.0;

        d.Ksaturated_rechr = 0.0;
        d.Ksaturated_lower = 0.0;
        d.Ksaturated_ground_water = 0.0;
        d.Ksaturated_organic = 0.0;

    }

};


void soil_module::set_ET_params(mesh_elem& face, soil_module::data& d)
{
    
    d.ground_cover_type = face->soil_attribute<int>("soil_groundcover_ET"_s);
    std::string type = face->soil_attribute<std::string>("soil_type"_s,"soils");

    // I'm about to do something very evil. 
    // CRHMs inconsistent soil typing is responsible
    // likely the result of empirical techniques not 
    // always classifying soils using the same exact categories
    // soil_ET, from CHRM SoiX, only takes sand, loam, or clay
    // physics/Soil.h takes 11 types, I use those types to infer the type here


    // NOTE: This soil_type is not the same as the soil_type used for other parameters like porosity. This is purely used for a computation in the ET calculation in the soil.
    int soil_type = 100;
    bool sandy = compare_substring(type,"sand");
    bool loamy = compare_substring(type,"loam");
    bool clayy = compare_substring(type,"clay");

        
    if (sandy && !loamy && !clayy) // mainly sand 
    {
        soil_type = 1;
    }
    else if (sandy && loamy) // loamy sand or sandy loam
    {
        if (type.substr(0,4) == "loam") //loamy sand
        {
            soil_type = 1;
        }
        else //sandy loam
            soil_type = 2;
    }
    else if (sandy && clayy) // sandy clay
        soil_type = 3;
    else if (clayy && loamy) // loam with clay (clay-y loam)
        soil_type = 2;
    else if (clayy) // mainly clay
        soil_type = 3;
    else if (loamy) // mainly loam
        soil_type = 2;
    
    // detaul is 100, this is handled by soil_ET by treating the soil as mainly organic
    d.soil_type_rechr = soil_type;
    d.soil_type_lower = soil_type; //Same for now 

};

int soil_module::compare_substring(std::string& type, std::string sub)
{
    return type.find(sub) != std::string::npos; // find returns npos
};


void soil_module::initial_soil_conditions(mesh_elem& face, soil_module::data& d)
{
    // TODO actually connect to stuff
    // requires MESHER or at least data for one station
    if (face->has_soil())
    {
        d.soil_storage = face->soil_attribute<double>("soil_storage"_s);
        d.soil_rechr_storage = face->soil_attribute<double>("soil_rechr_storage"_s);
        //d.thaw_fraction_rechr = face->soil_attribute<double>("thaw_fraction_rechr"_s);
        //d.thaw_fraction_lower = face->soil_attribute<double>("thaw_fraction_lower"_s);
        d.detention_snow_init = face->soil_attribute<double>("detention_snow_init"_s);
        d.detention_organic_init = face->soil_attribute<double>("detention_organic_init"_s);
        d.depression_storage = face->soil_attribute<double>("depression_storage"_s);
        d.ground_water_storage = face->soil_attribute<double>("ground_water_storage"_s);
        
    }
    else
    {
        d.soil_storage = 0.0;
        d.soil_rechr_storage = 0.0;
        d.thaw_fraction_rechr = 0.0;
        d.thaw_fraction_lower = 0.0;
        d.detention_snow_init = 0.0;
        d.detention_organic_init = 0.0;
        d.depression_storage = 0.0;
        d.ground_water_storage = 0.0;
        d.pore_size_dist = 1.0; // Set to 1 be default becuase there is a divide by zero with this value. It probably will actually be skipped so not a real worry.
        d.local_slope = 0.0;
    }

};

//bool soil_module::data::is_lake(soil_ET_DTO& DTO)
//{
//    try 
//    {
//        // TODO resolve this bug
//        soil_module::data& d = dynamic_cast<soil_module::data&>(DTO);
//        //bool temp = d.local_module->is_water(*d.my_face);
//        return false;//d.local_module->is_water(*d.my_face);
//    } catch (const std::bad_cast& e) {
//        SPDLOG_DEBUG("bad cast");
//        return false;
//    }
//};

int soil_module::data::get_dt()
{
    if (this->local_module)
        return this->local_module->global_param->dt();
    
    CHM_THROW_EXCEPTION(module_error,"local_module not set in soil_module::data");
    
};

bool soil_module::data::get_new_day()
{
    if (this->first_day)
    {
        this->first_day = false;
        return true;
    }
    else if (this->local_module)
        return this->local_module->is_new_day();
        
    CHM_THROW_EXCEPTION(module_error,"local_module not set in soil_module::data");
};

void soil_module::set_local_module(soil_module::data& d)
{
    d.local_module = this;
    //d.local_module = new soil_module::data::my_module(*this);
};
    // TODO this is just a copy of is_water and this is a bad practice but currently the is_water function is not accessible by the data class. Fix: create a separate object taht module_base inherits that contains these functions. face_info will also inherit these functions.
 //   soil_module::data& d = static_cast<soil_module::data&>(DTO);
 //   bool is = false;

 //   if(d.my_face->has_parameter("landcover"_s))
 //   {
 //       int LC = face->parameter("landcover"_s);
 //       is = global_param->parameters.get<bool>("landcover." + std::to_string(LC) + ".is_water",false);
 //   }
 //   return is
      
//}; 
