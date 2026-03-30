#include "katabatic_routing_glacier.hpp"
#include "Atmosphere.h"
#include "Glacier.hpp"
#include "PhysConst.h"
#include "melt_routing_glacier.hpp"

// data and view constructors
REGISTER_MODULE_CPP(katabatic_routing_glacier);

katabatic_routing_glacier::katabatic_routing_glacier(config_file cfg)
	: module_base("katabatic_routing_glacier", parallel::data, cfg)
{
	depends("Pa");
	depends("t");
	depends("rh");
    depends("vapour_pressure_surface");
	// vapour pressure surface??
	depends("t_lapse_rate");	
	depends("snowmelt_int");
	depends("swe");
	// TODO uncomment, only commented for CRHM coupling
    //depends("iswr_net");
	//depends("iswr_subcanopy");
	//depends("ilwr_subcanopy");
    // TODO remove, only for CRHM coupling 
    // also should check if this exists in CHM
    depends("T_rain");

	// All submodule outputs
	provides("glacier_water_equivalent");
	provides("total_depth");
	provides("firnmelt");
	provides("icemelt");
	provides("latent_heat");
	provides("sensible_heat");
	provides("snowmelt_delayed");
	provides("firnmelt_delayed");
	provides("icemelt_delayed");
	provides("total_delayed");

};
katabatic_routing_glacier::~katabatic_routing_glacier() {};
// module constructor
// module init
// module run

void katabatic_routing_glacier::init(mesh& domain)
{
	for (size_t i = 0; i < domain->size_local_faces(); ++i)
	{
		auto face = domain->face(i);
		auto& p_glacier = glacier.get_params();
		auto& p_routing = routing.get_params();
		auto& d = face->make_module_data<data>(ID,face,global_param,&cfg,&p_glacier,&p_routing);
		d.ice_emissivity = cfg.get<double>("ice_emissivity");
		d.firn_emissivity = cfg.get<double>("firn_emissivity");
	};

	{
		using namespace Glacier;
		auto& p = glacier.get_params();

		auto densify_version = cfg.get<std::string>("densify_version","Linear");
		if (densify_version == "Linear")
			p.densify_version = DensifyVersion::Linear;
		else if (densify_version == "HerronLangway")
			p.densify_version = DensifyVersion::HerronLangway;
		else 
		{
			p.densify_version = DensifyVersion::HerronLangway;
			SPDLOG_DEBUG("Unknown densify version specified. Using HerronLangway as the default");
		}

		switch (p.densify_version)
		{
			case DensifyVersion::Linear:
				p.small_increment = cfg.get<double>("small_increment",25.0);
				p.big_increment = cfg.get<double>("big_increment",50.0);
				break;
			case DensifyVersion::HerronLangway:
				break;
		};

        p.critical_density = cfg.get<double>("critical_density",550.0);
        p.firn_to_ice_density = cfg.get<double>("firn_to_ice_density",830.0);
		
		p.thermal_factor = cfg.get<double>("thermal_factor",0.95);
		p.seconds_per_step = global_param->dt();
	}
	{
		auto& p = katabatic.get_params();
		p.prandtl = cfg.get<double>("Pr",5.0);
		p.k = cfg.get<double>("k",4e-4);
		p.k2 = cfg.get<double>("k2",1.0);
		p.seconds_per_step = global_param->dt();
	}
	{
		auto& p = routing.get_params();
		p.chain_length = cfg.get<size_t>("routing_lag",10u);
		p.seconds_per_step = global_param->dt();
	}
};

double katabatic_routing_glacier::rain_sun_energy(const mesh_elem& face)
{
    using namespace PhysConst;
    constexpr auto M_PER_MM = 1 / 1000.0;
    auto Qsun = (*face)["Qnsn_Var"_s];
    auto Qrain = Cw() * water_reference_density() * (*face)["rainfall_int"] * M_PER_MM
        * ((*face)["T_rain"_s]  - 0.0) / global_param->dt();

    return Qsun + Qrain;
};

double katabatic_routing_glacier::rain_sun_energy(const mesh_elem& face,data& d)
{
	//auto iswr = (*face)["iswr_subcanopy"_s];
	//auto ilwr = (*face)["ilwr_subcanopy"_s];
	//const auto& firn = d.glacier_state.firn;
	//const auto& ice = d.glacier_state.ice;
	//auto emissivity = 0.0;
	//if (d.swe().value > 0.0)
	//	return 0.0;
	//else if (firn.water_equivalent().value)
	//	emissivity = d.firn_emissivity;
	//else if (ice.water_equivalent().value)
	//	emissivity = d.ice_emissivity;
	//else
	//	return 0.0;
	//
	//auto olwr = PhysConst::sbc() * emissivity
	//	* std::pow(d.air_temperature().value,4.0);
	//auto oswr = (*face)["glacier_albedo"_s] * iswr;
	//auto Qsun = (ilwr - olwr) + (iswr - oswr);
	//auto Qrain = (*face)["Qrain"_s];
    auto Qsun = 0.0;
    auto Qrain = 0.0;
	return Qsun + Qrain;
};

void katabatic_routing_glacier::run(mesh_elem& face)
{
	auto& d = face->get_module_data<data>(ID);

	{
		// run every step, make sure the outputs are accumulating
		// somewhere
		katabatic_view data_k(d,face);
		katabatic.execute(data_k);
		const auto& cache = d.get_cache();
        // TODO this version of rain_run_energy is for CRHM coupling, do not include
		d.total_energy += rain_sun_energy(face) + cache->latent_heat +
			cache->sensible_heat;
		// only run once per day
		if (is_new_day())
		{
			glacier_view data_g(d,face);
			glacier.execute(data_g);

			// TODO reset latent/sensible heat here?
			// could get_cache
			// accumulate total_energy on public d member
			// get d.melt_energy() return this value
			// reset total_energy here
			d.total_energy = 0.0;
		}

		// run every step
		routing_view data_r(d,face);
		routing.execute(data_r);
		//Implicitly output to face via deconstructors
		//of *_view objects
	}

	// *_view objects out of scope, safe to reset
	d.reset_cache();

};

bool katabatic_routing_glacier::is_new_day()
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

katabatic_routing_glacier::katabatic_view::~katabatic_view()
{	
	auto& cache = d.get_cache();

	(*face)["latent_heat"_s] = cache->latent_heat;	
	(*face)["sensible_heat"_s] = cache->sensible_heat;
};

katabatic_routing_glacier::routing_view::~routing_view()
{
	auto& cache = d.get_cache();

	(*face)["snowmelt_delayed"_s] = cache->snowmelt_delayed;
	(*face)["firnmelt_delayed"_s] = cache->firnmelt_delayed;
	(*face)["icemelt_delayed"_s] = cache->icemelt_delayed;
	(*face)["total_delayed"_s] = cache->snowmelt_delayed + cache->firnmelt_delayed + cache->icemelt_delayed;
};

katabatic_routing_glacier::glacier_view::~glacier_view()
{
	auto& cache = d.get_cache();
	(*face)["glacier_water_equivalent"] = d.glacier_state.total_water_equiv().value;
	(*face)["total_depth"] = d.glacier_state.total_depth().value;
	(*face)["firnmelt"] = cache->firnmelt;
	(*face)["icemelt"] = cache->icemelt;
};

katabatic_routing_glacier::data::data(const mesh_elem& face, std::shared_ptr<global> global_param, const config_file* cfg,
		const Glacier::Params* p_g, const GlacierRouting::Params* p_r) 
	: data_base(face,global_param,cfg), routing_state(p_r), glacier_state(p_g) {};

Units::Kelvin katabatic_routing_glacier::data::glacier_temperature()
{
    if (!cfg)
    {
        std::string err = std::format("config file is null in katabatic_routing_glacier::data::glacier_temperature()");
        CHM_THROW_EXCEPTION(module_error,err);
    };
	static const double T = 273.15;//cfg_.get<double>("glacier_temp"_s,273.15);

	return Units::Kelvin{T};
};
Units::Kelvin katabatic_routing_glacier::data::air_temperature()
{
	update_value( [this]() -> auto& { return cache_->air_temperature; },
			[this]() { return (*face)["t"_s]; } );

	return Units::Kelvin{cache_->air_temperature + 273.15};
};
Units::Pa katabatic_routing_glacier::data::air_pressure()
{
	update_value( [this]() -> auto& { return cache_->air_pressure; },
			[this]() { return (*face)["Pa"_s]; } );

	return Units::Pa{cache_->air_pressure};
};
Units::Pa katabatic_routing_glacier::data::vapour_pressure()
{
	update_value( [this]() -> auto& { return cache_->vapour_pressure; },
			[this]() { 
                auto rh = (*face)["rh"_s];
                auto output = rh * Atmosphere::saturatedVapourPressure(air_temperature().value);
                return output; }
                );

	return Units::Pa{cache_->vapour_pressure};
};
Units::Pa katabatic_routing_glacier::data::vapour_pressure_surface()
{
	update_value( [this]() -> auto& { return cache_->vapour_pressure_surface; },
			[this]() { return (*face)["vapour_pressure_surface"_s]; } );

	return Units::Pa{cache_->vapour_pressure_surface};
};
Units::LapseRateSI katabatic_routing_glacier::data::lapse_rate()
{
	update_value( [this]() -> auto& { return cache_->lapse_rate; },
			[this]() { return (*face)["t_lapse_rate"_s]; } );

	return Units::LapseRateSI{cache_->lapse_rate};
};
void katabatic_routing_glacier::data::latent_heat(const double v)
{
	set_output( [this]() -> auto& { return cache_->latent_heat; },
			v);
};
void katabatic_routing_glacier::data::sensible_heat(const double v)
{
	set_output( [this]() -> auto& { return cache_->sensible_heat; },
			v);
};

const Units::Milimetres katabatic_routing_glacier::data::snowmelt()
{
	update_value( [this]() -> auto& { return cache_->snowmelt; },
			[this]() { return (*face)["snowmelt_int"_s]; } );

	return Units::Milimetres{cache_->snowmelt};
};
const Units::Milimetres katabatic_routing_glacier::data::firnmelt()
{
	if (!cache_ || std::isnan(cache_->firnmelt))
		return Units::Milimetres{0.0};
	
	return Units::Milimetres{cache_->firnmelt};
};
const Units::Milimetres katabatic_routing_glacier::data::icemelt()
{
	if (!cache_ || std::isnan(cache_->icemelt))
		return Units::Milimetres{0.0};
	
	return Units::Milimetres{cache_->icemelt};
};
void katabatic_routing_glacier::data::snowmelt_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->snowmelt_delayed; },
			v);
};
void katabatic_routing_glacier::data::firnmelt_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->firnmelt_delayed; },
			v);

};
void katabatic_routing_glacier::data::icemelt_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->icemelt_delayed; },
			v);
};
void katabatic_routing_glacier::data::total_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->total_delayed; },
			v);
};

const Units::Milimetres katabatic_routing_glacier::data::swe()
{
	update_value( [this]() -> auto& { return cache_->swe; },
			[this]() { return (*face)["swe"_s]; } );

	return Units::Milimetres{cache_->swe};
};
const Units::Watts_per_m2 katabatic_routing_glacier::data::melt_energy()
{
	auto iswr = (*face)["iswr_subcanopy"_s];
	auto ilwr = (*face)["ilwr_subcanopy"_s];
	const auto& firn = glacier_state.firn;
	const auto& ice = glacier_state.ice;
	auto emissivity = 0.0;
	if (swe().value > 0.0)
		return Units::Watts_per_m2{0.0};
	else if (firn.water_equivalent().value)
		emissivity = firn_emissivity;
	else if (ice.water_equivalent().value)
		emissivity = ice_emissivity;
	else
		return Units::Watts_per_m2{0.0};
	

	auto olwr = PhysConst::sbc() * emissivity
		* std::pow(air_temperature().value,4.0);
	auto oswr = (*face)["glacier_albedo"_s] * iswr;
	auto Qsun = (ilwr - olwr) + (iswr - oswr);
	auto Qrain = (*face)["Qrain"_s];
	// Means that katabatic_melt_energy didn't run
	if (!cache_)
		return Units::Watts_per_m2{Qsun + Qrain};

	auto Q_sensible = 0.0;	
	auto Q_latent = 0.0;

	if ( !std::isnan(cache_->latent_heat) )
		Q_latent = cache_->latent_heat;

	if ( !std::isnan(cache_->sensible_heat) )
		Q_sensible = cache_->sensible_heat;

	return Units::Watts_per_m2{Qsun + Qrain + Q_sensible + Q_latent};
};
bool katabatic_routing_glacier::data::update_now()
{
	static const size_t day = cfg->get<size_t>("change_day"_s,300u);
	return day == global_param->day(); 
};

void katabatic_routing_glacier::data::glacier_water_equivalent(const double out)
{
	set_output( [this]() -> auto& { return cache_->glacier_water_equivalent; },
			out);
};
void katabatic_routing_glacier::data::total_depth(const double out)
{
	set_output( [this]() -> auto& { return cache_->total_depth; },
			out);
};
void katabatic_routing_glacier::data::firn_melt(const double out)
{
	set_output( [this]() -> auto& { return cache_->firnmelt; },
			out);
};
void katabatic_routing_glacier::data::ice_melt(const double out)
{
	set_output( [this]() -> auto& { return cache_->icemelt; },
			out);
};


