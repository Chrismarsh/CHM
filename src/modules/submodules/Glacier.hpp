#include "base_step.hpp"
#include <concepts>

namespace Glacier
{
    struct State;
    class katabatic_melt_energy
    {
    public:
        double get(State&);
    };

    template<class T>
    concept GlacierData = requires(T& t)
    {
        { t.snowmelt() } -> std::floating_point;
        { t.rainfall() } -> std::floating_point;
        { t.swe() } -> std::floating_point;
        { t.albedo() } -> std::floating_point;
				
      	{ t.swe(std::declval<double>()) } -> std::same_as<void>;
    		{ t.total_water_equivalent(std::declval<double>()) } -> std::same_as<void>;
        { t.total_depth(std::declval<double>()) } -> std::same_as<void>;
        { t.total_melt(std::declval<double>()) } -> std::same_as<void>;
    };

    template<GlacierData data>
    class Model : public base_step<Model<data>,data>
    {
    public:
        void execute_impl(data& d);
    private:
        katabatic_melt_energy melt_energy;
        //ice_manager ice;
        //firn_manager firn;
        //melt_router routing;
        //glacier_depths depths;
    };
};
