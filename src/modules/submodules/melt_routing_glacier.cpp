#include "melt_routing_glacier.hpp"
#include <stdexcept>
#include <format>

namespace GlacierRouting
{
	LinearReservoir::LinearReservoir(const double k, const double dt)
		: c0(get_c0(k,dt)), c1(get_c1(k,dt)) {
			last_input = 0.0;
			last_output = 0.0;
		};
    static void verifyStability(const double k, const double dt)
    {
        if ( k > 0.0 && k <= dt / 2 )
        {
            std::string err = std::format("In GlacierRouting, numerical stability conditions is not satisfied."
                    "k = {} s, and dt = {} s.",k,dt);
            throw std::runtime_error(err);
        }
    };
	const double LinearReservoir::get_c0(const double k, const size_t dt)
	{
        verifyStability(k,dt);
		return dt / (2 * k + dt);
	};
	const double LinearReservoir::get_c1(const double k, const size_t dt)
	{
		return ( 2*k - dt ) / ( 2*k + dt );
	};

	const double LinearReservoir::step(const double input) {
		auto result = c0 * ( input + last_input )
			+ c1 * last_output;
		
		last_input = input;
		last_output = result;

		return result;
	};

	GlacierReservoir::GlacierReservoir(const Params* _p) : p(_p)
	{
		reservoir_chain = [&]{
			std::vector<LinearReservoir> result;
			result.reserve(p->chain_length);
			for (size_t i = 0; i < p->chain_length-1; ++i)
			{
				// First chain_length - 1 elements are linear reservoirs but with k=0
				// This means that Q_J+1 = I_J+1
				result.emplace_back(0.0,p->seconds_per_step);
			}
			result.emplace_back(p->k,p->seconds_per_step);
			return result;
		}();
		
		old_outputs.resize(p->chain_length - 1, 0.0);
	};

	const double GlacierReservoir::step(const double input) {
		// Process reservoirs in sequence: each reservoir consumes the previous
		// reservoir's output from the last timestep (old_outputs) and produces
		// a new output that becomes the input for the next reservoir for the next timestep.
		// The number of old_outputs is one less than reservoirs because the first
		// reservoir takes the external input, not a previous output.
		
		auto old_output_upchain = old_outputs.begin();
		auto current_reservoir = reservoir_chain.begin();
		
		double output = current_reservoir->step(input);
		++current_reservoir;

		// Note: Extra ++ on reservoir iterator moves it one ahead from the old_outputs iterator
		// This is intentional and syncs up the right element of old_outputs with reservoir_chain
		// as described in the first comment.
		for (; current_reservoir != reservoir_chain.end() && old_output_upchain != old_outputs.end();)
		{
			auto output_temp = current_reservoir->step(*old_output_upchain);
			*old_output_upchain = output;
			output = output_temp;
			++old_output_upchain;
			++current_reservoir;
		};

		return output;
	};

	State::State(const Params* p) : snow_route(p), firn_route(p), ice_route(p) {};
	// Delete copy constructors

};
