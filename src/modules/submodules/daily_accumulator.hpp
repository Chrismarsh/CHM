#pragma once
#include "base_step.hpp"

template<typename T>
concept accumulator_data = requires(T& t) {
    // Must have is_new_day() const member function returning bool
    { t.is_new_day() } -> std::same_as<bool>;
    
    // Must have steps_per_day() const member function convertible to double
    { t.steps_per_day() } -> std::convertible_to<double>;
};

template<accumulator_data data>
class daily_accumulator : public base_step<data>
{
	const double* _target_var = nullptr;
	double _accumulator = 0.0;
	double _mean_value = 0.0;
public:
	// overloaded constructors
	// normal construction is the first line
	// second constructor if you want to initialize the mean_value (starting in the middle of a day).
	explicit daily_accumulator(data& _d) : base_step<data>(_d) {};
	explicit daily_accumulator(data& _d, const double value) : base_step<data>(_d), _mean_value(value) {};

	void bind_to_var(double& var)
	{
        if (!_target_var)
        {
		    _target_var = &var;
            return;
        }
        throw std::logic_error("_targer_var already bound"); 
	};

	const double& get_last_mean() const
	{
		return _mean_value;
	};

	void execute() override final
	{
		if (!_target_var) 
		{
			throw std::logic_error("_target_var should be bound before calling execute()");
            return;
		};
        
        if (this->d.is_new_day())
		{
			_mean_value = _accumulator / this->d.steps_per_day();
			_accumulator = 0.0;
		};
		
        _accumulator += *_target_var;

	};
};
