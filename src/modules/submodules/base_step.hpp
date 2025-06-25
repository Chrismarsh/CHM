#pragma once
#include <stdexcept>

#define THROW_NULL_POINTER_EXCEPTION() \
	throw std::runtime_error( \
			std::string("Null pointer at ") + __File__ + ":" + std::to_string(__Line__) \
			)


template<class data>
class base_step
{
public:
	explicit base_step(data& _d) : d(_d) {};
	virtual ~base_step() = default;

	virtual void execute() = 0;
protected:
	data& d;
};


