#pragma once
#include "submodule_base.hpp"

class I_K_estimate : public submodule_base
{
public:
    virtual ~I_K_estimate() = default;
    virtual void run(void) = 0;
};
