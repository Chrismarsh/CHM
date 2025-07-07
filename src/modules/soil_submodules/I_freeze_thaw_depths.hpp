#pragma once
#include "submodule_base.hpp"

class I_freeze_thaw_depths : public submodule_base
{
public:
    virtual ~I_freeze_thaw_depths() = default;
    virtual void run(void) = 0;

    
    double get_freeze_depth() { return freeze_front_depth;};

    double get_thaw_depth() { return thaw_front_depth; };

    double get_first_front_depth() { return first_front_depth; };

protected:
    double thaw_front_depth = 0.0;
    double freeze_front_depth = 0.0;
    double first_front_depth = 0.0;

};
