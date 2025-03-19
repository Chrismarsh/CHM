#pragma once

#include <cmath>
#include <iostream>
#include "logger.hpp"
#include <algorithm>

class submodule_base
{
public:
    virtual ~submodule_base() = default;

    virtual void run() = 0;
    
};

