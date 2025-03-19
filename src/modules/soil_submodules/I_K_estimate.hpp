#pragma once
class I_K_estimate
{
public:
    virtual ~I_K_estimate() = default;
    virtual void run(void) = 0;
};
