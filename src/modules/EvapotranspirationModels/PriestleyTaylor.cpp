#include "PriestleyTaylor.hpp"

PriestleyTaylor::PriestleyTaylor(const double& alpha_const, const double& Cp) : alpha(alpha_const), heat_capacity_air(Cp)
{

}

PriestleyTaylor::~PriestleyTaylor()
{
    // Do nothing
}

void PriestleyTaylor::CalcEvapT(var_base& basevar, model_output& output)
{
    const PT_vars& var = static_cast<const PT_vars&>(basevar);

    double Q = var.all_wave_net * (1 - Frac_to_ground);

    output.ET = alpha * delta(var.air_temperature) * Q / (delta(var.air_temperature) + gamma(var.P_atm,var.air_temperature,heat_capacity_air) );
}
