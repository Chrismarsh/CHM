#include "EvapotranspirationModels/evapbase.hpp"

double evapT_base::delta(const double& t) // Slope of sat vap p vs t, kPa/DEGREE_CELSIUS
{
  if (t > 0.0)
    return(2504.0*exp(17.27 * t/(t+237.3)) / pow(t+237.3,2));
  else
    return(3549.0*exp( 21.88 * t/(t+265.5)) / pow(t+265.5,2));
}

double evapT_base::lambda(const double& t) // Latent heat of vaporization (J/kg)
{
    // Equation 7-8 Dingman (2002) Second Edition
    return (2.501 - 0.002361 * t) * 1e6; // original is MegaJoules/kg, 1e6 returns it to joules/kg
}

double evapT_base::gamma(const double& Pa, const double& t, const double& c_a) // Psychrometric constant (kPa/DEGREE_CELSIUS)
{
   // Equation 7-13 Dingman Second Edition 2002
   return c_a * Pa / (0.622 * lambda(t)); // lambda (J/kg)
}
