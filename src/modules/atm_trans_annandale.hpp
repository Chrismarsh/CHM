#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "daily.hpp"

#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

/**
* \addtogroup modules
* @{
* \class atm_trans_annanadale
* \brief Calculates daily average atmospheric transmittance
*
* Calculates daily averaged atmospheric transmittance from Annandale (2002)
* Annandale, J., N. Jovanovic, N. Benadé, and R. Allen. “Software for Missing Data Error Analysis of Penman-Monteith Reference Evapotranspiration.” Irrigation Science 21, no. 2 (2002): 57–67. doi:10.1007/s002710100047.
*
* Depends:
* -  Air temperature "t" [degC]
*
* Provides:
* -
* - Atmospheric transmittance "atm_trans" [-]
*/
class atm_trans_annandale : public module_base
{
public:
    atm_trans_annandale(std::string ID);
    ~atm_trans_annandale();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);


};

/**
@}
*/