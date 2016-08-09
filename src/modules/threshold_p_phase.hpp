#pragma once

#include <boost/shared_ptr.hpp>

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"


#include <string>

/**
* \addtogroup modules
* @{
* \class threshold_p_phase : public module_base
{
* \brief Precip phase from air temp
*
* Calculates phase from air temperature. Defaults to 2degC, looks for the config parameter "threshold_temperature".
* Depends:
* - Air temperature (t)
*
* Provides:
* - Atmospheric transmittance "cloud_frac" [-]
*/
class threshold_p_phase : public module_base
{
public:
    threshold_p_phase(config_file cfg);

    ~threshold_p_phase();

    virtual void run(mesh_elem &face, boost::shared_ptr <global> global_param);

private:

    // the air temperature threshold above which the precip phase is liquid
    double t_thresh;

};