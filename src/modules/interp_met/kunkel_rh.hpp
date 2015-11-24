#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"

/**
* \addtogroup modules
* @{
* \class kunkel_rh
* \brief RH interpolation
*
* Interpolates RH with elevation. Uses eqn 15 in Kunkel to directly interplate RH.
*
* Depends from met:
* - Relative Humidity 'rh' [%]
*
* Provides:
* - Relative Humidity 'rh' [%]
 *
 * References:
 * - Kunkel, K. E. (1989). Simple procedures for extrapolation of humidity variables in the mountainous western United States. Journal of Climate, 2(7), 656–669. Retrieved from http://ams.allenpress.com/perlserv/?request=get-abstract&amp;doi=10.1175/1520-0442(1989)002<0656:SPFEOH>2.0.CO;2
*/
class kunkel_rh : public module_base
{
public:
    kunkel_rh();

    ~kunkel_rh();

    virtual void run(mesh_elem &elem, boost::shared_ptr <global> global_param);


};



/**
@}
*/