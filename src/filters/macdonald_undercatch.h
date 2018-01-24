#pragma once

#include "filter_base.h"
#include <math.h>


/**
 * \addtogroup filters
 * @{
 * \class macdonald_undercatch
 * \brief Computes undercatch correction
 *
 * Undercatch correction for a Alter shielded Geonor and tipping bucket via Macdonald, et al. 2007
 *
 * Depends:
 * - p [mm]
 * - u [m/s]
 *
 * References:
 * - Macdonald, J., & Pomeroy, J. (2007). Gauge Undercatch of Two Common Snowfall Gauges in a Prairie Environment. Proceedings of the 64th Eastern Snow Conference, St. John‘s, Canada., 119–126.
 * */
class macdonald_undercatch : public filter_base
{
private:
    std::string var;
public:
    macdonald_undercatch();
    ~macdonald_undercatch();
    void init(boost::shared_ptr<station> station);
    void process(boost::shared_ptr<station> station);
};