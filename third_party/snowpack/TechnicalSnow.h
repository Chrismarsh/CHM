/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
 */
/*  This file is part of Snowpack.
    Snowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TECHNICALSNOW_H
#define TECHNICALSNOW_H

#include <snowpack/DataClasses.h>

/**
 * @brief Implementation of snow grooming
 * @details
 * 
 * This module relies on the following configuration keys, all in the [TechSnow] section:
 *  - SNOW_GROOMING: if set to true, enables this module (default: false);
 *  - GROOMING_WEEK_START: ISO week number when to start grooming (default: 40)
 *  - GROOMING_WEEK_END: ISO week number when to stop grooming (default: 17)
 *  - GROOMING_HOUR: at what time should grooming be performed?  (default: 21 hour)
 *  - GROOMING_DEPTH_START: how much snow must be on the ground to start grooming [m]  (default: 0.4);
 *  - GROOMING_DEPTH_IMPACT: maximum depth of snow impacted by grooming [m]  (default: 0.4);
 * 
 * @author Mathias Bavay, Pirmin Ebner and others
 * @ingroup postprocessing
 */

class TechSnow {
	public:
		TechSnow(const SnowpackConfig& cfg);

		bool prepare(const mio::Date& current_date) const;
		void preparation(SnowStation& Xdata) const;
		static void productionPpt(const CurrentMeteo& Mdata, const double& cumu_precip, double &Tw, double &rho_hn, double &delta_cH, double &theta_w);
		
	private:
		double grooming_week_start, grooming_week_end;
		double grooming_hour;
		double min_depth, max_depth; //minimum depth of snow for grooming, maximum depth affected by grooming
};
#endif
