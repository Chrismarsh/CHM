/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ALLSKYGENERATOR_H
#define ALLSKYGENERATOR_H

#include <meteoio/dataGenerators/GeneratorAlgorithms.h>
#include <meteoio/dataGenerators/TauCLDGenerator.h>
#include <meteoio/meteoLaws/Sun.h>

#include <map>
#include <utility>

namespace mio {

/**
 * @class AllSkyLWGenerator
 * @brief ILWR all sky parametrization
 * Using air temperature (TA) and relative humidity (RH) and optionnally cloud transmissivity (TAU_CLD),
 * this offers the choice of several all-sky parametrizations:
 *  - OMSTEDT -- from Omstedt, <i>"A coupled one-dimensional sea ice-ocean model applied to a semi-enclosed basin"</i>,
 * Tellus, <b>42 A</b>, 568-582, 1990, DOI:10.1034/j.1600-0870.1990.t01-3-00007.
 *  - KONZELMANN -- from Konzelmann et al., <i>"Parameterization of global and longwave incoming radiation
 * for the Greenland Ice Sheet."</i> Global and Planetary change <b>9.1</b> (1994): 143-164.
 *  - UNSWORTH -- from Unsworth and Monteith, <i>"Long-wave radiation at the ground"</i>,
 * Q. J. R. Meteorolo. Soc., Vol. 101, 1975, pp 13-24 coupled with a clear sky emissivity following (Dilley, 1998).
 *  - CRAWFORD -- from Crawford and Duchon, <i>"An Improved Parametrization for Estimating Effective Atmospheric Emissivity for Use in Calculating Daytime
 * Downwelling Longwave Radiation"</i>, Journal of Applied Meteorology, <b>38</b>, 1999, pp 474-480
 *
 * If no cloud transmissivity is provided in the data, it is calculated from the solar index (ratio of measured iswr to potential iswr, therefore using
 * the current location (lat, lon, altitude) and ISWR to parametrize the cloud cover). This relies on (Kasten and Czeplak, 1980)
 * except for Crawford that provides its own parametrization.
 * The last evaluation of cloud transmissivity is used all along during the times when no ISWR is available if such ratio
 * is not too old (ie. no more than 1 day old).
 * If only RSWR is measured, the measured snow height is used to determine if there is snow on the ground or not.
 * In case of snow, a snow albedo of 0.85 is used while in the abscence of snow, a grass albedo of 0.23 is used
 * in order to compute ISWR from RSWR.
 * Finally, it is recommended to also use a clear sky generator (declared after this one)
 * for the case of no available short wave measurement (by declaring the ClearSky generator \em after AllSky).
 * @code
 * ILWR::generators = allsky_LW
 * ILWR::allsky_lw = Omstedt
 * @endcode
 *
 */
class AllSkyLWGenerator : public GeneratorAlgorithm {
	public:
		AllSkyLWGenerator(const std::vector<std::string>& vecArgs, const std::string& i_algo)
		               : GeneratorAlgorithm(vecArgs, i_algo), model(OMSTEDT), sun(), clf_model(TauCLDGenerator::KASTEN),
		                 last_cloudiness() { parse_args(vecArgs); }
		bool generate(const size_t& param, MeteoData& md);
		bool create(const size_t& param, std::vector<MeteoData>& vecMeteo);
	private:
		void parse_args(const std::vector<std::string>& vecArgs);

		typedef enum PARAMETRIZATION {
			OMSTEDT,
			KONZELMANN,
			UNSWORTH,
			CRAWFORD
		} parametrization;
		parametrization model;

		SunObject sun;
		TauCLDGenerator::clf_parametrization clf_model;
		std::map< std::string, std::pair<double, double> > last_cloudiness; //as < station_hash, <julian_gmt, cloudiness> >
};

} //end namespace mio

#endif
