/***********************************************************************************/
/*  Copyright 2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>
#include <meteoio/meteoStats/libinterpol2D.h>
#include <meteoio/MathOptim.h>

#include <meteoio/spatialInterpolations/ALSScaleAlgorithm.h>
#include <meteoio/spatialInterpolations/AvgAlgorithm.h>
#include <meteoio/spatialInterpolations/AvgLapseAlgorithm.h>
#include <meteoio/spatialInterpolations/ConstAlgorithm.h>
#include <meteoio/spatialInterpolations/IDWAlgorithm.h>
#include <meteoio/spatialInterpolations/IDWLapseAlgorithm.h>
#include <meteoio/spatialInterpolations/IDWLapseLocalAlgorithm.h>
#include <meteoio/spatialInterpolations/ILWREpsAlgorithm.h>
#include <meteoio/spatialInterpolations/ListonWindAlgorithm.h>
#include <meteoio/spatialInterpolations/NoneAlgorithm.h>
#include <meteoio/spatialInterpolations/ODKrigAlgorithm.h>
#include <meteoio/spatialInterpolations/ODKrigLapseAlgorithm.h>
#include <meteoio/spatialInterpolations/PPhaseAlgorithm.h>
#include <meteoio/spatialInterpolations/RHListonAlgorithm.h>
#include <meteoio/spatialInterpolations/RyanWindAlgorithm.h>
#include <meteoio/spatialInterpolations/SnowPsumAlgorithm.h>
#include <meteoio/spatialInterpolations/StdPressAlgorithm.h>
#include <meteoio/spatialInterpolations/SwRadAlgorithm.h>
#include <meteoio/spatialInterpolations/template.h>
#include <meteoio/spatialInterpolations/UserAlgorithm.h>
#include <meteoio/spatialInterpolations/WinstralAlgorithm.h>
#include <meteoio/spatialInterpolations/WinstralListonAlgorithm.h>

#include <sstream>
#include <vector>
#include <algorithm>

namespace mio {

/**
 * @page interpol2d Spatial interpolations
 * Using the vectors of MeteoData and StationData as filled by the IOInterface::readMeteoData call
 * as well as a grid of elevations (DEM, stored as a DEMObject), it is possible to get spatially
 * interpolated parameters.
 *
 * First, an interpolation method has to be selected for each variable which needs interpolation. Then the class computes
 * the interpolation for each 2D grid point, combining the inputs provided by the available data sources.
 * Any parameter of MeteoData can be interpolated, using the names given by \ref meteoparam. One has to keep
 * in mind that the interpolations are time-independent: each interpolation is done at a given time step and no
 * memory of (eventual) previous time steps is kept. This means that all parameters and variables that are
 * automatically calculated get recalculated anew for each time step.
 *
 * @section interpol2D_section Spatial interpolations section
 * Practically, the user
 * has to specify in his configuration file (typically io.ini), for each parameter to be interpolated, which
 * spatial interpolations algorithms should be considered, in the [Interpolations2D] section. This is provided as a space separated list of keywords
 * (one per interpolation algorithm). Please notice that some algorithms may require extra arguments.
 * Then, each algorithm will be evaluated (through the use of its rating method) and receive a grade (that might
 * depend on the number of available data, the quality of the data, etc). The algorithm that receives the higher
 * score within the user list, will be used for interpolating the selected variable at the given timestep. This means that at another
 * timestep, the same parameter might get interpolated by a different algorithm.
 * An example of such section is given below:
 * @code
 * [Interpolations2D]
 * TA::algorithms = IDW_LAPSE AVG_LAPSE
 * TA::avg_lapse = -0.008
 *
 * RH::algorithms = RH IDW_LAPSE AVG_LAPSE AVG
 *
 * PSUM::algorithms = PSUM_SNOW IDW_LAPSE AVG_LAPSE AVG CST
 * PSUM::psum_snow = avg_lapse
 * PSUM::avg_lapse = 0.0005 frac
 * PSUM::cst        = 0
 *
 * VW::algorithms = IDW_LAPSE AVG_LAPSE
 *
 * P::algorithms = STD_PRESS
 * @endcode
 *
 * @section interpol2D_keywords Available algorithms
 * The keywords defining the algorithms are the following:
 * - NONE: returns a nodata filled grid (see NoneAlgorithm)
 * - STD_PRESS: standard atmospheric pressure as a function of the elevation of each cell (see StandardPressureAlgorithm)
 * - CST: constant value in each cell (see ConstAlgorithm)
 * - AVG: average of the measurements in each cell (see AvgAlgorithm)
 * - AVG_LAPSE: constant value reprojected to the elevation of the cell (see AvgLapseRateAlgorithm)
 * - IDW: Inverse Distance Weighting averaging (see IDWAlgorithm)
 * - IDW_LAPSE: Inverse Distance Weighting averaging with reprojection to the elevation of the cell (see IDWLapseAlgorithm)
 * - LIDW_LAPSE: IDW_LAPSE restricted to a local scale (n neighbor stations, see LocalIDWLapseAlgorithm)
 * - LISTON_RH: the dew point temperatures are interpolated using IDW_LAPSE, then reconverted locally to relative humidity (see RHListonAlgorithm)
 * - ILWR_EPS: the incoming long wave radiation is converted to emissivity and then interpolated (see ILWREpsAlgorithm)
 * - SWRAD: The atmospheric attenuation and splitting coefficients are evaluated and used to compute the short wave radiation with topographic shading (see SWRadInterpolation)
 * - LISTON_WIND: the wind field (VW and DW) is interpolated using IDW_LAPSE and then altered depending on the local curvature and slope (taken from the DEM, see ListonWindAlgorithm)
 * - RYAN: the wind direction is interpolated using IDW and then altered depending on the local slope (see RyanAlgorithm)
 * - WINSTRAL: the solid precipitation is redistributed by wind according to (Winstral, 2002) (see WinstralAlgorithm)
 * - PSUM_SNOW: precipitation interpolation according to (Magnusson, 2011) (see SnowPSUMInterpolation)
 * - PPHASE: precipitation phase parametrization performed at each cell (see PPHASEInterpolation)
 * - ODKRIG: ordinary kriging (see OrdinaryKrigingAlgorithm)
 * - ODKRIG_LAPSE: ordinary kriging with lapse rate (see LapseOrdinaryKrigingAlgorithm)
 * - USER: user provided grids to be read from disk (if available, see USERInterpolation)
 * - ALS_SCALING: scaling from Airborn Laser Scan data (see ALS_Interpolation)
 *
 * @section interpol2D_trends Altitudinal trends
 * Several algorithms use elevation trends, all of them relying on the same principles: the lapse rates are recomputed at each time steps
 * (see section \ref interpol2D_lapse), all stations' data are detrended with this lapse rate, the residuals are spatially interpolated
 * with the algorithm as configured by the user and finally, the values at each cell are retrended (ie the lapse rates are re-applied
 * using the cell's elevation).
 *
 * @subsection interpol2D_lapse Lapse rates
 * The altitudinal trends are currently modelled as a linear relation. The slope of this linear relation can
 * sometimes be provided by the end user (through his io.ini configuration file), otherwise it is computed from the data.
 * In order to bring slightly more robustness, if the correlation between the input data and the computed linear regression
 * is not good enought (below 0.7, as defined in Interpol2D::LinRegression), the same regression will get re-calculated
 * with one point less (cycling throught all the points). The best result (ie: highest correlation coefficient) will be
 * kept. If the final correlation coefficient is less than 0.7, a warning is displayed.
 *
 * @section interpol2D_dev_use Developer usage
 * From the developer's point of view, all that has to be done is instantiate an IOManager object and call its
 * IOManager::interpolate method.
 * @code
 * 	Config cfg("io.ini");
 * 	IOManager io(cfg);
 *
 * 	//reading the dem (necessary for several spatial interpolations algoritms)
 * 	DEMObject dem;
 * 	io.readDEM(dem);
 *
 * 	//performing spatial interpolations
 * 	Grid2DObject param;
 * 	io.interpolate(date, dem, MeteoData::TA, param);
 *
 * @endcode
 *
 * @section interpol2D_biblio Bibliography
 * The interpolation algorithms have been inspired by the following papers:
 * - <i>"A Meteorological Distribution System for High-Resolution Terrestrial Modeling (MicroMet)"</i>, Liston and Elder, Journal of Hydrometeorology <b>7</b> (2006), 217-234.
 * - <i>"Simulating wind ﬁelds and snow redistribution using terrain-based parameters to model snow accumulation and melt over a semi-arid mountain catchment"</i>, Adam Winstral and Danny Marks, Hydrological Processes <b>16</b> (2002), 3585– 3603. DOI: 10.1002/hyp.1238
 * - <i>"Quantitative evaluation of different hydrological modelling approaches in a partly glacierized Swiss watershed"</i>, Jan Magnusson, Daniel Farinotti, Tobias Jonas and Mathias Bavay, Hydrological Processes, 2010, under review.
 * - <i>"Modelling runoff from highly glacierized alpine catchments in a changing climate"</i>, Matthias Huss, Daniel Farinotti, Andreas Bauder and Martin Funk, Hydrological Processes, <b>22</b>, 3888-3902, 2008.
 * - <i>"Geostatistics for Natural Resources Evaluation"</i>, Pierre Goovaerts, Oxford University Press, Applied Geostatistics Series, 1997, 483 p., ISBN 0-19-511538-4
 * - <i>"Statistics for spatial data"</i>, Noel A. C. Cressie, John Wiley & Sons, revised edition, 1993, 900 p.
 *
 * @author Mathias Bavay
 * @date   2010-04-12
 */

InterpolationAlgorithm* AlgorithmFactory::getAlgorithm(const std::string& i_algoname,
                                                       Meteo2DInterpolator& i_mi,
                                                       const std::vector<std::string>& i_vecArgs, TimeSeriesManager& tsm, GridsManager& gdm)
{
	const std::string algoname( IOUtils::strToUpper(i_algoname) );

	if (algoname == "NONE") {// return a nodata grid
		return new NoneAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "STD_PRESS") {// standard air pressure interpolation
		return new StandardPressureAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "CST") {// constant fill
		return new ConstAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "AVG") {// average fill
		return new AvgAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "AVG_LAPSE") {// average fill with an elevation lapse rate
		return new AvgLapseRateAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "IDW") {// Inverse Distance Weighting fill
		return new IDWAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "IDW_LAPSE") {// Inverse Distance Weighting with an elevation lapse rate fill
		return new IDWLapseAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "LIDW_LAPSE") {// Inverse Distance Weighting with an elevation lapse rate fill, restricted to a local scale
		return new LocalIDWLapseAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "LISTON_RH") {// relative humidity interpolation
		return new RHListonAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "ILWR_EPS") {// long wave radiation interpolation
		return new ILWREpsAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "LISTON_WIND") {// wind velocity interpolation (using a heuristic terrain effect)
		return new ListonWindAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "RYAN") {// RYAN wind direction
		return new RyanAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "WINSTRAL") {// Winstral wind exposure factor
		return new WinstralAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "WINSTRAL++") {// Winstral/Liston wind exposure factor
		return new WinstralListonAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "ODKRIG") {// ordinary kriging
		return new OrdinaryKrigingAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "ODKRIG_LAPSE") {// ordinary kriging with lapse rate
		return new LapseOrdinaryKrigingAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "USER") {// read user provided grid
		return new USERInterpolation(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "ALS_SCALING") {// scale from ALS grid
		return new ALS_Interpolation(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "PPHASE") {// precipitation phase parametrization
		return new PPHASEInterpolation(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "PSUM_SNOW") {// precipitation interpolation according to (Magnusson, 2010)
		return new SnowPSUMInterpolation(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "SWRAD") {// terrain shadding interpolation
		return new SWRadInterpolation(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "RH") {//HACK deprecated
		throw IOException("The 'RH' interpolation algorithm has been renamed as 'LISTON_RH' for consistency" , AT);
	} else {
		throw IOException("The interpolation algorithm '"+algoname+"' is not implemented" , AT);
	}
}

size_t InterpolationAlgorithm::getData(const Date& i_date, const MeteoData::Parameters& i_param,
                                       std::vector<double>& o_vecData)
{
	tsmanager.getMeteoData(i_date, vecMeteo);
	o_vecData.clear();
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		const double& val = vecMeteo[ii](i_param);
		if (val != IOUtils::nodata) {
			o_vecData.push_back( val );
		}
	}

	return o_vecData.size();
}

size_t InterpolationAlgorithm::getData(const Date& i_date, const MeteoData::Parameters& i_param,
                                       std::vector<double>& o_vecData, std::vector<StationData>& o_vecMeta)
{
	tsmanager.getMeteoData(i_date, vecMeteo);
	o_vecData.clear();
	o_vecMeta.clear();
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		const double& val = vecMeteo[ii](i_param);
		if (val != IOUtils::nodata){
			o_vecData.push_back( val );
			o_vecMeta.push_back( vecMeteo[ii].meta );
		}
	}

	return o_vecData.size();
}

size_t InterpolationAlgorithm::getStationAltitudes(const std::vector<StationData>& i_vecMeta,
                                                   std::vector<double>& o_vecData)
{
	o_vecData.clear();
	for (size_t ii=0; ii<i_vecMeta.size(); ii++){
		const double& alt = i_vecMeta[ii].position.getAltitude();
		if (alt != IOUtils::nodata) {
			o_vecData.push_back( alt );
		}
	}

	return o_vecData.size();
}

/**
 * @brief Return an information string about the interpolation process
 * @return string containing some information (algorithm used, number of stations)
*/
std::string InterpolationAlgorithm::getInfo() const
{
	std::ostringstream os;
	os << algo << ", " << nrOfMeasurments << " station";
	if (nrOfMeasurments!=1) os << "s"; //add plural mark

	const std::string tmp( info.str() );
	if (!tmp.empty()) os << ", " << tmp;

	return os.str();
}

/**
 * @brief Read the interpolation arguments and compute the trend accordingly
 *
 * @param vecAltitudes altitudes sorted similarly as the data in vecDat
 * @param vecDat data for the interpolated parameter
 * @param trend object containing the fitted trend to be used for detrending/retrending
*/
void InterpolationAlgorithm::getTrend(const std::vector<double>& vecAltitudes, const std::vector<double>& vecDat, Fit1D &trend) const
{
	bool status;
	if (vecArgs.empty()) {
		trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
		status = trend.fit();
	} else if (vecArgs.size() == 1) {
		double lapse_rate;
		IOUtils::convertString(lapse_rate, vecArgs.front());
		trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
		trend.setLapseRate(lapse_rate);
		status = trend.fit();
	} else if (vecArgs.size() == 2) {
		const std::string extraArg( vecArgs[1]);
		if (extraArg=="soft") { //soft
			trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
			status = trend.fit();
			if (!status) {
				double lapse_rate;
				IOUtils::convertString(lapse_rate, vecArgs[0]);
				trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
				trend.setLapseRate(lapse_rate);
				status = trend.fit();
			}
		} else if (extraArg=="frac") {
			double lapse_rate;
			IOUtils::convertString(lapse_rate, vecArgs[0]);
			trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
			const double avgData = Interpol1D::arithmeticMean(vecDat);
			trend.setLapseRate(lapse_rate*avgData);
			status = trend.fit();
			if (lapse_rate*avgData==0.) trend.setInfo(trend.getInfo() + " (null average input for frac lapse rate)");
		} else {
			throw InvalidArgumentException("Unknown argument \""+extraArg+"\" supplied for the "+algo+" algorithm", AT);
		}
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
	}

	if (!status)
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param) + ": " + trend.getInfo(), AT);
}


void InterpolationAlgorithm::detrend(const Fit1D& trend, const std::vector<double>& vecAltitudes, std::vector<double> &vecDat, const double& min_alt, const double& max_alt)
{
	if (vecDat.size() != vecAltitudes.size()) {
		std::ostringstream ss;
		ss << "Number of station data (" << vecDat.size() << ") and number of elevations (" << vecAltitudes.size() << ") don't match!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	for (size_t ii=0; ii<vecAltitudes.size(); ii++) {
		const double altitude = std::min( std::max(vecAltitudes[ii], min_alt), max_alt );
		double &val = vecDat[ii];
		if (val!=IOUtils::nodata)
			val -= trend( altitude );
	}
}

void InterpolationAlgorithm::retrend(const DEMObject& dem, const Fit1D& trend, Grid2DObject &grid, const double& min_alt, const double& max_alt)
{
	const size_t nxy = grid.size();
	if (dem.size() != nxy) {
		std::ostringstream ss;
		ss << "Dem size (" << dem.grid2D.getNx() << "," << dem.grid2D.getNy() << ") and";
		ss << "grid size (" << grid.getNx() << "," << grid.getNy() << ") don't match!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	for (size_t ii=0; ii<nxy; ii++) {
		const double altitude = std::min( std::max(dem(ii), min_alt), max_alt );
		double &val = grid(ii);
		if (val!=IOUtils::nodata)
			val += trend.f( altitude );
	}
}

//this interpolates VW, DW by converting to u,v and then doing IDW_LAPSE before reconverting to VW, DW
void InterpolationAlgorithm::simpleWindInterpolate(const DEMObject& dem, const std::vector<double>& vecDataVW, const std::vector<double>& vecDataDW, Grid2DObject &VW, Grid2DObject &DW)
{
	if (vecDataVW.size() != vecDataDW.size())
		throw InvalidArgumentException("VW and DW vectors should have the same size!", AT);

	//compute U,v
	std::vector<double> Ve, Vn;
	for (size_t ii=0; ii<vecDataVW.size(); ii++) {
		Ve.push_back( vecDataVW[ii]*sin(vecDataDW[ii]*Cst::to_rad) );
		Vn.push_back( vecDataVW[ii]*cos(vecDataDW[ii]*Cst::to_rad) );
	}

	//spatially interpolate U,V
	std::vector<double> vecAltitudes;
	getStationAltitudes(vecMeta, vecAltitudes);
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating wind", AT);

	if (vecDataVW.size()>=4) { //at least for points to perform detrending
		Fit1D trend;

		getTrend(vecAltitudes, Ve, trend);
		info << trend.getInfo();
		detrend(trend, vecAltitudes, Ve);
		Interpol2D::IDW(Ve, vecMeta, dem, VW);
		retrend(dem, trend, VW);

		getTrend(vecAltitudes, Vn, trend);
		info << trend.getInfo();
		detrend(trend, vecAltitudes, Vn);
		Interpol2D::IDW(Vn, vecMeta, dem, DW);
		retrend(dem, trend, DW);
	} else {
		Interpol2D::IDW(Ve, vecMeta, dem, VW);
		Interpol2D::IDW(Vn, vecMeta, dem, DW);
	}

	//recompute VW, DW in each cell
	const size_t nrCells = VW.size();
	for (size_t ii=0; ii<nrCells; ii++) {
		const double ve = VW(ii);
		const double vn = DW(ii);

		if (ve!=IOUtils::nodata && vn!=IOUtils::nodata) {
			VW(ii) = Optim::fastSqrt_Q3(ve*ve + vn*vn);
			DW(ii) = fmod( atan2(ve,vn) * Cst::to_deg + 360., 360.);
		}
	}
}

} //namespace
