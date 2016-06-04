/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/ResamplingAlgorithms.h>
#include <meteoio/MathOptim.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/meteoLaws/Sun.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/meteoFilters/ProcPSUMDistribute.h> //for the precipitation distribution

#include <cmath>
#include <algorithm>

using namespace std;

namespace mio {

ResamplingAlgorithms* ResamplingAlgorithmsFactory::getAlgorithm(const std::string& i_algoname, const std::string& parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs)
{
	const std::string algoname(IOUtils::strToUpper(i_algoname));

	if (algoname == "NONE" || algoname == "NO") {
		return new NoResampling(algoname, parname, dflt_window_size, vecArgs);
	} else if (algoname == "LINEAR") {
		return new LinearResampling(algoname, parname, dflt_window_size, vecArgs);
	} else if (algoname == "NEAREST") {
		return new NearestNeighbour(algoname, parname, dflt_window_size, vecArgs);
	} else if (algoname == "ACCUMULATE") {
		return new Accumulate(algoname, parname, dflt_window_size, vecArgs);
	} else if (algoname == "DAILY_SOLAR") {
		return new Daily_solar(algoname, parname, dflt_window_size, vecArgs);
	} else if (algoname == "DAILY_AVG") {
		return new DailyAverage(algoname, parname, dflt_window_size, vecArgs);
	} else {
		throw IOException("The resampling algorithm '"+algoname+"' is not implemented" , AT);
	}
}

//compute the partial accumulation at the left of curr_date within a sampling interval
double ResamplingAlgorithms::partialAccumulateAtLeft(const std::vector<MeteoData>& vecM, const size_t& paramindex,
                                                     const size_t& pos, const Date& curr_date)
{
	const size_t end = pos+1;
	if (end>=vecM.size()) return IOUtils::nodata; //reaching the end of the input vector

	const double valend = vecM[end](paramindex);
	if (valend == IOUtils::nodata) return IOUtils::nodata;

	const double jul1 = vecM[pos].date.getJulian(true);
	const double jul2 = vecM[end].date.getJulian(true);

	const double left_accumulation = linearInterpolation(jul1, 0., jul2, valend, curr_date.getJulian(true));

	return left_accumulation;
}

//compute the partial accumulation at the right of curr_date within a sampling interval
double ResamplingAlgorithms::partialAccumulateAtRight(const std::vector<MeteoData>& vecM, const size_t& paramindex,
                                                      const size_t& pos, const Date& curr_date)
{
	const size_t end = pos+1;
	if(end>=vecM.size()) return IOUtils::nodata; //reaching the end of the input vector

	const double valend = vecM[end](paramindex);
	if (valend == IOUtils::nodata) return IOUtils::nodata;

	const double jul1 = vecM[pos].date.getJulian(true);
	const double jul2 = vecM[end].date.getJulian(true);

	const double left_accumulation = linearInterpolation(jul1, 0., jul2, valend, curr_date.getJulian(true));

	return valend - left_accumulation;
}

/**
 * @brief This function returns the last and next valid points around a given position
 * @param pos current position (index)
 * @param paramindex meteo parameter to use
 * @param vecM vector of MeteoData
 * @param resampling_date date to resample
 * @param window_size size of the search window
 * @param indexP1 index of point before the current position (IOUtils::npos if none could be found)
 * @param indexP2 index of point after the current position (IOUtils::npos if none could be found)
 */
void ResamplingAlgorithms::getNearestValidPts(const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
                                              const double& window_size, size_t& indexP1, size_t& indexP2) //HACK
{
	indexP1=IOUtils::npos;
	indexP2=IOUtils::npos;

	const Date dateStart = resampling_date - window_size;
	for (size_t ii=pos; ii-- >0; ) {
		if (vecM[ii].date < dateStart) break;
		if (vecM[ii](paramindex) != IOUtils::nodata) {
			indexP1 = ii;
			break;
		}
	}

	//make sure the search window remains window_size
	const Date dateEnd = (indexP1 != IOUtils::npos)? vecM[indexP1].date+window_size : resampling_date+window_size;
	for (size_t ii=pos; ii<vecM.size(); ++ii) {
		if (vecM[ii].date > dateEnd) break;
		if (vecM[ii](paramindex) != IOUtils::nodata) {
			indexP2 = ii;
			break;
		}
	}
}

/**
 * @brief This function solves the equation y = ax + b for two given points and returns y for a given x
 * @param x1 x-coordinate of first point
 * @param y1 y-coordinate of first point
 * @param x2 x-coordinate of second point
 * @param y2 y-coordinate of second point
 * @param x x-coordinate of desired point
 * @return y-coordinate of desired point
 */
double ResamplingAlgorithms::linearInterpolation(const double& x1, const double& y1,
                                       const double& x2, const double& y2, const double& x)
{
	if (x1 == x2)
		throw IOException("Attempted division by zero", AT);

	//Solving y = ax + b
	const double a = (y2 - y1) / (x2 - x1);
	const double b = y2 - a*x2;

	return (a*x + b);
}

/**
 * @brief For a given date, find the start of the day, considering that for midnight we return the day before!
 * (as is necessary for daily averages, sums, etc that can be provided at midnight for the day before)
 * @param resampling_date current date
 * @return start of the day or start of the day before in case of midnight
 */
Date ResamplingAlgorithms::getDailyStart(const Date& resampling_date)
{
	Date Start( resampling_date );
	Start.rnd(24*3600, Date::DOWN);
	if (Start==resampling_date) //if resampling_date=midnight GMT, the rounding lands on the exact same date
		Start -= 1.;
	
	return Start;
}

/**
 * @brief Find a unique value in a given time interval. 
 * This is useful for retrieving a unique daily average, daily sum, etc
 * @param vecM meteo data time series for the station
 * @param paramindex index of the meteo parameter to process
 * @param pos index of the current time step
 * @param intervalStart start of the interval within which we will look for a valid value (often, start of the day)
 * @param intervalEnd end of the interval within which we will look for a valid value (often, end of the day)
 * @return index of the parameter's valid value
 */
size_t ResamplingAlgorithms::getDailyValue(const std::vector<MeteoData>& vecM, const size_t& paramindex, size_t pos, const Date& intervalStart, const Date& intervalEnd)
{
	size_t indexP1=IOUtils::npos;
	size_t indexP2=IOUtils::npos;
	
	//if pos was not properly pre-positioned, do it
	if (vecM[pos].date<intervalStart) {
		for (; pos<vecM.size(); pos++)
			if (vecM[pos].date>=intervalStart) break;
	}
	if (vecM[pos].date>intervalEnd) {
		for (size_t ii = pos; ii-- >0; ) {
			pos=ii; //because ii gets corrupted at the final iteration if going all the way down
			if (vecM[pos].date<=intervalEnd) break;
		}
	}

	//look for daily sum before the current point
	for (size_t ii=pos; ii-- >0; ) {
		if (vecM[ii].date < intervalStart) break;
		if (vecM[ii](paramindex) != IOUtils::nodata) {
			indexP1 = ii;
			break;
		}
	}

	//look for daily sum after the current point
	for (size_t ii=pos; ii<vecM.size(); ++ii) {
		if (vecM[ii].date > intervalEnd) break;
		if (vecM[ii](paramindex) != IOUtils::nodata) {
			indexP2 = ii;
			break;
		}
	}
	
	if (indexP1!=IOUtils::npos && indexP2!=IOUtils::npos) {
		//if the data was provided at 00:00, a sum for each day has been found
		//we only keep the later one as being the sum of the past day
		int hour1, min1;
		vecM[indexP1].date.getTime(hour1, min1);
		if (hour1==0 && min1==0)
			indexP1=IOUtils::npos;
		else { //otherwise, this means multiple daily sums have been found for the same day
			const string msg = "More than one daily sum of solar radiation found between "+intervalStart.toString(Date::ISO)+" and "+intervalEnd.toString(Date::ISO);
			throw IOException(msg, AT);
		}
	}

	if (indexP1!=IOUtils::npos) return indexP1;
	if (indexP2!=IOUtils::npos) return indexP2;
	return IOUtils::npos;
}

/**********************************************************************************
 * The following functions are implementations of different resampling algorithms *
 **********************************************************************************/

NoResampling::NoResampling(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs)
             : ResamplingAlgorithms(i_algoname, i_parname, dflt_window_size, vecArgs)
{
	if (!vecArgs.empty()) //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments for \""+i_parname+"::"+i_algoname+"\"", AT);
}

std::string NoResampling::toString() const
{
	ostringstream ss;
	ss << right << setw(10) << parname << "::"  << left << setw(15) << algo;
	ss << "[ ]";
	return ss.str();
}

void NoResampling::resample(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                            const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (index >= vecM.size())
		throw IOException("The index of the element to be resampled is out of bounds", AT);

	if (position == ResamplingAlgorithms::exact_match) {
		const double value = vecM[index](paramindex);
		if (value != IOUtils::nodata) {
			md(paramindex) = value; //propagate value
		}
	}

	return;
}


/**********************************************************************************/
NearestNeighbour::NearestNeighbour(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs)
                 : ResamplingAlgorithms(i_algoname, i_parname, dflt_window_size, vecArgs), extrapolate(false)
{
	const size_t nr_args = vecArgs.size();
	if (nr_args==0) return;
	if (nr_args==1) {
		if (vecArgs[0]=="extrapolate")
			extrapolate=true;
		else {
			IOUtils::convertString(window_size, vecArgs[0]);
			window_size /= 86400.; //user uses seconds, internally julian day is used
		}
	} else if (nr_args==2) {
		IOUtils::convertString(window_size, vecArgs[0]);
		window_size /= 86400.; //user uses seconds, internally julian day is used
		if (vecArgs[1]=="extrapolate")
			extrapolate=true;
		else
			throw InvalidArgumentException("Invalid argument \""+vecArgs[1]+"\" for \""+i_parname+"::"+i_algoname+"\"", AT);
	} else {
		throw InvalidArgumentException("Wrong number of arguments for \""+i_parname+"::"+i_algoname+"\"", AT);
	}
}

std::string NearestNeighbour::toString() const
{
	ostringstream ss;
	ss << right << setw(10) << parname << "::"  << left << setw(15) << algo;
	ss << "[ window_size=" << window_size << " extrapolate=" << boolalpha << extrapolate << noboolalpha << " ]";
	return ss.str();
}

void NearestNeighbour::resample(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                                const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (index >= vecM.size())
		throw IOException("The index of the element to be resampled is out of bounds", AT);

	if (position == ResamplingAlgorithms::exact_match) {
		const double value = vecM[index](paramindex);
		if (value != IOUtils::nodata) {
			md(paramindex) = value; //propagate value
			return;
		}
	}

	//if we are at the very beginning or end of vecM and !extrapolate, then there's nothing to do
	if (((!extrapolate) && (position == ResamplingAlgorithms::end))
	    || ((!extrapolate) && (position == ResamplingAlgorithms::begin)))
		return;

	const Date resampling_date = md.date;
	size_t indexP1=IOUtils::npos, indexP2=IOUtils::npos;
	getNearestValidPts(index, paramindex, vecM, resampling_date, window_size, indexP1, indexP2);
	const bool foundP1=(indexP1!=IOUtils::npos), foundP2=(indexP2!=IOUtils::npos);

	//Try to find the nearest neighbour, if there are two equally distant, then return the arithmetic mean
	if (foundP1 && foundP2) { //standard behavior
		const Duration diff1 = resampling_date - vecM[indexP1].date; //calculate time interval to element at index
		const Duration diff2 = vecM[indexP2].date - resampling_date; //calculate time interval to element at index
		const double val1 = vecM[indexP1](paramindex);
		const double val2 = vecM[indexP2](paramindex);

		if (IOUtils::checkEpsilonEquality(diff1.getJulian(true), diff2.getJulian(true), 0.1/1440.)) { //within 6 seconds
			md(paramindex) = Interpol1D::weightedMean(val1, val2, 0.5);
		} else if (diff1 < diff2) {
			md(paramindex) = val1;
		} else if (diff1 > diff2) {
			md(paramindex) = val2;
		}
	} else if (extrapolate) {
		if (foundP1 && !foundP2) { //nearest neighbour on found after index 'index'
			md(paramindex) = vecM[indexP1](paramindex);
		} else if (!foundP1 && foundP2) { //nearest neighbour on found before index 'index'
			md(paramindex) = vecM[indexP2](paramindex);
		} else { // no nearest neighbour with a value different from IOUtils::nodata
			return;
		}
	}
}


/**********************************************************************************/
LinearResampling::LinearResampling(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs)
                 : ResamplingAlgorithms(i_algoname, i_parname, dflt_window_size, vecArgs), extrapolate(false)
{
	const size_t nr_args = vecArgs.size();
	if (nr_args==0) return;
	if (nr_args==1) {
		if (vecArgs[0]=="extrapolate")
			extrapolate=true;
		else {
			IOUtils::convertString(window_size, vecArgs[0]);
			window_size /= 86400.; //user uses seconds, internally julian day is used
		}
	} else if (nr_args==2) {
		IOUtils::convertString(window_size, vecArgs[0]);
		window_size /= 86400.; //user uses seconds, internally julian day is used
		if (vecArgs[1]=="extrapolate")
			extrapolate=true;
		else
			throw InvalidArgumentException("Invalid argument \""+vecArgs[1]+"\" for \""+i_parname+"::"+i_algoname+"\"", AT);
	} else {
		throw InvalidArgumentException("Wrong number of arguments for \""+i_parname+"::"+i_algoname+"\"", AT);
	}
}

std::string LinearResampling::toString() const
{
	ostringstream ss;
	ss << right << setw(10) << parname << "::"  << left << setw(15) << algo;
	ss << "[ window_size=" << window_size << " extrapolate=" << boolalpha << extrapolate << noboolalpha << " ]";
	return ss.str();
}

void LinearResampling::resample(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                                const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (index >= vecM.size())
		throw IOException("The index of the element to be resampled is out of bounds", AT);

	if (position == ResamplingAlgorithms::exact_match) {
		const double value = vecM[index](paramindex);
		if (value != IOUtils::nodata) {
			md(paramindex) = value; //propagate value
			return;
		}
	}

	//if we are at the very beginning or end of vecM and !extrapolate, then there's nothing to do
	if (((!extrapolate) && (position == ResamplingAlgorithms::end))
	    || ((!extrapolate) && (position == ResamplingAlgorithms::begin)))
		return;

	const Date resampling_date = md.date;
	size_t indexP1=IOUtils::npos, indexP2=IOUtils::npos;
	getNearestValidPts(index, paramindex, vecM, resampling_date, window_size, indexP1, indexP2);
	bool foundP1=(indexP1!=IOUtils::npos), foundP2=(indexP2!=IOUtils::npos);

	//do nothing if we can't interpolate, and extrapolation is not explicitly activated
	if ((!extrapolate) && ((!foundP1) || (!foundP2)))
		return;
	//do nothing if not at least one value different from IOUtils::nodata has been found
	if (!foundP1 && !foundP2)
		return;

	//At this point we either have a valid indexP1 or indexP2 and we can at least try to extrapolate
	if (!foundP1 && foundP2) { //only nodata values found before index, try looking after indexP2
		for (size_t ii=indexP2+1; ii<vecM.size(); ii++) {
			if (vecM[ii](paramindex) != IOUtils::nodata) {
				indexP1 = ii;
				foundP1 = true;
				break;
			}
		}
	} else if (foundP1 && !foundP2) { //only nodata found after index, try looking before indexP1
		for (size_t ii=indexP1; (ii--) > 0; ){ 
			if (vecM[ii](paramindex) != IOUtils::nodata){
				indexP2=ii;
				foundP2 = true;
				break;
			}
		}
	}

	if (!foundP1 || !foundP2) //now at least two points need to be present
		return;

	//At this point indexP1 and indexP2 point to values that are different from IOUtils::nodata
	const double val1 = vecM[indexP1](paramindex);
	const double jul1 = vecM[indexP1].date.getJulian(true);
	const double val2 = vecM[indexP2](paramindex);
	const double jul2 = vecM[indexP2].date.getJulian(true);

	md(paramindex) = linearInterpolation(jul1, val1, jul2, val2, resampling_date.getJulian(true));
}


/**********************************************************************************/
Accumulate::Accumulate(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs)
           : ResamplingAlgorithms(i_algoname, i_parname, dflt_window_size, vecArgs),
             accumulate_period(IOUtils::nodata), strict(false)
{
	const size_t nr_args = vecArgs.size();
	if (nr_args<1 || nr_args>2)
		throw InvalidArgumentException("Please at least provide accumulation period (in seconds) for \""+i_parname+"::"+i_algoname+"\"", AT);

	bool period_read = false;
	for (size_t ii=0; ii<nr_args; ii++) {
		if (IOUtils::isNumeric(vecArgs[ii])) {
			if (period_read==true)
				throw InvalidArgumentException("Two arguments "+i_algoname+" resampling has been deprecated! Please use the \"PSUM_Distribute\" Processing Element instead!", AT);

			IOUtils::convertString(accumulate_period, vecArgs[ii]);
			accumulate_period /= 86400.; //user uses seconds, internally julian day is used
			if (accumulate_period<=0.) {
				std::ostringstream ss;
				ss << "Invalid accumulation period (" << accumulate_period << ") for \"" << i_parname << "::" << i_algoname << "\"";
				throw InvalidArgumentException(ss.str(), AT);
			}
			period_read = true;
		} else if (vecArgs[ii]=="strict" && !strict) {
			if (strict) //do not set strict more than once!
				throw InvalidArgumentException("Do not provide \"strict\" more than once for \""+i_parname+"::"+i_algoname+"\"", AT);
			strict = true;
		} else 
			throw InvalidArgumentException("Invalid argument \""+vecArgs[ii]+"\" for \""+i_parname+"::"+i_algoname+"\"", AT);
	}
}

std::string Accumulate::toString() const
{
	ostringstream ss;
	ss << right << setw(10) << parname << "::"  << left << setw(15) << algo;
	ss << "[ period=" << accumulate_period << " strict=" << boolalpha << strict << noboolalpha << " ]";
	return ss.str();
}

//find the index just before the start of accumulation period
size_t Accumulate::findStartOfPeriod(const std::vector<MeteoData>& vecM, const size_t& index, const Date& dateStart)
{
	size_t start_idx = IOUtils::npos;
	for (size_t idx=index; idx--> 0; ) {
		const Date curr_date = vecM[idx].date;
		if (curr_date <= dateStart) {
			start_idx = idx;
			break;
		}
	}

	return start_idx;
}

double Accumulate::easySampling(const std::vector<MeteoData>& vecM, const size_t& paramindex, const size_t& /*index*/, const size_t& start_idx, const Date& dateStart, const Date& resampling_date) const
{//to keep in mind: start_idx is last index <= dateStart and index is first index >= resampling_date
	double sum = IOUtils::nodata;
	const double start_val = partialAccumulateAtLeft(vecM, paramindex, start_idx, dateStart);
	const double end_val = partialAccumulateAtLeft(vecM, paramindex, start_idx, resampling_date);

	if (start_val!=IOUtils::nodata && end_val!=IOUtils::nodata)
		sum = end_val - start_val;

	return sum;
}

double Accumulate::complexSampling(const std::vector<MeteoData>& vecM, const size_t& paramindex, const size_t& index, const size_t& start_idx, const Date& dateStart, const Date& resampling_date) const
{//to keep in mind: start_idx is last index <= dateStart and index is first index >= resampling_date
	double sum = IOUtils::nodata;
	//resample begining point, in the [start_idx ; start_idx+1] interval
	const double start_value = partialAccumulateAtRight(vecM, paramindex, start_idx, dateStart);
	if (start_value!=IOUtils::nodata)
		sum=start_value;
	else if (strict) return IOUtils::nodata;

	//sum all whole periods AFTER the begining point, in the [start_idx+2 ; index-1] interval
	for (size_t idx=(start_idx+2); idx<index; idx++) {
		const double curr_value = vecM[idx](paramindex);
		if (curr_value!=IOUtils::nodata) {
			if (sum!=IOUtils::nodata) sum += curr_value;
			else sum = curr_value;
		} else if (strict) return IOUtils::nodata;
	}

	//resample end point, in the [index-1 ; index] interval
	const double end_val = partialAccumulateAtLeft(vecM, paramindex, index-1, resampling_date);
	if (end_val!=IOUtils::nodata) {
		if (sum!=IOUtils::nodata) sum += end_val;
		else sum = end_val;
	} else if (strict) return IOUtils::nodata;

	return sum;
}

//index is the first element AFTER the resampling_date
void Accumulate::resample(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                          const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (index >= vecM.size())
		throw IOException("The index of the element to be resampled is out of bounds", AT);
	if (position==ResamplingAlgorithms::begin || position==ResamplingAlgorithms::end)
		return;

	md(paramindex) = IOUtils::nodata;
	const Date resampling_date = md.date;

	const Date dateStart(resampling_date.getJulian() - accumulate_period, resampling_date.getTimeZone());
	const size_t start_idx = findStartOfPeriod(vecM, index, dateStart);
	if (start_idx==IOUtils::npos) {//No acceptable starting point found
		cerr << "[W] Could not accumulate " << vecM.at(0).getNameForParameter(paramindex) << ": ";
		cerr << "not enough data for accumulation period at date " << resampling_date.toString(Date::ISO) << "\n";
		return;
	}

	if ((index - start_idx) <= 1) {//easy upsampling when start & stop are in the same input time step
		//upsampling (for example, generate 15min values from hourly data)
		const double sum = easySampling(vecM, paramindex, index, start_idx, dateStart, resampling_date);
		md(paramindex) = sum; //if resampling was unsuccesful, sum==IOUtils::nodata
	} else {
		//downsampling (for example, generate daily values from hourly data)
		//and upsampling when resampled period falls accross a measurement timestamp
		const double sum = complexSampling(vecM, paramindex, index, start_idx, dateStart, resampling_date);
		md(paramindex) = sum; //if resampling was unsuccesful, sum==IOUtils::nodata
	}
}


/**********************************************************************************/
const double Daily_solar::soil_albedo = .23; //grass
const double Daily_solar::snow_albedo = .85; //snow
const double Daily_solar::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo
const size_t Daily_solar::samples_per_day = 24*3; //every 20 minutes

Daily_solar::Daily_solar(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs)
            : ResamplingAlgorithms(i_algoname, i_parname, dflt_window_size, vecArgs), radiation(), station_index(), dateStart(), dateEnd(), loss_factor()
{
	const size_t nr_args = vecArgs.size();
	if (nr_args>0) {
		throw InvalidArgumentException("Too many arguments for \""+i_parname+"::"+i_algoname+"\"", AT);
	}
}

std::string Daily_solar::toString() const
{
	ostringstream ss;
	ss << right << setw(10) << parname << "::"  << left << setw(15) << algo;
	return ss.str();
}

//compute the daily sum of solar radiation as well as fill a vector containing the solar radiation at a regular sampling rate
double Daily_solar::compRadiation(const double& lat, const double& lon, const double& alt, const double& HS, const size_t& stat_idx)
{
	const double P = Atmosphere::stdAirPressure(alt);
	double albedo = 0.5;
	if (HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
		albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;
	const double TA=274.98, RH=0.666; //the reduced precipitable water will get an average value

	SunObject sun(lat, lon, alt);
	double sum = 0.;
	size_t index=0;
	for (Date date(dateStart[stat_idx]); date<dateEnd[stat_idx]; date += 1./double(samples_per_day)) {
		//compute potential solar radiation at this time step
		sun.setDate(date.getJulian(), date.getTimeZone());
		sun.calculateRadiation(TA, RH, P, albedo);
		double toa, direct, diffuse;
		sun.getHorizontalRadiation(toa, direct, diffuse);
		const double global_h = direct+diffuse;

		//compute the integral by a simple triangle method
		sum += global_h * static_cast<double>(24*3600/samples_per_day);

		//store the radiation for later reuse
		radiation[stat_idx][index++] = global_h;
	}

	return sum;
}

//interpolate the solar radiation from the vector containing regularly sampled solar radiation
double Daily_solar::getSolarInterpol(const Date& resampling_date, const size_t& stat_idx) const
{
	const double in_day = (resampling_date.getJulian() - dateStart[stat_idx].getJulian()) / (dateEnd[stat_idx].getJulian() - dateStart[stat_idx].getJulian());
	const size_t vec_index = static_cast<size_t>( Optim::floor(in_day*samples_per_day) ); //so the sample will be between vec_index and vec_index+1
	const double weight = in_day*static_cast<double>(samples_per_day) - static_cast<double>(vec_index);

	return Interpol1D::weightedMean(radiation[stat_idx][vec_index], radiation[stat_idx][vec_index+1], weight);
}

//a new, previously unknown station has been found, allocate the memory
size_t Daily_solar::getStationIndex(const std::string& key)
{
	const size_t nr_stations = station_index.size();
	for (size_t ii=0; ii<nr_stations; ++ii) {
		if (station_index[ii]==key)
			return ii;
	}

	radiation.push_back( vector<double>(samples_per_day, 0.) );
	loss_factor.push_back( 0. );

	Date null_date(0., 0.);
	dateStart.push_back( null_date );
	dateEnd.push_back( null_date );

	station_index.push_back( key );
	return nr_stations; //the new index is the old nr_stations
}

void Daily_solar::resample(const size_t& index, const ResamplingPosition& /*position*/, const size_t& paramindex,
                           const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (index >= vecM.size())
		throw IOException("The index of the element to be resampled is out of bounds", AT);

	if (paramindex!=MeteoData::ISWR && paramindex!=MeteoData::RSWR)
		throw IOException("This method only applies to short wave radiation! (either ISWR or RSWR)", AT);

	const double lat = md.meta.position.getLat();
	const double lon = md.meta.position.getLon();
	const double alt = md.meta.position.getAltitude();
	if (lat==IOUtils::nodata || lon==IOUtils::nodata || alt==IOUtils::nodata) return;
	const double HS = md(MeteoData::HS);

	//get station index
	const size_t stat_idx = getStationIndex(md.meta.stationID);

	//has the radiation already been calculated for this day and station?
	if (md.date<dateStart[stat_idx] || md.date>=dateEnd[stat_idx]) {
		dateStart[stat_idx] = getDailyStart(md.date);
		dateEnd[stat_idx] = dateStart[stat_idx]+1.;
		//const size_t indexP = getNearestValidPt(vecM, paramindex, stat_idx, index);
		const size_t indexP = getDailyValue(vecM, paramindex, index, dateStart[stat_idx], dateEnd[stat_idx]);
		if (indexP==IOUtils::npos) { //no daily sum found for the current day
			loss_factor[stat_idx] = IOUtils::nodata;
			return;
		}

		const double daily_sum = compRadiation(lat, lon, alt, HS, stat_idx);
		loss_factor[stat_idx] = (daily_sum>0)? vecM[indexP](paramindex) / daily_sum : 0.; //in case of polar night...
	}

	if (loss_factor[stat_idx]==IOUtils::nodata) //the station could not be calculated for this day
		return;

	//interpolate radiation for this timestep and write it out
	const double rad = getSolarInterpol(md.date, stat_idx);
	if (paramindex==MeteoData::ISWR) {
		md(paramindex) = loss_factor[stat_idx] * rad;
	} else {
		double albedo = 0.5;
		if (HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
			albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;
		md(paramindex) = loss_factor[stat_idx] * rad * albedo;
	}

	md.setResampled(true);
}


/**********************************************************************************/
DailyAverage::DailyAverage(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs)
                 : ResamplingAlgorithms(i_algoname, i_parname, dflt_window_size, vecArgs), range(IOUtils::nodata), phase(0.)
{
	const size_t nr_args = vecArgs.size();
	
	if (nr_args==1) {
		IOUtils::convertString(range, vecArgs[0]);
	} else if (nr_args==2) {
		IOUtils::convertString(range, vecArgs[0]);
		IOUtils::convertString(phase, vecArgs[1]);
		phase *= -1; //shift the minimum *later* in the day
	} else if (nr_args>2) {
		throw InvalidArgumentException("Wrong number of arguments for \""+i_parname+"::"+i_algoname+"\"", AT);
	}
}

std::string DailyAverage::toString() const
{
	ostringstream ss;
	ss << right << setw(10) << parname << "::"  << left << setw(15) << algo;
	ss << "[ window_size=" << window_size;
	if (range!=IOUtils::nodata) 
		ss << " range=" << range;
	if (phase!=0.)
		ss << " phase=" << phase;
	ss << " ]";
	return ss.str();
}

double DailyAverage::getValue(const std::vector<MeteoData>& vecM, const size_t& paramindex, const size_t& index, const Date& dayStart, const double& frac_day) const
{
	const Date dayEnd = dayStart+1.;
	const size_t indexP = getDailyValue(vecM, paramindex, index, dayStart, dayEnd);
	if (indexP==IOUtils::npos)
		return IOUtils::nodata;
	
	const double avg = vecM[indexP](paramindex); //from getDailyValue we are guaranteed to have a value
	
	//search for min or max or both
	const string param_name = vecM[indexP].getNameForParameter(paramindex);
	const string min_name = param_name+"_MIN";
	const string max_name = param_name+"_MAX";
	const double min = (vecM[indexP].param_exists( min_name ))?  vecM[indexP](min_name) : IOUtils::nodata;
	const double max = (vecM[indexP].param_exists( max_name ))?  vecM[indexP](max_name) : IOUtils::nodata;
	double A=range; //fallback: user provided range
	
	if (min==IOUtils::nodata && max==IOUtils::nodata) {
		if (A==IOUtils::nodata) {
			ostringstream ss;
			ss << "Missing parameters for " << algo << " on " << dayStart.toString(Date::ISO_DATE) << ". Please consider providing a fallback range argument!";
			throw InvalidArgumentException(ss.str(), AT);
		}
	} else if (min!=IOUtils::nodata && max==IOUtils::nodata) {
		A = (avg -  min);
	} else if (min==IOUtils::nodata && max!=IOUtils::nodata) {
		A =(max - avg);
	} else
		throw InvalidArgumentException("Providing both AVG, MIN and MAX for \'"+algo+"\' is not supported!", AT);

	return A * sin( 2.*Cst::PI * (frac_day-.25+phase) ) + avg;
}

void DailyAverage::resample(const size_t& index, const ResamplingPosition& /*position*/, const size_t& paramindex,
                                const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (index >= vecM.size())
		throw IOException("The index of the element to be resampled is out of bounds", AT);

	//define various time parameters
	const double julian = md.date.getJulian();
	double frac_day = (julian+.5 - Optim::intPart(julian+.5)); // julian day starts at noon and we want the fraction of the day, from 0 to 1!
	if (frac_day==0.) frac_day = 1.; //because we considere here that midnight belongs to the day before
	
	const Date dayStart = getDailyStart(md.date);
	const double val1 = getValue(vecM, paramindex, index, dayStart, frac_day);
	if (val1==IOUtils::nodata)
		return; //we could also decide to use val0 & val2 
	
	//get the other values (day before, day after)
	const double val0 = getValue(vecM, paramindex, index, dayStart-1., frac_day);
	const double val2 = getValue(vecM, paramindex, index, dayStart+1., frac_day);
	
	double val = val1;
	if (frac_day!=0.5) { //at noon, we purely take val1
		const double w1 = 1./Optim::pow2( abs(.5-frac_day) );
		double norm = w1;
		val *= w1;
		
		if (val0!=IOUtils::nodata) {
			const double w0 = 1./Optim::pow2( frac_day+0.5 );
			val += w0*val0;
			norm += w0;
		}
		if (val2!=IOUtils::nodata) {
			const double w2 = 1./Optim::pow2( 1.5-frac_day );
			val += w2*val2;
			norm += w2;
		}
		val /= norm;
	}
	
	md(paramindex) = val;
	//very basic range checking
	//HACK this means that this should be implemented as a data creator so it could be filtered afterward
	//HACK or we could one more pass of filtering *after* the resampling
	if (paramindex==MeteoData::RH) {
		if (md(paramindex)<0.01) md(paramindex)=0.01;
		if (md(paramindex)>1.) md(paramindex)=1.;
	} else if (paramindex==MeteoData::TA ||  paramindex==MeteoData::VW ||  paramindex==MeteoData::VW_MAX || paramindex==MeteoData::ISWR || paramindex==MeteoData::RSWR || paramindex==MeteoData::ILWR || paramindex==MeteoData::TSG || paramindex==MeteoData::TSS || paramindex==MeteoData::TAU_CLD) {
		if (md(paramindex)<0.)
			md(paramindex)=0.;
	}
		
	md.setResampled(true);
}

} //namespace

