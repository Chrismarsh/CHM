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
#include "SMETIO.h"
#include <meteoio/IOUtils.h>

using namespace std;

namespace mio {
/**
 * @page smetio SMET
 * @section smetio_format Format
 * The Station meteo data files is a station centered, ascii file format that has been designed with flexibility and ease of use in mind. Please refer to its
 * <a href="../SMET_specifications.pdf">official format specification</a> for more information (including the list of standard parameters: TA, TSS, TSG,
 * RH, VW, DW, ISWR, OSWR, ILWR, OLWR, PINT, PSUM, HS). For PINT, it is assumed that the intensity (in mm/h) is valid for the whole period between the actual
 * time step and the previous one.
 *
 * This plugin can also provide Points Of Interest, given as a SMET file containing either latitude/longitude/altitude or easting/northing/altitude. For the latter,
 * the header must contain the epsg code (see example below).
 *
 * Non-standard parameters can also be given, such as extra snow temperatures. These parameters will then take the name that has been given in "fields", converted to uppercase.
 * It is usually a good idea to number these parameters, such as TS1, TS2, TS3 for a serie of temperatures at various positions.
 *
 * @section smetio_units Units
 * All units are MKSA, the only exception being the precipitations that are in mm/h or mm/{time step}. It is however possible to use  multipliers and offsets (but they must be specified in the file header).
 *
 * @section smetio_keywords Keywords
 * This plugin uses the following keywords:
 * - STATION#: input filename (in METEOPATH). As many meteofiles as needed may be specified
 * - METEOPATH: meteo files directory where to read/write the meteofiles; [Input] and [Output] sections
 * - METEOPARAM: output file format options (ASCII or BINARY that might be followed by GZIP)
 * - POIFILE: a path+file name to the a file containing grid coordinates of Points Of Interest (for special outputs)
 *
 * Example:
 * @code
 * [Input]
 * METEO = SMET
 * METEOPATH = ./input
 * STATION1 = uppper_station.smet
 * STATION2 = lower_station.smet
 * STATION3 = outlet_station.smet
 * [Output]
 * METEOPATH = ./output
 * METEOPARAM = ASCII GZIP
 * @endcode
 *
 * Below is an example of Points Of Interest input:
 * @code
 * SMET 1.1 ASCII
 * [HEADER]
 * station_id	= my_pts
 * epsg	= 21781
 * nodata	= -999
 * fields = easting northing altitude
 * [DATA]
 * 832781 187588 2115
 * 635954 80358 2428
 * @endcode
 *
 */

const char* SMETIO::dflt_extension = ".smet";

SMETIO::SMETIO(const std::string& configfile)
        : cfg(configfile),
          coordin(), coordinparam(), coordout(), coordoutparam(),
          vec_smet_reader(), vecFiles(), outpath(), in_dflt_TZ(0.), out_dflt_TZ(0.),
          plugin_nodata(IOUtils::nodata), nr_stations(0), outputIsAscii(true), outputIsGzipped(false)
{
	parseInputOutputSection();
}

SMETIO::SMETIO(const Config& cfgreader)
        : cfg(cfgreader),
          coordin(), coordinparam(), coordout(), coordoutparam(),
          vec_smet_reader(), vecFiles(), outpath(), in_dflt_TZ(0.), out_dflt_TZ(0.),
          plugin_nodata(IOUtils::nodata), nr_stations(0), outputIsAscii(true), outputIsGzipped(false)
{
	parseInputOutputSection();
}

SMETIO::~SMETIO() throw()
{

}

void SMETIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*_name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SMETIO::read2DGrid(Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SMETIO::readDEM(DEMObject& /*dem_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SMETIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SMETIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SMETIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{//HACK: It should support coordinates in the data, ie: it should use the given date! (and TZ)
	vecStation.clear();
	vecStation.reserve(nr_stations);

	//Now loop through all requested stations, open the respective files and parse them
	for (size_t ii=0; ii<vec_smet_reader.size(); ii++){
		StationData sd;
		smet::SMETReader& myreader = vec_smet_reader[ii];

		read_meta_data(myreader, sd);
		vecStation.push_back(sd);
	}
}

void SMETIO::parseInputOutputSection()
{
	//default timezones
	in_dflt_TZ = out_dflt_TZ = IOUtils::nodata;
	cfg.getValue("TIME_ZONE","Input",in_dflt_TZ,IOUtils::nothrow);
	cfg.getValue("TIME_ZONE","Output",out_dflt_TZ,IOUtils::nothrow);

	// Parse the [Input] and [Output] sections within Config object cfg
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);

	//Parse input section: extract number of files to read and store filenames in vecFiles
	std::string inpath, in_meteo;
	cfg.getValue("METEO", "Input", in_meteo, IOUtils::nothrow);
	if (in_meteo == "SMET") { //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("METEOPATH", "Input", inpath);
		std::vector<std::string> vecFilenames;
		cfg.getValues("STATION", "INPUT", vecFilenames);

		for (size_t ii=0; ii<vecFilenames.size(); ii++) {
			const string filename = vecFilenames[ii];
			const string extension = IOUtils::getExtension(filename);
			const std::string file_and_path = (!extension.empty())? inpath+"/"+filename : inpath+"/"+filename+dflt_extension;

			if (!IOUtils::validFileAndPath(file_and_path)) //Check whether filename is valid
				throw InvalidFileNameException(file_and_path, AT);
			vecFiles.push_back(file_and_path);
			vec_smet_reader.push_back(smet::SMETReader(file_and_path));
		}
	}

	//Parse output section: extract info on whether to write ASCII or BINARY format, gzipped or not
	outpath.clear();
	outputIsAscii = true;
	outputIsGzipped = false;

	vector<string> vecArgs;
	cfg.getValue("METEOPATH", "Output", outpath, IOUtils::nothrow);
	cfg.getValue("METEOPARAM", "Output", vecArgs, IOUtils::nothrow); //"ASCII|BINARY GZIP"

	if (outpath.empty())
		return;

	if (vecArgs.empty())
		vecArgs.push_back("ASCII");

	if (vecArgs.size() > 2)
		throw InvalidFormatException("Too many values for key METEOPARAM", AT);

	if (vecArgs[0] == "BINARY")
		outputIsAscii = false;
	else if (vecArgs[0] == "ASCII")
		outputIsAscii = true;
	else
		throw InvalidFormatException("The first value for key METEOPARAM may only be ASCII or BINARY", AT);

	if (vecArgs.size() == 2){
		if (vecArgs[1] != "GZIP")
			throw InvalidFormatException("The second value for key METEOPARAM may only be GZIP", AT);

		outputIsGzipped = true;
	}

}

void SMETIO::identify_fields(const std::vector<std::string>& fields, std::vector<size_t>& indexes,
                             bool& julian_present, MeteoData& md)
{
	/*
	 * This function associates a parameter index for MeteoData objects with the
	 * lineup of field types in a SMET header. The following SMET fields are treated
	 * exceptionally:
	 * - julian, associated with IOUtils::npos
	 * - latitude, associated with IOUtils::npos-1
	 * - longitude, associated with IOUtils::npos-2
	 * - easting, associated with IOUtils::npos-3
	 * - norhting, associated with IOUtils::npos-4
	 * - altitude, associated with IOUtils::npos-5
	 * If a paramter is unknown in the fields section, then it is added as separate field to MeteoData
	 */
	for (size_t ii=0; ii<fields.size(); ii++){
		const string& key = fields[ii];

		if(md.param_exists(key)) {
			indexes.push_back(md.getParameterIndex(key));
			continue;
		}

		//specific key mapping
		if (key == "OSWR") {
			indexes.push_back(md.getParameterIndex("RSWR"));
		} else if (key == "OLWR") {
			md.addParameter("OLWR");
			indexes.push_back(md.getParameterIndex("OLWR"));
		} else if (key == "PINT") { //in mm/h
			md.addParameter("PINT");
			indexes.push_back(md.getParameterIndex("PINT"));
		} else if (key == "julian") {
			julian_present = true;
			indexes.push_back(IOUtils::npos);
		} else if (key == "latitude") {
			indexes.push_back(IOUtils::npos-1);
		} else if (key == "longitude") {
			indexes.push_back(IOUtils::npos-2);
		} else if (key == "easting") {
			indexes.push_back(IOUtils::npos-3);
		} else if (key == "northing") {
			indexes.push_back(IOUtils::npos-4);
		} else if (key == "altitude") {
			indexes.push_back(IOUtils::npos-5);
		} else {
			//this is an extra parameter, we convert to uppercase
			const std::string extra_param = IOUtils::strToUpper(key);
			md.addParameter(extra_param);
			indexes.push_back(md.getParameterIndex(extra_param));
		}
	}
}

void SMETIO::read_meta_data(const smet::SMETReader& myreader, StationData& meta)
{
	/*
	 * This function reads in the header data provided by a SMETReader object.
	 * SMETReader objects read all the header info upon construction and can subsequently
	 * be queried for that info
	 */
	const double nodata_value = myreader.get_header_doublevalue("nodata");

	meta.position.setProj(coordin, coordinparam); //set the default projection from config file
	if (myreader.location_in_header(smet::WGS84)){
		const double lat = myreader.get_header_doublevalue("latitude");
		const double lon = myreader.get_header_doublevalue("longitude");
		const double alt = myreader.get_header_doublevalue("altitude");
		meta.position.setLatLon(lat, lon, alt);
	}

	if (myreader.location_in_header(smet::EPSG)){
		const double east  = myreader.get_header_doublevalue("easting");
		const double north = myreader.get_header_doublevalue("northing");
		const double alt   = myreader.get_header_doublevalue("altitude");
		const short int epsg  = (short int)(floor(myreader.get_header_doublevalue("epsg") + 0.1));
		meta.position.setEPSG(epsg); //this needs to be set before calling setXY(...)
		meta.position.setXY(east, north, alt);
	}

	meta.stationID = myreader.get_header_value("station_id");
	meta.stationName = myreader.get_header_value("station_name");

	const bool data_epsg = myreader.location_in_data(smet::EPSG);
	if (data_epsg){
		const double d_epsg = myreader.get_header_doublevalue("epsg");
		const short int epsg = (d_epsg != nodata_value)? (short int)(floor(d_epsg + 0.1)): IOUtils::snodata;
		meta.position.setEPSG(epsg);
	}
}

void SMETIO::copy_data(const smet::SMETReader& myreader,
                       const std::vector<std::string>& timestamps,
                       const std::vector<double>& mydata, std::vector<MeteoData>& vecMeteo)
{
	/*
	 * This function parses the data read from a SMETReader object, a vector<double>,
	 * and copies the values into their respective places in the MeteoData structure
	 * Meta data, whether in header or in data is also handled
	 */
	const string myfields = myreader.get_header_value("fields");
	vector<string> fields;
	IOUtils::readLineToVec(myfields, fields);

	bool julian_present = false;
	MeteoData md;
	vector<size_t> indexes;
	identify_fields(fields, indexes, julian_present, md);

	if ((timestamps.empty()) && (!julian_present)) return; //nothing to do

	const bool pint_present = md.param_exists("PINT");
	const bool data_wgs84 = myreader.location_in_data(smet::WGS84);
	const bool data_epsg = myreader.location_in_data(smet::EPSG);

	read_meta_data(myreader, md.meta);

	const double nodata_value = myreader.get_header_doublevalue("nodata");
	double current_timezone = myreader.get_header_doublevalue("tz");
	if (current_timezone == nodata_value)
		current_timezone = in_dflt_TZ;
	const bool timestamp_present = myreader.contains_timestamp();

	const size_t nr_of_fields = indexes.size();
	const size_t nr_of_lines = mydata.size() / nr_of_fields;

	double lat=IOUtils::nodata, lon=IOUtils::nodata, east=IOUtils::nodata, north=IOUtils::nodata, alt=IOUtils::nodata;
	size_t current_index = 0; //index to vec_data
	double previous_ts = IOUtils::nodata;
	for (size_t ii = 0; ii<nr_of_lines; ii++){
		vecMeteo.push_back(md);
		MeteoData& tmp_md = vecMeteo.back();

		if (timestamp_present)
			IOUtils::convertString(tmp_md.date, timestamps[ii], current_timezone);

		//Copy data points
		for (size_t jj=0; jj<nr_of_fields; jj++){
			const double& current_data = mydata[current_index];
			if (indexes[jj] >= IOUtils::npos-5){ //the special fields have high indexes
				if (indexes[jj] == IOUtils::npos){
					if (!timestamp_present){
						if (current_data != nodata_value)
							tmp_md.date.setDate(current_data, current_timezone);
					}
				} else if (indexes[jj] == IOUtils::npos-1){
					lat = current_data;
				} else if (indexes[jj] == IOUtils::npos-2){
					lon = current_data;
				} else if (indexes[jj] == IOUtils::npos-3){
					east = current_data;
				} else if (indexes[jj] == IOUtils::npos-4){
					north = current_data;
				} else if (indexes[jj] == IOUtils::npos-5){
					alt = current_data;
				}
			} else {
				if (current_data == nodata_value)
					tmp_md(indexes[jj]) = IOUtils::nodata;
				else
					tmp_md(indexes[jj]) = current_data;
			}

			if (data_epsg)
				tmp_md.meta.position.setXY(east, north, alt);

			if (data_wgs84)
				tmp_md.meta.position.setXY(lat, lon, alt);

			current_index++;
		}

		if ((pint_present) && (tmp_md(MeteoData::PSUM) == IOUtils::nodata)) {
			const double pint = tmp_md("PINT");
			if (pint==0.) {
				tmp_md(MeteoData::PSUM) = 0.;
			} else if (previous_ts!=IOUtils::nodata) {
				const double acc_period = (tmp_md.date.getJulian() - previous_ts) * 24.; //in hours
				tmp_md(MeteoData::PSUM) = pint * acc_period;
			}
		}

		previous_ts = tmp_md.date.getJulian();
	}
}

void SMETIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                           std::vector< std::vector<MeteoData> >& vecMeteo,
                           const size_t& stationindex)
{
	//Make sure that vecMeteo have the correct dimension and stationindex is valid
	size_t startindex=0, endindex=vecFiles.size();
	if (stationindex != (size_t)IOUtils::npos){ //HACK do we really still need stationindex??
		if ((stationindex < vecFiles.size()) || (stationindex < vecMeteo.size())){
			startindex = stationindex;
			endindex = stationindex+1;
		} else {
			throw IndexOutOfBoundsException("Invalid stationindex", AT);
		}

		vecMeteo[stationindex].clear();
	} else {
		vecMeteo.clear();
		vecMeteo = vector< vector<MeteoData> >(vecFiles.size());
		vecMeteo.reserve(nr_stations);
	}

	//Now loop through all requested stations, open the respective files and parse them
	for (size_t ii=startindex; ii<endindex; ii++){
		const string& filename = vecFiles.at(ii); //filename of current station

		if (!IOUtils::fileExists(filename))
			throw FileNotFoundException(filename, AT);

		smet::SMETReader& myreader = vec_smet_reader.at(ii);
		myreader.convert_to_MKSA(true); // we want converted values for MeteoIO

		vector<double> mydata; //sequentially store all data in the smet file
		vector<string> mytimestamps;

		if (myreader.contains_timestamp()){
			myreader.read(dateStart.toString(Date::ISO), dateEnd.toString(Date::ISO), mytimestamps, mydata);
		} else {
			myreader.read(dateStart.getJulian(), dateEnd.getJulian(), mydata);
		}

		copy_data(myreader, mytimestamps, mydata, vecMeteo[ii]);
	}
}

void SMETIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string&)
{
	//Loop through all stations
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		//1. check consistency of station data position -> write location in header or data section
		StationData sd;
		sd.position.setProj(coordout, coordoutparam);
		const bool isConsistent = checkConsistency(vecMeteo.at(ii), sd);

		if (sd.stationID.empty()){
			ostringstream ss;
			ss << "Station" << ii+1;
			sd.stationID = ss.str();
		}

		const string filename = outpath + "/" + sd.stationID + ".smet";
		if (!IOUtils::validFileAndPath(filename)) //Check whether filename is valid
			throw InvalidFileNameException(filename, AT);

		//2. check which meteo parameter fields are actually in use
		const size_t nr_of_parameters = getNrOfParameters(sd.stationID, vecMeteo[ii]);
		vector<bool> vecParamInUse = vector<bool>(nr_of_parameters, false);
		vector<string> vecColumnName = vector<string>(nr_of_parameters, "NULL");
		double timezone = IOUtils::nodata; //time zone of the data
		checkForUsedParameters(vecMeteo[ii], nr_of_parameters, timezone, vecParamInUse, vecColumnName);
		if(out_dflt_TZ != IOUtils::nodata) timezone=out_dflt_TZ; //if the user set an output time zone, all will be converted to it

		try {
			const smet::SMETType type = (outputIsAscii)? smet::ASCII : smet::BINARY;

			smet::SMETWriter mywriter(filename, type, outputIsGzipped);
			generateHeaderInfo(sd, outputIsAscii, isConsistent, timezone,
                               nr_of_parameters, vecParamInUse, vecColumnName, mywriter);

			vector<string> vec_timestamp;
			vector<double> vec_data;
			for (size_t jj=0; jj<vecMeteo[ii].size(); jj++) {
				if (outputIsAscii){
					if (out_dflt_TZ != IOUtils::nodata) { //user-specified time zone
						Date tmp_date(vecMeteo[ii][jj].date);
						tmp_date.setTimeZone(out_dflt_TZ);
						vec_timestamp.push_back(tmp_date.toString(Date::ISO));
					} else {
						vec_timestamp.push_back(vecMeteo[ii][jj].date.toString(Date::ISO));
					}
				} else {
					double julian;
					if(out_dflt_TZ!=IOUtils::nodata) {
						Date tmp_date(vecMeteo[ii][jj].date);
						tmp_date.setTimeZone(out_dflt_TZ);
						julian = tmp_date.getJulian();
					} else {
						julian = vecMeteo[ii][jj].date.getJulian();
					}
					vec_data.push_back(julian);
				}

				if (!isConsistent) { //Meta data changes
					vec_data.push_back(vecMeteo[ii][jj].meta.position.getLat());
					vec_data.push_back(vecMeteo[ii][jj].meta.position.getLon());
					vec_data.push_back(vecMeteo[ii][jj].meta.position.getAltitude());
				}

				for (size_t kk=0; kk<nr_of_parameters; kk++) {
					if (vecParamInUse[kk])
						vec_data.push_back(vecMeteo[ii][jj](kk)); //add data value
				}
			}

			if (outputIsAscii) mywriter.write(vec_timestamp, vec_data);
			else mywriter.write(vec_data);

		} catch(exception&) {
			throw;
		}
	}
}

void SMETIO::generateHeaderInfo(const StationData& sd, const bool& i_outputIsAscii, const bool& isConsistent,
                                const double& timezone, const size_t& nr_of_parameters,
                                const std::vector<bool>& vecParamInUse, const std::vector<std::string>& vecColumnName,
                                smet::SMETWriter& mywriter)
{
	/**
	 * This procedure sets all relevant information for the header in the SMETWriter object mywriter
	 * The following key/value pairs are set for the header:
	 * - station_id, station_name (if present)
	 * - nodata (set to IOUtils::nodata)
	 * - fields (depending on ASCII/BINARY format and whether the meta data is part of the header or data)
	 * - timezone
	 * - meta data (lat/lon/alt or east/north/alt/epsg if not part of data section)
	 */
	ostringstream ss;

	mywriter.set_header_value("station_id", sd.stationID);
	if (!sd.stationName.empty())
		mywriter.set_header_value("station_name", sd.stationName);
	mywriter.set_header_value("nodata", IOUtils::nodata);

	vector<int> myprecision, mywidth; //set meaningful precision/width for each column

	if (i_outputIsAscii) {
		ss << "timestamp";
	} else {
		ss << "julian";
		myprecision.push_back(8);
		mywidth.push_back(16);
	}

	if (isConsistent) {
		mywriter.set_header_value("latitude", sd.position.getLat());
		mywriter.set_header_value("longitude", sd.position.getLon());
		mywriter.set_header_value("easting", sd.position.getEasting());
		mywriter.set_header_value("northing", sd.position.getNorthing());
		mywriter.set_header_value("altitude", sd.position.getAltitude());
		mywriter.set_header_value("epsg", (double)sd.position.getEPSG());

		if ((timezone != IOUtils::nodata) && (timezone != 0.0))
			mywriter.set_header_value("tz", timezone);
	} else {
		ss << " latitude longitude altitude";
		myprecision.push_back(8); //for latitude
		mywidth.push_back(11);    //for latitude
		myprecision.push_back(8); //for longitude
		mywidth.push_back(11);    //for longitude
		myprecision.push_back(1); //for altitude
		mywidth.push_back(7);     //for altitude
	}

	//Add all other used parameters
	int tmpwidth, tmpprecision;
	for (size_t ll=0; ll<nr_of_parameters; ll++) {
		if (vecParamInUse[ll]) {
			string column = vecColumnName.at(ll);
			if (column == "RSWR") column = "OSWR";
			ss << " " << column;

			getFormatting(ll, tmpprecision, tmpwidth);
			myprecision.push_back(tmpprecision);
			mywidth.push_back(tmpwidth);
		}
	}

	mywriter.set_header_value("fields", ss.str());
	mywriter.set_width(mywidth);
	mywriter.set_precision(myprecision);
}

void SMETIO::getFormatting(const size_t& param, int& prec, int& width)
{
	/**
	 * When writing a SMET file, different meteo parameters require a different
	 * format with regard to precision and width when printing.
	 * This procedure sets the precision and width for each known parameter and
	 * defaults to a width of 8 and precision of 3 digits for each unknown parameter.
	 */
	if ((param == MeteoData::TA) || (param == MeteoData::TSS) || (param == MeteoData::TSG)){
		prec = 2;
		width = 8;
	} else if ((param == MeteoData::VW) || (param == MeteoData::VW_MAX)){
		prec = 1;
		width = 6;
	} else if (param == MeteoData::DW){
		prec = 0;
		width = 5;
	} else if ((param == MeteoData::ISWR) || (param == MeteoData::RSWR) || (param == MeteoData::ILWR)){
		prec = 0;
		width = 6;
	} else if (param == MeteoData::PSUM){
		prec = 3;
		width = 6;
	} else if (param == MeteoData::PSUM_PH){
		prec = 3;
		width = 4;
	} else if (param == MeteoData::HS){
		prec = 3;
		width = 8;
	} else if (param == MeteoData::RH){
		prec = 3;
		width = 7;
	} else {
		prec = 3;
		width = 8;
	}
}

size_t SMETIO::getNrOfParameters(const std::string& stationname, const std::vector<MeteoData>& vecMeteo)
{
	/**
	 * This function loops through all MeteoData objects present in vecMeteo and returns the
	 * number of meteo parameters that the MeteoData objects have. If there is an inconsistency
	 * in the number of meteo parameters in use within the vector of MeteoData then a warning
	 * is printed and MeteoData::nrOfParameters is returned, thus all additional meteo parameters
	 * that might be in use are ignored.
	 */

	size_t actual_nr_of_parameters = IOUtils::npos;
	bool has_one_element = false;

	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		has_one_element = true;
		const size_t current_size = vecMeteo[ii].getNrOfParameters();

		if (actual_nr_of_parameters == IOUtils::npos){
			actual_nr_of_parameters = current_size;
		} else if (actual_nr_of_parameters != current_size){
			//There is an inconsistency in the fields, print out a warning and proceed
			cerr << "[w] While writing SMET file: Inconsistency in number of meteo "
				<< "parameters for station " << stationname << endl;
			actual_nr_of_parameters = MeteoData::nrOfParameters;
			break;
		}
	}

	if (!has_one_element)
		return MeteoData::nrOfParameters;

	return actual_nr_of_parameters;
}

void SMETIO::checkForUsedParameters(const std::vector<MeteoData>& vecMeteo, const size_t& nr_parameters, double& timezone,
                                    std::vector<bool>& vecParamInUse, std::vector<std::string>& vecColumnName)
{
	/**
	 * This procedure loops through all MeteoData objects present in vecMeteo and finds out which
	 * meteo parameters are actually in use, i. e. have at least one value that differs from IOUtils::nodata.
	 * If a parameter is in use, then vecParamInUse[index_of_parameter] is set to true and the column
	 * name is set in vecColumnName[index_of_parameter]
	 */
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		for (size_t jj=0; jj<nr_parameters; jj++) {
			if (!vecParamInUse[jj]) {
				if (vecMeteo[ii](jj) != IOUtils::nodata) {
					vecParamInUse[jj] = true;
					vecColumnName.at(jj) = vecMeteo[ii].getNameForParameter(jj);
				}
			}
		}
	}

	if (!vecMeteo.empty())
		timezone = vecMeteo[0].date.getTimeZone();
}

bool SMETIO::checkConsistency(const std::vector<MeteoData>& vecMeteo, StationData& sd)
{
	/**
	 * This function checks whether all the MeteoData elements in vecMeteo are consistent
	 * regarding their meta data (position information, station name). If they are consistent
	 * true is returned, otherwise false
	 */

	if (!vecMeteo.empty()) //to get the station data even when in bug 87 conditions
		sd = vecMeteo[0].meta;

	for (size_t ii=1; ii<vecMeteo.size(); ii++){
		const Coords& p1 = vecMeteo[ii-1].meta.position;
		const Coords& p2 = vecMeteo[ii].meta.position;
		if (p1 != p2) {
			//we don't mind if p1==nodata or p2==nodata
			if(p1.isNodata()==false && p2.isNodata()==false) return false;
		}
	}

	return true;
}

void SMETIO::readPOI(std::vector<Coords>& pts)
{
	const std::string filename = cfg.get("POIFILE", "Input");
	if (!IOUtils::fileExists(filename)) {
		throw FileNotFoundException(filename, AT);
	}

	smet::SMETReader myreader(filename);
	vector<double> vec_data;
	myreader.read(vec_data);
	const size_t nr_fields = myreader.get_nr_of_fields();
	const int epsg = myreader.get_header_intvalue("epsg");
	const double smet_nodata = myreader.get_header_doublevalue("nodata");

	if (myreader.location_in_data(smet::WGS84)==true) {
		size_t lat_fd=IOUtils::unodata, lon_fd=IOUtils::unodata;
		size_t alt_fd=IOUtils::unodata;
		for (size_t ii=0; ii<nr_fields; ii++) {
			const string tmp = myreader.get_field_name(ii);
			if (tmp=="latitude") lat_fd=ii;
			if (tmp=="longitude") lon_fd=ii;
			if (tmp=="altitude") alt_fd=ii;
		}
		for (size_t ii=0; ii<vec_data.size(); ii+=nr_fields) {
			Coords point;
			point.setLatLon(vec_data[ii+lat_fd], vec_data[ii+lon_fd], vec_data[ii+alt_fd]);
			pts.push_back(point);
		}
	} else if (myreader.location_in_data(smet::EPSG)==true) {
		if (epsg==(int)floor(smet_nodata + 0.1))
			throw InvalidFormatException("In file \""+filename+"\", missing EPSG code in header!", AT);

		size_t east_fd=IOUtils::unodata, north_fd=IOUtils::unodata;
		size_t alt_fd=IOUtils::unodata;
		for (size_t ii=0; ii<nr_fields; ii++) {
			const string tmp = myreader.get_field_name(ii);
			if (tmp=="easting") east_fd=ii;
			if (tmp=="northing") north_fd=ii;
			if (tmp=="altitude") alt_fd=ii;
		}
		if ((east_fd == IOUtils::unodata) || (north_fd == IOUtils::unodata) || (alt_fd == IOUtils::unodata))
			throw InvalidFormatException("File \""+filename+"\" does not contain all data fields necessary for EPSG coordinates", AT);

		for (size_t ii=0; ii<vec_data.size(); ii+=nr_fields) {
			Coords point;
			point.setEPSG(epsg);
			point.setXY(vec_data[ii+east_fd], vec_data[ii+north_fd], vec_data[ii+alt_fd]);
			pts.push_back(point);
		}
	} else {
		throw InvalidFormatException("File \""+filename+"\" does not contain expected location information in DATA section!", AT);
	}
}

void SMETIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SMETIO::write2DGrid(const Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

} //namespace
