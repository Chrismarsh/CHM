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

#include <snowpack/plugins/SmetIO.h>

using namespace std;
using namespace mio;

/**
 * @page smet SMET
 * @section smet_format Format
 * This plugin reads the SMET files as specified in the 
 * <a href="https://models.slf.ch/p/meteoio">MeteoIO</a> pre-processing library documentation (under 
 * <i>"Available plugins and usage"</i>, then <i>"smet"</i>).
 *
 * @section layers_data Layers data
 * The snow/soil layers file has the structure described below:
 * - the SMET signature (to identify the file as SMET as well as the format version)
 * - a Header section containing the metadata for the location
 * - a Data section containing description of the layers (if any). Please note that the layers are given from the bottom to the top.
 *
 * The following points are important to remember:
 * - the station_id field is often used to generate output file names
 * - you don't have to provide lat/lon, you can provide easiting/northing instead alongside an epsg code (see MeteoIO documentation)
 * - the ProfileDate will be used as starting date for the simulation. Therefore, make sure you have meteorological data from this point on!
 * - the number of soil and snow layers <b>must</b> be right!
 * - timestamps follow the ISO format, temperatures are given in Kelvin, thicknesses in m, fractional volumes are between 0 and 1 (and the total sum <b>must</b> be exactly one), densities are in kg/m<SUP>3</SUP> (see the definition of the fields in the table below)
 *
 * <center><table border="0">
 * <caption>initial snow profile fields description</caption>
 * <tr><td>
 * <table border="1">
 * <tr><th>Field</th><th>Description</th></tr>
 * <tr><th>timestamp</th><td>ISO formatted time</td></tr>
 * <tr><th>Layer_Thick</th><td>layer thickness [mm]</td></tr>
 * <tr><th>T</th><td>layer temperature [K]</td></tr>
 * <tr><th>Vol_Frac_I</th><td>fractional ice volume [0-1]</td></tr>
 * <tr><th>Vol_Frac_W</th><td>fractional water volume [0-1]</td></tr>
 * <tr><th>Vol_Frac_V</th><td>fractional voids volume [0-1]</td></tr>
 * <tr><th>Vol_Frac_S</th><td>fractional soil volume [0-1]</td></tr>
 * <tr><th>Rho_S</th><td>soil density [kg/m3]</td></tr>
 * <tr><th>Conduc_S</th><td>soil thermal conductivity [w/(mK)]</td></tr>
 * <tr><th>HeatCapac_S</th><td>soil thermal capacity [J/K]</td></tr>
 * </table></td><td><table border="1">
 * <tr><th>Field</th><th>Description</th></tr>
 * <tr><th>rg</th><td>grain radius [mm]</td></tr>
 * <tr><th>rb</th><td>bond radius [mm]</td></tr>
 * <tr><th>dd</th><td>dendricity [0-1]</td></tr>
 * <tr><th>sp</th><td>spericity [0-1]</td></tr>
 * <tr><th>mk</th><td>marker</td></tr>
 * <tr><th>mass_hoar</th><td>mass of surface hoar []</td></tr>
 * <tr><th>ne</th><td>number of elements</td></tr>
 * <tr><th>CDot</th><td> </td></tr>
 * <tr><th>metamo</th><td> </td></tr>
 * <tr><th> <br></th><td> </td></tr>
 * </table></td></tr>
 * </table></center>
 *
 * Usually, simulations are started at a point in time when no snow is on the ground, therefore not requiring the definition of snow layers. An example is given below with one snow layer:
 * @code
 * SMET 1.1 ASCII
 * [HEADER]
 * station_id       = DAV2
 * station_name     = Davos:Baerentaelli
 * latitude         = 46.701
 * longitude        = 9.82
 * altitude         = 2560
 * nodata           = -999
 * tz               = 1
 * source           = WSL-Institute for Snow and Avalanche Research SLF; CFierz, 2011-10
 * ProfileDate      = 2009-10-01T00:00
 * HS_Last          = 0.0000
 * SlopeAngle       = 38.0
 * SlopeAzi         = 0.0
 * nSoilLayerData   = 0
 * nSnowLayerData   = 1
 * SoilAlbedo       = 0.20
 * BareSoil_z0      = 0.200
 * CanopyHeight     = 0.00
 * CanopyLeafAreaIndex     = 0.00
 * CanopyDirectThroughfall = 1.00
 * WindScalingFactor       = 1.19
 * ErosionLevel     = 0
 * TimeCountDeltaHS = 0.000000
 * fields           = timestamp Layer_Thick  T  Vol_Frac_I  Vol_Frac_W  Vol_Frac_V  Vol_Frac_S Rho_S Conduc_S HeatCapac_S  rg  rb  dd  sp  mk mass_hoar ne CDot metamo
 * [DATA]
 * 2009-09-19T02:30 0.003399 273.15 0.579671 0.068490 0.351839 0.000000 0.0 0.0 0.0 1.432384 1.028390 0.000000 1.000000 22 0.000000 1 0.000000 0.000000
 * @endcode
 *
 * @section hazard_data Hazard data
 * The hazards file contain the temporal history of various parameters that are relevant for avalanche warning (such as three hours new
 * snow fall, etc). If such file is not provided, the internal data structures for such data will be initialized to zero (which is
 * what you usually want when starting a simulation before the start of the snow season). The hazards file has the following structure:
 * @code
 * SMET 1.1 ASCII
 * [HEADER]
 * station_id       = DAV2
 * station_name     = Davos:Baerentaelli
 * latitude         = 46.701
 * longitude        = 9.82
 * altitude         = 2560
 * nodata           = -999
 * tz               = 1
 * ProfileDate      = 2012-06-11T17:30
 * fields           = timestamp SurfaceHoarIndex DriftIndex ThreeHourNewSnow TwentyFourHourNewSnow
 * [DATA]
 * 2010-06-08T18:00       -999       -999   0.000000   0.000000
 * 2010-06-08T18:30       -999       -999   0.000000   0.000000
 * 2010-06-08T19:00       -999       -999   0.000000   0.000000
 * 2010-06-08T19:30       -999       -999   0.000000   0.000000
 * 2010-06-08T20:00       -999       -999   0.000000   0.000000
 * 2010-06-08T20:30       -999       -999   0.000000   0.000000
 * 2010-06-08T21:00       -999       -999   0.000000   0.000000
 * ...
 * 2010-06-11T17:30       -999       -999   0.000000   0.000000
 * @endcode
 * 
 * As can be seen in this example, the various indices as well as the snow statistics are given every half an hour in reverse chronological order until
 * the profile date. 
 */
SmetIO::SmetIO(const SnowpackConfig& cfg, const RunInfo& run_info)
        : outpath(), o_snowpath(), snowpath(), experiment(), inpath(), i_snowpath(), sw_mode(),
          info(run_info),
          in_dflt_TZ(), useSoilLayers(false), perp_to_slope(false)
{
	cfg.getValue("TIME_ZONE", "Input", in_dflt_TZ);
	cfg.getValue("SNP_SOIL", "Snowpack", useSoilLayers);
	cfg.getValue("SW_MODE", "Snowpack", sw_mode);
	cfg.getValue("PERP_TO_SLOPE", "SnowpackAdvanced", perp_to_slope);

	cfg.getValue("EXPERIMENT", "Output", experiment);
	cfg.getValue("METEOPATH", "Output", outpath, IOUtils::nothrow);
	cfg.getValue("SNOWPATH", "Output", snowpath, IOUtils::nothrow);
	if (!snowpath.empty()) {
		o_snowpath = snowpath;
	} else {
		o_snowpath = outpath;
	}

	cfg.getValue("METEOPATH", "Input", inpath, IOUtils::nothrow);
	snowpath = string();
	cfg.getValue("SNOWPATH", "Input", snowpath, IOUtils::nothrow);
	if (!snowpath.empty()) {
		i_snowpath = snowpath;
	} else {
		i_snowpath = inpath;
	}
}

/**
 * @brief This routine checks if the specified snow cover data exists
 * @param i_snowfile file containing the initial state of the snowpack
 * @param stationID
 * @return true if the file exists
 */
bool SmetIO::snowCoverExists(const std::string& i_snowfile, const std::string& /*stationID*/) const
{
	string snofilename = getFilenamePrefix(i_snowfile, i_snowpath, false);

	if (snofilename.rfind(".sno") == string::npos) {
		snofilename += ".sno";
	}

	return IOUtils::fileExists(snofilename);
}

/**
 * @brief This routine reads the status of the snow cover at program start
 * @version 11.02
 * @param i_snowfile file containing the initial state of the snowpack
 * @param stationID
 * @param SSdata
 * @param Zdata
 */
void SmetIO::readSnowCover(const std::string& i_snowfile, const std::string& stationID,
                           SN_SNOWSOIL_DATA& SSdata, ZwischenData& Zdata)
{
	string snofilename = getFilenamePrefix(i_snowfile, i_snowpath, false);
	string hazfilename(snofilename);

	if (snofilename.rfind(".sno") == string::npos) {
		snofilename += ".sno";
		hazfilename += ".haz";
	} else {
		hazfilename.replace(hazfilename.rfind(".sno"), 4, ".haz");
	}

	const Date sno_date = read_snosmet(snofilename, stationID, SSdata);
	if (IOUtils::fileExists(hazfilename)) {
		const Date haz_date = read_hazsmet(hazfilename, Zdata);
		if (haz_date != sno_date)
			throw IOException("Inconsistent ProfileDate in files: " + snofilename + " and " + hazfilename, AT);
	} else {
		//prn_msg(__FILE__, __LINE__, "wrn", Date(), "Hazard file %s does not exist. Initialize Zdata to zero.", hazfilename.c_str());
		Zdata.reset();
	}
}

mio::Date SmetIO::read_hazsmet(const std::string& hazfilename, ZwischenData& Zdata)
{
	/*
	 * Read HAZ SMET file, parse header and fill Zdata with values from the [DATA] section
	 * The header is not parsed apart from the ProfileDate which is needed for a consistency
	 * check with the corresponding SNO SMET file
	 */
	smet::SMETReader haz_reader(hazfilename);
	const unsigned char nr_datalines = 144;

	Date profile_date;
	IOUtils::convertString(profile_date, haz_reader.get_header_value("ProfileDate"),  SmetIO::in_dflt_TZ);
	if (profile_date.isUndef())
			throw ("Invalid ProfileDate in file \""+hazfilename+"\"", AT);

	vector<string> vec_timestamp;
	vector<double> vec_data;
	haz_reader.read(vec_timestamp, vec_data);

	if (vec_timestamp.size() != nr_datalines)
		throw InvalidFormatException("There need to be 144 data lines in " + haz_reader.get_filename(), AT);

	size_t current_index = 0;
	for (size_t ii=0; ii<nr_datalines; ii++) {
		const size_t index = nr_datalines - 1 - ii;
		if (ii>=96){
			Zdata.hoar24[index]  = vec_data[current_index++];
			Zdata.drift24[index] = vec_data[current_index++];
		} else {
			current_index += 2;
		}

		Zdata.hn3[index]  = vec_data[current_index++];
		Zdata.hn24[index] = vec_data[current_index++];
	}

	return profile_date;
}

mio::Date SmetIO::read_snosmet(const std::string& snofilename, const std::string& stationID, SN_SNOWSOIL_DATA& SSdata)
{
	/*
	 * Read SNO SMET file, parse header and fill SSdata with values from the [DATA] section
	*/
	smet::SMETReader sno_reader(snofilename);
	Date profile_date = read_snosmet_header(sno_reader, stationID, SSdata);
	if (profile_date.isUndef())
		throw ("Invalid ProfileDate in file \""+snofilename+"\"", AT);
	profile_date.rnd(1.);

	//Read actual data
	vector<string> vec_timestamp;
	vector<double> vec_data;
	sno_reader.read(vec_timestamp, vec_data);

	if (SSdata.nLayers > 0)
		SSdata.Ldata.resize(SSdata.nLayers, LayerData());

	if (vec_timestamp.size() != SSdata.nLayers)
		throw InvalidFormatException("Xdata: Layers expected != layers read in " + sno_reader.get_filename(), AT);

	const size_t nr_of_fields = sno_reader.get_nr_of_fields();
	const size_t nr_of_solutes = (nr_of_fields - 18) / 4;

	if (SnowStation::number_of_solutes != nr_of_solutes)
		throw InvalidFormatException("Mismatch in number_of_solutes and fields in " + sno_reader.get_filename(), AT);

	//copy data to respective variables in SSdata
	size_t current_index = 0;
	for (size_t ll=0; ll<SSdata.nLayers; ll++) {
		//firstly deal with date
		IOUtils::convertString(SSdata.Ldata[ll].depositionDate, vec_timestamp[ll],  in_dflt_TZ);
		SSdata.Ldata[ll].depositionDate.rnd(1.);

		if (SSdata.Ldata[ll].depositionDate > SSdata.profileDate) {
			prn_msg(__FILE__, __LINE__, "err", Date(),
				   "Layer %d from bottom is younger (%lf) than ProfileDate (%lf) !!!",
				   ll+1, SSdata.Ldata[ll].depositionDate.getJulian(), SSdata.profileDate.getJulian());
			throw IOException("Cannot generate Xdata from file " + sno_reader.get_filename(), AT);
		}

		//secondly with the actual data
		SSdata.Ldata[ll].hl = vec_data[current_index++];
		SSdata.Ldata[ll].tl = vec_data[current_index++];
		SSdata.Ldata[ll].phiIce = vec_data[current_index++];
		SSdata.Ldata[ll].phiWater = vec_data[current_index++];
		SSdata.Ldata[ll].phiVoids = vec_data[current_index++];
		SSdata.Ldata[ll].phiSoil = vec_data[current_index++];

		if (SSdata.Ldata[ll].tl < 100.) {
			SSdata.Ldata[ll].tl = IOUtils::C_TO_K(SSdata.Ldata[ll].tl);
		}
		SSdata.Ldata[ll].SoilRho = vec_data[current_index++];
		SSdata.Ldata[ll].SoilK = vec_data[current_index++];
		SSdata.Ldata[ll].SoilC = vec_data[current_index++];
		SSdata.Ldata[ll].rg = vec_data[current_index++];
		SSdata.Ldata[ll].rb = vec_data[current_index++];
		SSdata.Ldata[ll].dd = vec_data[current_index++];
		SSdata.Ldata[ll].sp = vec_data[current_index++];
		SSdata.Ldata[ll].mk = static_cast<unsigned short int>(vec_data[current_index++]+.5); //int
		SSdata.Ldata[ll].hr = vec_data[current_index++];
		SSdata.Ldata[ll].ne = static_cast<unsigned int>(vec_data[current_index++]+.5); //int

		if (SSdata.Ldata[ll].rg>0. && SSdata.Ldata[ll].rb >= SSdata.Ldata[ll].rg) {
			//HACK To avoid surprises in lwsn_ConcaveNeckRadius()
			SSdata.Ldata[ll].rb = Metamorphism::max_grain_bond_ratio * SSdata.Ldata[ll].rg;
			prn_msg(__FILE__, __LINE__, "wrn", Date(), "Layer %d from bottom: bond radius rb/rg larger than Metamorphism::max_grain_bond_ratio=%lf (rb=%lf mm, rg=%lf mm)! Reset to Metamorphism::max_grain_bond_ratio", ll+1, Metamorphism::max_grain_bond_ratio, SSdata.Ldata[ll].rb, SSdata.Ldata[ll].rg);
		}

		SSdata.Ldata[ll].CDot = vec_data[current_index++];
		SSdata.Ldata[ll].metamo = vec_data[current_index++];

		for (size_t ii=0; ii<SnowStation::number_of_solutes; ii++) {
			SSdata.Ldata[ll].cIce[ii] = vec_data[current_index++];
			SSdata.Ldata[ll].cWater[ii] = vec_data[current_index++];
			SSdata.Ldata[ll].cVoids[ii] = vec_data[current_index++];
			SSdata.Ldata[ll].cSoil[ii] = vec_data[current_index++];
		}
	} //for loop over layers

	SSdata.nN = 1;
	SSdata.Height = 0.;
	for (size_t ll = 0; ll < SSdata.nLayers; ll++) {
		SSdata.nN    += SSdata.Ldata[ll].ne;
		SSdata.Height += SSdata.Ldata[ll].hl;
	}

	return profile_date;
}

mio::Date SmetIO::read_snosmet_header(const smet::SMETReader& sno_reader, const std::string& stationID,
                                      SN_SNOWSOIL_DATA& SSdata)
{
	/*
	 * Read values for certain header keys (integer and double values) and perform
	 * consistency checks upon them.
	 */
	const string station_name = sno_reader.get_header_value("station_name");
	IOUtils::convertString(SSdata.profileDate, sno_reader.get_header_value("ProfileDate"),  in_dflt_TZ);

	SSdata.HS_last = get_doubleval(sno_reader, "HS_Last");

	mio::Coords tmppos;
	const double alt = get_doubleval(sno_reader, "altitude");
	if (keyExists(sno_reader, "latitude")) {
		const double lat = get_doubleval(sno_reader, "latitude");
		const double lon = get_doubleval(sno_reader, "longitude");
		tmppos.setLatLon(lat, lon, alt);
	} else {
		const double easting = get_doubleval(sno_reader, "easting");
		const double northing = get_doubleval(sno_reader, "northing");
		const int epsg = get_intval(sno_reader, "epsg");
		tmppos.setEPSG(epsg);
		tmppos.setXY(easting, northing, alt);
	}

	const double slope_angle = get_doubleval(sno_reader, "SlopeAngle");
	const double azi = get_doubleval(sno_reader, "SlopeAzi");
	SSdata.meta.setStationData(tmppos, stationID, station_name);
	SSdata.meta.setSlope(slope_angle, azi);

	// Check consistency with radiation switch
	if ((sw_mode == "BOTH") && perp_to_slope && (SSdata.meta.getSlopeAngle() > Constants::min_slope_angle)) {
		prn_msg(__FILE__, __LINE__, "wrn", Date(),
		        "You want to use measured albedo in a slope steeper than 3 deg  with PERP_TO_SLOPE set!");
		throw IOException("Do not generate Xdata from file " + sno_reader.get_filename(), AT);
	}

	const int nSoilLayerData = get_intval(sno_reader, "nSoilLayerData");
	if (nSoilLayerData < 0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "'nSoilLayerData' < 0 !!!");
		throw InvalidFormatException("Cannot generate Xdata from file " + sno_reader.get_filename(), AT);
	} else if (useSoilLayers && (nSoilLayerData < 1)) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "useSoilLayers set but 'nSoilLayerData' < 1 !!!");
		throw InvalidFormatException("Cannot generate Xdata from file " + sno_reader.get_filename(), AT);
	} else if (!useSoilLayers && (nSoilLayerData > 0)) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "useSoilLayers not set but 'nSoilLayerData' > 0 !!!");
		throw InvalidFormatException("Cannot generate Xdata from file " + sno_reader.get_filename(), AT);
	}
	SSdata.nLayers = (unsigned int) nSoilLayerData;

	const int nSnowLayerData = get_intval(sno_reader, "nSnowLayerData");
	if (nSnowLayerData < 0) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "'nSnowLayerData' < 0  !!!");
		throw InvalidFormatException("Cannot generate Xdata from file " + sno_reader.get_filename(), AT);
	}
	SSdata.nLayers += (unsigned int)nSnowLayerData;

	SSdata.SoilAlb = get_doubleval(sno_reader, "SoilAlbedo");
	SSdata.BareSoil_z0 = get_doubleval(sno_reader, "BareSoil_z0");
	if (SSdata.BareSoil_z0 == 0.) {
		prn_msg(__FILE__, __LINE__, "wrn", Date(), "'BareSoil_z0'=0 from %s, reset to 0.02",
			   sno_reader.get_filename().c_str());
		SSdata.BareSoil_z0 = 0.02;
	}
	if (SSdata.HS_last > 0.05) {
		SSdata.Albedo = 0.9;
	} else {
		SSdata.Albedo = SSdata.SoilAlb;
	}

	SSdata.Canopy_Height = get_doubleval(sno_reader, "CanopyHeight");
	SSdata.Canopy_LAI = get_doubleval(sno_reader, "CanopyLeafAreaIndex");
	SSdata.Canopy_Direct_Throughfall = get_doubleval(sno_reader, "CanopyDirectThroughfall");
	SSdata.WindScalingFactor = get_doubleval(sno_reader, "WindScalingFactor");

	SSdata.ErosionLevel = get_intval(sno_reader, "ErosionLevel");
	SSdata.TimeCountDeltaHS = get_doubleval(sno_reader, "TimeCountDeltaHS");

	return SSdata.profileDate;
}

bool SmetIO::keyExists(const smet::SMETReader& reader, const std::string& key) const
{
	const double nodata = reader.get_header_doublevalue("nodata");
	const double value = reader.get_header_doublevalue(key);

	return value!=nodata;
}

double SmetIO::get_doubleval(const smet::SMETReader& reader, const std::string& key) const
{
	// Retrieve a double value from a SMETReader object header and make sure it exists.
	// If the header key does not exist or the value is not set throw an exception
	const double nodata = reader.get_header_doublevalue("nodata");
	const double value = reader.get_header_doublevalue(key);

	if (value == nodata){
		const string msg = "Missing key '" + key + "'";
		prn_msg(__FILE__, __LINE__, "err", Date(), msg.c_str());
		throw InvalidFormatException("Cannot generate Xdata from file " + reader.get_filename(), AT);
	}

	return value;
}

int SmetIO::get_intval(const smet::SMETReader& reader, const std::string& key) const
{
	// Retrieve an integer value from a SMETReader object header and make sure it exists.
	// If the header key does not exist or the value is not set throw an exception
	const double nodata = reader.get_header_doublevalue("nodata");
	const int inodata = static_cast<int>( floor(nodata + .1) );

	const int value = reader.get_header_intvalue(key);
	if (value == inodata){
		const string msg = "Missing key '" + key + "'";
		prn_msg(__FILE__, __LINE__, "err", Date(), msg.c_str());
		throw InvalidFormatException("Cannot generate Xdata from file " + reader.get_filename(), AT);
	}

	return value;
}

/**
 * @brief This routine writes the status of the snow cover at program termination and at specified backup times
 * @version 11.02
 * @param date current
 * @param Xdata
 * @param Zdata
 * @param forbackup dump Xdata on the go
 */
void SmetIO::writeSnowCover(const mio::Date& date, const SnowStation& Xdata,
                            const ZwischenData& Zdata, const bool& forbackup)
{
	string snofilename = getFilenamePrefix(Xdata.meta.getStationID().c_str(), o_snowpath) + ".sno";
	string hazfilename = getFilenamePrefix(Xdata.meta.getStationID().c_str(), o_snowpath) + ".haz";

	if (forbackup){
		stringstream ss;
		ss << (int)(date.getJulian() + 0.5);
		snofilename += ss.str();
		hazfilename += ss.str();
	}

	writeSnoFile(snofilename, date, Xdata, Zdata);
	writeHazFile(hazfilename, date, Xdata, Zdata);
}

void SmetIO::writeHazFile(const std::string& hazfilename, const mio::Date& date, const SnowStation& Xdata,
                          const ZwischenData& Zdata)
{
	/*
	 * This procedure creates a SMETWriter object, sets its header and copies all required
	 * data and timestamps into vec_timestamp and vec_data (copied from Zdata).
	 * The SMETWriter object finally writes out the HAZ SMET file
	 */
	vector<string> vec_timestamp;
	vector<double> vec_data;

	smet::SMETWriter haz_writer(hazfilename);
	setBasicHeader(Xdata, "timestamp SurfaceHoarIndex DriftIndex ThreeHourNewSnow TwentyFourHourNewSnow", haz_writer);
	haz_writer.set_header_value("ProfileDate", date.toString(Date::ISO));

	haz_writer.set_width(vector<int>(4,10));
	haz_writer.set_precision(vector<int>(4,6));

	Date hrs72(date - Date(3.0,0.0));
	const Duration half_hour(1.0/48.0,0.0);

	for (size_t ii=0; ii<144; ii++){
		//hn3 and hn24 have 144 fields, the newest values have index 0
		//hoar24 and drift 24 have 48 fields, the newest values have index 0
		hrs72 += half_hour;
		vec_timestamp.push_back(hrs72.toString(Date::ISO));
		if (ii >= 96){
			const size_t index = 143-ii;
			vec_data.push_back(Zdata.hoar24[index]);  //Print out the hoar hazard data info
			vec_data.push_back(Zdata.drift24[index]); //Print out the drift hazard data info
		} else {
			vec_data.push_back(IOUtils::nodata);
			vec_data.push_back(IOUtils::nodata);
		}

		vec_data.push_back(Zdata.hn3[143-ii]);  //Print out the 3 hour new snowfall hazard data info
		vec_data.push_back(Zdata.hn24[143-ii]); //Print out the 24 hour new snowfall hazard data info
	}
	haz_writer.write(vec_timestamp, vec_data);
}

void SmetIO::writeSnoFile(const std::string& snofilename, const mio::Date& date, const SnowStation& Xdata,
                          const ZwischenData& /*Zdata*/)
{
	/*
	 * This procedure creates a SMETWriter object, sets its header and copies all required
	 * data and timestamps into vec_timestamp and vec_data (from Xdata).
	 * The SMETWriter object finally writes out the SNO SMET file
	 */
	smet::SMETWriter sno_writer(snofilename);
	stringstream ss;
	ss << "timestamp Layer_Thick  T  Vol_Frac_I  Vol_Frac_W  Vol_Frac_V  Vol_Frac_S Rho_S " //8
	   << "Conduc_S HeatCapac_S  rg  rb  dd  sp  mk mass_hoar ne CDot metamo";
	for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
		ss << " cIce cWater cAir  cSoil";
	}

	setBasicHeader(Xdata, ss.str(), sno_writer);
	setSnoSmetHeader(Xdata, date, sno_writer);

	vector<string> vec_timestamp;
	vector<double> vec_data;
	vector<int> vec_width, vec_precision;
	setFormatting(Xdata.number_of_solutes, vec_width, vec_precision);
	sno_writer.set_width(vec_width);
	sno_writer.set_precision(vec_precision);

	//Fill vec_data with values
	const vector<ElementData>& EMS = Xdata.Edata;
	for (size_t e = 0; e < Xdata.getNumberOfElements(); e++) {
		vec_timestamp.push_back(EMS[e].depositionDate.toString(Date::ISO));

		vec_data.push_back(EMS[e].L);
		vec_data.push_back(Xdata.Ndata[e+1].T);
		vec_data.push_back(EMS[e].theta[ICE]);
		vec_data.push_back(EMS[e].theta[WATER]);
		vec_data.push_back(EMS[e].theta[AIR]);
		vec_data.push_back(EMS[e].theta[SOIL]);
		vec_data.push_back(EMS[e].soil[SOIL_RHO]);
		vec_data.push_back(EMS[e].soil[SOIL_K]);
		vec_data.push_back(EMS[e].soil[SOIL_C]);
		vec_data.push_back(EMS[e].rg);
		vec_data.push_back(EMS[e].rb);
		vec_data.push_back(EMS[e].dd);
		vec_data.push_back(EMS[e].sp);
		vec_data.push_back(static_cast<double>(EMS[e].mk));

		vec_data.push_back(Xdata.Ndata[e+1].hoar);
		vec_data.push_back(1.);
		vec_data.push_back(EMS[e].CDot);
		vec_data.push_back(EMS[e].metamo);

		for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
			vec_data.push_back(EMS[e].conc(ICE,ii));
			vec_data.push_back(EMS[e].conc(WATER,ii));
			vec_data.push_back(EMS[e].conc(AIR,ii));
			vec_data.push_back(EMS[e].conc(SOIL,ii));
		}
	}

	sno_writer.write(vec_timestamp, vec_data);
}

void SmetIO::setBasicHeader(const SnowStation& Xdata, const std::string& fields, smet::SMETWriter& smet_writer)
{
	/*
	 * Set the basic, mandatory header key/value pairs for a SMET file
	 */
	smet_writer.set_header_value("station_id", Xdata.meta.getStationID());
	smet_writer.set_header_value("station_name", Xdata.meta.getStationName());
	smet_writer.set_header_value("nodata", IOUtils::nodata);
	smet_writer.set_header_value("fields", fields);

	// Latitude, Longitude, Altitude
	smet_writer.set_header_value("latitude", Xdata.meta.position.getLat());
	smet_writer.set_header_value("longitude", Xdata.meta.position.getLon());
	smet_writer.set_header_value("altitude", Xdata.meta.position.getAltitude());
}

void SmetIO::setSnoSmetHeader(const SnowStation& Xdata, const Date& date, smet::SMETWriter& smet_writer)
{
	/*
	 * The non-compulsory header key/value pairs for SNO files are set in this procedure
	 * The sequence in which they are handed to smet_writer will be preserved when
	 * the SMETWriter actually writes the header
	 */
	stringstream ss; //we use the stringstream to produce strings in desired format

	smet_writer.set_header_value("ProfileDate", date.toString(Date::ISO));

	// Last checked calculated snow depth used for albedo control of next run
	ss.str(""); ss << fixed << setprecision(6) << (Xdata.cH - Xdata.Ground);
	smet_writer.set_header_value("HS_Last", ss.str());

	// Latitude, Longitude, Altitude NOTE:redundant?, Slope Angle, Slope Azimut
	smet_writer.set_header_value("latitude", Xdata.meta.position.getLat());
	smet_writer.set_header_value("longitude", Xdata.meta.position.getLon());
	smet_writer.set_header_value("altitude", Xdata.meta.position.getAltitude());
	ss.str(""); ss << fixed << setprecision(2) << Xdata.meta.getSlopeAngle();
	smet_writer.set_header_value("SlopeAngle", ss.str());
	ss.str(""); ss << fixed << setprecision(2) << Xdata.meta.getAzimuth();
	smet_writer.set_header_value("SlopeAzi", ss.str());

	// Number of Soil Layer Data; in case of no soil used to store the erosion level
	smet_writer.set_header_value("nSoilLayerData", static_cast<double>(Xdata.SoilNode));
	// Number of Snow Layer Data
	smet_writer.set_header_value("nSnowLayerData", static_cast<double>(Xdata.getNumberOfElements() - Xdata.SoilNode));

	// Ground Characteristics (introduced June 2006)
	ss.str(""); ss << fixed << setprecision(2) << Xdata.SoilAlb;
	smet_writer.set_header_value("SoilAlbedo", ss.str());
	ss.str(""); ss << fixed << setprecision(3) << Xdata.BareSoil_z0;
	smet_writer.set_header_value("BareSoil_z0", ss.str());

	// Canopy Characteristics
	ss.str(""); ss << fixed << setprecision(2) << Xdata.Cdata.height;
	smet_writer.set_header_value("CanopyHeight", ss.str());
	ss.str(""); ss << fixed << setprecision(6) << Xdata.Cdata.lai;
	smet_writer.set_header_value("CanopyLeafAreaIndex", ss.str());
	ss.str(""); ss << fixed << setprecision(2) << Xdata.Cdata.direct_throughfall;
	smet_writer.set_header_value("CanopyDirectThroughfall", ss.str());

	// Additional parameters
	ss.str(""); ss << fixed << setprecision(2) << Xdata.WindScalingFactor;
	smet_writer.set_header_value("WindScalingFactor", ss.str());
	smet_writer.set_header_value("ErosionLevel", static_cast<double>(Xdata.ErosionLevel));
	ss.str(""); ss << fixed << setprecision(6) << Xdata.TimeCountDeltaHS;
	smet_writer.set_header_value("TimeCountDeltaHS", ss.str());
}

void SmetIO::setFormatting(const size_t& nr_solutes,
                           std::vector<int>& vec_width, std::vector<int>&  vec_precision)
{
	/*
	 * When writing a SNOW SMET file each written parameter may have a different
	 * column width and precision. This procedure sets the vectors vec_width
	 * and vec_precision which are subsequently handed to a SMETWriter object
	 * It is paramount that the number of fields (not counting the timestamp)
	 * in the SMET file corresponds to the number of elements in both
	 * vec_width and vec_precision
	 */
	vec_width.clear();
	vec_precision.clear();

	vec_width.push_back(12); vec_precision.push_back(6); //EMS[e].L
	vec_width.push_back(12); vec_precision.push_back(6); //Xdata.Ndata[e+1].T
	vec_width.push_back(12); vec_precision.push_back(6); //EMS[e].theta[ICE]
	vec_width.push_back(12); vec_precision.push_back(6); //EMS[e].theta[WATER]
	vec_width.push_back(12); vec_precision.push_back(6); //EMS[e].theta[AIR]
	vec_width.push_back(12); vec_precision.push_back(6); //EMS[e].theta[SOIL]
	vec_width.push_back(9); vec_precision.push_back(1);  //EMS[e].soil[SOIL_RHO]
	vec_width.push_back(9); vec_precision.push_back(1);  //EMS[e].soil[SOIL_K]
	vec_width.push_back(12); vec_precision.push_back(1);  //EMS[e].soil[SOIL_C]
	vec_width.push_back(11); vec_precision.push_back(6);  //EMS[e].rg
	vec_width.push_back(10); vec_precision.push_back(6);  //EMS[e].rb
	vec_width.push_back(10); vec_precision.push_back(6);  //EMS[e].dd
	vec_width.push_back(10); vec_precision.push_back(6);  //EMS[e].sp
	vec_width.push_back(6); vec_precision.push_back(0);  //EMS[e].mk

	vec_width.push_back(13); vec_precision.push_back(6); //Xdata.Ndata[e+1].hoar
	vec_width.push_back(4); vec_precision.push_back(0);  //ne
	vec_width.push_back(15); vec_precision.push_back(6); //EMS[e].CDot
	vec_width.push_back(15); vec_precision.push_back(6); //EMS[e].metamo

	for (size_t ii = 0; ii < nr_solutes; ii++) {
		vec_width.push_back(17); vec_precision.push_back(6); //EMS[e].conc(ICE,ii)
		vec_width.push_back(18); vec_precision.push_back(7); //EMS[e].conc(WATER,ii)
		vec_width.push_back(18); vec_precision.push_back(7); //EMS[e].conc(AIR,ii)
		vec_width.push_back(18); vec_precision.push_back(7); //EMS[e].conc(SOIL,ii)
	}
}


void SmetIO::writeTimeSeries(const SnowStation& /*Xdata*/, const SurfaceFluxes& /*Sdata*/, const CurrentMeteo& /*Mdata*/,
                               const ProcessDat& /*Hdata*/, const double /*wind_trans24*/)
{
	throw IOException("Nothing implemented here!", AT);
}

void SmetIO::writeProfile(const mio::Date& /*date*/, const SnowStation& /*Xdata*/)
{
	throw IOException("Nothing implemented here!", AT);
}

bool SmetIO::writeHazardData(const std::string& /*stationID*/, const std::vector<ProcessDat>& /*Hdata*/,
                             const std::vector<ProcessInd>& /*Hdata_ind*/, const size_t& /*num*/)
{
	throw IOException("Nothing implemented here!", AT);
}

// complete filename_prefix
std::string SmetIO::getFilenamePrefix(const std::string& fnam, const std::string& path, const bool addexp) const
{
	//TODO: read only once (in constructor)
	string filename_prefix = path + "/" + fnam;

	if (addexp && (experiment != "NO_EXP")) //NOTE usually, experiment == NO_EXP in operational mode
		filename_prefix += "_" + experiment;

	return filename_prefix;
}
