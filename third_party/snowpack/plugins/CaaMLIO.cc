/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of Snowpack.
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
#include <snowpack/plugins/CaaMLIO.h>
#include <snowpack/Utils.h>
//#include <meteoio/meteolaws/Atmosphere.h>

#include <sstream>
#include <fstream>
#include <iostream>

#include <sys/time.h>

using namespace std;
using namespace mio;

/**
 * @page caaml CAAML
 * @section caaml_format Format
 * This plugin reads the CAAML files as generated according <A HREF="http://caaml.org/">CAAML V6.0.3</A>'s
 * <A HREF="http://caaml.org/Schemas/SnowProfileIACS/v6.0.3">specification</A>. In order to validate
 * a CAAML file, download the <A HREF="http://caaml.org/Schemas/SnowProfileIACS/v6.0.3/CAAMLv6_SnowProfileIACS.xsd">xsd</A>
 * file and use either an online XML validator or an offline tool (such as Notepad++'s XML tool).
 *
 * @section caaml_keywords Keywords
 * This plugin uses the following keywords (all specified in the [Input] section):
 * - COORDSYS:  input coordinate system (see Coords)
 * - SNOW:     specify CAAML to read in a caaml file
 * - SNOWPATH: string containing the path to the caaml files to be read
 * - SNOWFILE: specify the caaml file to read the data from
 * - XML_ENCODING: force the input file encoding, overriding the file's own encoding declaration (optional, see \ref caaml_encoding "XML encoding" below)
 * - CAAML_MAX_ELEMENT_THICKNESS: if set (and non-zero), the thickness of the elements will be set to this value, otherwise each element will correspond to
 *                                one stratigraphic layer. Recommendation: set this value to 0.002 (= 2 mm)
 * - CAAML_WRITEOUT_AS_READIN: if set to true, a caaml will be written just after reading in, to check if the reading of the caaml was correct.
 *
 * @section caaml_reading Reading a caaml-file
 * Data which is important for a snowpack-simulation but usually not given in a caaml-file (like dendricity for example), can be included in a caaml-file
 * as snowpack-custom-data with the prefix "snp". A caaml-file written out by this plugin will contain this data. However, if this data is not available,
 * the corresponding values will be estimated or set to default values:
 * - snowpack-custom-data for the whole snowpack (like Albedo, WindScalingFactor,...) will be set to default values.
 * - layer-custom-data (like dendricity, sphericity, maker,...) will be estimated from primary grain form, grain size and temperature.
 * - the formation time of a layer will be estimated depending on snow type and layer location.
 *
 * The liquid water content (lwc) of a layer will be read in from the lwc-profile. If there is no lwc-profile in the caaml-file, the lwc will be estimated
 * from the wetness-codes (D, M, W,...) given in the stratigraphic profile.
 *
 * Besides the lwc-profile the density- and temperature-profile are read in from the caaml-file. The lwc-profile is optional, but the temperature- and
 * density-profile have to be given in the caaml-file.
 *
 * <b>Consistency checks</b>
 *
 * When reading in a caaml-file, the data will be checked for consistency. If there is an inconsistency a warning will be printed and values will be adjusted,
 * if possible. Otherwise an exception is thrown.
 *
 * Warnings:
 * - total snow height given in <caaml:snowPackCond> is different from the sum of the layer-thicknesses. The snow height given in <caaml:snowPackCond> will be
 *   ignored.
 * - grain size of a surface-hoar-layer differs from the surface-hoar-layer-thickness. Adjustment: Grain size will be set to layer thickness.
 * - temperature of a layer is above 0&deg;C. Adjustment: Set temperature to 0&deg;C.
 * - liquid water content (lwc) is greater than 0 and temperature is below 0&deg;C. Adjustment: Set lwc to 0.
 * - grain form is "FC" and grain size is above 1.6 mm. Adjustment: nothing.
 * - grain form is "DH" and grain size is below 1.4 mm. Adjustment: nothing.
 *
 * Exceptions:
 * - grain size is 0 and grain form is not "IF"
 * - slope angle is > 0&deg; and no azimuth is given. (If the slope angle is not given, the slope angle is set to 0&deg;.)
 * - missing data / wrong camml-version ( != 6.0.3)
 * - wrong syntax in caaml-file (pointy brackets, matching quotes,...)
 *
 * @subsection caaml_encoding XML encoding
 * Each XML document should specify its encoding. However this information might sometimes be missing or even worse, be false. This makes the XML document non-compliant.
 * Normally, CAAML reads the file encoding in the file itself. If this does not work (one of the two cases given above), it is possible to force the
 * encoding of the input file by using the "XML_ENCODING" option. This option takes one of the following values
 * ("LE" stands for "Little Endian" and "BE" for "Big Endian"):
 *  UTF-8, UTF-16, UTF-16-LE, UTF-16-BE, UTF-32, UTF-32-LE, UTF-32-BE, LATIN1, ISO-8859-1, WCHAR
 *
 * @section caaml_writing Writing a caaml-file
 * This is an example of a caaml-file written with this plugin:
 * @code
 *<?xml version="1.0" encoding="UTF-8"?>
 * <caaml:SnowProfile gml:id="SLF_WFJ2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:gml="http://www.opengis.net/gml" xmlns:caaml="http://caaml.org/Schemas/SnowProfileIACS/v6.0.3" xmlns:snp="http://www.slf.ch/snowpack/1.0" xmlns:slf="http://www.slf.ch/snowprofile/1.0">
 *   <caaml:timeRef>
 *     <caaml:recordTime>
 *       <caaml:TimeInstant>
 *         <caaml:timePosition>2008-05-16T09:00:00+02:00</caaml:timePosition>
 *       </caaml:TimeInstant>
 *     </caaml:recordTime>
 *     <caaml:dateTimeReport>2018-12-18T17:38:41+01:00</caaml:dateTimeReport>
 *   </caaml:timeRef>
 *   <caaml:srcRef>
 *     <caaml:Operation gml:id="OPERATION_ID">
 *       <caaml:name>SNOWPACK v3.45 compiled on Dec 10 2018, 15:45:05. user: theile</caaml:name>
 *     </caaml:Operation>
 *   </caaml:srcRef>
 *...
 * @endcode
 * Furthermore some data about the profile location and following profiles will be written:
 * - stratigraphic profile
 * - temperature profile
 * - density profile
 * - lwc profile
 * - specific surface area profile
 * - strength profile
 *
 * @section caaml_example Example
 * @code
 * [Input]
 * COORDSYS	= CH1903
 * SNOW	= CAAML
 * SNOWPATH	= ./input/snowCAAMLdata
 * SNOWFILE	= 5WJ_20120229.caaml
 * CAAML_MAX_ELEMENT_THICKNESS = 0.002
 * @endcode
 *
 */

/**
 * @class DataResampler
 * @brief Resample data.
 * The data (double) is given as a x-vector xVec and y-Vector yVec. The new sampling can be defined by an interval dx or by an arbitrary
 * x-vector xVecGoal. Two methods are implemented one for point data and one for range data.
 * For point data the new values are interpolated and extrapolated.
 * For range data the new range-values are averaged over the old ranges.
 * The size of a range (=thickness of a layer) is defined by the difference: xVec[i+1]-xVec[i] and for the last range: xMax-xVec[i]
 * Assumptions:
 * - xVec[i] < xVec[i+1]
 * - xVecGoal[i] < xVecGoal[i+1]
 *
 * Additional assumptions for the range data:
 * - xVec[0] >= 0
 * - xVecGoal[0] >= 0
 * - the last value of xVec and xVecGoal has to be smaller (not equal!) than xMax
 *
 * Examples:
 * Point-data example with rotation:
 * Input:                    xVec:{10, 50, 184}  yVec:{273.15, 263.15, 253.15} xVecGoal:{0, 84} xMax=184
 * Rotation result:          xVec:{0, 134, 174} yVec:{253.15, 263.15, 273.15}
 * Final (resampled) result: xVec:{0,84}        yVec:{253.15, 259.42}
 *
 * Range-data example with rotation:
 * Input:                    xVec:{0, 10, 50}   yVec:{300,400,400} xVecGoal:{0,84} xMax=184
 * Rotation result: 		 xVec:{0, 134, 174} yVec:{400,400,300}
 * Final (resampled) result: xVec:{0,84}        yVec:{400,390}
 *
 * @author Thiemo Theile
 * @date   2018
 */
class DataResampler {
	public:
		DataResampler(const std::vector<double>& xVecIn, const std::vector<double>& yVecIn, const double dxIn, const double xMaxIn,
					  const bool isRangeMeasurement, const bool changeDirectionIn);
		DataResampler(const std::vector<double>& xVecIn, const std::vector<double>& yVecIn, const std::vector<double>& xVecGoalIn,
					  const double xMaxIn, const bool isRangeMeasurement, const bool changeDirectionIn);
		std::vector<double> getResampledXVec();
		std::vector<double> getResampledYVec();
		void printProfile(const std::string message);
		bool resample();

	private:
		bool checkIfDataIsValid();
		void resamplePointValues();
		void resampleRangeValues();
		void createXVecGoalFromDx(const double dx);
		void changeDirectionOfXAxis();

		std::vector<double> xVec;
		std::vector<double> yVec;
		std::vector<double> xVecGoal;
		const double xMax;
		const bool useMethodForRangeValues; //if true use rangeValueMethod, otherwise use pointValueMethod
		const bool changeDirection;         //if true change the direction of the x-Axis before resampling
};

DataResampler::DataResampler(const std::vector<double>& xVecIn, const std::vector<double>& yVecIn, const double dxIn, const double xMaxIn,
							const bool isRangeMeasurement, const bool changeDirectionIn)
		: xVec(xVecIn),yVec(yVecIn),xVecGoal(),xMax(xMaxIn),useMethodForRangeValues(isRangeMeasurement),
		  changeDirection(changeDirectionIn)
{
	createXVecGoalFromDx(dxIn);
}

DataResampler::DataResampler(const std::vector<double>& xVecIn, const std::vector<double>& yVecIn, const std::vector<double>& xVecGoalIn,
							 double xMaxIn, const bool isRangeMeasurement, const bool changeDirectionIn)
		: xVec(xVecIn),yVec(yVecIn),xVecGoal(xVecGoalIn),xMax(xMaxIn),
		  useMethodForRangeValues(isRangeMeasurement), changeDirection(changeDirectionIn)
{
}

void DataResampler::printProfile(const std::string message="profile:")
{
	std::cout << message << std::endl;
	std::cout << "Depth:    Value:   " << std::endl;
	if(xVec.size() == yVec.size()){
		for (size_t ii=0;ii<xVec.size();ii++){
			std::cout << xVec[ii] << "   " << yVec[ii] << std::endl;
		}
	}
	std::cout << std::endl;
}

std::vector<double> DataResampler::getResampledXVec()
{
	return xVec;
}

std::vector<double> DataResampler::getResampledYVec()
{
	return yVec;
}

bool DataResampler::checkIfDataIsValid()
{
	//check if the assumptions for the data are fulfilled
	if(xVec.size() != yVec.size()){
		return false;
	}
	if(xVec.size() < 1){
		return false;
	}
	if(xVecGoal.size() <1){
		return false;
	}
	for(size_t ii=0; ii<xVec.size()-1; ii++){
		if(xVec[ii] >= xVec[ii+1]){
			return false;
		}
	}
	for(size_t ii=0; ii<xVecGoal.size()-1; ii++){
		if(xVecGoal[ii] >= xVecGoal[ii+1]){
			return false;
		}
	}
	if(useMethodForRangeValues){
		if(xVecGoal[xVecGoal.size()-1] >= xMax){
			return false;
		}
		if(xVec[xVec.size()-1] >= xMax){
			return false;
		}
		if(xVecGoal[0] < 0){
			return false;
		}
		if(xVec[0] < 0){
			return false;
		}
	}
	return true;
}

bool DataResampler::resample()
{
	const bool dataIsValid = checkIfDataIsValid();
	if (dataIsValid==false) return false;

	//printProfile("original profile");
	if( changeDirection ){
		changeDirectionOfXAxis();
		//printProfile("rotated profile");
	}
	if(useMethodForRangeValues){
		resampleRangeValues();
	}else{
		resamplePointValues();
	}
	//printProfile("resampled profile");
	return true;
}

void DataResampler::resamplePointValues()
{
	std::vector<double> xVecResampled;
	std::vector<double> yVecResampled;
	size_t index=1;
	for (size_t indexGoal=0; indexGoal < xVecGoal.size(); indexGoal++){
		const double x=xVecGoal[indexGoal];
		//1. search for index so that xVec[index-1]<x<=xVec[index]
		while (xVec[index] < x && index < xVec.size()-1){
			index++;
		}
		//2. interpolate new y-value linearly
		//Example: x=5 xVec[index-1]=4 xVec[index]=6 yVec[index-1]=20 yVec[index]=10
		//         => y = 20 +(10-20)/(6-4)*(5-4) = 15
		const double y=yVec[index-1]+(yVec[index]-yVec[index-1])/(xVec[index]-xVec[index-1])*(x-xVec[index-1]);
		xVecResampled.push_back(x);
		yVecResampled.push_back(y);
	}
	xVec=xVecResampled;
	yVec=yVecResampled;
}

void DataResampler::resampleRangeValues()
{
	//append some values to the vectors for correct treatment of the last value:
	xVec.push_back(xMax);
	yVec.push_back(yVec[yVec.size()-1]);
	xVecGoal.push_back(xMax);
	std::vector<double> xVecResampled;
	std::vector<double> yVecResampled;

	size_t index=1;
	for (size_t indexGoal=0; indexGoal < xVecGoal.size()-1; indexGoal++){
		const double x=xVecGoal[indexGoal];
		const double dx=xVecGoal[indexGoal+1]-x;
		//1. search for index so that x<xVec[index]<xVec[index+k]<x+dx
		while (xVec[index]<x && index < xVec.size()-1){
			index++;
		}
		size_t k=0;
		while (xVec[index+k]<x+dx && index+k < xVec.size()-1){
			k++;
		}
		//2. calculate the value for the new range x:x+dx by weighted averaging of xVec[index:index+k]
		//example: x=0, dx=10, xVec=[0,5,15,...], yVec=[10,20,100,...] will give y = 10*(5-0)/10 + 20*(10-5)/10 = 15
		double y=0;
		double xiiPrevious=x;
		for (size_t ii=0; ii<=k; ii++){
			double xii = xVec[index+ii];
			double yi = yVec[index-1+ii];
			if(xii>x+dx){
				xii=x+dx;
			}
			double weight = (xii-xiiPrevious)/dx;
			xiiPrevious=xii;
			y=y+weight*yi;
		}
		xVecResampled.push_back(x);
		yVecResampled.push_back(y);
	}
	xVec=xVecResampled;
	yVec=yVecResampled;
}

/**
 * @brief xVecGoal with a regular sampling (dx) is created
 * Example: for dx=50, xMax=184: xVecGoal={0,50,100,150}
 */
void DataResampler::createXVecGoalFromDx(const double dx)
{
	xVecGoal.clear();
	double x=0;
	while(x < xMax){
		xVecGoal.push_back(x);
		x=x+dx;
	}
}

/**
 * @brief The direction of the x-axis is changed and the origin is shifted to xMax.
 * Example: xVec={0,40,90} xMax=100 is changed to:
 *			     {10,60,100} for point-data
 *				 {0,10,60} for range-data
 */
void DataResampler::changeDirectionOfXAxis()
{
	std::vector<double> dxVec;
	std::vector<double> xVecTemp;
	std::vector<double> yVecTemp;

	//calculate thicknesses of the layers:
	for (size_t ii=0; ii<xVec.size()-1; ii++){
		double dz = xVec[ii+1]-xVec[ii];
		dxVec.push_back(dz);
	}
	dxVec.push_back(xMax-xVec[xVec.size()-1]);
	//change x-direction:
	for (size_t ii=0; ii<xVec.size(); ii++){
		double thickness = 0;
		if (useMethodForRangeValues){
			thickness = dxVec[ii];
		}
		double zTemp = xMax - xVec[ii] - thickness;
		xVecTemp.push_back(zTemp);
		yVecTemp.push_back(yVec[ii]);
	}
	//reverse the order:
	for (size_t ii=0; ii<xVec.size(); ii++){
		xVec[xVec.size()-ii-1] = xVecTemp[ii];
		yVec[xVec.size()-ii-1] = yVecTemp[ii];
	}
}


//Define namespaces and abbreviations
const char* CaaMLIO::xml_ns_caaml = (const char*) "http://caaml.org/Schemas/SnowProfileIACS/v6.0.3";
const char* CaaMLIO::xml_ns_abrev_caaml = (const char*) "caaml";
const char* CaaMLIO::xml_ns_gml = (const char*) "http://www.opengis.net/gml";
const char* CaaMLIO::xml_ns_abrev_gml = (const char*) "gml";
const char* CaaMLIO::xml_ns_xsi = (const char*) "http://www.w3.org/2001/XMLSchema-instance";
const char* CaaMLIO::xml_ns_abrev_xsi = (const char*) "xsi";
const char* CaaMLIO::xml_ns_slf = (const char*) "http://www.slf.ch/snowprofile/1.0";
const char* CaaMLIO::xml_ns_abrev_slf = (const char*) "slf";
const char* CaaMLIO::xml_ns_snp = (const char*) "http://www.slf.ch/snowpack/1.0";
const char* CaaMLIO::xml_ns_abrev_snp = (const char*) "snp";
// const std::string xml_schemaLocation_snp = "http://www.slf.ch/snowpack/snowpack.xsd";
const std::string namespaceCAAML = "caaml";
const std::string namespaceSNP = "snp";
const std::string namespaceGML = "gml";
//Define paths in xml-file
//const std::string CaaMLIO::TimeData_xpath = "/caaml:SnowProfile/caaml:validTime";
const std::string CaaMLIO::TimeData_xpath = "/caaml:SnowProfile/caaml:timeRef/caaml:recordTime/caaml:TimeInstant/caaml:timePosition";

const std::string CaaMLIO::StationMetaData_xpath = "/caaml:SnowProfile/caaml:locRef";

const std::string CaaMLIO::SnowData_xpath = "/caaml:SnowProfile/caaml:snowProfileResultsOf/caaml:SnowProfileMeasurements";

CaaMLIO::CaaMLIO(const SnowpackConfig& cfg, const RunInfo& run_info)
           : info(run_info), i_snowpath(), o_snowpath(), experiment(), i_max_element_thickness(IOUtils::nodata),
             caaml_writeout_as_readin(false), haz_write(true), in_tz(), inDoc(),inEncoding(),hoarDensitySurf(0),grainForms()
{
	init(cfg);
	grainForms.clear();
}

/**
 * @brief This routine initializes the caamlIO (filenames, XML-encoding...)
 * @param cfg
 */
void CaaMLIO::init(const SnowpackConfig& cfg)
{
	std::string tmpstr;

	cfg.getValue("TIME_ZONE", "Input", in_tz);

	cfg.getValue("METEOPATH", "Input", tmpstr, IOUtils::nothrow);
	cfg.getValue("SNOWPATH", "Input", i_snowpath, IOUtils::nothrow);
	if (i_snowpath.empty())
		i_snowpath = tmpstr;
	cfg.getValue("CAAML_MAX_ELEMENT_THICKNESS", "Input", i_max_element_thickness, IOUtils::nothrow);
	if (i_max_element_thickness==0.) i_max_element_thickness=IOUtils::nodata; //so inishell can expose it
	cfg.getValue("CAAML_WRITEOUT_AS_READIN", "Input", caaml_writeout_as_readin, IOUtils::nothrow);

	cfg.getValue("EXPERIMENT", "Output", experiment);
	cfg.getValue("METEOPATH", "Output", tmpstr, IOUtils::nothrow);
	cfg.getValue("SNOWPATH", "Output", o_snowpath, IOUtils::nothrow);
	if (o_snowpath.empty())
		o_snowpath = tmpstr;
	cfg.getValue("HAZ_WRITE", "Output", haz_write, IOUtils::nothrow);

	//input encoding forcing, inherited from CosmoXMLIO
	inEncoding=pugi::encoding_auto;
	tmpstr.clear();
	cfg.getValue("XML_ENCODING", "INPUT", tmpstr, IOUtils::nothrow);
	if (!tmpstr.empty()) {
		if (tmpstr=="UTF-8") inEncoding=pugi::encoding_utf8;
		else if (tmpstr=="UTF-16-LE") inEncoding=pugi::encoding_utf16_le;
		else if (tmpstr=="UTF-16-BE") inEncoding=pugi::encoding_utf16_be;
		else if (tmpstr=="UTF-16") inEncoding=pugi::encoding_utf16;
		else if (tmpstr=="LATIN1") inEncoding=pugi::encoding_latin1;
		else if (tmpstr=="ISO-8859-1") inEncoding=pugi::encoding_latin1;
		else if (tmpstr=="UTF-32") inEncoding=pugi::encoding_utf32;
		else if (tmpstr=="UTF-32-LE") inEncoding=pugi::encoding_utf32_le;
		else if (tmpstr=="UTF-32-BE") inEncoding=pugi::encoding_utf32_be;
		else if (tmpstr=="WCHAR") inEncoding=pugi::encoding_wchar;
		else
			throw InvalidArgumentException("Encoding \""+tmpstr+"\" is not supported!", AT);
	}
	//get density for surface hoar
	cfg.getValue("HOAR_DENSITY_SURF", "SnowpackAdvanced", hoarDensitySurf); // Density of SH at surface node (kg m-3)
}

/**
 * @brief This routine opens the caaml-file and checks the syntax of the caaml-file.
 * @param in_snowfile filename (with path) of the caaml-file
 */
void CaaMLIO::openIn_CAAML(const std::string& in_snowfile)
{
	const pugi::xml_parse_result result = inDoc.load_file(in_snowfile.c_str(),pugi::parse_default,inEncoding);
	if (!result){
		throw IOException("CAAML [" + in_snowfile + "] parsed with errors. Error description: " + result.description(), AT);
	}
}


/**
 * @brief This routine checks if the specified snow cover data exists
 * @param i_snowfile file containing the initial state of the snowpack
 * @param stationID
 * @return true if the file exists
 */
bool CaaMLIO::snowCoverExists(const std::string& i_snowfile, const std::string& /*stationID*/) const
{
	std::string snofilename( getFilenamePrefix(i_snowfile, i_snowpath, false) );

	if (snofilename.rfind(".caaml") == string::npos) {
		snofilename += ".caaml";
	}

	return FileUtils::fileExists(snofilename);
}

/**
 * @brief This routine reads the status of the snow cover at program start
 * @param i_snowfile file containing the initial state of the snowpack
 * @param stationID
 * @param SSdata
 * @param Zdata
 * @param read_salinity
 */
void CaaMLIO::readSnowCover(const std::string& i_snowfile, const std::string& stationID,
                            SN_SNOWSOIL_DATA& SSdata, ZwischenData& Zdata, const bool& /*read_salinity*/)
{
	std::string snofilename( getFilenamePrefix(i_snowfile, i_snowpath, false) );
	std::string hazfilename(snofilename);

	if (snofilename.rfind(".caaml") == string::npos) {
		snofilename += ".caaml";
		hazfilename += ".haz";
	} else {
		hazfilename.replace(hazfilename.rfind(".caaml"), 6, ".haz");
	}

	read_snocaaml(snofilename, stationID, SSdata);
	Zdata.reset();

	if (caaml_writeout_as_readin){
		SnowStation Xdata(true, false);
		Xdata.initialize(SSdata,0);
		mio::Date now;
		now.setFromSys();
		writeSnowCover(now,Xdata,Zdata,true);
	}
}

/**
 * @brief This routine returns the filename_prefix
 * @param fnam e.g. "20180516-5wj_84934_en-v6"
 * @param path e.g. "input"
 * @param addexp if true the experiment-string will be added to the filename
 * @return the filename-prefix e.g. "input/20180516-5wj_84934_en-v6"
 */
std::string CaaMLIO::getFilenamePrefix(const std::string& fnam, const std::string& path, const bool addexp) const
{
	//TODO: read only once (in constructor)
	std::string filename_prefix( path + "/" + fnam );

	if (addexp && (experiment != "NO_EXP"))
		filename_prefix += "_" + experiment;

	return filename_prefix;
}

/**
 * @brief This routine reads the caaml file into SSdata
 * @param in_snowFilename filename of the caaml file
 * @param stationID
 * @param SSdata data structure including all important station parameters as well as LayerData
 * @return true if file was read successfully
 */
bool CaaMLIO::read_snocaaml(const std::string& in_snowFilename, const std::string& stationID, SN_SNOWSOIL_DATA& SSdata)
{
	// Read CAAML snow profile file
	openIn_CAAML(in_snowFilename);

	//Read profile date
	SSdata.profileDate = xmlGetDate();

	//Read station metadata
	SSdata.meta = xmlGetStationData(stationID);

	//Snow-Soil properties: set to default if not available in file
	setCustomSnowSoil(SSdata);

	//Read layers
	xmlReadLayerData(SSdata);

	if (SSdata.nLayers>0) {
		const bool directionTopDown = getLayersDir(); //Read profile direction
		//Read quantity profiles (density, temperature and lwc)
		//temperature data is needed. if not available an exception is thrown:
		getAndSetProfile("/caaml:tempProfile/caaml:Obs","caaml:snowTemp",directionTopDown,false,SSdata.Ldata);
		//density data is needed. if not available an exception is thrown:
		getAndSetProfile("/caaml:densityProfile/caaml:Layer","caaml:density",directionTopDown,true,SSdata.Ldata);
		//lwc data is optional. if not available, no problem:
		getAndSetProfile("/caaml:lwcProfile/caaml:Layer","caaml:lwc",directionTopDown,true,SSdata.Ldata);

		adjustToSlopeAngle(SSdata); // #cmb
		
		checkAllDataForConsistencyAndSetMissingValues(SSdata);
	}
	//checkWhatWasReadIn(SSdata);
	return true;
}


/**
 * @brief This routine adjust the snow pack thickness/height to the slope angle (internal snow profile is always vertical to ground)
 * @param SSdata data structure including all important station parameters as well as LayerData
 * #cmb (needs to be called before checkAllDataForConsistency...)
 */
void CaaMLIO::adjustToSlopeAngle(SN_SNOWSOIL_DATA& SSdata)
{
	const double cos_sl = cos(SSdata.meta.getSlopeAngle()*mio::Cst::to_rad);
	SSdata.Height = SSdata.Height * cos_sl;
	SSdata.HS_last = SSdata.HS_last * cos_sl;
	for (size_t ii=0; ii<SSdata.nLayers; ii++) {
		SSdata.Ldata[ii].hl = SSdata.Ldata[ii].hl * cos_sl;
		if(i_max_element_thickness != IOUtils::nodata){
			SSdata.Ldata[ii].ne = (size_t) ceil(SSdata.Ldata[ii].hl/i_max_element_thickness);
		}else{
			SSdata.Ldata[ii].ne=1;
		}
	}
}


/**
 * @brief This routine reads the date of the snowpack from the caaml-file
 */
Date CaaMLIO::xmlGetDate()
{
	Date date;
	pugi::xml_node dateNode =  inDoc.first_element_by_path(TimeData_xpath.c_str());
	const std::string dateStr( (char*)dateNode.child_value() );
	if (dateStr==""){
		throw NoDataException("No data found for '"+TimeData_xpath+"'. Please check the version of the caaml-file. Only caaml-version 6 is supported. ", AT);
	}
	IOUtils::convertString(date, dateStr, in_tz);
	return date;
}

/**
 * @brief This routine reads the station data from the caaml-file
 * @param stationID
 */
StationData CaaMLIO::xmlGetStationData(const std::string& stationID)
{
	double x=IOUtils::nodata, y=IOUtils::nodata, z=IOUtils::nodata;
	double slopeAngle=IOUtils::nodata, azimuth=IOUtils::nodata;

	pugi::xml_node stationNode =  inDoc.first_element_by_path(StationMetaData_xpath.c_str());
	if( stationNode.empty() ){
		throw NoDataException("No data found for '"+StationMetaData_xpath+"'", AT);
	}

	const std::string stationName = std::string( (const char*)stationNode.child("caaml:name").child_value() );
	sscanf(stationNode.child("caaml:validElevation").child("caaml:ElevationPosition").child("caaml:position").child_value(),"%lf",&z);
	const char* azimuthStr = stationNode.child("caaml:validAspect").child("caaml:AspectPosition").child("caaml:position").child_value();
	azimuth = IOUtils::bearing(azimuthStr);  //convert aspect-string to a bearing (0-360) if given as a char (N, NW, S,...)
	if(azimuth == IOUtils::nodata){
		sscanf(azimuthStr,"%lf",&azimuth);   //if aspect was given as a number, just convert aspect-string to double
	}
	sscanf(stationNode.child("caaml:validAspect").child("caaml:AspectPosition").child("caaml:position").child_value(),"%lf",&azimuth);
	sscanf(stationNode.child("caaml:validSlopeAngle").child("caaml:SlopeAnglePosition").child("caaml:position").child_value(),"%lf",&slopeAngle);
	sscanf(stationNode.child("caaml:pointLocation").child("gml:Point").child("gml:pos").child_value(),"%lf %lf",&x,&y);

	Coords tmppos;
	tmppos.setLatLon(x, y, z);
	StationData metatmp(tmppos, stationID, stationName);
	metatmp.setSlope(slopeAngle, azimuth);
	return metatmp;
}

/**
 * @brief This routine reads custom-snow-soil-data from the caaml-file into Xdata.
 * @param Xdata data structure including all important station parameters as well as LayerData
 */
void CaaMLIO::setCustomSnowSoil(SN_SNOWSOIL_DATA& Xdata)
{
	const std::string xpath( "/caaml:metaData/caaml:customData/snp" );
	//if (xmlDoesPathExist(SnowData_xpath+xpath+":Albedo") == false)
		//std::cout << "There is no snowpack-custom-data in the caaml-file. Setting all values (albedo, erosion level,...) to default..." << std::endl;
	Xdata.Albedo = xmlReadValueFromPath(xpath,"Albedo",0.6);
	Xdata.SoilAlb = xmlReadValueFromPath(xpath,"SoilAlb",0.2);
	Xdata.BareSoil_z0 = xmlReadValueFromPath(xpath,"BareSoil_z0",0.02);
	Xdata.Canopy_Height = xmlReadValueFromPath(xpath,"CanopyHeight",0.);
	Xdata.Canopy_LAI = xmlReadValueFromPath(xpath,"CanopyLAI",0.);
	Xdata.Canopy_BasalArea = xmlReadValueFromPath(xpath,"CanopyBasalArea",0.);
	Xdata.Canopy_Direct_Throughfall = xmlReadValueFromPath(xpath,"CanopyDirectThroughfall",1.);
	Xdata.WindScalingFactor = xmlReadValueFromPath(xpath,"WindScalingFactor",1.);
	Xdata.ErosionLevel = xmlReadValueFromPath(xpath,"ErosionLevel",0);
	Xdata.TimeCountDeltaHS = xmlReadValueFromPath(xpath,"TimeCountDeltaHS",0.);
}

/**
 * @brief Direction in which the layers are stored in the caaml-file. "top down" or "bottom up".
 *        Snowpack needs the "bottom up" direction.
 * @return true if the direction is "bottom up"
 */
bool CaaMLIO::getLayersDir()
{
	const std::string direction( xmlReadAttributeFromPath(SnowData_xpath,"dir") );
	return (direction!="bottom up");
}

/**
 * @brief Read one stratigraphic layer from the caaml-file
 * @param nodeLayer xml-node pointing to <caaml:Layer> in the caaml-file
 * @return data of one layer which was read in from the caaml-file
 */
LayerData CaaMLIO::xmlGetLayer(pugi::xml_node nodeLayer, std::string& grainFormCode)
{
	grainFormCode="";
	LayerData Layer;
	for (pugi::xml_node node = nodeLayer.first_child(); node; node = node.next_sibling()){
		const std::string fieldName( node.name() );
		if(fieldName == "caaml:grainSize"){
			const std::string unitMeasured( node.attribute("uom").as_string() );
			if (xmlReadValueFromNode(node.child("caaml:Components").child("caaml:avg"),"caaml:avg",Layer.rg,"mm",unitMeasured,0.5)){
				Layer.rb = Layer.rg/4.; //this value will be replaced if there is this value in the customData...
			}else if(xmlReadValueFromNode(node.child("caaml:Components").child("caaml:avgMax"),"caaml:avgMax",Layer.rg,"mm",unitMeasured,0.5)){
				Layer.rb = Layer.rg/4.; //this value will be replaced if there is this value in the customData...
			}
		}
		if (xmlReadValueFromNode(node,"caaml:thickness",Layer.hl,"m")){
			if(i_max_element_thickness != IOUtils::nodata){
				Layer.ne = (size_t) ceil(Layer.hl/i_max_element_thickness);
			}else{
				Layer.ne=1;
			}
		}
		if (fieldName == "caaml:wetness"){ //this value will be replaced if a lwc-profile is in the caaml-file!!!
			Layer.phiWater = lwc_codeToVal((char*) node.child_value());
		}
		if (fieldName == "caaml:hardness"){
			//const double hardness = hardness_codeToVal((char*) node.child_value());
			//std::cout << "hardness (we have no member (of Layer) we could set this value to!!!): " << hardness << " " << std::endl;
		}
		if (fieldName == "caaml:grainFormPrimary"){
			grainFormCode = std::string( (char*)node.child_value() );
			grainShape_codeToVal(grainFormCode, Layer.sp, Layer.dd, Layer.mk); // these values will be replaced if there are values in the custom-snowpack-data
		}
		if (fieldName == "caaml:validFormationTime"){
			const std::string date_str( (char*) node.child("caaml:TimeInstant").child("caaml:timePosition").child_value() );
			Date date;
			IOUtils::convertString(date, date_str, in_tz);
			Layer.depositionDate = date;
		}
	}

	//read custom-snowpack-data. do this at the end to be sure to use these values (if available)
	for (pugi::xml_node node = nodeLayer.child("caaml:metaData").child("caaml:customData").first_child(); node; node = node.next_sibling())
	{
		xmlReadValueFromNode(node,"snp:dendricity",Layer.dd);
		xmlReadValueFromNode(node,"snp:sphericity",Layer.sp);
		xmlReadValueFromNode(node,"snp:bondSize",Layer.rb);
		xmlReadValueFromNode(node,"snp:phiSoil",Layer.phiSoil);
		xmlReadValueFromNode(node,"snp:SurfaceHoarMass",Layer.hr);
		xmlReadValueFromNode(node,"snp:StressRate",Layer.CDot);
		xmlReadValueFromNode(node,"snp:Metamorphism",Layer.metamo);
		const std::string fieldName( node.name() );
		if (fieldName == "snp:marker") {
			sscanf((const char*) node.child_value(),"%hu",&Layer.mk);
		}
	}
	return Layer;
}

/**
 * @brief Read the stratigraphic layers from the caaml-file to SSdata
 * @param SSdata
 */
void CaaMLIO::xmlReadLayerData(SN_SNOWSOIL_DATA& SSdata)
{
	const bool directionTopDown = getLayersDir(); //Read profile direction

	std::string path = SnowData_xpath+"/caaml:snowPackCond/caaml:hS/caaml:Components/caaml:height";
	const pugi::xml_node nodeHS =  inDoc.first_element_by_path( (char*) path.c_str() );
	SSdata.HS_last=IOUtils::nodata;
	xmlReadValueFromNode(nodeHS,"caaml:height",SSdata.HS_last,"m");

	size_t nLayers=0;
	path = SnowData_xpath+"/caaml:stratProfile/caaml:Layer";
	const pugi::xml_node nodeFirstLayer =  inDoc.first_element_by_path( (char*) path.c_str() );
	if(nodeFirstLayer.empty()){
		throw NoDataException("Layer data not found in caaml-file. Expected path: '"+path+"'", AT);
	}
	pugi::xml_node nodeLastLayer;
	for (pugi::xml_node node = nodeFirstLayer; node; node = node.next_sibling("caaml:Layer")){
		nLayers++;
		nodeLastLayer=node;
	}
	SSdata.nLayers=nLayers;
	SSdata.Ldata.resize(SSdata.nLayers, LayerData());
	grainForms.clear();
	std::string grainForm="";
	if (SSdata.nLayers>0) {
		if (!directionTopDown) {
			size_t ii=0;
			for (pugi::xml_node node = nodeFirstLayer; node; node = node.next_sibling("caaml:Layer")){
				SSdata.Ldata[ii] = xmlGetLayer(node,grainForm);
				grainForms.push_back(grainForm);
				ii++;
			}
		}else{
			size_t ii=0;
			for (pugi::xml_node node = nodeLastLayer; node; node = node.previous_sibling("caaml:Layer")){
				SSdata.Ldata[ii] = xmlGetLayer(node,grainForm);
				grainForms.push_back(grainForm);
				ii++;
			}
		}
		//Estimate depostion dates (in case it was not in the caaml-file):
		estimateValidFormationTimesIfNotSetYet(SSdata.Ldata,SSdata.profileDate);
	}
}

/**
 * @brief Read profile data (density, temperature or lwc) from the caaml-file.
 *        The density- and temperature-profile have to be in the caaml-file. Otherwise an
 *        exception is thrown. The LWC-profile is optional and ignored if not available.
 * @param path path to the profile data. e.g. "/caaml:densityProfile/caaml:Layer"
 * @param name name of the xml-element to read. e.g. "caaml:snowTemp"
 * @param zVec output: depth-values in meter
 * @param valVec output: profile-values (valVec has the same length as zVec)
 * @return true if the profile was read in successfully
 */
bool CaaMLIO::xmlGetProfile(const std::string path, const std::string name, std::vector<double>& zVec, std::vector<double>& valVec)
{
	//check if data exists:
	if( !xmlDoesPathExist(SnowData_xpath+path)){
		if(name == "caaml:lwc"){
			//std::cout << "No lwc-profile in the caaml-file. lwc will be estimated from the wetness in the stratigraphic profile. " << std::endl;
			return false;
		}else{
			throw NoDataException("Invalid path for: '"+name+"'-profile. path: '"+path+"'", AT);
		}
	}

	std::string profilePath = SnowData_xpath+path;
	const pugi::xml_node nodeFirstProfile =  inDoc.first_element_by_path( (char*) profilePath.c_str() );
	zVec.resize(0);
	valVec.resize(0);
	for (pugi::xml_node nodeProfile = nodeFirstProfile; nodeProfile; nodeProfile = nodeProfile.next_sibling()){
		for (pugi::xml_node node = nodeProfile.first_child(); node; node = node.next_sibling()){
			const std::string fieldName( node.name() );
			if (fieldName == name){
					double val=0;
				sscanf((const char*) node.child_value(), "%lf", &val);
				if (name=="caaml:snowTemp"){
					val = unitConversion(val,(char*)node.attribute("uom").as_string(),(char*)"K");
				}
				valVec.push_back(val);
			}
			if (fieldName == "caaml:depth" || fieldName == "caaml:depthTop") {
				double z=0;
				sscanf((const char*) node.child_value(), "%lf", &z);
				z = unitConversion(z,(char*)node.attribute("uom").as_string(),(char*)"m");
				zVec.push_back(z);
			}
		}
	}
	return true;
}

/**
 * @brief Read profile data (density, temperature or lwc) from the caaml-file, rotate the direction
 *        of the data to "top up" and resample the data to the sampling of the stratigraphic data.
 * @param path path to the profile data. e.g. "/caaml:densityProfile/caaml:Layer"
 * @param name name of the xml-element to read. e.g. "caaml:snowTemp"
 * @param directionTopDown true if the direction is top down
 * @param isRangeMeasurement true if the data is a range measurement (like density or lwc),
 *                           false if the data is a point measurement (like temperature)
 * @param Layers output
 */
void CaaMLIO::getAndSetProfile(const std::string path, const std::string name,
							   const bool directionTopDown, const bool isRangeMeasurement,std::vector<LayerData>& Layers)
{
	std::vector<double> zVec;
	std::vector<double> valVec;
	//read the data from the Caaml-file:
	bool profileExists = xmlGetProfile(path, name, zVec, valVec);
	if (profileExists){
		//direction has to be bottomUp:
		bool changeDirection=false;
		if (directionTopDown) {
			changeDirection=true;
		}
		//create a vector of the z-positions from the layerData. the profile will be resampled to these layers.
		std::vector<double> zVecLayers;
		double z=0.;
		for (size_t ii=0; ii<Layers.size(); ii++) {
			if(isRangeMeasurement){
				zVecLayers.push_back(z);
				z += Layers[ii].hl;
			}else{
				z += Layers[ii].hl;
				zVecLayers.push_back(z);
			}
		}
		//resample the profile data
		DataResampler resampler(zVec,valVec,zVecLayers,z,isRangeMeasurement,changeDirection);
		if (resampler.resample()==false){
			throw IOException("Something wrong with "+name+"-profile data in CaaML-file", AT);
		}
		valVec = resampler.getResampledYVec();

		//Set temperature, density and lwc from the profiles
		for (size_t ii=0; ii<Layers.size(); ii++) { //loop over the number of layers
			if (name=="caaml:snowTemp"){
				Layers[ii].tl = valVec[ii];
			}
			if (name=="caaml:density"){
				Layers[ii].phiIce = valVec[ii]/Constants::density_ice;
			}
			if (name=="caaml:lwc"){
				Layers[ii].phiWater = valVec[ii]/100.0;  //change from % to relative value
			}
		}
	}
}

/**
 * @brief If the deposition-date-value for a layer is not defined in the caaml-file, the deposition-date
 *        will be estimated from profile-date and the snow-type.
 * @param Layers output
 * @param profileDate date of the snow-profile
 */
void CaaMLIO::estimateValidFormationTimesIfNotSetYet(std::vector<LayerData> &Layers, const Date profileDate)
{
	for (size_t ii=0; ii<Layers.size(); ii++) {
		if(Layers[ii].depositionDate.isUndef()){
			const unsigned int snowType = ElementData::snowType(Layers[ii].dd,Layers[ii].sp,Layers[ii].rg,Layers[ii].mk,Layers[ii].phiWater,
																ElementData::snowResidualWaterContent(Layers[ii].phiIce));
			const unsigned int a = (unsigned int) (snowType/100.);
			if (ii==0) {
				if (a==6) {
					Layers[ii].depositionDate = profileDate;
				} else if (a==0 || a==1) {
					Layers[ii].depositionDate = profileDate-(Date)1./2440638.;
				} else {
					Layers[ii].depositionDate = profileDate-(Date)2./2440638.;
				}
			} else {
				if ((a==0 || a==1) && (Layers[ii-1].depositionDate > profileDate-(Date)2./2440638.)) {
					Layers[ii].depositionDate = profileDate-(Date)1./2440638.;
				} else {
					Layers[ii].depositionDate = profileDate-(Date)2./2440638.;
				}
			}
		}
	}
}

/**
 * @brief Check station- and layer-data for consistency. If possible set reasonable values, otherwise throw exceptions or write warnings.
 *        Furthermore determine some missing values (total number of elements (nodes), height and phiVoids).
 * @param SSdata
 */
void CaaMLIO::checkAllDataForConsistencyAndSetMissingValues( SN_SNOWSOIL_DATA& SSdata )
{
	//check station data:
	double azimuth = SSdata.meta.getAzimuth();
	double slopeAngle = SSdata.meta.getSlopeAngle();
	if(slopeAngle ==IOUtils::nodata){
		slopeAngle=0;
	}
	if(azimuth==IOUtils::nodata){
		if(slopeAngle==0){
			azimuth=0;
		}
		if(slopeAngle > 0){
			throw NoDataException("No data found for 'caaml:validAspect'. If the slope-angle is >0 degree, we also need the azimuth. ", AT);
		}
	}
	SSdata.meta.setSlope(slopeAngle,azimuth);

	//check layer data:
	for (size_t ii = 0; ii < SSdata.nLayers; ii++) {
		std::string grainFormCode = grainForms[ii];
		if (grainFormCode=="MF"){
			if (SSdata.Ldata[ii].tl < Constants::meltfreeze_tk-0.05){
				SSdata.Ldata[ii].mk = (unsigned short int) (SSdata.Ldata[ii].mk + 10);
			}
		}
		if (grainFormCode=="SH" && ii==SSdata.nLayers-1){ //set parameters for surface hoar (only at the surface, not buried)
			SSdata.Ldata[ii].phiWater=0;
			SSdata.Ldata[ii].phiIce = hoarDensitySurf/Constants::density_ice;
			double grainRadius = M_TO_MM(SSdata.Ldata[ii].hl/2.0);
			if(grainRadius != SSdata.Ldata[ii].rg){
				std::cout << "WARNING! Inconsistent input data in caaml-file: Grain size for surface-hoar-layer should be about the same value as the "
				          << "surface-hoar-layer-thickness (" << grainRadius*2 << " mm). Adjusting grain size from: "
				          << SSdata.Ldata[ii].rg*2 << " mm to " << grainRadius*2 << " mm." << std::endl;
				SSdata.Ldata[ii].rg = grainRadius;
			}
			SSdata.Ldata[ii].rb = SSdata.Ldata[ii].rg/3.;
			SSdata.Ldata[ii].dd = 0.;
			SSdata.Ldata[ii].sp = 0.;
			SSdata.Ldata[ii].mk = 3;
		}
		if (SSdata.Ldata[ii].rg == 0.) {
			if (grainFormCode=="IF") {
				SSdata.Ldata[ii].rg = 3./2.;
				SSdata.Ldata[ii].rb = 3./8.;
			} else throw IOException("Grain size missing for a non-ice layer!", AT);
		}
		// set temperature to 0°C if warmer than 0°C:
		if (SSdata.Ldata[ii].tl > Constants::meltfreeze_tk){
			std::cout << "WARNING! Inconsistent input data in caaml-file: Temperature in layer " << ii << ": " << SSdata.Ldata[ii].tl-Constants::meltfreeze_tk
			          << " degree Celsius. Temperature above 0 degree Celsius not possible! Setting temperature to 0 degree Celsius." << std::endl;
			SSdata.Ldata[ii].tl = Constants::meltfreeze_tk;
		}
		// set lwc to 0 if colder than 0°C:
		if (SSdata.Ldata[ii].tl < Constants::meltfreeze_tk && SSdata.Ldata[ii].phiWater > 0 ){
			std::cout << "WARNING! Inconsistent input data in caaml-file: LWC: " << SSdata.Ldata[ii].phiWater
			          << " temperature: " << SSdata.Ldata[ii].tl << " in layer: " << ii << std::endl;
			std::cout << "Setting lwc to 0! Since liquid water is not possible in snow below 0 degree Celsius." << std::endl;
			//throw IOException("LWC > 0 but temperature below 0°C!", AT);
			SSdata.Ldata[ii].phiWater = 0;
		}
		if (grainFormCode=="FC" && SSdata.Ldata[ii].rg>0.8){
			std::cout << "WARNING! Inconsistent input data in caaml-file: Grain shape 'FC' and grain size > 1.5mm! "
			          << "Faceted crystals should be smaller than 1.5mm. " << std::endl;
		}
		if (grainFormCode=="DH" && SSdata.Ldata[ii].rg<0.7){
			std::cout << "WARNING! Inconsistent input data in caaml-file: Grain shape 'DH' and grain size < 1.5mm! "
			          << "Depth hoar crystals should be larger than 1.5mm. " << std::endl;
		}
	}
	//Compute total number of elements (nodes), height and phiVoids
	SSdata.nN = 1;
	SSdata.Height = 0.;
	for (size_t ii = 0; ii < SSdata.nLayers; ii++) {
		SSdata.nN += SSdata.Ldata[ii].ne;
		SSdata.Height += SSdata.Ldata[ii].hl;
		SSdata.Ldata[ii].phiVoids = 1. - SSdata.Ldata[ii].phiSoil - SSdata.Ldata[ii].phiWater - SSdata.Ldata[ii].phiIce;
	}

	if (SSdata.HS_last != IOUtils::nodata){
		if(SSdata.HS_last > SSdata.Height+0.0005 || SSdata.HS_last < SSdata.Height-0.0005){
			std::cout << "WARNING! Inconsistent input data in caaml-file: Snow-height given in caaml:snowPackCond (" << SSdata.HS_last
			          << " m) is different from the sum of the layer-thicknesses (" << SSdata.Height
			          << " m)! For the simulation the sum of the layer-thicknesses will be used." << std::endl;
		}
	}
	SSdata.HS_last = SSdata.Height;
}

/**
 * @brief This routine writes the status of the snow cover at program termination and at specified backup times
 * @param date current
 * @param Xdata
 * @param Zdata
 * @param forbackup dump Xdata on the go
 */
void CaaMLIO::writeSnowCover(const Date& date, const SnowStation& Xdata,
                             const ZwischenData& Zdata, const bool& forbackup)
{
	const std::string bak = (forbackup)? "_" + date.toString(mio::Date::NUM) : "";
	std::string snofilename( getFilenamePrefix(Xdata.meta.getStationID().c_str(), o_snowpath) + bak + ".caaml" );
	std::string hazfilename( getFilenamePrefix(Xdata.meta.getStationID().c_str(), o_snowpath) + bak + ".haz" );

	writeSnowFile(snofilename, date, Xdata);
	if (haz_write) SmetIO::writeHazFile(hazfilename, date, Xdata, Zdata);
}

/**
 * @brief This routine writes the status of the snow cover to a caaml-file.
 * @param snofilename
 * @param date date of the snow-profile
 * @param Xdata data structure containing all the data describing a snow-profile (which will be written to
 *              the caaml-file)
 */
void CaaMLIO::writeSnowFile(const std::string& snofilename, const Date& date, const SnowStation& Xdata)
{
	pugi::xml_document doc;
	// Generate XML declaration
	pugi::xml_node declarationNode = doc.append_child(pugi::node_declaration);
	declarationNode.append_attribute("version")    = "1.0";
	declarationNode.append_attribute("encoding")   = "UTF-8";
	// A valid XML doc must contain a single root node of any name
	pugi::xml_node root = doc.append_child( (namespaceCAAML+":SnowProfile").c_str() );
	root.append_attribute("gml:id") = ("SLF_"+Xdata.meta.stationID).c_str();
	root.append_attribute(("xmlns:"+string((const char*)xml_ns_abrev_xsi)).c_str()) = string((const char*)xml_ns_xsi).c_str();
	root.append_attribute(("xmlns:"+string((const char*)xml_ns_abrev_gml)).c_str()) = string((const char*)xml_ns_gml).c_str();
	root.append_attribute(("xmlns:"+string((const char*)xml_ns_abrev_caaml)).c_str()) = string((const char*)xml_ns_caaml).c_str();
	root.append_attribute("xmlns:snp") = string((const char*)xml_ns_snp).c_str();
	root.append_attribute(("xmlns:"+string((const char*)xml_ns_abrev_slf)).c_str()) = string((const char*)xml_ns_slf).c_str();

	// Write profile date
	pugi::xml_node dateNode = root.append_child( (namespaceCAAML+":timeRef").c_str() );
	dateNode.append_child( (namespaceCAAML+":recordTime").c_str() )
	        .append_child( (namespaceCAAML+":TimeInstant").c_str() )
	        .append_child( (namespaceCAAML+":timePosition").c_str() )
	        .append_child(pugi::node_pcdata).set_value( date.toString(mio::Date::ISO_TZ).c_str() );
	mio::Date now;
	now.setFromSys();
	dateNode.append_child( (namespaceCAAML+":dateTimeReport").c_str() )
	        .append_child(pugi::node_pcdata).set_value( now.toString(mio::Date::ISO_TZ).c_str() );

	//Write srcRef
	pugi::xml_node srcNode = root.append_child( (namespaceCAAML+":srcRef").c_str() )
	                             .append_child( (namespaceCAAML+":Operation").c_str() );
	srcNode.append_attribute("gml:id") = "OPERATION_ID";
	const std::string snowpackString = "SNOWPACK v"+info.version+" compiled on "+info.compilation_date+". user: "+info.user;
	xmlWriteElement(srcNode,(namespaceCAAML+":name").c_str(),snowpackString.c_str(),"","");

	// Write station data locRef
	writeStationData(root,Xdata);

	// Write stratigraphic profile
	pugi::xml_node stratNode = root.append_child( (namespaceCAAML+":snowProfileResultsOf").c_str() )
	                               .append_child( (namespaceCAAML+":SnowProfileMeasurements").c_str() );
	stratNode.append_attribute("dir") = "top down";

	//Write custom snow/soil data
	pugi::xml_node snowSoilNode = stratNode.append_child( (namespaceCAAML+":metaData").c_str() )
	                                       .append_child( (namespaceCAAML+":customData").c_str() );
	writeCustomSnowSoil(snowSoilNode,Xdata);

	// Write profile depth
	sprintf(valueStr,"%.4f",100.*Xdata.cH/Xdata.cos_sl); // cmb
	xmlWriteElement(stratNode,(namespaceCAAML+":profileDepth").c_str(),valueStr,"uom","cm");

	//Write height of snow and Snow Water Equivalent (SWE)
	pugi::xml_node tempNode = stratNode.append_child( (namespaceCAAML+":snowPackCond").c_str() )
	                                   .append_child( (namespaceCAAML+":hS").c_str() )
	                                   .append_child( (namespaceCAAML+":Components").c_str() );
	sprintf(valueStr,"%.4f",100.*(Xdata.cH - Xdata.Ground)/Xdata.cos_sl); // cmb
	xmlWriteElement(tempNode,(namespaceCAAML+":height").c_str(),valueStr,"uom","cm");
	sprintf(valueStr,"%.2f",Xdata.swe);
	xmlWriteElement(tempNode,(namespaceCAAML+":waterEquivalent").c_str(),valueStr,"uom","kgm-2");

	//Write layers and quantity profiles
	writeLayers(stratNode,Xdata);
	writeProfiles(stratNode,Xdata);

	// Save XML tree to file.
	// Remark: second optional param is indent string to be used;
	// default indentation is tab character.
	bool saveSucceeded = doc.save_file(snofilename.c_str(), PUGIXML_TEXT("  "), pugi::format_default, pugi::encoding_utf8);
	if(!saveSucceeded){
		std::cout << "Something went wrong with saving the caaml-file: " << snofilename << std::endl;
	}
}


/**
 * @brief Write the custom-snow-soil-data to the caaml-file.
 */
void CaaMLIO::writeCustomSnowSoil(pugi::xml_node& node, const SnowStation& Xdata)
{
	sprintf(valueStr,"%.4f",Xdata.Albedo);
	xmlWriteElement(node,(namespaceSNP+":Albedo").c_str(),valueStr,"","");
	sprintf(valueStr,"%.4f",Xdata.SoilAlb);
	xmlWriteElement(node,(namespaceSNP+":SoilAlb").c_str(),valueStr,"","");
	sprintf(valueStr,"%.4f",Xdata.BareSoil_z0);
	xmlWriteElement(node,(namespaceSNP+":BareSoil_z0").c_str(),valueStr,"uom","m");
	sprintf(valueStr,"%.4f",Xdata.Cdata.height);
	xmlWriteElement(node,(namespaceSNP+":CanopyHeight").c_str(),valueStr,"uom","m");
	sprintf(valueStr,"%.4f",Xdata.Cdata.lai);
	xmlWriteElement(node,(namespaceSNP+":CanopyLAI").c_str(),valueStr,"","");
	sprintf(valueStr,"%.4f",Xdata.Cdata.BasalArea);
	xmlWriteElement(node,(namespaceSNP+":CanopyBasalArea").c_str(),valueStr,"","");
	sprintf(valueStr,"%.4f",Xdata.Cdata.throughfall);
	xmlWriteElement(node,(namespaceSNP+":CanopyDirectThroughfall").c_str(),valueStr,"","");
	sprintf(valueStr,"%.4f",Xdata.WindScalingFactor);
	xmlWriteElement(node,(namespaceSNP+":WindScalingFactor").c_str(),valueStr,"","");
	sprintf(valueStr,"%d",static_cast<unsigned int>(Xdata.ErosionLevel));
	xmlWriteElement(node,(namespaceSNP+":ErosionLevel").c_str(),valueStr,"","");
	sprintf(valueStr,"%.4f",Xdata.TimeCountDeltaHS);
	xmlWriteElement(node,(namespaceSNP+":TimeCountDeltaHS").c_str(),valueStr,"","");
}

/**
 * @brief Write the stratigraphic-layer-data to the caaml-file.
 */
void CaaMLIO::writeLayers(pugi::xml_node& node, const SnowStation& Xdata)
{
	pugi::xml_node stratNode = node.append_child( (namespaceCAAML+":stratProfile").c_str() );
	stratNode.append_child( (namespaceCAAML+":stratMetaData").c_str() );

	if (!Xdata.Edata.empty()) {
		for (size_t ii = Xdata.Edata.size(); ii-->0;) {
			const bool snowLayer = (ii >= Xdata.SoilNode);
			pugi::xml_node layerNode = stratNode.append_child( (namespaceCAAML+":Layer").c_str() );
			// Write custom layer data
			pugi::xml_node customNode = layerNode.append_child( (namespaceCAAML+":metaData").c_str() )
			                                     .append_child( (namespaceCAAML+":customData").c_str() );
			writeCustomLayerData(customNode,Xdata.Edata[ii],Xdata.Ndata[ii]);

			// Write snow and soil layer data 
			sprintf(layerDepthTopStr,"%.4f",100.*(Xdata.cH - Xdata.Ndata[ii+1].z)/Xdata.cos_sl); // cmb
			xmlWriteElement(layerNode,(namespaceCAAML+":depthTop").c_str(),layerDepthTopStr,"uom","cm");
			sprintf(layerThicknessStr,"%.4f",100.*Xdata.Edata[ii].L/Xdata.cos_sl); // cmb
			xmlWriteElement(layerNode,(namespaceCAAML+":thickness").c_str(),layerThicknessStr,"uom","cm");

			if (snowLayer) {
				//const unsigned int snowType = ElementData::snowType(Xdata.Edata[ii].dd, Xdata.Edata[ii].sp, Xdata.Edata[ii].rg, Xdata.Edata[ii].mk, Xdata.Edata[ii].theta[WATER],  Xdata.Edata[ii].res_wat_cont);
				const unsigned int snowType = Xdata.Edata[ii].getSnowType();
				const unsigned int a = snowType/100;
				const unsigned int b = (snowType-100*a)/10;
				const unsigned int c = snowType-100*a-10*b;
				if (c != 2) {
					xmlWriteElement(layerNode,(namespaceCAAML+":grainFormPrimary").c_str(),grainShape_valToAbbrev(a).c_str(),"","");
				} else {
					xmlWriteElement(layerNode,(namespaceCAAML+":grainFormPrimary").c_str(),"MFcr","","");
				}
				xmlWriteElement(layerNode,(namespaceCAAML+":grainFormSecondary").c_str(),grainShape_valToAbbrev(b).c_str(),"","");

				pugi::xml_node grainSizeNode = layerNode.append_child( (namespaceCAAML+":grainSize").c_str() );
				grainSizeNode.append_attribute("uom") = "mm";
				pugi::xml_node compNode = grainSizeNode.append_child( (namespaceCAAML+":Components").c_str() );
				sprintf(layerValStr,"%.3f",2.*Xdata.Edata[ii].rg);
				xmlWriteElement(compNode,(namespaceCAAML+":avg").c_str(),layerValStr,"","");
			}
			pugi::xml_node timeNode = layerNode.append_child( (namespaceCAAML+":validFormationTime").c_str() )
			                                   .append_child( (namespaceCAAML+":TimeInstant").c_str() );
			xmlWriteElement(timeNode,(namespaceCAAML+":timePosition").c_str(),Xdata.Edata[ii].depositionDate.toString(mio::Date::ISO_TZ).c_str(),"","");
			if (snowLayer) {
				xmlWriteElement(layerNode,(namespaceCAAML+":hardness").c_str(),hardness_valToCode(Xdata.Edata[ii].hard).c_str(),"uom",""); //HACK: check values... seem always the same!
			}
			xmlWriteElement(layerNode,(namespaceCAAML+":wetness").c_str(),lwc_valToCode(Xdata.Edata[ii].theta[WATER]).c_str(),"uom","");
		}
	}
}

/**
 * @brief Write the custom-layer-data to the caaml-file.
 */
void CaaMLIO::writeCustomLayerData(pugi::xml_node& node, const ElementData& Edata, const NodeData& Ndata)
{
	sprintf(valueStr,"%.4f",Edata.theta[SOIL]);
	xmlWriteElement(node,(namespaceSNP+":phiSoil").c_str(),valueStr,"","");
	sprintf(valueStr,"%.4f",Edata.soil[2]);
	xmlWriteElement(node,(namespaceSNP+":SoilRho").c_str(),valueStr,"uom","kgm-3");
	sprintf(valueStr,"%.4f",Edata.soil[0]);
	xmlWriteElement(node,(namespaceSNP+":SoilK").c_str(),valueStr,"uom","Wm-1s-1");
	sprintf(valueStr,"%.4f",Edata.soil[1]);
	xmlWriteElement(node,(namespaceSNP+":SoilC").c_str(),valueStr,"uom","Jkg-1");
	sprintf(valueStr,"%.4f",Edata.rb);
	xmlWriteElement(node,(namespaceSNP+":bondSize").c_str(),valueStr,"","");
	sprintf(valueStr,"%.2f",Edata.dd);
	xmlWriteElement(node,(namespaceSNP+":dendricity").c_str(),valueStr,"","");
	sprintf(valueStr,"%.2f",Edata.sp);
	xmlWriteElement(node,(namespaceSNP+":sphericity").c_str(),valueStr,"","");
	sprintf(valueStr,"%4u",static_cast<int>(Edata.mk));
	xmlWriteElement(node,(namespaceSNP+":marker").c_str(),valueStr,"","");
	sprintf(valueStr,"%.4f",Ndata.hoar);
	xmlWriteElement(node,(namespaceSNP+":SurfaceHoarMass").c_str(),valueStr,"uom","kgm-2");
	xmlWriteElement(node,(namespaceSNP+":ne").c_str(),"1","","");
	sprintf(valueStr,"%.4f",Edata.CDot);
	xmlWriteElement(node,(namespaceSNP+":StressRate").c_str(),valueStr,"uom","Nm-2s-1");
	sprintf(valueStr,"%.4f",Edata.metamo);
	xmlWriteElement(node,(namespaceSNP+":Metamorphism").c_str(),valueStr,"","");
	sprintf(valueStr,"%5u",Edata.ID);
	xmlWriteElement(node,(namespaceSNP+":ID").c_str(),valueStr,"","");
}

/**
 * @brief Write the density-, temperature, lwc-, SSA- and strength-profile to the caaml-file.
 */
void CaaMLIO::writeProfiles(pugi::xml_node& node, const SnowStation& Xdata)
{
	// temperature profile
	pugi::xml_node tempNode = node.append_child( (namespaceCAAML+":tempProfile").c_str() );
	tempNode.append_child( (namespaceCAAML+":tempMetaData").c_str() );
	if (!Xdata.Ndata.empty()) {
		for (size_t ii = Xdata.Ndata.size(); ii-->0;) {
			pugi::xml_node obsNode = tempNode.append_child( (namespaceCAAML+":Obs").c_str() );
			sprintf(layerDepthTopStr,"%.4f",100*(Xdata.cH - Xdata.Ndata[ii].z)/Xdata.cos_sl); // cmb
			xmlWriteElement(obsNode,(namespaceCAAML+":depth").c_str(),layerDepthTopStr,"uom","cm");
			sprintf(valueStr,"%.3f",unitConversion(Xdata.Ndata[ii].T,(char*)"degK",(char*)"degC"));
			xmlWriteElement(obsNode,(namespaceCAAML+":snowTemp").c_str(),valueStr,"uom","degC");
		}
	}//end temperature profile

	// density profile;
	pugi::xml_node densityNode = node.append_child( (namespaceCAAML+":densityProfile").c_str() );
	pugi::xml_node metaNode1 = densityNode.append_child( (namespaceCAAML+":densityMetaData").c_str() );
	xmlWriteElement(metaNode1,(namespaceCAAML+":methodOfMeas").c_str(),"other","","");
	if (!Xdata.Edata.empty()) {
		for (size_t ii = Xdata.Edata.size(); ii-->0;) {
			pugi::xml_node layerNode = densityNode.append_child( (namespaceCAAML+":Layer").c_str() );
			sprintf(layerDepthTopStr,"%.4f",100*(Xdata.cH - Xdata.Ndata[ii+1].z)/Xdata.cos_sl); // cmb
			xmlWriteElement(layerNode,(namespaceCAAML+":depthTop").c_str(),layerDepthTopStr,"uom","cm");
			sprintf(layerThicknessStr,"%.4f",100*Xdata.Edata[ii].L/Xdata.cos_sl); // cmb
			xmlWriteElement(layerNode,(namespaceCAAML+":thickness").c_str(),layerThicknessStr,"uom","cm");
			sprintf(valueStr,"%.2f",Xdata.Edata[ii].Rho);
			xmlWriteElement(layerNode,(namespaceCAAML+":density").c_str(),valueStr,"uom","kgm-3");
		}
	}//end density profile


	// lwc profile;
	pugi::xml_node lwcNode = node.append_child( (namespaceCAAML+":lwcProfile").c_str() );
	pugi::xml_node metaNode2 = lwcNode.append_child( (namespaceCAAML+":lwcMetaData").c_str() );
	xmlWriteElement(metaNode2,(namespaceCAAML+":methodOfMeas").c_str(),"other","","");
	if (!Xdata.Edata.empty()) {
		for (size_t ii = Xdata.Edata.size(); ii-->0;) {
			pugi::xml_node layerNode = lwcNode.append_child( (namespaceCAAML+":Layer").c_str() );
			sprintf(layerDepthTopStr,"%.4f",100*(Xdata.cH - Xdata.Ndata[ii+1].z)/Xdata.cos_sl); // cmb
			xmlWriteElement(layerNode,(namespaceCAAML+":depthTop").c_str(),layerDepthTopStr,"uom","cm");
			sprintf(layerThicknessStr,"%.4f",100*Xdata.Edata[ii].L/Xdata.cos_sl); // cmb
			xmlWriteElement(layerNode,(namespaceCAAML+":thickness").c_str(),layerThicknessStr,"uom","cm");
			sprintf(valueStr,"%.2f",Xdata.Edata[ii].theta[2]*100);
			xmlWriteElement(layerNode,(namespaceCAAML+":lwc").c_str(),valueStr,"uom","% by Vol");
		}
	}//end lwc profile


	// specSurfAreaProfile profile;
	pugi::xml_node ssaNode = node.append_child( (namespaceCAAML+":specSurfAreaProfile").c_str() );
	pugi::xml_node metaNode3 = ssaNode.append_child( (namespaceCAAML+":specSurfAreaMetaData").c_str() );
	xmlWriteElement(metaNode3,(namespaceCAAML+":methodOfMeas").c_str(),"other","","");
	if (!Xdata.Edata.empty()) {
		for (size_t ii = Xdata.Edata.size(); ii-->0;) {
			pugi::xml_node layerNode = ssaNode.append_child( (namespaceCAAML+":Layer").c_str() );
			sprintf(layerDepthTopStr,"%.4f",100*(Xdata.cH - Xdata.Ndata[ii+1].z)/Xdata.cos_sl); // cmb
			xmlWriteElement(layerNode,(namespaceCAAML+":depthTop").c_str(),layerDepthTopStr,"uom","cm");
			sprintf(layerThicknessStr,"%.4f",100*Xdata.Edata[ii].L/Xdata.cos_sl); // cmb
			xmlWriteElement(layerNode,(namespaceCAAML+":thickness").c_str(),layerThicknessStr,"uom","cm");
			double ogs = Xdata.Edata[ii].ogs;
			double rho = Xdata.Edata[ii].Rho;
			double ssa = 6.0/(rho*ogs/1000.0); //conversion from optical grain size to ssa
			sprintf(valueStr,"%.2f",ssa);
			xmlWriteElement(layerNode,(namespaceCAAML+":specSurfArea").c_str(),valueStr,"uom","m2kg-1");
		}
	}//end specSurfAreaProfile

	// strengthProfile;
	pugi::xml_node strengthNode = node.append_child( (namespaceCAAML+":strengthProfile").c_str() );
	pugi::xml_node metaNode4 = strengthNode.append_child( (namespaceCAAML+":strengthMetaData").c_str() );
	xmlWriteElement(metaNode4,(namespaceCAAML+":strengthType").c_str(),"shear","","");
	xmlWriteElement(metaNode4,(namespaceCAAML+":methodOfMeas").c_str(),"other","","");
	if (!Xdata.Edata.empty()) {
		for (size_t ii = Xdata.Edata.size(); ii-->0;) {
			pugi::xml_node layerNode = strengthNode.append_child( (namespaceCAAML+":Layer").c_str() );
			sprintf(layerDepthTopStr,"%.4f",100*(Xdata.cH - Xdata.Ndata[ii+1].z)/Xdata.cos_sl); // cmb
			xmlWriteElement(layerNode,(namespaceCAAML+":depthTop").c_str(),layerDepthTopStr,"uom","cm");
			sprintf(layerThicknessStr,"%.4f",100*Xdata.Edata[ii].L/Xdata.cos_sl); // cmb
			xmlWriteElement(layerNode,(namespaceCAAML+":thickness").c_str(),layerThicknessStr,"uom","cm");
			sprintf(valueStr,"%.2f",(Xdata.Edata[ii].s_strength)*1000); //conversion from kPa to Nm-2: *1000
			xmlWriteElement(layerNode,(namespaceCAAML+":strengthValue").c_str(),valueStr,"uom","Nm-2");
		}
	}//end strengthProfile
}

/**
 * @brief Write the station-data to the caaml-file.
 */
void CaaMLIO::writeStationData(pugi::xml_node& root, const SnowStation& Xdata)
{
	pugi::xml_node node = root.append_child( (namespaceCAAML+":locRef").c_str() );
	node.append_attribute("gml:id") = ("SLF_"+Xdata.meta.stationID+"_1").c_str();
	xmlWriteElement(node,(namespaceCAAML+":name").c_str(),(const char*) Xdata.meta.stationName.c_str(),"","");
	xmlWriteElement(node,(namespaceCAAML+":obsPointSubType").c_str(),"","","");

	pugi::xml_node elevNode = node.append_child( (namespaceCAAML+":validElevation").c_str() )
	                              .append_child( (namespaceCAAML+":ElevationPosition").c_str() );
	elevNode.append_attribute("uom") = "m";
	char elevStr[5];
	sprintf(elevStr,"%.0f",Xdata.meta.position.getAltitude());
	xmlWriteElement(elevNode,(namespaceCAAML+":position").c_str(),elevStr,"","");

	pugi::xml_node aspectNode = node.append_child( (namespaceCAAML+":validAspect").c_str() )
	                                .append_child( (namespaceCAAML+":AspectPosition").c_str() );
	xmlWriteElement(aspectNode,(namespaceCAAML+":position").c_str(),IOUtils::bearing(Xdata.meta.getAzimuth()).c_str(),"","");

	if(Xdata.meta.getSlopeAngle()!=IOUtils::nodata){
		pugi::xml_node slopeNode = node.append_child( (namespaceCAAML+":validSlopeAngle").c_str() )
		                               .append_child( (namespaceCAAML+":SlopeAnglePosition").c_str() );
		slopeNode.append_attribute("uom") = "deg";
		char slopeStr[4];
		sprintf(slopeStr,"%.0f",Xdata.meta.getSlopeAngle());
		xmlWriteElement(slopeNode,(namespaceCAAML+":position").c_str(),slopeStr,"","");
	}
	pugi::xml_node pointNode = node.append_child( (namespaceCAAML+":pointLocation").c_str() )
								   .append_child( (namespaceGML+":Point").c_str() );
	pointNode.append_attribute("gml:id") = ("SLF_"+Xdata.meta.stationID+"_2").c_str();
	pointNode.append_attribute("srsName") = "urn:ogc:def:crs:OGC:1.3:CRS84";
	pointNode.append_attribute("srsDimension") = "2";
	char posStr[30];
	sprintf(posStr,"%f %f",Xdata.meta.position.getLon(),Xdata.meta.position.getLat());
	xmlWriteElement(pointNode,"gml:pos",posStr,"","");
}

void CaaMLIO::writeTimeSeries(const SnowStation& /*Xdata*/, const SurfaceFluxes& /*Sdata*/, const CurrentMeteo& /*Mdata*/,
                              const ProcessDat& /*Hdata*/, const double /*wind_trans24*/)
{
	throw IOException("Nothing implemented here!", AT);
}

void CaaMLIO::writeProfile(const Date& /*date*/, const SnowStation& /*Xdata*/)
{
	throw IOException("Nothing implemented here!", AT);
}

bool CaaMLIO::writeHazardData(const std::string& /*stationID*/, const std::vector<ProcessDat>& /*Hdata*/,
                             const std::vector<ProcessInd>& /*Hdata_ind*/, const size_t& /*num*/)
{
	throw IOException("Nothing implemented here!", AT);
}

/**
 * @brief Convert from liquid water content code to value
 * @author Adrien Gaudard
 * @param code Liquid water content code (one character)
 * return Liquid water content value (fraction)
 */
double CaaMLIO::lwc_codeToVal(const char* code)
{
	if (!strcmp(code,"D")) return 0.;
	if (!strcmp(code,"M")) return 0.01;
	if (!strcmp(code,"W")) return 0.03;
	if (!strcmp(code,"V")) return 0.08;
	if (!strcmp(code,"S")) return 0.15;

	throw IOException("Unrecognized liquid water content code.", AT);
}

/**
 * @brief Convert from liquid water content value to code
 * @author Adrien Gaudard
 * @param val Liquid water content value (fraction)
 * return Liquid water content code (one character)
 */
std::string CaaMLIO::lwc_valToCode(const double val)
{
	if (val == 0.00) return "D";
	if (val <  0.03) return "M";
	if (val <  0.08) return "W";
	if (val <  0.15) return "V";
	if (val <  1.00) return "S";

	throw IOException("Invalid liquid water content value.", AT);
}

/**
 * @brief Convert from hardness code to value
 * @author Adrien Gaudard
 * @param code Hardness code
 * return Hardness value (1 to 6)
 */
double CaaMLIO::hardness_codeToVal(char* code)
{
	const std::string codeString(code);
	if( codeString == "n/a")
		return IOUtils::nodata;

	double val = 0.;
	unsigned int n = 0;
	char* c[2];
	c[0] = strtok(code,"-");
	c[1] = strtok(NULL,"-");

	for (size_t i=0; i<2; i++) {
		if (c[i]) {
			n++;
			if (!strcmp(c[i],"F")) {
				val += 1.;
			} else if (!strcmp(c[i],"4F")) {
				val += 2.;
			} else if (!strcmp(c[i],"1F")) {
				val += 3.;
			} else if (!strcmp(c[i],"P")) {
				val += 4.;
			} else if (!strcmp(c[i],"K")) {
				val += 5.;
			} else if (!strcmp(c[i],"I")) {
				val += 6.;
			} else {
				throw IOException("Unrecognized hardness code.", AT);
			}
		}
	}
	return val/n;
}

/**
 * @brief Convert from hardness value to code
 * @author Adrien Gaudard
 * @param val Hardness value (1 to 6)
 * return Hardness code
 */
std::string CaaMLIO::hardness_valToCode(const double val)
{
	if (val <= 1.0) return "F";
	if (val <= 1.5) return "F-4F";
	if (val <= 2.0) return "4F";
	if (val <= 2.5) return "4F-1F";
	if (val <= 3.0) return "1F";
	if (val <= 3.5) return "1F-P";
	if (val <= 4.0) return "P";
	if (val <= 4.5) return "P-K";
	if (val <= 5.0) return "K";
	if (val <= 5.5) return "K-I";
	if (val <= 6.0) return "I";
	if (val == IOUtils::nodata) return "n/a";
	std::cout<< "Hardness value: " << val << std::endl;
	throw IOException("Unrecognized hardness value: ", AT);
}

/**
 * @brief Convert from grain shape code to values (sphericity, dendricity, marker)
 * @author Adrien Gaudard
 * @param[in] code Grain shape code
 * @param[out] sp sphericity
 * @param[out] dd dendricity
 * @param[out] mk micro-structure marker
 */
void CaaMLIO::grainShape_codeToVal(const std::string& code, double &sp, double &dd, unsigned short int &mk)
{
	//first check some special grain shapes (with 4 letters):
	if (code=="PPgp") {
		sp = 1.; dd = 0.; mk = 4;
	} else {
		//otherwise take only the first two letters of the code
		const std::string code2Letters = code.substr(0,2);
		if (code2Letters=="PP") {
			sp = 0.5; dd = 1.; mk = 0;
		} else if (code2Letters=="DF") {
			sp = 0.5; dd = 0.5; mk = 0;
		} else if (code2Letters=="RG") {
			sp = 1.; dd = 0.; mk = 2; //why marker=2?
		} else if (code2Letters=="FC") {
			sp = 0.2; dd = 0.; mk = 1;
		} else if (code2Letters=="DH") { //FC or DH is distinguished by grain size. A message is given if there is DH with grain size < 1.5mm. or if FC with grain size > 1.5 mm!
			sp = 0.2; dd = 0.; mk = 1;
		} else if (code2Letters=="SH") {
			sp = 0.; dd = 0.; mk = 3;
		} else if (code2Letters=="MF") {
			sp = 1.; dd = 0.; mk = 12;
		} else if (code2Letters=="IF") {
			sp = 1.; dd = 0.; mk = 7;
		} else if (code2Letters=="MM") {
			sp = 1.; dd = 0.; mk = 2; //which values make sense here???
		} else {
			throw IOException("Unrecognized grain shape code.", AT);
		}
	}
}


/**
 * @brief Convert from grain shape value to code
 * @author Adrien Gaudard
 * @param var Grain shape value
 * return Grain shape code
 */
std::string CaaMLIO::grainShape_valToAbbrev(const unsigned int var)
{
	if (var == 0) return "PPgp";
	if (var == 1) return "PP";
	if (var == 2) return "DF";
	if (var == 3) return "RG";
	if (var == 4) return "FC";
	if (var == 5) return "DH";
	if (var == 6) return "SH";
	if (var == 7) return "MF";
	if (var == 8) return "IF";
	if (var == 9) return "FCxr";

	throw IOException("Unrecognized grain shape code.", AT);
}

/**
 * @brief Convert from grain shape values (sphericity, dendricity, marker) to two-character code
 * @author Adrien Gaudard
 * @param var Grain shape values (sphericity, dendricity, marker) (AMBIGUOUS)
 * return Grain shape two-character code (NOT ALL REPRESENTED)
 */
std::string CaaMLIO::grainShape_valToAbbrev_old(const double* var)
{
	const double sp = ((int)(var[0]*10+0.5))/10.;
	const double dd = ((int)(var[1]*10+0.5))/10.;
	const double mk = ((int)(var[2]*10+0.5))/10.;

	if (sp == 0.5 && dd == 1. && mk == 0.) return "PP";
	if (sp == 0.5 && dd == 0.5 && mk == 0.) return "DF";
	if (sp == 1. && dd == 0. && (mk == 2. || mk == 12.)) return "RG";
	if (sp == 0. && dd == 0. && mk == 1.) return "FC";
	if (sp == 0. && dd == 0. && mk == 1.) return "DH";
	if (sp == 0. && dd == 0. && mk == 1.) return "SH";
	if (sp == 1. && dd == 0. && mk == 2.) return "MF";
	if (sp == 1. && dd == 0. && mk == 2.) return "IF";

	throw IOException("Unrecognized set of grain shape values.", AT);
}

/**
 * @brief Add a child-xml-element with a value and optionally an attribute to an xml-element. Example:
 *        Current xml-element:                       <caaml:Layer>
 *        New child-element which will be added:         <caaml:depthTop uom="cm">1.91</caaml:depthTop>
 * @param node pointing to the xml-element where we want to add a child-element. e.g. "caaml:Layer"
 * @param name name of the new xml-element. e.g. "caaml:depthTop"
 * @param content value of the new xml-elemnt. e.g. "1.91"
 * @param att_name name of the attribute. e.g. "uom"
 * @return att_val value of the attribute. e.g. "cm"
 */
void CaaMLIO::xmlWriteElement(pugi::xml_node& node, const char* name, const char* content, const char* att_name, const char* att_val)
{
	pugi::xml_node child = node.append_child( name );
	child.append_child(pugi::node_pcdata).set_value( content );
	if (strcmp(att_name,"")) //ie: string not empty
		child.append_attribute( att_name ) = att_val;
}

/**
 * @brief Check if a certain path exists in an xml-file.
 * @param path path to check. e.g. "/caaml:SnowProfile/caaml:locRef"
 * return true if the path exists
 */
bool CaaMLIO::xmlDoesPathExist(const std::string& path)
{
	const pugi::xml_node node =  inDoc.first_element_by_path(path.c_str());
	if( node.empty() ){
		return false;
	}
	return true;
}

/**
 * @brief Read a double-value from a path.
 * @param path e.g. "/caaml:SnowProfile/caaml:snowProfileResultsOf/.../caaml:customData/snp"
 * @param property e.g. "Albedo"
 * @param dflt default-value to return if the path does not exist in the xml-file
 * return the value which was read in (or the default-value)
 */
double CaaMLIO::xmlReadValueFromPath (const string& xpath, const std::string& property, const double& dflt)
{
	double val = IOUtils::nodata;
	const std::string path( SnowData_xpath+xpath+":"+property );
	const pugi::xml_node node =  inDoc.first_element_by_path(path.c_str());
	const char* valStr = node.child_value();
	if (node.empty() ){
		val=dflt;
	}else{
		sscanf(valStr, "%lf", &val);
	}
	//std::cout << property << ": " << val << "." << std::endl;
	return val;
}

/**
 * @brief Read an int-value from a path.
 * @param path e.g. "/caaml:SnowProfile/caaml:snowProfileResultsOf/.../caaml:customData/snp"
 * @param property e.g. "Albedo"
 * @param dflt default-value to return if the path does not exist in the xml-file
 * return the value which was read in (or the default-value)
 */
int CaaMLIO::xmlReadValueFromPath (const string& xpath, const std::string& property, const int& dflt)
{
	int val = IOUtils::inodata;
	const std::string path( SnowData_xpath+xpath+":"+property );
	const pugi::xml_node node =  inDoc.first_element_by_path(path.c_str());
	const char* valStr = node.child_value();
	if (node.empty() ){
		val=dflt;
	}else{
		sscanf(valStr, "%d", &val);
	}
	//std::cout << property << ": " << val << "." << std::endl;
	return val;
}

/**
 * @brief Read an attribute from a path.
 * @param path e.g.  "/caaml:SnowProfile/caaml:locRef"
 * @param attributeName e.g. "gml:id"
 * return the value of the attribute which was read in e.g. "OPERATION_ID"
 */
std::string CaaMLIO::xmlReadAttributeFromPath (const string& path, const std::string& attributeName)
{
	const pugi::xml_node node =  inDoc.first_element_by_path(path.c_str());
	const std::string attrStr(node.attribute(attributeName.c_str()).as_string());
	//std::cout << attributeName << ": " << attrStr << "." << std::endl;
	return attrStr;
}

/**
 * @brief Read a value from an xml-element (node) and convert the unit if wished.
 * @param node pointing to the xml-element we want to read the value from
 * @param propertyName only read the value from the xml-element if its name is equal to this parameter
 * @param variableToSet output
 * @param unitOut wished unit of the variable
 * @param unitMeasured the (current) unit of the variable
 * @param factor this factor is applied to the variable
 * return true if value was read successfully
 */
bool CaaMLIO::xmlReadValueFromNode (const pugi::xml_node node, const std::string propertyName, double& variableToSet,
							        const std::string unitOut,const std::string unitMeasured, const double factor)
{
	const std::string fieldName( node.name() );
	if (fieldName == propertyName){
		double temp;
		sscanf((const char*) node.child_value(),"%lf",&temp);
		if (unitOut != ""){
			char* unitOfMeasurement = (char*)unitMeasured.c_str();
			if(unitMeasured==""){
				unitOfMeasurement = (char*) node.attribute("uom").as_string();
			}
			temp = unitConversion(temp,unitOfMeasurement,(char*) unitOut.c_str());
		}
		variableToSet = temp*factor;
		//std::cout << propertyName << ": " << variableToSet << unitOut << std::endl;
		return true;
	}
	return false;
}

/**
 * @brief Print all the important data which was read in from the caaml-file. Just to check if
 *        everything was read in correctly.
 * @param SSdata contains all the data which was read in from the caaml-file.
 */
bool CaaMLIO::checkWhatWasReadIn(SN_SNOWSOIL_DATA& SSdata)
{
	std::cout << "Summary of all snow-soil-data (from Caaml-file): " << std::endl;
	std::cout << SSdata.toString() << std::endl;
	return true;
}
