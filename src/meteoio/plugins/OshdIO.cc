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
#include <meteoio/plugins/OshdIO.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/FileUtils.h>

#include <matio.h>
#include <algorithm>

using namespace std;

namespace mio {
/**
 * @page oshd OshdIO
 * This plugin reads the meteorological forecast data downscaled for each of the Swiss meteorological networks IMIS/ANETZ stations
 * as preprocessed by the <A HREF="www.wsl.ch/fe/gebirgshydrologie/schnee_hydro/oshd/index_EN">Operational Snow-Hydrological Service</A> 
 * of the <A HREF="www.wsl.ch">WSL/SLF</A>. The data is written as Matlab
 * <A HREF="http://www.mathworks.com/help/pdf_doc/matlab/matfile_format.pdf">binary files (.mat)</A>, one per meteorological parameter and per timestep, 
 * available on an access-controlled server after each new <A HREF="www.cosmo-model.org/">COSMO</A> run. It therefore requires a third party 
 * library to read this file format: the Open Source <A HREF="https://sourceforge.net/projects/matio/">MatIO</A> library. This can be installed directly from
 * the default repositories under Linux or installed by downloading the proper package for Windows or OsX.
 * 
 * \note If non-ascii characters have been used and the file has been created under Windows, some of the strings might end up using the UTF-16 encoding.
 * This requires a recent version of libmatio (see <A HREF="https://github.com/tbeu/matio/issues/34">this issue</A>). Another option would be to 
 * add at the begining of the Matlab routine a call to *feature('DefaultCharacterSet', 'UTF8')* in order to switch from the current default (which can be read by the same call, 
 * ommitting the 'UTF8' option) to the (<A HREF="http://blog.omega-prime.co.uk/?p=150">partial</A>) UTF-8  encoding of Matlab.
 * 
 * @section oshd_data_structure Data structure
 * The files are named with the following schema: <i>{parameter}_{timestep}_{cosmo model version}_F_{run time}.mat</i> with the following possible values:
 *     + *parameter* is one of idfc, idrc, albd, ilwc, pair, prec, rcor, tcor, wcor, wdir;
 *     + *timestep* is written as purely numeric ISO with minute resolution;
 *     + *cosmo model version* could be any of cosmo7, cosmo2, cosmo1, cosmoE;
 *     + *run time* is the purely numeric ISO date and time of when COSMO produced the dataset.
 * 
 * The files have the following internal data structure (represented as "name {data type}"):
 * @verbatim
      stat {1x1 struct}
        ├── time {1x1 array of doubles}
        ├── data {1x623 array of doubles}
        ├── acro {1x623 array of arrays of char}
        ├── dunit {array of char}
        ├── type {array of char}
        └── name {array of char}
  @endverbatim
 * 
 * The stations' acronyms follow a fixed order but their coordinates must be provided in a separate file, given as *METAFILE* key (see below). This file
 * must have the following structure (the *x* and *y* coordinates being the CH1903 easting and northing, respectively): 
 * @verbatim
      statlist {1x1 struct}
        ├── acro {1x623 array of arrays of char}
        ├── name {1x623 array of arrays of char}
        ├── x {1x623 array of doubles}
        ├── y {1x623 array of doubles}
        └── z {1x623 array of doubles}
  @endverbatim
 *
 *
 * @section oshd_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - METEOPATH: directory containing all the data files with the proper file naming schema; [Input] section
 * - METEOPATH_RECURSIVE: should *meteopath* be searched recursively for files? (default: no); [Input] section
 * - STATION#: input stations' IDs (in METEOPATH). As many meteofiles as needed may be specified
 * - METAFILE: file containing the stations' IDs, names and location; [Input] section (either within METEOPATH if not path is 
 provided or within the provided path)
 * - OSHD_DEBUG: write out extra information to better show what is in the files
 *
 * @section oshd_example Example use
 * @code
 * [Input]
 * METEO = OSHD
 * METEOPATH = /local/LATEST_03h_RUN
 * METEOPATH_RECURSIVE = true
 * METAFILE  = STAT_LIST.mat ;another possibility could be /local/metadata/STAT_LIST.mat
 * STATION1  = ATT2
 * STATION2  = WFJ2
 * @endcode
 *
 */

/**********************************************************************************
 * Here we define some wrappers around libmatio. These are not declared as class methods in
 * order to avoid having to expose matio.h when including OshdIO.h
 **********************************************************************************/
void listFields(matvar_t *matvar)
{
	const unsigned int nrFields = Mat_VarGetNumberOfFields(matvar);
	char * const *fields = Mat_VarGetStructFieldnames(matvar);
	for (unsigned int ii=0; ii<nrFields; ii++) 
		printf("field[%d] = %s\n", ii, fields[ii]);
}

void printStructure(matvar_t *matvar)
{
	//Mat_VarPrint(field, 0);
	printf("name=%s class_type=%d data_type=%d rank=%d", matvar->name, matvar->class_type, matvar->data_type, matvar->rank);
	for (int ii=0; ii<matvar->rank; ii++) 
		printf("\tdims[%d]=%d", ii, (int)matvar->dims[ii]);
	printf("\n");
}


std::string readString(const std::string &filename, const std::string &fieldname, mat_t *matfp, matvar_t *matvar)
{
	matvar_t *field = Mat_VarGetStructFieldByName(matvar, fieldname.c_str(), 0);
	if (matvar==NULL) throw NotFoundException("could not read field '"+fieldname+"' in file '"+filename+"'", AT);
	if (Mat_VarReadDataAll(matfp, field)) 
		throw InvalidFormatException("could not read field '"+fieldname+"' in file '"+filename+"'", AT);
	if (field->class_type!=MAT_C_CHAR) throw InvalidFormatException("field '"+fieldname+"' in file '"+filename+"' is not a type string", AT);
	
	return std::string( static_cast<char*>(field->data) );
}

void readStringVector(const std::string &filename, const std::string &fieldname, mat_t *matfp, matvar_t *matvar, std::vector<std::string> &vecString)
{
	vecString.clear();
	
	matvar_t *field = Mat_VarGetStructFieldByName(matvar, fieldname.c_str(), 0);
	if (matvar==NULL) 	throw NotFoundException("could not read field '"+fieldname+"' in file '"+filename+"'", AT);
	if (Mat_VarReadDataAll(matfp, field)) 
		throw InvalidFormatException("could not read field '"+fieldname+"' in file '"+filename+"'", AT);
	
	if (field->class_type!=MAT_C_CELL) throw InvalidFormatException("field '"+fieldname+"' in file '"+filename+"' is not a cell type", AT);
	if (field->data_type!=MAT_T_CELL) throw InvalidFormatException("field '"+fieldname+"' in file '"+filename+"' is not a cell array data type", AT);
	
	matvar_t *cell = Mat_VarGetCell(matvar, 0);
	if (cell==NULL) throw InvalidFormatException("could not read data in field '"+fieldname+"' in file '"+filename+"'", AT);
	if (field->rank!=2) throw InvalidFormatException("invalid rank for field '"+fieldname+"' in file '"+filename+"'", AT);
	
	const size_t nrows = field->dims[0];
	const size_t ncols = field->dims[1];
	if (nrows!=1) throw InvalidFormatException("invalid nrows for field '"+fieldname+"' in file '"+filename+"'", AT);
	
	vecString.resize( ncols );
	for (size_t ii=0; ii<ncols; ii++) {
		cell = Mat_VarGetCell(field, static_cast<int>(ii));
		if (cell->rank!=2) throw InvalidFormatException("invalid cell rank in file '"+filename+"'", AT);
		if (cell->class_type!=MAT_C_CHAR) throw InvalidFormatException("field '"+fieldname+"' in file '"+filename+"' is not a type string", AT);
		vecString[ii] = static_cast<char*>(cell->data);
	}
}

void readDoubleVector(const std::string &filename, const std::string &fieldname, mat_t *matfp, matvar_t *matvar, std::vector<double> &vecData)
{
	vecData.clear();
	
	matvar_t *field = Mat_VarGetStructFieldByName(matvar, fieldname.c_str(), 0);
	if (matvar==NULL) 	throw NotFoundException("could not read field '"+fieldname+"' in file '"+filename+"'", AT);
	if (Mat_VarReadDataAll(matfp, field)) 
		throw InvalidFormatException("could not read field '"+fieldname+"' in file '"+filename+"'", AT);
	
	if (field->class_type!=MAT_C_DOUBLE) throw InvalidFormatException("field '"+fieldname+"' in file '"+filename+"' is not a double type", AT);
	if (field->rank!=2) throw InvalidFormatException("invalid rank for field '"+fieldname+"' in file '"+filename+"'", AT);
	
	const size_t nrows = field->dims[0];
	const size_t ncols = field->dims[1];
	if (nrows!=1) throw InvalidFormatException("invalid nrows for field '"+fieldname+"' in file '"+filename+"'", AT);
	
	const double* matData( static_cast<double*>( field->data ) );
	vecData.resize( ncols );
	for (size_t ii=0; ii<ncols; ii++) {
		vecData[ii] = matData[ii];
	}
}

void printFileStructure(const std::string& filename)
{
	mat_t *matfp = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
	if ( NULL == matfp ) throw AccessException(filename, AT);

	std::cout << "<" << FileUtils::getFilename( filename ) << ">\n";
	matvar_t *matvar;
	while ( (matvar = Mat_VarReadNextInfo(matfp)) != NULL ) {
		std::cout << "\t" << matvar->name << " [";
		for (int ii=0; ii<matvar->rank; ii++) {
			std::cout << (int)matvar->dims[ii];
			if (ii<(matvar->rank-1)) std::cout << "x";
		}
		std::cout << "]\n";

		const unsigned int nrFields = Mat_VarGetNumberOfFields(matvar);
		char * const *fields = Mat_VarGetStructFieldnames(matvar);
		for (unsigned int ii=0; ii<nrFields; ii++) {
			const std::string field_name( fields[ii] );
			matvar_t *field = Mat_VarGetStructFieldByName(matvar, field_name.c_str(), 0);
			const std::string prefix = (ii<(nrFields-1))? "├──" : "└──";
			std::cout << "\t" << prefix << field_name;
			if (field->class_type==MAT_C_CHAR)
				std::cout << " = \"" << readString(filename, field_name, matfp, matvar) << "\"";
			if (field->class_type==MAT_C_DOUBLE) {
				std::cout << " [";
				size_t count=1;
				for (int jj=0; jj<field->rank; jj++) {
					std::cout << field->dims[jj];
					if (jj<(field->rank-1)) std::cout << "x";
					count *= field->dims[jj];
				}
				std::cout << "]";
				if (count==1) {
					if (Mat_VarReadDataAll(matfp, field))
						throw InvalidFormatException("could not read field '"+field_name+"' in file '"+filename+"'", AT);
					const double val = static_cast<double*>(field->data)[0];
					if (field_name=="time") {
						Date timestep;
						timestep.setMatlabDate( val, OshdIO::in_dflt_TZ );
						std::cout << " = " << timestep.toString(Date::ISO_TZ);
					} else
						std::cout << " = " << val;
				}
			}
			if (field->class_type==MAT_C_CELL) {
				std::cout << " [";
				for (int jj=0; jj<field->rank; jj++) {
					std::cout << field->dims[jj];
					if (jj<(field->rank-1)) std::cout << "x";
				}
				std::cout << "]";
			}

			std::cout << "\n";
		}
	}
	std::cout << "</" << FileUtils::getFilename( filename ) << ">\n\n";
	Mat_VarFree(matvar);
	matvar = NULL;
	Mat_Close(matfp);
}

/**********************************************************************************
 * Now really implementing the OshdIO class
 **********************************************************************************/
const char* OshdIO::meteo_ext = ".mat";
const double OshdIO::in_dflt_TZ = 0.; //COSMO data is always GMT

OshdIO::OshdIO(const std::string& configfile) : cfg(configfile), cache_meteo_files(), vecMeta(), vecIDs(), params_map(), vecIdx(), 
               in_meteopath(), in_metafile(), debug(false)
{
	parseInputOutputSection();
}

OshdIO::OshdIO(const Config& cfgreader) : cfg(cfgreader), cache_meteo_files(), vecMeta(), vecIDs(), params_map(), vecIdx(), 
               in_meteopath(), in_metafile(), debug(false)
{
	parseInputOutputSection();
}

void OshdIO::parseInputOutputSection()
{
	cfg.getValue("OSHD_DEBUG", "INPUT", debug, IOUtils::nothrow);
	//IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	
	cfg.getValues("STATION", "INPUT", vecIDs);
	cfg.getValue("METEOPATH", "Input", in_meteopath);
	bool is_recursive = false;
	cfg.getValue("METEOPATH_RECURSIVE", "Input", is_recursive, IOUtils::nothrow);
	scanMeteoPath(in_meteopath, is_recursive, cache_meteo_files);
	
	cfg.getValue("METAFILE", "INPUT", in_metafile);
	if (FileUtils::getFilename(in_metafile) == in_metafile) { //ie there is no path in the provided filename
		in_metafile = in_meteopath + "/" + in_metafile;
	}
	
	//fill the params mapping vector
	params_map.push_back( std::make_pair(MeteoData::ILWR, "ilwc") );
	params_map.push_back( std::make_pair(MeteoData::P, "pair") );
	params_map.push_back( std::make_pair(MeteoData::PSUM, "prec") ); //in mm/ts
	params_map.push_back( std::make_pair(MeteoData::RH, "rcor") ); //old: rhum
	params_map.push_back( std::make_pair(MeteoData::TA, "tcor") ); //old:tair
	params_map.push_back( std::make_pair(MeteoData::VW, "wcor") ); //old: wind
	params_map.push_back( std::make_pair(MeteoData::DW, "wdir") );
}

void OshdIO::scanMeteoPath(const std::string& meteopath_in, const bool& is_recursive,  std::vector< struct file_index > &meteo_files)
{
	meteo_files.clear();

	const std::list<std::string> dirlist( FileUtils::readDirectory(meteopath_in, "prec", is_recursive) ); //we consider that if we have found one parameter, the others are also there

	std::map<std::string, size_t> mapIdx; //make sure each timestamp only appears once, ie remove duplicates
	for (std::list<std::string>::const_iterator it = dirlist.begin(); it != dirlist.end(); ++it) {
		const std::string file_and_path( *it );
		const std::string filename( FileUtils::getFilename(file_and_path) );

		//we need to split the file name into its components: parameter, date, run_date
		const std::string::size_type pos_param = filename.find('_');
		if (pos_param==string::npos) continue;
		const std::string::size_type date_start = filename.find_first_of("0123456789");
		if (date_start==string::npos) continue;
		const std::string::size_type date_end = filename.find('_', date_start);
		if (date_end==string::npos) continue;
		const std::string::size_type rundate_start = filename.rfind('_');
		if (rundate_start==string::npos) continue;
		const std::string::size_type rundate_end = filename.rfind('.');
		if (rundate_end==string::npos) continue;

		const std::string date_str( filename.substr(date_start, date_end-date_start) );
		const std::string run_date( filename.substr(rundate_start+1, rundate_end-rundate_start-1) );

		//do we already have an entry for this date?
		size_t idx = IOUtils::npos;
		const std::map<std::string, size_t>::const_iterator it_map = mapIdx.find( date_str );
		if (it_map!=mapIdx.end()) {
			idx = it_map->second;
			if (meteo_files[idx].run_date>run_date) continue;
		}

		//we don't have an entry or it is too old -> create new entry / replace existing one
		const std::string path( FileUtils::getPath(file_and_path) );
		Date date;
		IOUtils::convertString(date, date_str, in_dflt_TZ);
		const file_index elem(date, path, filename.substr(pos_param), run_date);
		if (idx==IOUtils::npos) {
			meteo_files.push_back( elem );
			mapIdx[ date_str ] = meteo_files.size()-1;
		} else {
			meteo_files[ idx] = elem;
			mapIdx[ date_str ] = idx;
		}
	}

	std::sort(meteo_files.begin(), meteo_files.end());
}

size_t OshdIO::getFileIdx(const Date& start_date) const
{
	if (cache_meteo_files.empty())
		throw InvalidArgumentException("No input files found or configured!", AT);

	//find which file we should open
	if (cache_meteo_files.size()==1) {
		return 0;
	} else {
		for (size_t idx=1; idx<cache_meteo_files.size(); idx++) {
			if (start_date>=cache_meteo_files[idx-1].date && start_date<cache_meteo_files[idx].date) {
				return --idx;
			}
		}

		//not found, we take the closest timestamp we have
		if (start_date<cache_meteo_files.front().date)
			return 0;
		else
			return cache_meteo_files.size()-1;
	}
}

void OshdIO::readStationData(const Date& /*date*/, std::vector<StationData>& vecStation)
{
	vecStation.clear();
	if (vecMeta.empty()) fillStationMeta();
	vecStation = vecMeta;
}

void OshdIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                             std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();
	size_t file_idx = getFileIdx( dateStart );
	Date station_date( cache_meteo_files[file_idx].date );
	if (station_date<dateStart || station_date>dateEnd) return; //the requested period is NOT in the available files

	const size_t nr_files = cache_meteo_files.size();
	const size_t nrIDs = vecIDs.size();
	
	if (vecMeta.empty()) fillStationMeta(); //this also fills vecIdx
	vecMeteo.resize( nrIDs );
	do {
		//create empty MeteoData for the current timestep
		for (size_t jj=0; jj<nrIDs; jj++) {
			const MeteoData md( station_date, vecMeta[jj] );
			vecMeteo[jj].push_back( md );
		}
		
		//read the data and fill vecMeteo
		const std::string file_suffix( cache_meteo_files[ file_idx ].file_suffix );
		const std::string path( in_meteopath + "/" + cache_meteo_files[ file_idx ].path );
		std::vector<double> vecData;
		for (size_t ii=0; ii<params_map.size(); ii++) {
			const MeteoData::Parameters param( params_map[ii].first );
			const std::string filename( path + "/" + params_map[ii].second + file_suffix );
			vecData.resize( nrIDs, IOUtils::nodata );
			readFromFile(filename, param, station_date, vecData);

			for (size_t jj=0; jj<nrIDs; jj++)
				vecMeteo[jj].back()( param ) =  vecData[jj];
		}
		
		readSWRad(station_date, path, file_suffix, nrIDs, vecMeteo); //the short wave radiation is a little different...
		readPPhase(station_date, path, file_suffix, nrIDs, vecMeteo); //the precipitation phase is a little different...

		file_idx++;
		station_date = ((file_idx)<nr_files)? cache_meteo_files[file_idx].date : dateEnd+1.;
	} while (file_idx<nr_files && station_date<=dateEnd);
}

void OshdIO::readSWRad(const Date& station_date, const std::string& path, const std::string& file_suffix, const size_t& nrIDs, std::vector< std::vector<MeteoData> >& vecMeteo) const
{
	std::vector<double> vecDir;
	vecDir.resize( nrIDs, IOUtils::nodata );
	const std::string filename_dir( path + "/" + "idrc" + file_suffix );
	readFromFile(filename_dir, MeteoData::ISWR, station_date, vecDir);
	
	std::vector<double> vecDiff;
	vecDiff.resize( nrIDs, IOUtils::nodata );
	const std::string filename_diff( path + "/" + "idfc" + file_suffix );
	readFromFile(filename_diff, MeteoData::ISWR, station_date, vecDiff);

	std::vector<double> vecAlbd;
	const std::string filename_albd( path + "/" + "albd" + file_suffix );
	if (FileUtils::fileExists(filename_albd)) {
		vecAlbd.resize( nrIDs, IOUtils::nodata );
		readFromFile(filename_albd, MeteoData::RSWR, station_date, vecAlbd); //We read ALBD and use it to build RSWR
	}
	
	const double albedo = !vecAlbd.empty();
	for (size_t jj=0; jj<nrIDs; jj++) {
		vecMeteo[jj].back()( MeteoData::ISWR ) =  vecDir[jj]+vecDiff[jj];
		if (albedo) vecMeteo[jj].back()( MeteoData::RSWR ) =  (vecDir[jj]+vecDiff[jj])*vecAlbd[jj];
	}
}

void OshdIO::readPPhase(const Date& station_date, const std::string& path, const std::string& file_suffix, const size_t& nrIDs, std::vector< std::vector<MeteoData> >& vecMeteo) const
{
	const std::string filename( path + "/" + "snfl" + file_suffix );

	if (FileUtils::fileExists(filename)) {
		static const double half_elevation_band = 50.;  //we consider that there are mixed precip in the elevation range snow_line ± half_elevation_band
		std::vector<double> vecSnowLine;
		vecSnowLine.resize( nrIDs, IOUtils::nodata );
		readFromFile(filename, MeteoData::PSUM_PH, station_date, vecSnowLine);

		for (size_t jj=0; jj<nrIDs; jj++) {
			const double altitude = vecMeteo[jj].front().meta.getAltitude();
			if (altitude>(vecSnowLine[jj]+half_elevation_band))
				vecMeteo[jj].back()( MeteoData::PSUM_PH ) = 0.;
			else if (altitude<(vecSnowLine[jj]-half_elevation_band))
				vecMeteo[jj].back()( MeteoData::PSUM_PH ) = 1.;
			else
				vecMeteo[jj].back()( MeteoData::PSUM_PH ) = .5;
		}
	}
}

void OshdIO::readFromFile(const std::string& filename, const MeteoData::Parameters& param, const Date& in_timestep, std::vector<double> &vecData) const
{
	if (debug) printFileStructure(filename);
	mat_t *matfp = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
	if ( NULL == matfp ) throw AccessException(filename, AT);

	//open the file and read some metadata
	matvar_t *matvar = Mat_VarReadInfo(matfp, "stat");
	if (matvar==NULL) throw NotFoundException("structure 'stat' not found in file '"+filename+"'", AT);
	if (matvar->class_type!=MAT_C_STRUCT) throw InvalidFormatException("The matlab file should contain 1 structure", AT);
	
	const std::string type( readString(filename, "type", matfp, matvar) );
	checkFieldType(param, type);

	//check that the timestep is as expected
	std::vector<double> vecTime;
	readDoubleVector(filename, "time", matfp, matvar, vecTime);
	if (vecTime.size()!=1) throw InvalidFormatException("one and only one time step must be present in the 'time' vector", AT);
	Date timestep;
	timestep.setMatlabDate( vecTime[0], in_dflt_TZ );
	if (in_timestep!=timestep) throw InvalidArgumentException("the in-file timestep and the filename time step don't match for for '"+filename+"'", AT);
	
	//check that each station is still at the same index, build the index cache if necessary
	std::vector<std::string> vecAcro;
	readStringVector(filename, "acro", matfp, matvar, vecAcro);
	const size_t nrIDs = vecIDs.size();
	for (size_t ii=0; ii<nrIDs; ii++) { //check that the IDs still match
		if (vecIDs[ii] != vecAcro[ vecIdx[ii] ])
			throw InvalidFormatException("station '"+vecIDs[ii]+"' is not listed in the same position as previously in file '"+filename+"'", AT);
	}
	
	//extract the data for the selected stations
	const std::string units( readString(filename, "dunit", matfp, matvar) );
	std::vector<double> vecRaw;
	readDoubleVector(filename, "data", matfp, matvar, vecRaw);
	if (vecAcro.size() != vecRaw.size()) throw InvalidFormatException("'acro' and 'data' arrays don't match in file '"+filename+"'", AT);
	for (size_t ii=0; ii<nrIDs; ii++)
		vecData[ii] = convertUnits( vecRaw[ vecIdx[ii] ], units, param);
	
	Mat_VarFree(matvar);
	Mat_Close(matfp);
}

void OshdIO::checkFieldType(const MeteoData::Parameters& param, const std::string& type)
{
	if (param==MeteoData::TA && type=="TA") return;
	if (param==MeteoData::RH && type=="RH") return;
	if (param==MeteoData::PSUM && type=="PREC") return;
	if (param==MeteoData::VW && type=="WS") return;
	if (param==MeteoData::DW && type=="WD") return;
	if (param==MeteoData::ILWR && type=="LWR") return;
	if (param==MeteoData::ISWR && type=="SWR") return;
	if (param==MeteoData::P && type=="other") return;
	if (param==MeteoData::PSUM_PH && type=="other") return;
	if (param==MeteoData::RSWR && type=="other") return; //this is in fact ALBD
	
	throw InvalidArgumentException("trying to read "+MeteoData::getParameterName(param)+" but found '"+type+"'", AT);
}

double OshdIO::convertUnits(const double& val, const std::string& units, const MeteoData::Parameters& param)
{
	if (units=="%") return val/100.;
	if (units=="cm") return val/100.;
	if (units=="mm") {
		if (param==MeteoData::PSUM) return val;
		else return val/1000.;
	}
	if (units=="\xB0\x43") return val+Cst::t_water_freezing_pt; //ISO-8859-1 hex for '°C'
// 	//usually skip these extra tests
// 	if (units=="\xB0") return val; //ISO-8859-1 hex for '°'
// 	if (units.empty()) return val;
// 	if (units=="Pa") return val;
// 	if (units=="W/m2") return val;
// 	if (units=="m/s") return val;
// 	else 
// 		throw IOException("Unknown units '"+units+"'", AT);
	
	return val;
}

void OshdIO::fillStationMeta()
{
	vecMeta.resize( vecIDs.size(), StationData() );
	mat_t *matfp = Mat_Open(in_metafile.c_str(), MAT_ACC_RDONLY);
	if ( NULL == matfp ) throw AccessException(in_metafile, AT);

	matvar_t *matvar = Mat_VarReadInfo(matfp, "statlist");
	if (matvar==NULL) throw NotFoundException("structure 'statlist' not found in file '"+in_metafile+"'", AT);
	
	std::vector<std::string> vecAcro;
	readStringVector(in_metafile, "acro", matfp, matvar, vecAcro);
	
	std::vector<std::string> vecNames;
	readStringVector(in_metafile, "name", matfp, matvar, vecNames);
	
	std::vector<double> easting, northing, altitude;
	readDoubleVector(in_metafile, "x", matfp, matvar, easting);
	readDoubleVector(in_metafile, "y", matfp, matvar, northing);
	readDoubleVector(in_metafile, "z", matfp, matvar, altitude);
	
	Mat_VarFree(matvar);
	Mat_Close(matfp);
	
	if (debug) {
		for (size_t ii=0; ii<vecAcro.size(); ii++) 
			std::cout << std::setw(8) << vecAcro[ii] << std::setw(40) << vecNames[ii] << std::setw(8) << easting[ii] << std::setw(8) << northing[ii] << std::setw(8) << altitude[ii] << "\n";
		std::cout << endl;
	}
	
	buildVecIdx(vecAcro);
	for (size_t ii=0; ii<vecIdx.size(); ii++) {
		const size_t idx = vecIdx[ii];
		Coords position("CH1903", "");
		position.setXY(easting[idx], northing[idx], altitude[idx]);
		std::string name( vecNames[idx] );

		//if the network name has been appended, remove it. We also remove spaces, just in case
		const size_t netz_pos = name.find(" (");
		if (netz_pos!=std::string::npos) name.erase(netz_pos);
		std::replace( name.begin(), name.end(), ' ', '_');

		const StationData sd(position, vecAcro[idx], name);
		vecMeta[ii] = sd;
	}
}

void OshdIO::buildVecIdx(const std::vector<std::string>& vecAcro)
{
	const size_t nrIDs = vecIDs.size();
	if (nrIDs==0)
		throw InvalidArgumentException("Please provide at least one station ID to read!", AT);
	vecIdx.resize( nrIDs, 0 );
	
	for (size_t ii=0; ii<nrIDs; ii++) {
		bool found = false;
		for (size_t jj=0; jj<vecAcro.size(); jj++) {
			if (vecIDs[ii]==vecAcro[jj]) {
				vecIdx[ii] = jj;
				found = true;
				break;
			}
		}
		if (!found) 
			throw NotFoundException("station ID '"+vecIDs[ii]+"' could not be found in the provided data", AT);
	}
}

} //namespace
