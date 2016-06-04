/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __NetCDFIO_H__
#define __NetCDFIO_H__

#include <meteoio/IOInterface.h>
#include <meteoio/Config.h>

#include <string>

namespace mio {

/**
 * @class NetCDFIO
 * @brief This plug-in allows reading and writing of NetCDF files formatted according to CNRM standard.
 *
 * @ingroup plugins
 * @author Thomas Egger
 * @date   2014-03-13
 */
class NetCDFIO : public IOInterface {
	public:
		NetCDFIO(const std::string& configfile);
		NetCDFIO(const NetCDFIO&);
		NetCDFIO(const Config& cfgreader);
		~NetCDFIO() throw();

		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo,
		                           const size_t& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
		                            const std::string& name="");

		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readPOI(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	private:
		typedef struct ATTRIBUTES {
			ATTRIBUTES() : var(), standard_name(), long_name(), units(), height(IOUtils::nodata) {};
			ATTRIBUTES(const std::string& str1, const std::string& str2, const std::string& str3, const std::string& str4, const double& hgt)
			                     : var(str1), standard_name(str2), long_name(str3), units(str4), height(hgt) {};
  
			std::string var;
			std::string standard_name;
			std::string long_name;
			std::string units;
			double height;
		} attributes;

		void initAttributesMap(const std::string& schema, std::map<MeteoGrids::Parameters, attributes> &attr);
		void setTimeTransform(const std::string& schema, double &time_offset, double &time_multiplier);
		void parseInputOutputSection();
		void check_consistency(const int& ncid, const Grid2DObject& grid, double*& lat_array, double*& lon_array,
		                       int& did_lat, int& did_lon, int& vid_lat, int& vid_lon);
		//void read2DGrid_meta(const std::string& filename, const std::string& varname, int &ncid, int &varid, std::vector<size_t> &dimlen, size_t &time_index, size_t &lat_index, size_t &lon_index, double &missing_value,double &lat, double &lon);
		bool read2DGrid_internal(Grid2DObject& grid_out, const std::string& filename, const MeteoGrids::Parameters& parameter, const Date& date=Date());
		bool read2DGrid_internal(Grid2DObject& grid_out, const std::string& full_name, const std::string& varname, const Date& date=Date(), const bool& isPrecip=false);
		void write2DGrid_internal(Grid2DObject grid_in, const std::string& filename, const attributes& attr, const Date& date=Date(), const bool& isPrecip=false);
		void add_attributes_for_variable(const int& ncid, const int& varid, const attributes& attr, const double& nodata_out);
		void getTimeTransform(const int& ncid, const int& varid, double &time_offset, double &time_multiplier) const;
		void create_latlon_dimensions(const int& ncid, const Grid2DObject& grid, int& did_lat, int& did_lon, int& vid_lat, int& vid_lon);
		void create_time_dimension(const int& ncid, int& did_time, int& vid_time);
		void readWind(const std::string& filename, const Date& date);

		// Private variables
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		static const double epsilon; //for numerical comparisons of double values
		static const std::string cf_time, cf_latitude, cf_longitude, cf_altitude;

		const Config cfg;
		std::map <MeteoGrids::Parameters, attributes> in_attributes, out_attributes;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		double in_dflt_TZ, out_dflt_TZ; //default time zones
		double in_time_offset, in_time_multiplier; //each schema defines its own time specification...
		bool dem_altimeter, in_strict, out_strict;
		std::vector<StationData> vecMetaData;
};

} //namespace
#endif
