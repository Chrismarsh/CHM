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
#ifndef __CNRMIO_H__
#define __CNRMIO_H__

#include <meteoio/IOInterface.h>
#include <meteoio/Config.h>

#include <string>

namespace mio {

/**
 * @class CNRMIO
 * @brief This plug-in allows reading and writing of NetCDF files formatted according to CNRM standard.
 *
 * @ingroup plugins
 * @author Thomas Egger
 * @date   2014-03-13
 */
class CNRMIO : public IOInterface {
	public:
		CNRMIO(const std::string& configfile);
		CNRMIO(const CNRMIO&);
		CNRMIO(const Config& cfgreader);
		~CNRMIO() throw();

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
		enum TimeUnit { seconds, hours, days };
		enum Naming { cf, cnrm, ecmwf };

		void parseInputOutputSection();
		void create_parameters(const int& ncid, const int& did_time, const int& did_points, const size_t& number_of_records,
		                       const size_t& number_of_stations, std::map<size_t, std::string>& map_param_name,
		                       std::map<std::string, double*>& map_data_2D, std::map<std::string, int>& varid);
		void create_meta_data(const int& ncid, const int& did, std::map<std::string, double*>& map_data_1D, std::map<std::string, int>& varid);
		void get_parameters(const std::vector< std::vector<MeteoData> >& vecMeteo, std::map<size_t, std::string>& map_param_name,
		                    std::map<std::string, double*>& map_data_1D, double*& dates);
		void get_parameters(const int& ncid, std::map<std::string, size_t>& map_parameters, MeteoData& meteo_data);
		size_t get_dates(const std::vector< std::vector<MeteoData> >& vecMeteo, double*& dates);
		void copy_data(const size_t& number_of_stations, const size_t& number_of_records, const std::vector< std::vector<MeteoData> >& vecMeteo,
                         const std::map<size_t, std::string>& map_param_name, std::map<std::string, double*>& map_data_2D);
		void copy_data(const int& ncid, const std::map<std::string, size_t>& map_parameters, const std::map<std::string, double*> map_data,
		               const size_t& number_of_stations, const size_t& number_of_records, std::vector< std::vector<MeteoData> >& vecMeteo);
		void readData(const int& ncid, const size_t& index_start, const std::vector<Date>& vec_date, const std::map<std::string, size_t>& map_parameters,
		              const MeteoData& meteo_data, std::vector< std::vector<MeteoData> >& vecMeteo);
		void readMetaData(const int& ncid, std::vector<StationData>& vecStation);
		void get_meta_data_ids(const int& ncid, std::map<std::string, int>& map_vid);
		std::string get_varname(const MeteoGrids::Parameters& parameter);
		void get_indices(const int& ncid, const Date& dateStart, const Date& dateEnd, size_t& indexStart, size_t& indexEnd, std::vector<Date>& vecDate);
		void calculate_offset(const std::string& units, CNRMIO::TimeUnit& time_unit, Date& offset);
		void check_consistency(const int& ncid, const Grid2DObject& grid, double*& lat_array, double*& lon_array,
		                       int& did_lat, int& did_lon, int& vid_lat, int& vid_lon);
		void read2DGrid_internal(Grid2DObject& grid_out, const std::string& full_name, const std::string& varname, const Date& date=Date());
		void write2DGrid_internal(const Grid2DObject& grid_in, const std::string& filename, const std::string& varname, const Date& date=Date());
		void add_attributes_for_variable(const int& ncid, const int& varid, const std::string& varname);
		void create_latlon_dimensions(const int& ncid, const Grid2DObject& grid, int& did_lat, int& did_lon, int& vid_lat, int& vid_lon);
		void create_time_dimension(const int& ncid, int& did_time, int& vid_time);

		// Private variables
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		static const double epsilon; //for numerical comparisons of double values
		static const std::string cf_time, cf_units, cf_days, cf_hours, cf_seconds, cf_latitude, cf_longitude, cf_altitude, cf_ta, cf_rh, cf_p;
		static const std::string cnrm_points, cnrm_latitude, cnrm_longitude, cnrm_altitude, cnrm_aspect, cnrm_slope, cnrm_ta, cnrm_rh, cnrm_vw, cnrm_dw, cnrm_qair;
		static const std::string cnrm_co2air, cnrm_theorsw, cnrm_neb, cnrm_psum, cnrm_snowf, cnrm_swr_direct, cnrm_swr_diffuse, cnrm_p, cnrm_ilwr, cnrm_timestep;
		static const std::string ecmwf_ta, ecmwf_p, ecmwf_iswr, ecmwf_ilwr, ecmwf_psum, ecmwf_u10, ecmwf_v10, ecmwf_td;

		static std::map<std::string, size_t> paramname; ///<Associate a name with meteo parameters in Parameters
		static std::map<std::string, std::string> map_name; ///Associate MeteoIO parameter names with CNRM parameter names
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map

		const Config cfg;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		double in_dflt_TZ, out_dflt_TZ;     //default time zones
		bool in_strict, out_strict;
		std::vector<StationData> vecMetaData;
};

} //namespace
#endif
