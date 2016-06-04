/***********************************************************************************/
/*  Copyright 2009 EPFL                                                            */
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
#ifndef __GEOTOPIO_H__
#define __GEOTOPIO_H__

#include <meteoio/Config.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/IOExceptions.h>

#include <string>
#include <sstream>
#include <iostream>

namespace mio {

/**
 * @class GeotopIO
 * @brief This class enables the access meteo data in legacy Geotop format
 *
 * @ingroup plugins
 * @author Thomas Egger
 * @date   2009-07-02
 */
class GeotopIO : public IOInterface {
	public:
		GeotopIO(const std::string& configfile);
		GeotopIO(const GeotopIO&);
		GeotopIO(const Config&);
		~GeotopIO() throw();

		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);

		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo,
		                           const size_t& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string& name="");

		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readPOI(std::vector<Coords>& pts);

		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	private:
		std::string getValueForKey(const std::string& line);
		void initParamNames(std::map<std::string, size_t>& mapParam);
		void readMetaData(const std::string& metafile);
		void identify_fields(const std::vector<std::string>& tmpvec, const std::string& filename,
		                     std::vector<size_t>& indices, MeteoData& md);
		void convertUnits(MeteoData& meteo);
		void convertUnitsBack(MeteoData& meteo);
		void cleanup() throw();
		void parseDate(const std::string& datestring, const std::string& fileandline, Date& date);
		void parseMetaData(const std::string& head, const std::string& datastr, std::vector<std::string>& tmpvec);

		const Config cfg;
		double in_tz, out_tz;
		size_t nr_of_stations;
		std::ifstream fin; //Input file streams
		std::ofstream fout; //Output file streams
		std::vector< std::map <Date, std::streampos> > vec_streampos; //in order to save file pointers
		std::vector<mio::StationData> vecStation;
		std::map<std::string, size_t> mapColumnNames;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		static const size_t sw_direct, sw_diffuse, cloud_factor;
};

} //end namespace

#endif
