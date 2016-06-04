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
#ifndef __A3DIO_H__
#define __A3DIO_H__

#include <meteoio/IOInterface.h>

#include <meteoio/dataClasses/Coords.h>
#include <meteoio/Config.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>

#include <string>
#include <vector>
#include <map>

#ifdef _MSC_VER
	#pragma warning(disable:4512) //we don't need any = operator!
#endif

namespace mio {

class A3DIO : public IOInterface {
	public:
		A3DIO(const std::string& configfile);
		A3DIO(const Config&);

		virtual void read2DGrid(Grid2DObject& dem_out, const std::string& name="");
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

		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	private:
		void read1DStation(StationData& sd);
		void read1DMeteo(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >&);
		void read2DStations(const Date& timestamp, std::vector<StationData>& vecStation);
		void read2DMeteo(std::vector< std::vector<MeteoData> >&);

		void constructMeteo2DFilenames(const Date& i_startDate, const Date& i_endDate, std::vector<std::string>& i_filenames);
		bool readMeteoDataLine(std::string& line, MeteoData& tmpdata, std::string filename);
		void convertUnits(MeteoData& meteo);
		void read2DMeteoData(const std::string&, const std::string&, std::map<std::string,size_t>& hashStations,
		                     std::vector< std::vector<MeteoData> >&, size_t& bufferindex);
		void read2DMeteoHeader(const std::string& filename, std::map<std::string, size_t>& hashStations,
		                       std::vector<StationData>&);
		size_t getNrOfStations(std::vector<std::string>& filenames,
		                       std::map<std::string, size_t>& hashStations);

		int create1DFile(const std::vector< std::vector<MeteoData> >& data);
		int writeHeader(std::ofstream &file, const std::vector< std::vector<MeteoData> >& stations, const std::string& parameter_name);
		void open2DFile(const std::vector< std::vector<MeteoData> >& stations,
		                const std::string& fileprefix, const std::string& label, const double& year,
		                std::ofstream& file);
		int write2DmeteoFile(const std::vector< std::vector<MeteoData> >& data, const unsigned int& parindex,
		                     const std::string& filename, const std::string& label);
		void write2DMeteo(const std::vector< std::vector<MeteoData> >& data);

		const Config cfg;
		double in_tz, out_tz; //timezones
		std::string meteo1d;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters

		static const double plugin_nodata; //plugin specific nodata value, e.g. -9999
};
} //end namespace

#endif
