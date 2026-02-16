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
#ifndef CAAMLIO_H
#define CAAMLIO_H

#include <meteoio/MeteoIO.h>
#include <snowpack/Constants.h>
#include <snowpack/Hazard.h>
#include <snowpack/plugins/SnowpackIOInterface.h>
#include <snowpack/plugins/SmetIO.h>

#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include <snowpack/plugins/pugixml/pugixml.hpp>
#pragma GCC diagnostic pop


/**
 * @class CaaMLIO
 * @brief Reading snow profile data in CAAML format.
 * Reads in CAAML snow profile data, the CAA-IACS international standard
 *
 * @author Charles Fierz (Mathias Bavay)
 * @date   2014
 */
class CaaMLIO : public SnowpackIOInterface {

	public:
		CaaMLIO(const SnowpackConfig& i_cfg, const RunInfo& run_info);
		CaaMLIO(const CaaMLIO&);

		virtual bool snowCoverExists(const std::string& i_snowfile, const std::string& stationID) const;

		virtual void readSnowCover(const std::string& i_snowfile, const std::string& stationID,
		                           SN_SNOWSOIL_DATA& SSdata, ZwischenData& Zdata, const bool& read_salinity);

		virtual void writeSnowCover(const mio::Date& date, const SnowStation& Xdata,
		                            const ZwischenData& Zdata, const bool& forbackup=false);

		virtual void writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata,
		                             const ProcessDat& Hdata, const double wind_trans24);

		virtual void writeProfile(const mio::Date& date, const SnowStation& Xdata);

		virtual bool writeHazardData(const std::string& stationID, const std::vector<ProcessDat>& Hdata,
		                             const std::vector<ProcessInd>& Hdata_ind, const size_t& num);

	private:
		void cleanup() throw();
		void init(const SnowpackConfig& cfg);

		void setBasicHeader(const SnowStation& Xdata, const std::string& fields, smet::SMETWriter& smet_writer) const;
		void setSnoSmetHeader(const SnowStation& Xdata, const SN_SNOWSOIL_DATA& SSdata, const mio::Date& date, smet::SMETWriter& smet_writer) const;
		void setFormatting(const size_t& nr_solutes, std::vector<int>& vec_width, std::vector<int>&  vec_precision) const;
		std::string getFilenamePrefix(const std::string& fnam, const std::string& path, const bool addexp=true) const;

		//functions for reading caaml-files:
		bool read_snocaaml(const std::string& snofilename, const std::string& stationID, SN_SNOWSOIL_DATA& SSdata);
		void openIn_CAAML(const std::string& in_snowfile);
		mio::Date xmlGetDate();
		mio::StationData xmlGetStationData(const std::string& stationID);
		void setCustomSnowSoil(SN_SNOWSOIL_DATA& Xdata);
		void xmlReadLayerData(SN_SNOWSOIL_DATA& SSdata);
		void adjustToSlopeAngle(SN_SNOWSOIL_DATA& SSdata); // #cmb
		LayerData xmlGetLayer(pugi::xml_node nodeLayer, std::string& grainFormCode);
		bool getLayersDir();
		void getAndSetProfile(const std::string path, const std::string name,const bool directionTopDown,
		                      const bool isRangeMeasurement,std::vector<LayerData>& Layers);
		bool xmlGetProfile(const std::string path, const std::string name, std::vector<double>& zVec, std::vector<double>& valVec);
		void estimateValidFormationTimesIfNotSetYet(std::vector<LayerData> &Layers, const mio::Date);
		void checkAllDataForConsistencyAndSetMissingValues( SN_SNOWSOIL_DATA& SSdata );
		bool checkWhatWasReadIn(SN_SNOWSOIL_DATA& SSdata);

		//functions for writing caaml-file:
		void writeSnowFile(const std::string& snofilename, const mio::Date& date, const SnowStation& Xdata);
		void writeCustomSnowSoil(pugi::xml_node& node, const SnowStation& Xdata);
		void writeLayers(pugi::xml_node& node, const SnowStation& Xdata);
		void writeCustomLayerData(pugi::xml_node& node, const ElementData& Edata, const NodeData& Ndata);
		void writeProfiles(pugi::xml_node& node, const SnowStation& Xdata);
		void writeStationData(pugi::xml_node& root, const SnowStation& Xdata);

		double lwc_codeToVal(const char* code);
		std::string lwc_valToCode(const double val);
		double hardness_codeToVal(char* code);
		std::string hardness_valToCode(const double code);
		void grainShape_codeToVal(const std::string& code, double &sp, double &dd, unsigned short int &mk);
		std::string grainShape_valToAbbrev(const unsigned int var);
		std::string grainShape_valToAbbrev_old(const double* var);

		//xml functions:
		double xmlReadValueFromPath(const std::string& xpath, const std::string& property, const double& dflt);
		int xmlReadValueFromPath(const std::string& xpath, const std::string& property, const int& dflt);
		void xmlWriteElement(pugi::xml_node& node, const char* name, const char* content, const char* att_name, const char* att_val);
		bool xmlDoesPathExist(const std::string& path);
		bool xmlReadValueFromNode(const pugi::xml_node node, const std::string propertyName, double& variableToSet,
		                          const std::string unitOut = "",const std::string unitMeasured = "", const double factor=1.0);
		std::string xmlReadAttributeFromPath (const std::string& path, const std::string& attributeName);

		const RunInfo info;
		std::string i_snowpath, o_snowpath, experiment;
		double i_max_element_thickness;
		bool caaml_writeout_as_readin, haz_write;
		/*static const*/ double in_tz; //plugin specific time zones

		pugi::xml_document inDoc;
		pugi::xml_encoding inEncoding;

		//charEncoding in_encoding;
		static const char *xml_ns_caaml, *xml_ns_abrev_caaml;
		static const char *xml_ns_gml, *xml_ns_abrev_gml;
		static const char *xml_ns_xsi, *xml_ns_abrev_xsi;
		static const char *xml_ns_slf, *xml_ns_abrev_slf;
		static const char *xml_ns_snp, *xml_ns_abrev_snp;
		static const std::string TimeData_xpath, StationMetaData_xpath, SnowData_xpath;

		char layerDepthTopStr[10], layerThicknessStr[10], layerValStr[10], valueStr[10];
		double hoarDensitySurf;
		std::vector<std::string> grainForms;
};

#endif //End of CAAMLIO.h
