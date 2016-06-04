/***********************************************************************************/
/*  Copyright 2009, 2010 WSL Institute for Snow and Avalanche Research   SLF-DAVOS */
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
#include <meteoio/plugins/ImisIO.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/MathOptim.h>

using namespace std;
using namespace oracle;
using namespace oracle::occi;
using namespace mio;

namespace mio {
/**
 * @page imis IMIS
 * @section imis_format Format
 * This plugin reads data directly from the IMIS network database (Oracle database).
 * It retrieves standard IMIS data as well as ENETZ and ANETZ data.
 *
 * @section imis_units Units
 * The units are assumed to be the following:
 * - temperatures in celsius
 * - relative humidity in %
 * - wind speed in m/s
 * - precipitations in mm/h
 * - radiation in W/m²
 *
 * @section imis_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: input coordinate system (see Coords) specified in the [Input] section
 * - COORDPARAM: extra input coordinates parameters (see Coords) specified in the [Input] section
 * - COORDSYS: output coordinate system (see Coords) specified in the [Output] section
 * - COORDPARAM: extra output coordinates parameters (see Coords) specified in the [Output] section
 * - DBNAME: name of the database to use (exemple: sdbo)
 * - DBUSER: user name to use when connecting to the database
 * - DBPASS: password to use when connecting to the database
 * - STATION#: station code for the given number #
 * - USEANETZ: use ANETZ stations to provide precipitations for normal IMIS stations. Almost each IMIS station is associated with one or two ANETZ stations and does a weighted average to get what should be its local precipitations if no local precipitation has been found (either nodata or 0).
 * - USE_IMIS_PSUM: if set to false (default), all IMIS precipitation will be deleted (since IMIS stations don't have heated rain gauges, their precipitation measurements are not good in winter conditions). If set to true, the precipitation measurements will be accepted from IMIS stations. In this case, it is strongly advised to apply the filter FilterUnheatedPSUM to detect snow melting in the rain gauge.
 * - USE_SNOWPACK_PSUM: if set to true, the SNOWPACK simulated Snow Water Equivalent from the database will be used to compute PSUM. Data gaps greater than 3 hours on SWE will lead to unchanged psum while all data that can properly be computed will <b>overwrite</b> psum. (default=false)
 *
 * It is possible to use both USE_IMIS_PSUM and USE_SNOWPACK_PSUM to create composite PSUM (from SNOWPACK in the snow season and from IMIS otherwise).
 * In such a case, as soon as SNOWPACK SWE > 0, all previous PSUM data will be deleted (ie those potentially coming from IMIS_PSUM).
 * But if there is no SNOWPACK data, the IMIS measurements will be kept.
 */

const double ImisIO::plugin_nodata = -999.; ///< plugin specific nodata value
const double ImisIO::in_tz = 1.; ///< All IMIS data is in gmt+1, that is UTC+1 (a quelques secondes près;-)

const string ImisIO::sqlQueryStationIDs = "SELECT station_name, drift_stat_abk, drift_stao_nr FROM station2.v_snow_drift_standort WHERE application_code='snowpack' AND station_code=:1"; ///< Wind drift station meta data

const string ImisIO::sqlQueryStationMetaData = "SELECT stao_name, stao_x, stao_y, stao_h FROM station2.standort WHERE stat_abk LIKE :1 AND stao_nr=:2"; ///< Snow station meta data

const string ImisIO::sqlQuerySensorDepths = "SELECT hts1_1, hts1_2, hts1_3 FROM station2.standort WHERE stat_abk LIKE :1 AND stao_nr=:2"; ///< Sensor depths at station

const string ImisIO::sqlQueryMeteoDataDrift = "SELECT TO_CHAR(a.datum, 'YYYY-MM-DD HH24:MI') AS thedate, a.ta, a.iswr, a.vw, a.dw, a.vw_max, a.rh, a.ilwr, a.hnw, a.tsg, a.tss, a.hs, a.rswr, b.vw AS vw_drift, b.dw AS dw_drift, a.ts1, a.ts2, a.ts3 FROM (SELECT * FROM ams.v_ams_raw WHERE stat_abk=:1 AND stao_nr=:2 AND datum>=:3 AND datum<=:4) a LEFT OUTER JOIN (SELECT case when to_char(datum,'MI')=40 then trunc(datum,'HH24')+0.5/24 else datum end as datum, vw, dw FROM ams.v_ams_raw WHERE stat_abk=:5 AND stao_nr=:6 AND datum>=:3 AND datum<=:4) b ON a.datum=b.datum ORDER BY thedate"; ///< C. Marty's Data query with wind drift station; gets wind from enet stations for imis snow station too! [2010-02-24]

const string ImisIO::sqlQueryMeteoData = "SELECT TO_CHAR(datum, 'YYYY-MM-DD HH24:MI') AS thedate, ta, iswr, vw, dw, vw_max, rh, ilwr, hnw, tsg, tss, hs, rswr, ts1, ts2, ts3 FROM ams.v_ams_raw WHERE stat_abk=:1 AND stao_nr=:2 AND datum>=:3 AND datum<=:4 ORDER BY thedate ASC"; ///< Data query without wind drift station

const string ImisIO::sqlQuerySWEData = "SELECT TO_CHAR(datum, 'YYYY-MM-DD HH24:MI') AS thedate, swe FROM snowpack.ams_pmod WHERE stat_abk=:1 AND stao_nr=:2 AND datum>=:3 AND datum<=:4 ORDER BY thedate ASC"; ///< Query SWE as calculated by SNOWPACK to feed into PSUM

std::map<std::string, AnetzData> ImisIO::mapAnetz;
const bool ImisIO::__init = ImisIO::initStaticData();

bool ImisIO::initStaticData()
{
	//Associate string with AnetzData
	//map[station ID] = (#stations, STA1, STA2, STA3, #coeffs, coeff1, coeff2, coeff3)
	mapAnetz["AMD2"] = AnetzData(2,"*GLA","*SAE","",3,1.2417929,0.548411708,-0.0692799);
	mapAnetz["ANV2"] = AnetzData(2,"*EVO","*MVE","",2,0.7920454,0.771111962,IOUtils::nodata);
	mapAnetz["ANV3"] = AnetzData(1,"*EVO","","",1,1.6468,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["ARO2"] = AnetzData(2,"*EVO","*GSB","",2,0.9692294,0.218384531,IOUtils::nodata);
	mapAnetz["ARO3"] = AnetzData(2,"*EVO","*ZER","",3,1.0748285,1.649860092,-0.0728015);
	mapAnetz["BED2"] = AnetzData(2,"*PIO","*ULR","",3,0.9934869,1.047586006,-0.05489259);
	mapAnetz["BED3"] = AnetzData(2,"*PIO","*ULR","",2,0.6999,0.4122,IOUtils::nodata);
	mapAnetz["BER2"] = AnetzData(2,"*ROB","*COV","",3,1.4454061,0.558775717,-0.05063568);
	mapAnetz["BER3"] = AnetzData(2,"*ROB","*COV","",2,0.378476,0.817976734,IOUtils::nodata);
	mapAnetz["BEV2"] = AnetzData(2,"*SAM","*COV","",3,1.8237643,0.853292298,-0.33642156);
	mapAnetz["BOG2"] = AnetzData(1,"*ROE","","",1,1.0795,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["BOR2"] = AnetzData(1,"*VIS","","",1,1.0662264,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["BOV2"] = AnetzData(2,"*GSB","*EVO","",2,0.3609309,0.934922978,IOUtils::nodata);
	mapAnetz["CAM2"] = AnetzData(2,"*PIO","*COM","",2,0.750536,0.426864157,IOUtils::nodata);
	mapAnetz["CHA2"] = AnetzData(2,"*AIG","*SIO","",2,0.7107216,0.99869915,IOUtils::nodata);
	mapAnetz["CON2"] = AnetzData(2,"*SIO","*MVE","",3,3.5344378,1.952708399,-0.74509918);
	mapAnetz["DAV2"] = AnetzData(2,"*WFJ","*DAV","",3,0.594108,1.091565634,-0.12150025);
	mapAnetz["DAV3"] = AnetzData(2,"*WFJ","*DAV","",3,0.9266618,0.815816241,-0.06248703);
	mapAnetz["DAV4"] = AnetzData(2,"*WFJ","*DAV","",3,0.9266618,0.815816241,-0.06248703);
	mapAnetz["DAV5"] = AnetzData(2,"*WFJ","*DAV","",3,0.9266618,0.815816241,-0.06248703);
	mapAnetz["DTR2"] = AnetzData(2,"*PIO","*COM","",2,0.0384,0.9731,IOUtils::nodata);
	mapAnetz["DVF2"] = AnetzData(1,"*WFJ","","",1,1,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["ELM2"] = AnetzData(1,"*GLA","","",1,1.4798048,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["ELS2"] = AnetzData(2,"*ABO","*INT","",3,1.0886792,0.568730457,-0.07758286);
	mapAnetz["FAE2"] = AnetzData(1,"*ABO","","",1,2.1132038,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["FIR2"] = AnetzData(2,"*INT","*GRH","",3,1.2416838,0.243226327,-0.02392287);
	mapAnetz["FIS2"] = AnetzData(1,"*ABO","","",1,1.1991,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["FNH2"] = AnetzData(2,"*AIG","*GSB","",2,1.3949428,0.297933922,IOUtils::nodata);
	mapAnetz["FOU2"] = AnetzData(1,"*GSB","","",1,0.8448844,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["FUL2"] = AnetzData(2,"*FEY","*AIG","",2,1.070156,0.587972864,IOUtils::nodata);
	mapAnetz["FUS2"] = AnetzData(1,"*PIO","","",1,1.3557753,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["GAD2"] = AnetzData(2,"*ENG","*GUE","",3,0.9764334,0.814293499,-0.07074082);
	mapAnetz["GAN2"] = AnetzData(2,"*ABO","*VIS","",2,0.520224,0.825813298,IOUtils::nodata);
	mapAnetz["GLA2"] = AnetzData(1,"*GLA","","",1,1.7186314,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["GOM2"] = AnetzData(2,"*ULR","*GRH","",2,0.4413,0.4235,IOUtils::nodata);
	mapAnetz["GOM3"] = AnetzData(2,"*ULR","*GRH","",2,0.3269755,0.62995601,IOUtils::nodata);
	mapAnetz["GUT2"] = AnetzData(2,"*GRH","*ENG","",2,0.3977985,0.463100458,IOUtils::nodata);
	mapAnetz["GUT3"] = AnetzData(2,"*GRH","*ENG","",2,0.3977985,0.463100458,IOUtils::nodata);
	mapAnetz["HTR2"] = AnetzData(2,"*HIR","*COM","",2,0.8668,0.5939,IOUtils::nodata);
	mapAnetz["HTR3"] = AnetzData(2,"*SBE","*COM","",2,1.3023275,-0.663411226,IOUtils::nodata);
	mapAnetz["ILI2"] = AnetzData(1,"*AIG","","",1,1.2341516,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["JUL2"] = AnetzData(2,"*COV","*SAM","",2,0.4900961,0.871078269,IOUtils::nodata);
	mapAnetz["KES2"] = AnetzData(2,"*SAM","*DAV","",2,0.847596,1.112635571,IOUtils::nodata);
	mapAnetz["KLO2"] = AnetzData(1,"*DAV","","",1,1.585,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["KLO3"] = AnetzData(2,"*DAV","*WFJ","",3,0.8352,0.9493,-0.0526);
	mapAnetz["LAU2"] = AnetzData(2,"*ABO","*SIO","",2,0.3037172,0.791695555,IOUtils::nodata);
	mapAnetz["LUK2"] = AnetzData(2,"*DIS","*PIO","",3,0.8593029,0.378261758,0.85930291);
	mapAnetz["MEI2"] = AnetzData(3,"*ENG","*GUE","*ALT",3,0.3882119,0.399244859,0.3298324);
	mapAnetz["MES2"] = AnetzData(2,"*HIR","*COM","",2,1.3552818,-0.393843912,IOUtils::nodata);
	mapAnetz["MUN2"] = AnetzData(1,"*VIS","","",1,0.8624804,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["NAR2"] = AnetzData(2,"*PIO","*COM","",3,0.4089981,0.873419792,-0.028464);
	mapAnetz["NEN2"] = AnetzData(2,"*SIO","*EVO","",3,0.9352699,1.312867984,-0.14543389);
	mapAnetz["OBM2"] = AnetzData(2,"*AIG","*MLS","",3,1.9413387,1.64250639,-0.37210579);
	mapAnetz["OBW2"] = AnetzData(2,"*GRH","*ULR","",3,0.2471352,1.219258485,-0.02153657);
	mapAnetz["OBW3"] = AnetzData(2,"*GRH","*ULR","",2,0.5274,0.4815,IOUtils::nodata);
	mapAnetz["OFE2"] = AnetzData(1,"*SCU","","",1,1.8758744,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["ORT2"] = AnetzData(1,"*GLA","","",1,1.6214,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["OTT2"] = AnetzData(1,"*ABO","","",1,1.3759903,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["PAR2"] = AnetzData(1,"*WFJ","","",1,1.6252986,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["PUZ2"] = AnetzData(2,"*DIS","*GUE","",2,0.9481811,0.1490937,IOUtils::nodata);
	mapAnetz["ROA2"] = AnetzData(2,"*INT","*NAP","",3,1.748338,0.574491521,-0.1670437);
	mapAnetz["SAA2"] = AnetzData(2,"*ZER","*VIS","",3,0.6316695,1.210149675,-0.11760175);
	mapAnetz["SAA3"] = AnetzData(1,"*VIS","","",1,1.2905,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["SCA2"] = AnetzData(2,"*ALT","*DIS","",2,0.8118627,0.360141586,IOUtils::nodata);
	mapAnetz["SCA3"] = AnetzData(2,"*ALT","*GLA","",2,0.4768725,0.819642544,IOUtils::nodata);
	mapAnetz["SCB2"] = AnetzData(2,"*ENG","*INT","",3,1.0535332,1.21234263,-0.1307221);
	mapAnetz["SCH2"] = AnetzData(1,"*INT","","",1,1.54557,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["SHE2"] = AnetzData(1,"*INT","","",1,1.1065938,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["SIM2"] = AnetzData(2,"*COM","*SBE","",2,0.6861131,0.296215066,IOUtils::nodata);
	mapAnetz["SLF2"] = AnetzData(1,"*WFJ","","",1,0.9585787,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["SMN2"] = AnetzData(1,"*SCU","","",1,0.6979953,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["SPN2"] = AnetzData(2,"*VIS","*ZER","",2,1.1049,1.4598,IOUtils::nodata);
	mapAnetz["SPN3"] = AnetzData(1,"*VIS","","",1,1.0244902,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["STH2"] = AnetzData(2,"*PLF","*ABO","",3,1.1252659,0.893324895,-0.13194965);
	mapAnetz["STN2"] = AnetzData(2,"*EVO","*MVE","",2,0.9042348,0.687519213,IOUtils::nodata);
	mapAnetz["TAM2"] = AnetzData(2,"*VAD","*GLA","",2,0.6304286,0.738150034,IOUtils::nodata);
	mapAnetz["TAM3"] = AnetzData(2,"*VAD","*GLA","",3,1.5515584,0.407868299,-0.0800763);
	mapAnetz["TRU2"] = AnetzData(2,"*MVE","*VIS","",2,1.1359,0.6577,IOUtils::nodata);
	mapAnetz["TUJ2"] = AnetzData(2,"*GUE","*DIS","",2,0.3636322,0.591777057,IOUtils::nodata);
	mapAnetz["TUJ3"] = AnetzData(2,"*GUE","*DIS","",2,0.4742,0.7791,IOUtils::nodata);
	mapAnetz["TUM2"] = AnetzData(1,"*DIS","","",1,1.752091,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["URS2"] = AnetzData(2,"*GUE","*GRH","",3,0.6847615,0.277707092,-0.03085219);
	mapAnetz["VAL2"] = AnetzData(2,"*PIO","*GUE","",3,1.2130704,0.508735389,-0.02905053);
	mapAnetz["VDS2"] = AnetzData(1,"*MVE","","",1,1.8282525,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["VIN2"] = AnetzData(1,"*SCU","","",1,0.8245,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["VLS2"] = AnetzData(2,"*DIS","*HIR","",2,0.5764952,0.613916765,IOUtils::nodata);
	mapAnetz["WFJ2"] = AnetzData(1,"*WFJ","","",1,1.,IOUtils::nodata,IOUtils::nodata);
	mapAnetz["ZER2"] = AnetzData(2,"*ZER","*EVO","",2,0.8707182,0.988158355,IOUtils::nodata);
	mapAnetz["ZER4"] = AnetzData(2,"*ZER","*EVO","",2,0.8707182,0.988158355,IOUtils::nodata);
	mapAnetz["ZNZ2"] = AnetzData(1,"*WFJ","","",1,0.9980525,IOUtils::nodata,IOUtils::nodata);

	return true;
}

void ImisIO::getDBParameters()
{
	cfg.getValue("DBNAME", "Input", oracleDBName_in);
	cfg.getValue("DBUSER", "Input", oracleUserName_in);
	cfg.getValue("DBPASS", "Input", oraclePassword_in);

	cfg.getValue("USEANETZ", "Input", useAnetz, IOUtils::nothrow);
	cfg.getValue("USE_IMIS_PSUM", "Input", use_imis_psum, IOUtils::nothrow);
	cfg.getValue("USE_SNOWPACK_PSUM", "Input", use_psum_snowpack, IOUtils::nothrow);
}

ImisIO::ImisIO(const std::string& configfile)
        : cfg(configfile), coordin(), coordinparam(), coordout(), coordoutparam(), vecStationMetaData(), mapDriftStation(),
          oracleUserName_in(), oraclePassword_in(), oracleDBName_in(), useAnetz(false), use_imis_psum(false), use_psum_snowpack(false)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getDBParameters();
}

ImisIO::ImisIO(const Config& cfgreader)
        : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam(), vecStationMetaData(), mapDriftStation(),
          oracleUserName_in(), oraclePassword_in(), oracleDBName_in(), useAnetz(false), use_imis_psum(false), use_psum_snowpack(false)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getDBParameters();
}

ImisIO::~ImisIO() throw()
{
	cleanup();
}

void ImisIO::read2DGrid(Grid2DObject&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::read2DGrid(Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readDEM(DEMObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readLanduse(Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readAssimilationData(const Date&, Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readPOI(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::write2DGrid(const Grid2DObject&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::write2DGrid(const Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::writeMeteoData(const std::vector< std::vector<MeteoData> >&,
                            const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::openDBConnection(oracle::occi::Environment*& env, oracle::occi::Connection*& conn)
{
	env  = Environment::createEnvironment();// static OCCI function
	conn = env->createConnection(oracleUserName_in, oraclePassword_in, oracleDBName_in);
}

void ImisIO::closeDBConnection(oracle::occi::Environment*& env, oracle::occi::Connection*& conn)
{
	try {
		if (conn != NULL)
			env->terminateConnection(conn);
		Environment::terminateEnvironment(env); // static OCCI function
	} catch (const exception&){
		Environment::terminateEnvironment(env); // static OCCI function
	}
}

void ImisIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{
	vecStation.clear();

	if (vecStationMetaData.empty()){//Imis station meta data cannot change between time steps
		Environment *env = NULL;
		Connection *conn = NULL;

		try {
			openDBConnection(env, conn);

			readStationMetaData(conn); //reads all the station meta data into the vecStationMetaData (member vector)

			closeDBConnection(env, conn);
		} catch (const exception& e){
			closeDBConnection(env, conn);
			throw IOException("Oracle Error when reading stations' metedata: " + string(e.what()), AT); //Translation of OCCI exception to IOException
		}
	}

	vecStation = vecStationMetaData; //vecStationMetaData is a global vector holding all meta data
}

/**
 * @brief A meta function that extracts all station names from the Config,
 *        parses them and retrieves all meta data from SDB
 */
void ImisIO::readStationMetaData(oracle::occi::Connection*& conn)
{
	vector<string> vecStationID;
	readStationIDs(vecStationID);


	Statement *stmt = conn->createStatement();
	for (size_t ii=0; ii<vecStationID.size(); ii++) {
		// Retrieve the station IDs - this only needs to be done once per instance
		string stat_abk, stao_nr, station_name;
		parseStationID(vecStationID[ii], stat_abk, stao_nr);
		vector<string> stnIDs;
		string drift_stat_abk, drift_stao_nr;
		getStationIDs(vecStationID[ii], sqlQueryStationIDs, stnIDs, stmt);
		IOUtils::convertString(station_name, stnIDs.at(0));
		IOUtils::convertString(drift_stat_abk, stnIDs.at(1));
		IOUtils::convertString(drift_stao_nr, stnIDs.at(2));
		const string drift_stationID = drift_stat_abk + drift_stao_nr;
		if (!drift_stationID.empty()) {
			mapDriftStation[vecStationID[ii]] = drift_stationID;
		} else {
			throw ConversionFailedException("Error! No drift station for station "+stat_abk+stao_nr, AT);
		}

		// Retrieve the station meta data - this only needs to be done once per instance
		vector<string> stationMetaData;
		string stao_name;
		getStationMetaData(stat_abk, stao_nr, sqlQueryStationMetaData, stationMetaData, stmt);
		double east, north, alt;
		IOUtils::convertString(stao_name, stationMetaData.at(0));
		IOUtils::convertString(east, stationMetaData.at(1), std::dec);
		IOUtils::convertString(north, stationMetaData.at(2), std::dec);
		IOUtils::convertString(alt, stationMetaData.at(3), std::dec);

		//obtain a valid station_name w/o spaces within
		if (station_name.empty()) {
			if (!stao_name.empty()) {
				station_name += vecStationID[ii] + ":" + stao_name;
			} else {
				station_name += vecStationID[ii];
			}
		} else {
			vector<string> tmpname;
			IOUtils::readLineToVec(station_name, tmpname, ' ');
			size_t jj=1;
			while (jj < tmpname.size()) {
				if (tmpname.at(jj) != "-") {
					tmpname.at(0) += "_" + tmpname.at(jj);
				} else {
					tmpname.at(0) += ":";
					if (jj < tmpname.size()-1)
						tmpname.at(0) += tmpname.at(++jj);
				}
				jj++;
			}
			station_name = tmpname.at(0);
		}
		Coords myCoord(coordin, coordinparam);
		myCoord.setXY(east, north, alt);
		vecStationMetaData.push_back(StationData(myCoord, vecStationID[ii], station_name));
	}
	conn->terminateStatement(stmt);
}

/**
 * @brief This function breaks up the station name into two components (a string and a number e.g. KLO2 -> "KLO","2")
 * @param stationID The full name of the station (e.g. "KLO2")
 * @param stName      The string part of the name  (e.g. "KLO")
 * @param stNumber    The integer part of the name (e.g. "2")
 */
void ImisIO::parseStationID(const std::string& stationID, std::string& stat_abk, std::string& stao_nr)
{
	stat_abk = stationID.substr(0, stationID.length()-1); //The station name: e.g. KLO
	stao_nr = stationID.substr(stationID.length()-1, 1); //The station number: e.g. 2
	if(!std::isdigit(stao_nr[0])) {
		//the station is one of these non-imis stations that don't contain a number...
		stat_abk = stationID;
		stao_nr = "0";
	}
}

/**
 * @brief This function extracts all info about the stations that are to be used from global Config object
 * @param vecStationID A vector that will hold all relevant stations as std::strings
 */
void ImisIO::readStationIDs(std::vector<std::string>& vecStationID)
{
	vecStationID.clear();
	cfg.getValues("STATION", "INPUT", vecStationID);

	if(vecStationID.empty()) {
		cerr << "\tNo stations specified for IMISIO... is this what you want?\n";
	}
}

void ImisIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                           std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& stationindex)
{
	Environment *env = NULL;
	Connection *conn = NULL;
	Statement *stmt = NULL;

	try {
		if (vecStationMetaData.empty()) {
			openDBConnection(env, conn);
			readStationMetaData(conn); //reads all the station meta data into the vecStationMetaData (member vector)
		}

		if (vecStationMetaData.empty()) { //if there are no stations -> return
			if ((env != NULL) || (conn != NULL)) closeDBConnection(env, conn);
			return;
		}

		size_t indexStart=0, indexEnd=vecStationMetaData.size();

		//The following part decides whether all the stations are rebuffered or just one station
		if (stationindex == IOUtils::npos) {
			vecMeteo.clear();
			vecMeteo.insert(vecMeteo.begin(), vecStationMetaData.size(), vector<MeteoData>());
		} else {
			if (stationindex < vecMeteo.size()) {
				indexStart = stationindex;
				indexEnd   = stationindex+1;
			} else {
				throw IndexOutOfBoundsException("You tried to access a stationindex in readMeteoData that is out of bounds", AT);
			}
		}

		if ((env == NULL) || (conn == NULL))
			openDBConnection(env, conn);
		stmt = conn->createStatement();

		for (size_t ii=indexStart; ii<indexEnd; ii++) { //loop through relevant stations
			readData(dateStart, dateEnd, vecMeteo, ii, vecStationMetaData, env, stmt);
		}

		if (useAnetz) { //Important: we don't care about the metadata for ANETZ stations
			vector<StationData> vecAnetzStation;       //holds the unique ANETZ stations that need to be read
			vector< vector<MeteoData> > vecMeteoAnetz; //holds the meteo data of the ANETZ stations
			map<string, size_t> mapAnetzNames;   //associates an ANETZ station with an index within vecMeteoAnetz

			findAnetzStations(indexStart, indexEnd, mapAnetzNames, vecAnetzStation);
			vecMeteoAnetz.insert(vecMeteoAnetz.begin(), vecAnetzStation.size(), vector<MeteoData>());

			//date_anetz_start/end must be changed to be a multiple of 6h before the original dateStart, dateEnd
			Date date_anetz_start = Date(floor(dateStart.getJulian(true) * 4.0) / 4.0, 0.);
			date_anetz_start.setTimeZone(in_tz);
			Date date_anetz_end   = Date(floor(dateEnd.getJulian(true) * 4.0) / 4.0, 0.);
			date_anetz_end.setTimeZone(in_tz);

			//read Anetz Data
			for (size_t ii=0; ii<vecAnetzStation.size(); ii++)
				readData(date_anetz_start, dateEnd, vecMeteoAnetz, ii, vecAnetzStation, env, stmt);

			//We got all the data, now calc psum for all ANETZ stations
			vector< vector<double> > vec_of_psums; //6 hour accumulations of psum
			calculatePsum(date_anetz_start, date_anetz_end, vecMeteoAnetz, vec_of_psums);

			for (size_t ii=indexStart; ii<indexEnd; ii++){ //loop through relevant stations
				const map<string,AnetzData>::const_iterator it = mapAnetz.find(vecStationMetaData.at(ii).getStationID());
				if (it != mapAnetz.end())
					assimilateAnetzData(date_anetz_start, it->second, vec_of_psums, mapAnetzNames, ii, vecMeteo);
			}
		}

		if(use_psum_snowpack) {
			for (size_t ii=indexStart; ii<indexEnd; ii++) { //loop through relevant stations
				readSWE(dateStart, dateEnd, vecMeteo, ii, vecStationMetaData, env, stmt);
			}
		}

		conn->terminateStatement(stmt);
		closeDBConnection(env, conn);
	} catch (const exception& e){
		closeDBConnection(env, conn);
		throw IOException("Oracle Error when reading stations' data: " + string(e.what()), AT); //Translation of OCCI exception to IOException
	}
}

void ImisIO::assimilateAnetzData(const Date& dateStart, const AnetzData& ad,
                                 const std::vector< std::vector<double> > vec_of_psums,
                                 const std::map<std::string, size_t>& mapAnetzNames, const size_t& stationindex,
                                 std::vector< std::vector<MeteoData> >& vecMeteo)
{
	//Do coefficient calculation (getPSUM) for every single station and data point
	vector<double> current_station_psum;
	getAnetzPSUM(ad, mapAnetzNames, vec_of_psums, current_station_psum);

	size_t counter = 0;
	Date current_slice_date = dateStart;
	current_slice_date.setTimeZone(in_tz);
	for (size_t jj=0; jj<vecMeteo[stationindex].size(); jj++){
		while (vecMeteo[stationindex][jj].date > (current_slice_date+0.2485)){
			counter++;
			const double julian = floor((current_slice_date.getJulian(true) +0.25001) * 4.0) / 4.0;
			current_slice_date = Date(julian, 0.);
			current_slice_date.setTimeZone(in_tz);
		}

		if (counter >= current_station_psum.size()) { break; } //should never happen

		double& psum = vecMeteo[stationindex][jj](MeteoData::PSUM);
		if ((psum == IOUtils::nodata) || (IOUtils::checkEpsilonEquality(psum, 0.0, 0.001))){
			//replace by psum if there is no own value measured
			psum = current_station_psum.at(counter);
		}
	}
}

void ImisIO::getAnetzPSUM(const AnetzData& ad, const std::map<std::string, size_t>& mapAnetzNames,
                         const std::vector< std::vector<double> >& vec_of_psums, std::vector<double>& psum)
{
	vector<size_t> vecIndex; //this vector will hold up to three indexes for the Anetz stations (position in vec_of_psums)
	for (size_t ii=0; ii<ad.nrOfAnetzStations; ii++){
		const map<string, size_t>::const_iterator it = mapAnetzNames.find(ad.anetzstations[ii]);
		vecIndex.push_back(it->second);
	}

	if (ad.nrOfAnetzStations == ad.nrOfCoefficients){
		//1, 2, or 3 ANETZ stations without interaction
		for (size_t kk=0; kk<vec_of_psums.at(vecIndex.at(0)).size(); kk++){
			double sum = 0.0;
			for (size_t ii=0; ii<ad.nrOfCoefficients; ii++){
				sum += ad.coeffs[ii] * vec_of_psums.at(vecIndex[ii])[kk];
			}
			psum.push_back(sum/12.0);
		}
	} else {
		if (ad.nrOfCoefficients != 3)
			throw IOException("Misconfiguration in ANETZ data", AT);

		// Exactly two ANETZ stations with one interaction term
		for (size_t kk=0; kk<vec_of_psums.at(vecIndex.at(0)).size(); kk++){
			double sum = 0.0;
			const double& psum0 = vec_of_psums.at(vecIndex.at(0))[kk];
			const double& psum1 = vec_of_psums.at(vecIndex.at(1))[kk];
			sum += ad.coeffs[0] * psum0;
			sum += ad.coeffs[1] * psum1;
			sum += ad.coeffs[2] * psum0 * psum1;

			psum.push_back(sum/12.0);
		}
	}
}

void ImisIO::calculatePsum(const Date& dateStart, const Date& dateEnd,
                           const std::vector< std::vector<MeteoData> >& vecMeteoAnetz,
                           std::vector< std::vector<double> >& vec_of_psums)
{
	const unsigned int nr_of_slices = (unsigned int)((dateEnd.getJulian(true) - dateStart.getJulian(true) + 0.00001) * 4.0) + 1;

	for (size_t ii=0; ii<vecMeteoAnetz.size(); ii++){
		double tmp_psum = 0.0;
		Date current_date = dateStart;
		current_date.setTimeZone(in_tz);

		vector<double> vec_current_station;
		size_t counter_of_elements = 0;
		for (size_t jj=0; jj<vecMeteoAnetz[ii].size(); jj++){
			const Date& anetzdate = vecMeteoAnetz[ii][jj].date;
			const double& psum = vecMeteoAnetz[ii][jj](MeteoData::PSUM);

			if ((current_date < anetzdate) && ((current_date+0.25) > anetzdate)){
				;
			} else {
				if ((counter_of_elements > 0) && (counter_of_elements < 6)) //this is mystical, but kind of a guess of the future
					tmp_psum = tmp_psum * 6.0 / (double)counter_of_elements;

				vec_current_station.push_back(tmp_psum);

				current_date += 0.25;
				tmp_psum = 0.0;
				counter_of_elements = 0;
			}

			if (psum != IOUtils::nodata){
				tmp_psum += psum;
				counter_of_elements++;
			}

		}

		if ((counter_of_elements > 0) && (counter_of_elements < 6)) //this is mystical, but kind of a guess of the future
			tmp_psum = tmp_psum*6./static_cast<double>(counter_of_elements);

		vec_current_station.push_back(tmp_psum);

		for (size_t jj=vec_current_station.size(); jj<nr_of_slices; jj++){ //To fill up the vector
			vec_current_station.push_back(0.0);
		}

		vec_of_psums.push_back(vec_current_station);
	}

	for (size_t ii=1; ii<vec_of_psums.size(); ii++){
		if (vec_of_psums[ii].size() != vec_of_psums[ii-1].size())
			throw IOException("Error while summing up the precipitation data for the ANETZ stations", AT);
	}
}

void ImisIO::findAnetzStations(const size_t& indexStart, const size_t& indexEnd,
                               std::map<std::string, size_t>& mapAnetzNames,
                               std::vector<StationData>& vecAnetzStation)
{
	set<string> uniqueStations;

	for (size_t ii=indexStart; ii<indexEnd; ii++){ //loop through stations
		const map<string, AnetzData>::const_iterator it = mapAnetz.find(vecStationMetaData.at(ii).getStationID());
		if (it != mapAnetz.end()){
			for (size_t jj=0; jj<it->second.nrOfAnetzStations; jj++){
				uniqueStations.insert(it->second.anetzstations[jj]);
			}
		}
	}

	size_t pp = 0;
	for (set<string>::const_iterator ii=uniqueStations.begin(); ii!=uniqueStations.end(); ii++){
		mapAnetzNames[*ii] = pp;
		pp++;

		StationData sd;
		sd.stationID = *ii;
		vecAnetzStation.push_back(sd);
	}
}

/**
 * @brief A meta function to read meteo data for one specific station (specified by the stationindex)
 * @param dateStart     The beginning of the interval to retrieve data for
 * @param dateEnd       The end of the interval to retrieve data for
 * @param vecMeteo      The vector that will hold all MeteoData for each station
 * @param stationindex  The index of the station as specified in the Config
 * @param vecStationIDs Vector of station IDs
 * @param env           Create Oracle environnment
 * @param conn          Create connection to SDB
 */
void ImisIO::readData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo,
                      const size_t& stationindex, const std::vector<StationData>& vecStationIDs,
                      oracle::occi::Environment*& env, oracle::occi::Statement*& stmt)
{
	vecMeteo.at(stationindex).clear();

	string stat_abk, stao_nr;
	vector< vector<string> > vecResult;

	// Moving back to the IMIS timezone (UTC+1)
	Date dateS(dateStart), dateE(dateEnd);
	dateS.setTimeZone(in_tz);
	dateE.setTimeZone(in_tz);

	//get data for one specific station
	std::vector<std::string> vecHTS1;
	parseStationID(vecStationIDs.at(stationindex).getStationID(), stat_abk, stao_nr);
	getSensorDepths(stat_abk, stao_nr, sqlQuerySensorDepths, vecHTS1, stmt);
	bool fullStation = getStationData(stat_abk, stao_nr, dateS, dateE, vecHTS1, vecResult, env, stmt);

	MeteoData tmpmd;
	tmpmd.meta = vecStationIDs.at(stationindex);
	for (size_t ii=0; ii<vecResult.size(); ii++){
		parseDataSet(vecResult[ii], tmpmd, fullStation);
		convertUnits(tmpmd);

		//For IMIS stations the psum value is a rate (kg m-2 h-1), therefore we need to
		//divide it by two to conjure the accumulated value for the half hour
		if (tmpmd.meta.stationID.length() > 0){
			if (tmpmd.meta.stationID[0] != '*') { //only consider IMIS stations (ie: not ANETZ)
				if(use_imis_psum==false) {
					tmpmd(MeteoData::PSUM) = IOUtils::nodata;
				} else {
					double& psum = tmpmd(MeteoData::PSUM);
					if(psum!=IOUtils::nodata) {
						psum /= 2.; //half hour accumulated value for IMIS stations only
					}
				}
			}
		}

		vecMeteo.at(stationindex).push_back(tmpmd); //Now insert tmpmd
	}
}

/**
 * @brief Read simulated SWE from the database, compute Delta(SWE) and use it as PSUM for one specific station (specified by the stationindex)
 * @param dateStart     The beginning of the interval to retrieve data for
 * @param dateEnd       The end of the interval to retrieve data for
 * @param vecMeteo      The vector that will hold all MeteoData for each station
 * @param stationindex  The index of the station as specified in the Config
 * @param vecStationIDs Vector of station IDs
 * @param env           Create Oracle environnment
 * @param conn          Create connection to SDB
 */
void ImisIO::readSWE(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo,
                      const size_t& stationindex, const std::vector<StationData>& vecStationIDs,
                      oracle::occi::Environment*& env, oracle::occi::Statement*& stmt)
{
	const double max_interval = 3./24.; //max hours between two SWE values
	const double swe_threshold = 1.5; //precip less than this are delayed until its sum gets greater
	const double eps_swe = 0.1; //very small variations on SWE are simply ignored

	// Moving back to the IMIS timezone (UTC+1)
	Date dateS(dateStart), dateE(dateEnd);
	dateS.setTimeZone(in_tz);
	dateE.setTimeZone(in_tz);

	//build stat_abk and stao_nr from station name
	string stat_abk, stao_nr;
	parseStationID(vecStationIDs.at(stationindex).getStationID(), stat_abk, stao_nr);

	const unsigned int max_row = static_cast<unsigned int>( Optim::ceil( (dateE.getJulian()-dateS.getJulian())*24.*2. ) ); //for prefetching

	//query
	try {
		stmt->setSQL(sqlQuerySWEData);
		stmt->setPrefetchRowCount(max_row);

		// construct the oracle specific Date object: year, month, day, hour, minutes
		int year, month, day, hour, minutes, seconds;
		dateS.getDate(year, month, day, hour, minutes, seconds);
		const occi::Date begindate(env, year, month, day, hour, minutes, seconds);
		dateE.getDate(year, month, day, hour, minutes, seconds);
		const occi::Date enddate(env, year, month, day, hour, minutes, seconds);
		stmt->setString(1, stat_abk); // set 1st variable's value (station name)
		stmt->setString(2, stao_nr);  // set 2nd variable's value (station number)
		stmt->setDate(3, begindate);  // set 3rd variable's value (begin date)
		stmt->setDate(4, enddate);    // set 4th variable's value (end date)

		ResultSet *rs = stmt->executeQuery(); // execute the statement stmt
		const vector<MetaData> cols = rs->getColumnListMetaData();

		double prev_swe = IOUtils::nodata;
		Date prev_date;
		size_t ii_serie = 0; //index in meteo time serie
		const size_t serie_len = vecMeteo[stationindex].size();
		double accumulator = 0.;

		while (rs->next() == true) { //loop over timesteps
			if (cols.size()!=2) {
				ostringstream ss;
				ss << "For station " << vecStationIDs.at(stationindex).getStationID() << ", ";
				ss << "snowpack SWE query returned " << cols.size() << " columns, while 2 were expected";
				throw UnknownValueException(ss.str(), AT);
			}
			Date curr_date;
			IOUtils::convertString(curr_date, rs->getString(1), 1.);
			double curr_swe;
			IOUtils::convertString(curr_swe, rs->getString(2));
			if (curr_swe==IOUtils::nodata || curr_swe<0.) continue;

			//looking for matching timestamp in the vecMeteo
			while (ii_serie<serie_len && vecMeteo[stationindex][ii_serie].date<curr_date) ii_serie++;
			if (ii_serie>=serie_len) return;


			if (prev_swe==IOUtils::nodata) {
				 //this looks like the first valid data point that we find
				prev_swe = curr_swe;
				prev_date = curr_date;
				continue;
			}

			if ((curr_date.getJulian()-prev_date.getJulian())<=max_interval || curr_swe==0.) {
				vecMeteo[stationindex][ii_serie](MeteoData::PSUM) = 0.;
				//data not too far apart, so we accept it for Delta SWE
				if (vecMeteo[stationindex][ii_serie].date==curr_date) {
					//we found the matching timestamp -> writing Delta(SWE) as psum
					const double new_psum_sum = curr_swe - prev_swe;
					if (new_psum_sum>eps_swe) {
						accumulator += new_psum_sum;
						if (accumulator>=swe_threshold) {
							vecMeteo[stationindex][ii_serie](MeteoData::PSUM) = accumulator;
							accumulator = 0.;
						}
					}
				}
				prev_swe = curr_swe;
				prev_date = curr_date;
			} else {
				//data points in SWE too far apart, we could not use it for psum but we reset our prev_swe to this new point
				prev_swe = curr_swe;
				prev_date = curr_date;
			}
		}

		stmt->closeResultSet(rs);
	} catch (const exception& e){
		throw IOException("Oracle Error when SWE data: " + string(e.what()), AT); //Translation of OCCI exception to IOException
	}
}


/**
 * @brief Puts the data that has been retrieved from the database into a MeteoData object
 * @param i_meteo a row of meteo data from the database (NOTE order important, matches SQL query, see also MeteoData.[cch])
 * @param md     the object to copy the data to
 * @param fullstation
 * 	- true if it is a combined snow_drift station (station2.v_snow_drift_standort)
 * 	- false if it is a "conventional" station, for example an ANETZ-station (station2.standort)
 */
void ImisIO::parseDataSet(const std::vector<std::string>& i_meteo, MeteoData& md, bool& fullStation)
{
	IOUtils::convertString(md.date, i_meteo.at(0), in_tz, dec);
	IOUtils::convertString(md(MeteoData::TA),     i_meteo.at(1),  std::dec);
	IOUtils::convertString(md(MeteoData::ISWR),   i_meteo.at(2),  std::dec);
	IOUtils::convertString(md(MeteoData::VW),     i_meteo.at(3),  std::dec);
	IOUtils::convertString(md(MeteoData::DW),     i_meteo.at(4),  std::dec);
	IOUtils::convertString(md(MeteoData::VW_MAX), i_meteo.at(5),  std::dec);
	IOUtils::convertString(md(MeteoData::RH),     i_meteo.at(6),  std::dec);
	IOUtils::convertString(md(MeteoData::ILWR),   i_meteo.at(7),  std::dec);
	IOUtils::convertString(md(MeteoData::PSUM),    i_meteo.at(8),  std::dec);
	IOUtils::convertString(md(MeteoData::TSG),    i_meteo.at(9),  std::dec);
	IOUtils::convertString(md(MeteoData::TSS),    i_meteo.at(10), std::dec);
	IOUtils::convertString(md(MeteoData::HS),     i_meteo.at(11), std::dec);
	IOUtils::convertString(md(MeteoData::RSWR),   i_meteo.at(12), std::dec);

	unsigned int ii = 13;
	if (fullStation) {
		if (!md.param_exists("VW_DRIFT")) md.addParameter("VW_DRIFT");
		IOUtils::convertString(md("VW_DRIFT"), i_meteo.at(ii++), std::dec);
		if (!md.param_exists("DW_DRIFT")) md.addParameter("DW_DRIFT");
		IOUtils::convertString(md("DW_DRIFT"), i_meteo.at(ii++), std::dec);
	}

	// additional snow station parameters
	if (!md.param_exists("TS1")) md.addParameter("TS1");
	IOUtils::convertString(md("TS1"), i_meteo.at(ii++), std::dec);
	if (!md.param_exists("TS2")) md.addParameter("TS2");
	IOUtils::convertString(md("TS2"), i_meteo.at(ii++), std::dec);
	if (!md.param_exists("TS3")) md.addParameter("TS3");
	IOUtils::convertString(md("TS3"), i_meteo.at(ii++), std::dec);
	if (fullStation) {
		if (!md.param_exists("HTS1")) md.addParameter("HTS1");
		IOUtils::convertString(md("HTS1"), i_meteo.at(ii++), std::dec);
		if (!md.param_exists("HTS2")) md.addParameter("HTS2");
		IOUtils::convertString(md("HTS2"), i_meteo.at(ii++), std::dec);
		if (!md.param_exists("HTS3")) md.addParameter("HTS3");
		IOUtils::convertString(md("HTS3"), i_meteo.at(ii++), std::dec);
	}
}

/**
 * @brief This function gets IDs from table station2.v_snow_drift_standort and fills vecStationIDs
 * @param station_code  a string key corresponding to stationID
 * @param vecStationIDs string vector in which data will be filled
 * @param conn          create connection to SDB
 * @param return number of columns retrieved
 */
size_t ImisIO::getStationIDs(const std::string& station_code, const std::string& sqlQuery,
                                   std::vector<std::string>& vecStationIDs,
                                   oracle::occi::Statement*& stmt)
{
	vecStationIDs.clear();

	try {
		stmt->setSQL(sqlQuery);
		stmt->setString(1, station_code); // set 1st variable's value

		ResultSet *rs = stmt->executeQuery();    // execute the statement stmt
		const vector<MetaData> cols = rs->getColumnListMetaData();

		while (rs->next() == true) {
			for (unsigned int ii=1; ii<=static_cast<unsigned int>(cols.size()); ii++) {
				vecStationIDs.push_back(rs->getString(ii));
			}
		}

		if (vecStationIDs.size() < 3) { //if the station has not been found
			string stat_abk, stao_nr;
			parseStationID(station_code, stat_abk, stao_nr);
			vecStationIDs.push_back(station_code);
			vecStationIDs.push_back(stat_abk);
			vecStationIDs.push_back(stao_nr);
		}

		stmt->closeResultSet(rs);
		return cols.size();
	} catch (const exception& e){
		throw IOException("Oracle Error when reading stations' id: " + string(e.what()), AT); //Translation of OCCI exception to IOException
	}
}

/**
 * @brief This function gets IDs from table station2.v_snow_drift_standort and fills vecStationIDs
 * @param stat_abk a string key of table station2
 * @param stao_nr  a string key of table station2
 * @param vecHTS1   vector of string to retieve sensor depths
 * @param conn     create connection to SDB
 * @param return number of columns retrieved
 */
size_t ImisIO::getSensorDepths(const std::string& stat_abk, const std::string& stao_nr,
                                     const std::string& sqlQuery, std::vector<std::string>& vecHTS1,
                                     oracle::occi::Statement*& stmt)
{
	vecHTS1.clear();

	try {
		stmt->setSQL(sqlQuery);
		stmt->setString(1, stat_abk); // set 1st variable's value
		stmt->setString(2, stao_nr);  // set 2nd variable's value

		ResultSet *rs = stmt->executeQuery();    // execute the statement stmt
		const vector<MetaData> cols = rs->getColumnListMetaData();

		while (rs->next() == true) {
			for (unsigned int ii=1; ii<=static_cast<unsigned int>(cols.size()); ii++) {
				vecHTS1.push_back(rs->getString(ii));
			}
		}

		stmt->closeResultSet(rs);
		return cols.size();
	} catch (const exception& e){
		throw IOException("Oracle Error when reading sensors' depths: " + string(e.what()), AT); //Translation of OCCI exception to IOException
	}
}

/**
 * @brief This function gets meta data from table station2.standort and fills vecStationMetaData.
 * This is also the moment to take the opportunity to check if the station really does exist.
 * @param stat_abk           a string key of table
 * @param stao_nr            a string key of table
 * @param sqlQuery           the query to execute
 * @param vecStationMetaData string vector in which data will be filled
 * @param conn               create connection to SDB
 * @param return             number of metadata read (ie. retrieve and not NULL)
 */
size_t ImisIO::getStationMetaData(const std::string& stat_abk, const std::string& stao_nr,
                                        const std::string& sqlQuery, std::vector<std::string>& vecMetaData,
                                        oracle::occi::Statement*& stmt)
{
	vecMetaData.clear();

	try {
		stmt->setSQL(sqlQuery);
		stmt->setString(1, stat_abk); // set 1st variable's value
		stmt->setString(2, stao_nr);  // set 2nd variable's value

		ResultSet *rs = stmt->executeQuery();    // execute the statement stmt
		const vector<MetaData> cols = rs->getColumnListMetaData();

		while (rs->next() == true) {
			for (unsigned int ii=1; ii<=static_cast<unsigned int>(cols.size()); ii++) {
				vecMetaData.push_back(rs->getString(ii));
			}
		}

		stmt->closeResultSet(rs);
	} catch (const exception& e){
		throw IOException("Oracle Error when reading stations' metadata: " + string(e.what()), AT); //Translation of OCCI exception to IOException
	}

	const size_t nr_metadata = vecMetaData.size();
	if(nr_metadata==0)
			throw NoAvailableDataException("Station " + stat_abk+stao_nr + " not found in the database", AT);
	if(nr_metadata<4)
			throw ConversionFailedException("Error while converting station meta data for station "+stat_abk+stao_nr, AT);
	return nr_metadata;
}

/**
 * @brief Gets data from ams.v_ams_raw which is a table of SDB and
 * retrieves the temperature sensor depths from station2.standort \n
 * Each record returned are vector of strings which are pushed back in vecMeteoData.
 * @param stat_abk :     a string key of ams.v_ams_raw
 * @param stao_nr :      a string key of ams.v_ams_raw
 * @param dateS :        begining of the recording date
 * @param dateE :        end of the recording date
 * @param vecMeteoData : a vector of vector of string in which data will be filled
 * @param return number of columns retrieved
 */
bool ImisIO::getStationData(const std::string& stat_abk, const std::string& stao_nr,
                            const Date& dateS, const Date& dateE,
                            const std::vector<std::string>& vecHTS1,
                            std::vector< std::vector<std::string> >& vecMeteoData,
                            oracle::occi::Environment*& env, oracle::occi::Statement*& stmt)
{
	vecMeteoData.clear();
	bool fullStation = true;
	const unsigned int max_row = static_cast<unsigned int>( Optim::ceil( (dateE.getJulian()-dateS.getJulian())*24.*2. ) ); //for prefetching
	try {
		const map<string, string>::const_iterator it = mapDriftStation.find(stat_abk+stao_nr);
		if (it != mapDriftStation.end()) {
			stmt->setSQL(sqlQueryMeteoDataDrift);
			string drift_stat_abk, drift_stao_nr;
			parseStationID(it->second, drift_stat_abk, drift_stao_nr);
			stmt->setString(5, drift_stat_abk);
			stmt->setString(6, drift_stao_nr);
		} else {
			stmt->setSQL(sqlQueryMeteoData);
			fullStation = false;
		}
		stmt->setPrefetchRowCount(max_row);

		// construct the oracle specific Date object: year, month, day, hour, minutes
		int year, month, day, hour, minutes, seconds;
		dateS.getDate(year, month, day, hour, minutes, seconds);
		const occi::Date begindate(env, year, month, day, hour, minutes, seconds);
		dateE.getDate(year, month, day, hour, minutes, seconds);
		const occi::Date enddate(env, year, month, day, hour, minutes, seconds);
		stmt->setString(1, stat_abk); // set 1st variable's value (station name)
		stmt->setString(2, stao_nr);  // set 2nd variable's value (station number)
		stmt->setDate(3, begindate);  // set 3rd variable's value (begin date)
		stmt->setDate(4, enddate);    // set 4th variable's value (end date)

		ResultSet *rs = stmt->executeQuery(); // execute the statement stmt
		const vector<MetaData> cols = rs->getColumnListMetaData();

		vector<string> vecData;
		while (rs->next() == true) {
			vecData.clear();
			for (unsigned int ii=1; ii<=static_cast<unsigned int>(cols.size()); ii++) {
				vecData.push_back(rs->getString(ii));
			}
			if (fullStation) {
				for (unsigned int ii=0; ii<static_cast<unsigned int>(vecHTS1.size()); ii++) {
					vecData.push_back(vecHTS1.at(ii));
				}
			}
			vecMeteoData.push_back(vecData);
		}

		stmt->closeResultSet(rs);
		return fullStation;
	} catch (const exception& e){
		throw IOException("Oracle Error when reading stations' data: " + string(e.what()), AT); //Translation of OCCI exception to IOException
	}
}

void ImisIO::convertSnowTemperature(MeteoData& meteo, const std::string& parameter)
{
	if (meteo.param_exists(parameter)) {
		const size_t idx = meteo.getParameterIndex(parameter);
		if(meteo(idx)!=IOUtils::nodata)
			meteo(idx) += Cst::t_water_freezing_pt; //C_TO_K
	}
}

void ImisIO::convertSensorDepth(MeteoData& meteo, const std::string& parameter)
{
	if (meteo.param_exists(parameter)) {
		const size_t idx = meteo.getParameterIndex(parameter);
		if(meteo(idx)!=IOUtils::nodata)
			meteo(idx) /= 100.; // centimetre to metre
	}
}

void ImisIO::convertUnits(MeteoData& meteo)
{
	meteo.standardizeNodata(plugin_nodata);

	//converts C to Kelvin, converts RH to [0,1]
	double& ta = meteo(MeteoData::TA);
	ta = IOUtils::C_TO_K(ta);

	double& tsg = meteo(MeteoData::TSG);
	tsg = IOUtils::C_TO_K(tsg);

	double& tss = meteo(MeteoData::TSS);
	tss = IOUtils::C_TO_K(tss);

	double& rh = meteo(MeteoData::RH);
	if (rh != IOUtils::nodata)
		rh /= 100.;

	double& hs = meteo(MeteoData::HS);
	if (hs != IOUtils::nodata)
		hs /= 100.0;

	//convert extra parameters (if present) //HACK TODO: find a dynamic way...
	convertSnowTemperature(meteo, "TS1");
	convertSnowTemperature(meteo, "TS2");
	convertSnowTemperature(meteo, "TS3");
	convertSensorDepth(meteo, "HTS1");
	convertSensorDepth(meteo, "HTS2");
	convertSensorDepth(meteo, "HTS3");
}

void ImisIO::cleanup() throw()
{
}

} //namespace
