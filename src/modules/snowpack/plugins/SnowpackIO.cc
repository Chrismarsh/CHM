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

#include <snowpack/plugins/SnowpackIO.h>

using namespace std;
using namespace mio;

SnowpackIO::SnowpackIO(const SnowpackConfig& cfg):
#ifdef IMISDBIO
	imisdbio(NULL),
#endif
#ifdef CAAMLIO
	caamlio(NULL),
#endif
	smetio(NULL), asciiio(NULL),
	input_snow_as_smet(false), output_snow_as_smet(false),
	input_snow_as_caaml(false), output_snow_as_caaml(false),
	input_snow_as_ascii(false), output_snow_as_ascii(false),
	output_prf_as_ascii(false), output_prf_as_caaml(false), output_prf_as_imis(false),
	output_ts_as_ascii(false), output_haz_as_imis(false)

{
	//Format of initial snow profile:
	//TODO: document Input::SNOW = SMET, CAAML, or SNOOLD
	const string in_snow = cfg.get("SNOW", "Input", IOUtils::nothrow);
	if (in_snow == "SNOOLD") {
		input_snow_as_ascii = true;
	} else if (in_snow == "CAAML") {
		input_snow_as_caaml = true;
	} else if (in_snow == "SMET") {
		input_snow_as_smet = true;
	} else
		throw InvalidArgumentException("Invalid input snow profile format '"+in_snow+"'. Please choose from SMET, CAAML, SNOOLD", AT);

	//Format of transitional and final snow profile(s):
	//TODO: document ouput::SNOW = SMET, CAAML, or SNOOLD
	const string out_snow = cfg.get("SNOW", "Output", IOUtils::nothrow);
	if (out_snow == "SNOOLD") {
		output_snow_as_ascii = true;
	} else if (out_snow == "CAAML") {
		output_snow_as_caaml = true;
	} else if (out_snow == "SMET") {
		output_snow_as_smet = true;
	} else
		throw InvalidArgumentException("Invalid output snow profile format '"+out_snow+"'. Please choose from SMET, CAAML, SNOOLD", AT);

	/* Profiles may also be dumped in up to 4 formats specified by the key PROFILE_FORMAT in section [Output]:
	 * PRO   : ASCII-format for visulization with SN_GUI, includes soil elements
	 * PRF   : Snow profiles in tabular ASCII-format, aggregated if AGGREGATE_PRF is set
	 * IMIS  : aggregated profiles for upload to the SLF database sdbo
	 */
	std::vector<string> vecProfileFmt = cfg.get("PROFILE_FORMAT", "Output", IOUtils::nothrow);
	if (vecProfileFmt.size() > 4) {
		throw InvalidArgumentException("The key PROFILE_FORMAT in section [Output] can have three values at most", AT);
	} else {
		for (size_t ii=0; ii<vecProfileFmt.size(); ii++) {
			if (vecProfileFmt[ii] == "PRO") {
				output_prf_as_ascii = true;
			} else if (vecProfileFmt[ii] == "PRF") {
				output_prf_as_ascii  = true;
			} else if (vecProfileFmt[ii] == "IMIS") {
				output_prf_as_imis  = true;
			} else {
				throw InvalidArgumentException("Key PROFILE_FORMAT in section [Output] takes only PRO, PRF or IMIS values", AT);
			}
		}
	}
	//Format of meteo time series:
	output_ts_as_ascii = cfg.get("TS_WRITE", "Output", IOUtils::nothrow);

	//set the "plugins" pointers
	RunInfo run_info;
	if (input_snow_as_smet || output_snow_as_smet) smetio = new SmetIO(cfg, run_info);
	if (input_snow_as_ascii || output_snow_as_ascii || output_prf_as_ascii || output_ts_as_ascii) asciiio = new AsciiIO(cfg, run_info);
#ifdef CAAMLIO
	if (input_snow_as_caaml || output_snow_as_caaml) caamlio = new CaaMLIO(cfg, run_info);
#endif
#ifdef IMISDBIO
	output_haz_as_imis = output_prf_as_imis;
	if (output_prf_as_imis || output_haz_as_imis) imisdbio = new ImisDBIO(cfg, run_info);
#endif
}

SnowpackIO::SnowpackIO(const SnowpackIO& source) :
#ifdef IMISDBIO
	imisdbio(source.imisdbio),
#endif
#ifdef CAAMLIO
	caamlio(source.caamlio),
#endif
	smetio(source.smetio), asciiio(source.asciiio),
	input_snow_as_smet(source.input_snow_as_smet), output_snow_as_smet(source.input_snow_as_smet),
	input_snow_as_caaml(source.input_snow_as_caaml), output_snow_as_caaml(source.output_snow_as_caaml),
	input_snow_as_ascii(source.input_snow_as_ascii), output_snow_as_ascii(source.output_snow_as_ascii),
	output_prf_as_ascii(source.output_prf_as_ascii), output_prf_as_caaml(source.output_prf_as_caaml), output_prf_as_imis(source.output_prf_as_imis),
	output_ts_as_ascii(source.output_ts_as_ascii), output_haz_as_imis(source.output_haz_as_imis)
{}

SnowpackIO::~SnowpackIO()
{
	if (smetio != NULL) delete smetio;
	if (asciiio != NULL) delete asciiio;
#ifdef CAAMLIO
	if (caamlio != NULL) delete caamlio;
#endif
#ifdef IMISDBIO
	if (imisdbio != NULL) delete imisdbio;
#endif
}

bool SnowpackIO::snowCoverExists(const std::string& i_snowfile, const std::string& stationID) const
{
	if (input_snow_as_ascii) {
		return asciiio->snowCoverExists(i_snowfile, stationID);
#ifdef CAAMLIO
	} else if (input_snow_as_caaml){
		return caamlio->snowCoverExists(i_snowfile, stationID);
#endif
	} else {
		return smetio->snowCoverExists(i_snowfile, stationID);
	}
}

void SnowpackIO::readSnowCover(const std::string& i_snowfile, const std::string& stationID,
                               SN_SNOWSOIL_DATA& SSdata, ZwischenData& Zdata)
{
	if (input_snow_as_ascii) {
		asciiio->readSnowCover(i_snowfile, stationID, SSdata, Zdata);
#ifdef CAAMLIO
	} else if (input_snow_as_caaml) {
		caamlio->readSnowCover(i_snowfile, stationID, SSdata, Zdata);
#endif
	} else {
		smetio->readSnowCover(i_snowfile, stationID, SSdata, Zdata);
	}
}

void SnowpackIO::writeSnowCover(const mio::Date& date, const SnowStation& Xdata,
                                const ZwischenData& Zdata, const bool& forbackup)
{
	if (output_snow_as_ascii) {
		asciiio->writeSnowCover(date, Xdata, Zdata, forbackup);
#ifdef CAAMLIO
	} else if (output_snow_as_caaml) {
		caamlio->writeSnowCover(date, Xdata, Zdata, forbackup);
#endif
	} else {
		smetio->writeSnowCover(date, Xdata, Zdata, forbackup);
	}
}

void SnowpackIO::writeTimeSeries(const SnowStation& Xdata, const SurfaceFluxes& Sdata, const CurrentMeteo& Mdata,
                                 const ProcessDat& Hdata, const double wind_trans24)
{
	if (output_ts_as_ascii)
		asciiio->writeTimeSeries(Xdata, Sdata, Mdata, Hdata, wind_trans24);
}

void SnowpackIO::writeProfile(const mio::Date& date, const SnowStation& Xdata)
{
	if (output_prf_as_ascii)
		asciiio->writeProfile(date, Xdata);

#ifdef CAAMLIO
	if (output_prf_as_caaml)
		caamlio->writeProfile(date, Xdata);
#endif

#ifdef IMISDBIO
	if (output_prf_as_imis)
		imisdbio->writeProfile(date, Xdata);
#endif
}

#ifdef IMISDBIO
bool SnowpackIO::writeHazardData(const std::string& stationID, const std::vector<ProcessDat>& Hdata,
                                 const std::vector<ProcessInd>& Hdata_ind, const size_t& num)
{
	if(output_haz_as_imis)
		return imisdbio->writeHazardData(stationID, Hdata, Hdata_ind, num);
	return false;
}
#else
bool SnowpackIO::writeHazardData(const std::string& /*stationID*/, const std::vector<ProcessDat>& /*Hdata*/,
                                 const std::vector<ProcessInd>& /*Hdata_ind*/, const size_t& /*num*/)
{
	return false;
}
#endif

SnowpackIO& SnowpackIO::operator=(const SnowpackIO& source)
{
	if(this != &source) {
#ifdef IMISDBIO
		imisdbio = source.imisdbio;
#endif
#ifdef CAAMLIO
		caamlio = source.caamlio;
#endif
		asciiio = source.asciiio;
		smetio = source.smetio;
		output_prf_as_ascii = source.output_prf_as_ascii;
		output_prf_as_caaml = source.output_prf_as_caaml;
		output_prf_as_imis = source.output_prf_as_imis;
		output_snow_as_caaml = source.output_snow_as_caaml;
		output_snow_as_smet = source.output_snow_as_smet;
		input_snow_as_caaml = source.input_snow_as_caaml;
		input_snow_as_smet = source.input_snow_as_smet;
		output_ts_as_ascii = source.output_ts_as_ascii;
		output_haz_as_imis = source.output_haz_as_imis;
	}
	return *this;
}

