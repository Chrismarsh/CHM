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
#ifndef PROCESSINGBLOCK_H
#define PROCESSINGBLOCK_H

#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/Config.h>
#include <vector>
#include <string>
#include <set>

#ifdef _MSC_VER
	#pragma warning(disable:4512) //we don't need any = operator!
#endif

namespace mio {

class ProcessingProperties {
	public:
		typedef enum PROC_STAGE { none, ///< never activate this block
		                     first, ///< activate at first stage
		                     second, ///< activate at second stage
		                     both ///< activate at both first and second stage
		                     //once ///< activate at stage one or two, but only once
		                   } proc_stage;

		ProcessingProperties() : time_before(0., 0.), time_after(0., 0.),
		                         points_before(0), points_after(0),
		                         stage(first) {}

		const std::string toString() const;

		Duration time_before;
		Duration time_after;

		size_t points_before;
		size_t points_after;

		proc_stage stage;
};

/**
 * @class  ProcessingBlock
 * @brief  An abstract class
 * @author Thomas Egger
 * @date   2011-01-02
 */
class ProcessingBlock {
	public:
		virtual ~ProcessingBlock();

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec) = 0;

		std::string getName() const;
		const ProcessingProperties& getProperties() const;
		const std::string toString() const;

	protected:
		ProcessingBlock(const std::string& name); ///< protected constructor only to be called by children

		void convert_args(const size_t& min_nargs, const size_t& max_nargs,
		                  const std::vector<std::string>& vec_args, std::vector<double>& dbl_args) const;
		static bool is_soft(std::vector<std::string>& vec_args);
		static void readCorrections(const std::string& filter, const std::string& filename, const char& c_type, const double& init, std::vector<double> &corrections);

		ProcessingProperties properties;
		const std::string block_name;

		static const double soil_albedo, snow_albedo, snow_thresh; ///< parametrize the albedo from HS
};

class BlockFactory {
	public:
		static ProcessingBlock* getBlock(const std::string& blockname, const std::vector<std::string>& vec_args, const Config& cfg);
};

} //end namespace

#endif
