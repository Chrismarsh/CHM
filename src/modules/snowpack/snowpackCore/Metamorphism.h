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
/**
 * @file Metamorphism.h
 */

#ifndef __METAMORPHISM_H__
#define __METAMORPHISM_H__

#include <snowpack/Constants.h>
#include <snowpack/DataClasses.h>
#include <snowpack/Utils.h>
#include <snowpack/Laws_sn.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <fcntl.h>
#include <meteoio/MeteoIO.h>
#include <map>
#include <string>

class SnowStation;
class Metamorphism;

typedef void (Metamorphism::*MetaModelFn)(const CurrentMeteo&, SnowStation&);
typedef double (Metamorphism::*MetaSpRateFn)(const ElementData&);

class Metamorphism {

	public:
		Metamorphism(const SnowpackConfig& i_cfg);

		void runMetamorphismModel(const CurrentMeteo& Mdata, SnowStation& Xdata) throw();

		static double csPoreArea(const ElementData& Edata);

		static double getCoordinationNumberN3(const double& Rho);

		static double ddRate(const ElementData& Edata);

		static const double mm_tg_dpdz, ba_g_fudge, sa_g_fudge, max_grain_growth, bond_size_stop;
		static const double max_grain_bond_ratio, wind_slab_enhance, wind_slab_vw, wind_slab_depth;

	private:
		double TGBondRate(const ElementData& Edata);

		double LatticeConstant0(const double& th_ice);

		double TGGrainRate(const ElementData& Edata, const double& Tbot, const double& Ttop,
		                   const double& gradTSub, const double& gradTSup);

		double ETBondRate(ElementData& Edata);
		double ETGrainRate(const ElementData& Edata);

		double PressureSintering(ElementData& Edata);

		void metamorphismDEFAULT(const CurrentMeteo& Mdata, SnowStation& Xdata);
		void metamorphismNIED(const CurrentMeteo& Mdata, SnowStation& Xdata);


		double spRateDEFAULT(const ElementData& Edata);
		double spRateNIED(const ElementData& Edata);

		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static containers
		static std::map<std::string, MetaModelFn> mapMetamorphismModel;
		static std::map<std::string, MetaSpRateFn> mapSpRate;

		std::string metamorphism_model;
		double sn_dt, new_snow_grain_size;
};

#endif //End of Metamorphism.h
