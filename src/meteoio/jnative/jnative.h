/* * Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
 * modular unstructured mesh based approach for hydrological modelling
 * Copyright (C) 2018 Christopher Marsh
 *
 * This file is part of Canadian Hydrological Model.
 *
 * Canadian Hydrological Model is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Canadian Hydrological Model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Canadian Hydrological Model.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

/*
 * jnative.h
 *
 *  Created on: 08.01.2010
 *      Author: perot
 */


#ifndef _Included_JNative
#define _Included_JNative

#if  defined(_METEOIO_JNI) ||  defined(_METEOIO_JNA)

#include <meteoio/MeteoIO.h>


void loadMeteoAndStationData(double* cMetadata, double* cData,
		const int nbStation, const int nbDataPerStation,
		const std::string& algorithm,
		const std::string metaCoordinateSystem, std::vector<mio::StationData>& vecStation, 
                std::vector<mio::MeteoData>& vecData, enum mio::MeteoData::Parameters& interpolation_type);

void fulfillDoubleArray(const mio::Grid2DObject&  p, const std::string& cellOrder, double* dest);


#endif // defined(_METEOIO_JNI) ||  defined(_METEOIO_JNA)


#ifdef _METEOIO_JNA

/**
 *
 * Originally, these methods are dedicated to be called from JAVA with JNA framework.
 *
 *
 */


double* executeInterpolationSubDem
  (char*, char*, char*,char*, double, double, double, double,double*, int, double*, int, char*, char*);


double* executeInterpolation
(char*, char*, char*, char*, double*, int, double*, int, char*, char*);


#endif //_METEOIO_JNA

#endif//_Included_JNative
