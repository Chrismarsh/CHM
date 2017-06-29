/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#ifndef EXPORTS_H
#define EXPORTS_H

// Generic helper definitions for shared library support
#if defined _WIN32 || defined __MINGW32__ || defined __CYGWIN__
	#define MIO_HELPER_DLL_IMPORT __declspec(dllimport)
	#define MIO_HELPER_DLL_EXPORT __declspec(dllexport)
	#define MIO_HELPER_DLL_LOCAL
#else
	#if defined __clang__ || __GNUC__ >= 4
		#define MIO_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
		#define MIO_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
		#define MIO_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
	#else
		#define MIO_HELPER_DLL_IMPORT
		#define MIO_HELPER_DLL_EXPORT
		#define MIO_HELPER_DLL_LOCAL
	#endif
#endif

// Now we use the generic helper definitions above to define
// MIO_API (for public API symbols) and MIO_LOCAL (non-api symbols)
#ifdef MIO_DLL //dynamic compile
	#ifdef MIO_DLL_EXPORTS //building as a dll
		#define MIO_API MIO_HELPER_DLL_EXPORT
	#else
		#define MIO_API MIO_HELPER_DLL_IMPORT
	#endif
	#define MIO_LOCAL MIO_HELPER_DLL_LOCAL
#else //static compile
	#define MIO_API
	#define MIO_LOCAL
#endif

#endif
