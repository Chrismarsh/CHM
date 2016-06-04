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
#include <cmath>
#include <cstdio>
#include <iomanip>

#include <meteoio/dataClasses/Coords.h>
#include <meteoio/MathOptim.h>
#include <meteoio/meteoLaws/Meteoconst.h> //for math constants

#ifdef PROJ4
	#include <proj_api.h>
#endif

using namespace std;

namespace mio {
 /**
 * @page coords Available coordinate systems
 * Geographic coordinates will be transparently and automatically converted to lat/lon and any other coordinate system that
 * the client program uses. However, in order to do so, the input coordinate system must be specified. In order to output
 * geolocalized data, the desired coordinate system must also be specified for the outputs (in the output section).
 * This is done through the use of the COORDIN and COORDPARAM keys (see the documentation for each plugin).
 *
 * There are two ways of supporting a given coordinate system: through the use of an adhoc implementation
 * (that becomes part of MeteoIO) or through the use of an external library, Proj4 [ref: http://trac.osgeo.org/proj/].
 * The current internal implementations are the following (given by their keyword):
 * - CH1903 or CH1903+ for coordinates in the Swiss Grid [ref: http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf]
 * - UTM for UTM coordinates (the zone must be specified in the parameters, for example 31T) [ref: http://www.oc.nps.edu/oc2902w/maps/utmups.pdf]
 * - UPS for Universal Polar Stereographic coordinates (the zone, either N or S, must be specified in the parameters). [ref: J. Hager, J. Behensky, B. Drew, <i>THE UNIVERSAL GRIDS: Universal Transverse Mercator (UTM) and Universal Polar Stereographic (UPS)</i>, 1989, Defense Mapping Agency, DMATM 8358.2]
 * - LOCAL for local coordinate system (using the horizontal and vertical distance from a reference point, see Coords::geo_distances for the available choice of distance algorithms)
 *
 * Such an example of use is the following:
 * @code
 * COORDSYS	= UTM
 * COORDPARAM	= 31T
 * @endcode
 *
 * On the other hand, when using the Proj4 library for handling the coordinate conversion, the EPSG codes of
 * the chosen projection must be specified (such codes can be found at http://spatialreference.org/ref/epsg/?page=1)
 * as illustrated below (21781 is the EPSG code for the CH1903 coordinate system. Such a code is 32767 at the maximum):
 * @code
 * COORDSYS	= PROJ4
 * COORDPARAM	= 21781
 * @endcode
 *
 */

//Adding a coordinate system is done by editing convert_to_WGS84 and convert_from_WGS84
//and implementing the XXX_to_WGS84 and WGS84_to_XXX methods.

const struct Coords::ELLIPSOID Coords::ellipsoids[6] = {
	{ 6378137.,	6356752.3142 }, ///< E_WGS84
	{ 6378137.,	6356752.3141 }, ///< E_GRS80
	{ 6377563.396,	6356256.909 }, ///< E_AIRY
	{ 6378388.,	6356911.946 }, ///< E_INTL1924
	{ 6378249.145,	6356514.86955 }, ///< E_CLARKE1880
	{ 6378160.,	6356774.719 } ///< E_GRS67
};

/**
* @brief Equality operator that checks that two coordinates objects represent the same 3D point.
* The comparison checks either lat/lon/alt or easting/northing/alt. If both objects have nodata coordinates,
* then they are equal (even if the internal projections might be set to different systems).
* @param[in] in Coord object to compare to
* @return true or false
*/
bool Coords::operator==(const Coords& in) const {
	//check on lat/lon
	if(latitude!=IOUtils::nodata && longitude!=IOUtils::nodata) {
		const bool comparison = ( IOUtils::checkEpsilonEquality(getLat(), in.getLat(), IOUtils::lat_epsilon) &&
		                          IOUtils::checkEpsilonEquality(getLon(), in.getLon(), IOUtils::lon_epsilon) &&
		                          IOUtils::checkEpsilonEquality(getAltitude(), in.getAltitude(), IOUtils::grid_epsilon) );
		return comparison;
	}
	//check on easting/northing
	if(easting!=IOUtils::nodata && northing!=IOUtils::nodata) {
		//in this case, it means that we don't know anything about the projection parameters
		//otherwise the lat/long would have been calculated. So EPSG should be nodata
		const bool comparison = ( IOUtils::checkEpsilonEquality(getEasting(), in.getEasting(), IOUtils::grid_epsilon) &&
		                          IOUtils::checkEpsilonEquality(getNorthing(), in.getNorthing(), IOUtils::grid_epsilon) &&
		                          IOUtils::checkEpsilonEquality(getAltitude(), in.getAltitude(), IOUtils::grid_epsilon) &&
		                          getEPSG()==in.getEPSG());
		return comparison;
	}
	if(grid_i!=IOUtils::nodata && grid_j!=IOUtils::nodata && grid_k!=IOUtils::nodata) {
		//only available information is grid indices
		const bool comparison = ( grid_i==in.grid_i && grid_j==in.grid_j && grid_k==in.grid_k );
		return comparison;
	}
	//every field is nodata... the objects can only be equal if both are nodata
	if(in.isNodata()==true) return true;
	else return false;
}

/**
* @brief Inequality operator that checks that lat/lon don't match
* @param[in] in Coord object to compare to
* @return true or false
*/
bool Coords::operator!=(const Coords& in) const {
	return !(*this==in);
}

Coords& Coords::operator=(const Coords& source) {
	if(this != &source) {
		altitude = source.altitude;
		latitude = source.latitude;
		longitude = source.longitude;
		easting = source.easting;
		northing = source.northing;
		ref_latitude = source.ref_latitude;
		ref_longitude = source.ref_longitude;
		distance_algo = source.distance_algo;
		coordsystem = source.coordsystem;
		coordparam = source.coordparam;
		grid_i = source.grid_i;
		grid_j = source.grid_j;
		grid_k = source.grid_k;
	}
	return *this;
}

bool Coords::isNodata() const {
	if( latitude==IOUtils::nodata && longitude==IOUtils::nodata &&
	    easting==IOUtils::nodata && northing==IOUtils::nodata &&
	    altitude==IOUtils::nodata &&
	    grid_i==IOUtils::nodata && grid_j==IOUtils::nodata && grid_k==IOUtils::nodata) {
		return true;
	}
	return false;
}

///< move the point by the specified distance (in m) along easting and northing
void Coords::moveByXY(const double& x_displacement, const double& y_displacement) {
	setXY(easting+x_displacement, northing+y_displacement, altitude, true);
}

///< move the point by the specified bearing and distance (in m)
void Coords::moveByBearing(const double& i_bearing, const double& i_distance) {
	double new_lat, new_lon;

	switch(distance_algo) {
		case GEO_COSINE:
			cosineInverse(latitude, longitude, i_distance, i_bearing, new_lat, new_lon);
			break;
		case GEO_VINCENTY:
			VincentyInverse(latitude, longitude, i_distance, i_bearing, new_lat, new_lon);
			break;
	}

	setLatLon(new_lat, new_lon, altitude, true);
}

/**
* @brief Simple merge strategy.
* If some fields of the first argument are empty, they will be filled by the matching field from the
* second argument.
* @param coord1 first Coords to merge, highest priority
* @param coord2 second Coords to merge, lowest priority
* @return new Coords object
*/
Coords Coords::merge(const Coords& coord1, const Coords& coord2) {
	Coords tmp(coord1);
	tmp.merge(coord2);
	return tmp;
}

/**
* @brief Simple merge strategy.
* If some fields of the current object are empty, they will be filled by the matching field from the
* provided argument.
* @param coord2 extra Coords to merge, lowest priority
*/
void Coords::merge(const Coords& coord2) {
	if(altitude==IOUtils::nodata) altitude=coord2.altitude;
	if(latitude==IOUtils::nodata) latitude=coord2.latitude;
	if(longitude==IOUtils::nodata) longitude=coord2.longitude;
	if(easting==IOUtils::nodata) easting=coord2.easting;
	if(northing==IOUtils::nodata) northing=coord2.northing;

	if(grid_i==IOUtils::nodata) grid_i=coord2.grid_i;
	if(grid_j==IOUtils::nodata) grid_j=coord2.grid_j;
	if(grid_k==IOUtils::nodata) grid_k=coord2.grid_k;

	if(ref_latitude==IOUtils::nodata) ref_latitude=coord2.ref_latitude;
	if(ref_longitude==IOUtils::nodata) ref_longitude=coord2.ref_longitude;

	if(coordsystem=="NULL") coordsystem=coord2.coordsystem;
	if(coordparam=="NULL") coordparam=coord2.coordparam;

	if(distance_algo==IOUtils::nodata) distance_algo=coord2.distance_algo;

	//in LOCAL projection, the check for the existence of the ref point will be done in the projection functions
	if(coordsystem=="LOCAL" && !coordparam.empty()) {
		parseLatLon(coordparam, ref_latitude, ref_longitude);
	}
	if(latitude!=IOUtils::nodata && coordsystem!="NULL") {
		convert_from_WGS84(latitude, longitude, easting, northing);
	}
	if(latitude==IOUtils::nodata && coordsystem!="NULL") {
		convert_to_WGS84(easting, northing, latitude, longitude);
	}
}

/**
* @brief Print the content of the Coords object (usefull for debugging)
* The Coords is bound by "<Coords>" and "</Coords>" on separate lines
*/
const std::string Coords::toString() const {
	std::ostringstream os;
	os << "<Coords>\n";
	os << "Altitude\t" << altitude << "\n";
	os << "Lat/Long\t" << printLatLon() << "\n";
	std::streamsize p = os.precision();
	os << "Lat/Long\t" <<std::fixed << std::setprecision(6) << "(" << getLat() << " , " << getLon() << ")" << "\n";
	os << "X/Y_coords\t" << std::fixed << std::setprecision(0) << "(" << getEasting() << " , " << getNorthing() << ")" << "\n";
	os << std::resetiosflags(std::ios_base::fixed|std::ios_base::floatfield);
	os.precision(p);
	os << "I/J_indices\t" << "(" << getGridI() << " , " << getGridJ() << ")" << "\n";
	os << "Projection\t" << coordsystem << " " << coordparam << "\n";
	os << "EPSG\t\t" << getEPSG() << "\n";
	os << "</Coords>\n";
	return os.str();
}

std::iostream& operator<<(std::iostream& os, const Coords& coord) {
	os.write(reinterpret_cast<const char*>(&coord.ref_latitude), sizeof(coord.ref_latitude));
	os.write(reinterpret_cast<const char*>(&coord.ref_longitude), sizeof(coord.ref_longitude));
	os.write(reinterpret_cast<const char*>(&coord.altitude), sizeof(coord.altitude));
	os.write(reinterpret_cast<const char*>(&coord.latitude), sizeof(coord.latitude));
	os.write(reinterpret_cast<const char*>(&coord.longitude), sizeof(coord.longitude));
	os.write(reinterpret_cast<const char*>(&coord.easting), sizeof(coord.easting));
	os.write(reinterpret_cast<const char*>(&coord.northing), sizeof(coord.northing));
	os.write(reinterpret_cast<const char*>(&coord.grid_i), sizeof(coord.grid_i));
	os.write(reinterpret_cast<const char*>(&coord.grid_j), sizeof(coord.grid_j));
	os.write(reinterpret_cast<const char*>(&coord.grid_k), sizeof(coord.grid_k));

	const size_t s_coordsystem = coord.coordsystem.size();
	os.write(reinterpret_cast<const char*>(&s_coordsystem), sizeof(size_t));
	os.write(reinterpret_cast<const char*>(&coord.coordsystem[0]), s_coordsystem*sizeof(coord.coordsystem[0]));
	const size_t s_coordparam = coord.coordparam.size();
	os.write(reinterpret_cast<const char*>(&s_coordparam), sizeof(size_t));
	os.write(reinterpret_cast<const char*>(&coord.coordparam[0]), s_coordparam*sizeof(coord.coordparam[0]));

	os.write(reinterpret_cast<const char*>(&coord.distance_algo), sizeof(coord.distance_algo));
	return os;
}

std::iostream& operator>>(std::iostream& is, Coords& coord) {
	is.read(reinterpret_cast<char*>(&coord.ref_latitude), sizeof(coord.ref_latitude));
	is.read(reinterpret_cast<char*>(&coord.ref_longitude), sizeof(coord.ref_longitude));
	is.read(reinterpret_cast<char*>(&coord.altitude), sizeof(coord.altitude));
	is.read(reinterpret_cast<char*>(&coord.latitude), sizeof(coord.latitude));
	is.read(reinterpret_cast<char*>(&coord.longitude), sizeof(coord.longitude));
	is.read(reinterpret_cast<char*>(&coord.easting), sizeof(coord.easting));
	is.read(reinterpret_cast<char*>(&coord.northing), sizeof(coord.northing));
	is.read(reinterpret_cast<char*>(&coord.grid_i), sizeof(coord.grid_i));
	is.read(reinterpret_cast<char*>(&coord.grid_j), sizeof(coord.grid_j));
	is.read(reinterpret_cast<char*>(&coord.grid_k), sizeof(coord.grid_k));

	size_t s_coordsystem, s_coordparam;
	is.read(reinterpret_cast<char*>(&s_coordsystem), sizeof(size_t));
	coord.coordsystem.resize(s_coordsystem);
	is.read(reinterpret_cast<char*>(&coord.coordsystem[0]), s_coordsystem*sizeof(coord.coordsystem[0]));
	is.read(reinterpret_cast<char*>(&s_coordparam), sizeof(size_t));
	coord.coordparam.resize(s_coordparam);
	is.read(reinterpret_cast<char*>(&coord.coordparam[0]), s_coordparam*sizeof(coord.coordparam[0]));

	is.read(reinterpret_cast<char*>(&coord.distance_algo), sizeof(coord.distance_algo));
	return is;
}

/**
* @brief Default constructor
* This constructor builds a dummy object that performs no conversions but can be used for comparison
* purpose. This is more or less the equivalent of NULL for a pointer...
*/
Coords::Coords() : ref_latitude(IOUtils::nodata), ref_longitude(IOUtils::nodata),
                   altitude(IOUtils::nodata), latitude(IOUtils::nodata), longitude(IOUtils::nodata),
                   easting(IOUtils::nodata), northing(IOUtils::nodata),
                   grid_i(IOUtils::inodata), grid_j(IOUtils::inodata), grid_k(IOUtils::inodata),
                   coordsystem("NULL"), coordparam("NULL"), distance_algo(GEO_COSINE) {}

/**
* @brief Regular constructor: usually, this is the constructor to use
* @param[in] in_coordinatesystem string identifying the coordinate system to use
* @param[in] in_parameters string giving some additional parameters for the projection (optional)
*
* See setProj() for a full description of these strings
*/
Coords::Coords(const std::string& in_coordinatesystem, const std::string& in_parameters) :
                   ref_latitude(IOUtils::nodata), ref_longitude(IOUtils::nodata),
                   altitude(IOUtils::nodata), latitude(IOUtils::nodata), longitude(IOUtils::nodata),
                   easting(IOUtils::nodata), northing(IOUtils::nodata),
                   grid_i(IOUtils::inodata), grid_j(IOUtils::inodata), grid_k(IOUtils::inodata),
                   coordsystem(in_coordinatesystem), coordparam(in_parameters), distance_algo(GEO_COSINE)
{
	setProj(in_coordinatesystem, in_parameters);
}

/**
* @brief Local projection onstructor: this constructor is only suitable for building a local projection.
* Such a projection defines easting and northing as the distance (in meters) to a reference point
* which coordinates have to be provided here.
* @param[in] in_lat_ref latitude of the reference point
* @param[in] in_long_ref longitude of the reference point
*/
Coords::Coords(const double& in_lat_ref, const double& in_long_ref) :
                   ref_latitude(in_lat_ref), ref_longitude(in_long_ref),
                   altitude(IOUtils::nodata), latitude(IOUtils::nodata), longitude(IOUtils::nodata),
                   easting(IOUtils::nodata), northing(IOUtils::nodata),
                   grid_i(IOUtils::inodata), grid_j(IOUtils::inodata), grid_k(IOUtils::inodata),
                   coordsystem("LOCAL"), coordparam(), distance_algo(GEO_COSINE)
{
	setLocalRef(in_lat_ref, in_long_ref);
	setProj("LOCAL", "");
}

Coords::Coords(const Coords& c) : ref_latitude(c.ref_latitude), ref_longitude(c.ref_longitude),
                   altitude(c.altitude), latitude(c.latitude), longitude(c.longitude),
                   easting(c.easting), northing(c.northing),
                   grid_i(c.grid_i), grid_j(c.grid_j), grid_k(c.grid_k),
                   coordsystem(c.coordsystem), coordparam(c.coordparam), distance_algo(c.distance_algo) {}

/**
* @brief Local projection onstructor: this constructor is only suitable for building a local projection.
* Such a projection defines easting and northing as the distance (in meters) to a reference point
* which coordinates have to be provided here.
* @param[in] in_coordinatesystem string identifying the coordinate system to use
* @param[in] in_parameters string giving some additional parameters for the projection (empty string if not applicable)
* @param[in] coord_spec coordinate specification
*
* The coordinate specification is given as either: "easting northing epsg" or "lat lon".
*/
Coords::Coords(const std::string& in_coordinatesystem, const std::string& in_parameters, const std::string& coord_spec)
       : ref_latitude(IOUtils::nodata), ref_longitude(IOUtils::nodata),
         altitude(IOUtils::nodata), latitude(IOUtils::nodata), longitude(IOUtils::nodata),
         easting(IOUtils::nodata), northing(IOUtils::nodata),
         grid_i(IOUtils::inodata), grid_j(IOUtils::inodata), grid_k(IOUtils::inodata),
         coordsystem(in_coordinatesystem), coordparam(in_parameters), distance_algo(GEO_COSINE)
{
	std::istringstream iss(coord_spec);
	double coord1=IOUtils::nodata, coord2=IOUtils::nodata;
	int epsg=IOUtils::inodata;

	iss >> std::skipws >> coord1;
	iss >> std::skipws >> coord2;
	iss >> std::skipws >> epsg;

	if(coord1!=IOUtils::nodata && coord2!=IOUtils::nodata && epsg!=IOUtils::inodata) {
		setEPSG(epsg);
		setXY(coord1, coord2, IOUtils::nodata);
		setProj(in_coordinatesystem, in_parameters);
	} else if(coord1!=IOUtils::nodata && coord2!=IOUtils::nodata) {
		setLatLon(coord1, coord2, IOUtils::nodata);
	}
}

/**
* @brief Returns the East coordinate in the configured projection system
* @return easting
*/
double Coords::getEasting() const {
	return easting;
}

/**
* @brief Returns the North coordinate in the configured projection system
* @return northing
*/
double Coords::getNorthing() const {
	return northing;
}

/**
* @brief Returns the Latitude in the configured projection system
* @return latitude
*/
double Coords::getLat() const {
	return latitude;
}

/**
* @brief Returns the Latitude in the configured projection system
* @return longitude
*/
double Coords::getLon() const {
	return longitude;
}

/**
* @brief Returns the Altitude. This is currently independent of the configured projection system
* @return altitude
*/
double Coords::getAltitude() const {
	return altitude;
}

/**
* @brief Returns the grid index along the X axis
* @return grid index i
*/
int Coords::getGridI() const {
	return grid_i;
}

/**
* @brief Returns the grid index along the Y axis
* @return grid index j
*/
int Coords::getGridJ() const {
	return grid_j;
}

/**
* @brief Returns the grid index along the Z axis
* @return grid index k
*/
int Coords::getGridK() const {
	return grid_k;
}

/**
* @brief Returns the projection parameters
* @param[out] proj_type projection type
* @param[out] proj_args optional arguments
*/
void Coords::getProj(std::string& proj_type, std::string& proj_args) const {
	proj_type = coordsystem;
	if(coordsystem=="LOCAL") {
		std::ostringstream dms;
		dms << "(" << decimal_to_dms(ref_latitude) << " , " << decimal_to_dms(ref_longitude) << ")";
		proj_args=dms.str();
	} else {
		proj_args = coordparam;
	}
}

/**
* @brief Print a nicely formatted lat/lon in degrees, minutes, seconds
* @return lat/lon
*/
std::string Coords::printLatLon() const {
	std::ostringstream dms;
	dms << "(" << decimal_to_dms(latitude) << " , " << decimal_to_dms(longitude) << ")";

	return dms.str();
}

/**
* @brief Set latitude and longitude
* The automatic update of the easting/northing can be turned off so that
* both lat/lon and east/north coordinates can be provided in order to thereafter check the
* coordinates by calling check().
* @param[in] in_coordinates string containing the lat/lon to read
* @param[in] in_altitude altitude to set (optional)
* @param[in] in_update should the easting/northing be updated? (default=true)
*/
void Coords::setLatLon(const std::string& in_coordinates, const double in_altitude, const bool in_update) {
	double lat, lon;
	parseLatLon(in_coordinates, lat, lon);
	setLatLon(lat, lon, in_altitude, in_update);
}

/**
* @brief Set latitude and longitude
* The automatic update of the easting/northing can be turned off so that
* both lat/lon and east/north coordinates can be provided in order to thereafter check the
* coordinates by calling check().
* @param[in] in_latitude latitude to set
* @param[in] in_longitude longitude to set
* @param[in] in_altitude altitude to set
* @param[in] in_update should the easting/northing be updated? (default=true)
*/
void Coords::setLatLon(const double in_latitude, const double in_longitude, const double in_altitude, const bool in_update) {
	latitude = in_latitude;
	longitude = in_longitude;
	if(in_altitude!=IOUtils::nodata) {
		altitude = in_altitude;
	}
	if(coordsystem!="NULL" && in_update==true) {
		convert_from_WGS84(latitude, longitude, easting, northing);
	}
	grid_i = grid_j = grid_k = IOUtils::inodata;
}

/**
* @brief Set easting and northing
* The automatic update of the latitude/longitude can be turned off so that
* both lat/lon and east/north coordinates can be provided in order to thereafter check the
* coordinates by calling check().
* @param[in] in_easting easting to set
* @param[in] in_northing northing to set
* @param[in] in_altitude altitude to set
* @param[in] in_update should the easting/northing be updated? (default=true)
*/
void Coords::setXY(const double in_easting, const double in_northing, const double in_altitude, const bool in_update) {
	easting = in_easting;
	northing = in_northing;
	if(in_altitude!=IOUtils::nodata) {
		altitude = in_altitude;
	}
	if(coordsystem!="NULL" && in_update==true) {
		convert_to_WGS84(easting, northing, latitude, longitude);
	}
	grid_i = grid_j = grid_k = IOUtils::inodata;
}

/**
* @brief Set grid indices
* This index represent the position in a cartesian grid. It can not be automatically matched with
* a set of geographic coordinates because it needs the information about the said grid.
* Therefore, the coordinate object needs to be given to a grid object that will either set (i,j) or
* (lat,lon)/(easting,northing) as well as the grid's matching altitude if it is not already set.
* Any subsequent change of either (lat,lon) or (easting,northing) will reset these indexes to IOUtils::inodata.
* By default, setting (i,j) will invalidate (ie: delete) ALL geographic coordinates for the object (since we can not
* convert from grid indices to/from geographic coordinates in the current object by lack of information).
* Finally, the given indices are <b>NOT checked for validity</b>: such check must be done
* by calling Grid2DObject::gridify or Grid3DObject::gridify .
*
* @note To make it short: <b>use this method with caution!!</b>
* @param[in] in_grid_i grid index along the X direction
* @param[in] in_grid_j grid index along the Y direction
* @param[in] in_grid_k grid index along the Z direction
* @param[in] in_invalidate should the geographic coordinates be invalidated? (default=true, this flag should be false ONLY when calling from Grid2/3DObject)
*/
void Coords::setGridIndex(const int in_grid_i, const int in_grid_j, const int in_grid_k, const bool in_invalidate) {
	grid_i = in_grid_i;
	grid_j = in_grid_j;
	grid_k = in_grid_k;
	if(in_invalidate==true) {
		latitude = IOUtils::nodata;
		longitude = IOUtils::nodata;
		easting = IOUtils::nodata;
		northing = IOUtils::nodata;
		altitude = IOUtils::nodata;
	}
}

/**
* @brief Set altitude at a given value.
* If the i,j,k indices were set, reset them to inodata,
* except if specified otherwise with in_update=false.
* @param[in] in_altitude altitude above sea level, in meters
* @param[in] in_update should the indices be (if necessary) recalculated? (default=true)
*/
void Coords::setAltitude(const double in_altitude, const bool in_update) {
	altitude = in_altitude;
	if(in_update==true) {
		grid_i = grid_j = grid_k = IOUtils::inodata;
	}
}

/**
* @brief Set projection to use
* This projection will be used for converting between lat/lon and East/North
* @param[in] in_coordinatesystem string identifying the coordinate system to use
* @param[in] in_parameters string giving some additional parameters for the projection (optional)
*
*  \anchor Coordinate_types
* The coordinate system can be any of the following:
* - CH1903 for coordinates in the Swiss Grid [ref: http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf]
* - UTM for UTM coordinates (the zone must be specified in the parameters, for example 31T) [ref: http://www.oc.nps.edu/oc2902w/maps/utmups.pdf]
* - PROJ4 for coordinate conversion relying on the Proj4 library [ref: http://trac.osgeo.org/proj/]
* - LOCAL for local coordinate system (using the horizontal and vertical distance from a reference point, see Coords::geo_distances for the available choice of distance algorithms)
*/
void Coords::setProj(const std::string& in_coordinatesystem, const std::string& in_parameters) {
	//the latitude/longitude had not been calculated, so we do it first in order to have our reference
	//before further conversions (usage scenario: giving a x,y and then converting to anyother x,y in another system
	if ((coordsystem != "NULL") && ((latitude==IOUtils::nodata) || (longitude==IOUtils::nodata))) {
		convert_to_WGS84(easting, northing, latitude, longitude);
	}

	if(in_coordinatesystem.empty()) {
		coordsystem = std::string("NULL");
	} else {
		coordsystem = in_coordinatesystem;
	}
	coordparam  = in_parameters;
	if(coordsystem=="LOCAL" && !coordparam.empty()) {
		parseLatLon(coordparam, ref_latitude, ref_longitude);
	}

	//since lat/long is our reference, we refresh x,y (only if lat/lon exist)
	if(latitude!=IOUtils::nodata && longitude!=IOUtils::nodata) {
		convert_from_WGS84(latitude, longitude, easting, northing);
	}
	//if we only had x/y but not even a coord system, we could not compute lat/long. We now do it
	if( (latitude==IOUtils::nodata || longitude==IOUtils::nodata) &&
	    (easting!=IOUtils::nodata && northing!=IOUtils::nodata) &&
	    (coordsystem != "NULL") ) {
		convert_to_WGS84(easting, northing, latitude, longitude);
	}
}

/**
* @brief Set the local projection reference coordinates
* This projection will be used for converting between lat/lon and East/North
* @param[in] in_ref_latitude latitude of the local origin
* @param[in] in_ref_longitude longitude of the local origin
*/
void Coords::setLocalRef(const double in_ref_latitude, const double in_ref_longitude) {
	if(in_ref_latitude==IOUtils::nodata || in_ref_longitude==IOUtils::nodata) {
		throw InvalidArgumentException("For LOCAL projection, please provide both reference latitude and longitude!", AT);
	}
	ref_latitude = in_ref_latitude;
	ref_longitude = in_ref_longitude;
	if(coordsystem=="LOCAL") {
		convert_from_WGS84(latitude, longitude, easting, northing);
	}
}

/**
* @brief Set the local projection reference coordinates
* This projection will be used for converting between lat/lon and East/North
* @param[in] in_coordparam string containing the (lat,lon) of the local origin
*/
void Coords::setLocalRef(const std::string in_coordparam) {
	coordparam = in_coordparam;
	parseLatLon(coordparam, ref_latitude, ref_longitude);
	if(coordsystem=="LOCAL") {
		convert_from_WGS84(latitude, longitude, easting, northing);
	}
}

/**
* @brief Set the algorithm to use to compute distances
* Various algorithm exist that offer various precision/complexity tradeoffs.
* @param[in] in_algo enum giving the algorithm to be used (see documentation for geo_distances)
*/
void Coords::setDistances(const geo_distances in_algo) {
	distance_algo = in_algo;
	if(coordsystem=="LOCAL") {
		convert_from_WGS84(latitude, longitude, easting, northing);
	}
}

/**
* @brief Check consistency of coordinates
* When both latitude/longitude and easting/northing are given, there is
* a risk of inconsistency between these two sets of coordinates.
* This method checks that enough information is available (ie: at least one set
* of coordinates is present) and if more than one is present, that it is consistent (within 5 meters)
* It throws and exception if something is not right.
*/
void Coords::check()
{
	//calculate/check coordinates if necessary
	if(coordsystem=="LOCAL" && (ref_latitude==IOUtils::nodata || ref_longitude==IOUtils::nodata)) {
		throw InvalidArgumentException("please define a reference point for LOCAL coordinate system", AT);
	}

	if(latitude==IOUtils::nodata || longitude==IOUtils::nodata) {
		if(easting==IOUtils::nodata || northing==IOUtils::nodata) {
			throw InvalidArgumentException("missing positional parameters (easting,northing) or (lat,long) for coordinate", AT);
		}
		convert_to_WGS84(easting, northing, latitude, longitude);
	} else {
		if(easting==IOUtils::nodata || northing==IOUtils::nodata) {
			convert_from_WGS84(latitude, longitude, easting, northing);
		} else {
			double tmp_lat, tmp_lon;
			convert_to_WGS84(easting, northing, tmp_lat, tmp_lon);

			if(!IOUtils::checkEpsilonEquality(latitude, tmp_lat, IOUtils::lat_epsilon) || !IOUtils::checkEpsilonEquality(longitude, tmp_lon, IOUtils::lon_epsilon)) {
				throw InvalidArgumentException("Latitude/longitude and xllcorner/yllcorner don't match for coordinate", AT);
			}
		}
	}
}

/**
* @brief Calculate the distance between two points
* @param destination destination coordinate
* @return distance in meters
*/
double Coords::distance(const Coords& destination) const {
	double dist, bearing;
	distance(destination, dist, bearing);
	return dist;
}

/**
* @brief Check if two Coords object are using the same projection
* @param target coordinate to compare to
* @return true or false
*/
bool Coords::isSameProj(const Coords& target) const {
	if(coordsystem=="LOCAL") {
		return ( target.coordsystem=="LOCAL" && ref_latitude==target.ref_latitude && ref_longitude==target.ref_longitude );
	} else {
		return ( coordsystem==target.coordsystem && coordparam==target.coordparam );
	}
}

/**
* @brief Copy the projection parameters of another Coords object
* @param source source object to copy the projection from
* @param i_update should the necessary coordinates be updated? (default=true)
*/
void Coords::copyProj(const Coords& source, const bool i_update) {
	if(!isSameProj(source)) {
		//we only do a copy if we are not already using the same projection
		if(source.coordsystem=="LOCAL") {
			coordsystem="LOCAL";
			coordparam=source.coordparam;
			ref_latitude=source.ref_latitude;
			ref_longitude=source.ref_longitude;
		} else {
			coordsystem=source.coordsystem;
			coordparam=source.coordparam;
		}

		if(i_update==true) {
			if((latitude!=IOUtils::nodata) && (longitude!=IOUtils::nodata)) {
				convert_from_WGS84(latitude, longitude, easting, northing);
			} else {
				convert_to_WGS84(easting, northing, latitude, longitude);
			}
		}
	}
}

/**
* @brief returns the epsg code of the current projection
* @return epsg code
*/
short int Coords::getEPSG() const {
	if(coordsystem=="CH1903") return 21781;
	if(coordsystem=="CH1903+") return 2056;
	if(coordsystem=="UTM") {
		//UTM Zone information
		short int zoneNumber;
		char zoneLetter;
		parseUTMZone(coordparam, zoneLetter, zoneNumber);
		if(zoneLetter >= 'N') {
			//northern hemisphere. We KNOW it will fit in a short int
			return (static_cast<short int>(32600+zoneNumber));
		} else {
			//southern hemisphere. We KNOW it will fit in a short int
			return (static_cast<short int>(32700+zoneNumber));
		}
	}
	if(coordsystem=="UPS") {
		//UPS Zone
		if(coordparam == "S") {
			//southern hemisphere
			return (32761);
		} else {
			//northern hemisphere
			return (32661);
		}
	}
	if(coordsystem=="PROJ4") {
		const int tmp = atoi(coordparam.c_str());
		if(tmp<0 || tmp>32767) {
			std::ostringstream ss;
			ss << "Invalid EPSG code argument: " << tmp << ". It should be between 0 and 32767! (please check EPSG registry)";
			throw InvalidArgumentException(ss.str(), AT);
		}
		return static_cast<short>(tmp);
	}

	//all others have no associated EPSG code
	return -1;
}

/**
* @brief set the current projection to a given EPSG-defined projection
* @param epsg epsg code
*/
void Coords::setEPSG(const int epsg) {
//TODO: get rid of the zone letter. This is not part of the standard and redundant (and messy)
	bool found=false;
	std::string coord_sys, coord_param;

	if(epsg<0 || epsg>32767) {
		std::ostringstream ss;
		ss << "Invalid epsg code " << epsg << " (it should be between 0 and 32767)!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	if(!found && (epsg==21781)) {
		coord_sys.assign("CH1903");
		coord_param.clear();
		found=true;
	}
	if(!found && (epsg==2056)) {
		coord_sys.assign("CH1903+");
		coord_param.clear();
		found=true;
	}
	if(!found && (epsg>=32601) && (epsg<=32660)) {
		//northern hemisphere
		coord_sys.assign("UTM");
		const int zoneNumber = epsg-32600;
		std::ostringstream osstream;
		osstream << zoneNumber << "N";
		coord_param=osstream.str();
		found=true;
	}
	if(!found && (epsg>=32701) && (epsg<=32760)) {
		//southern hemisphere
		coord_sys.assign("UTM");
		const int zoneNumber = epsg-32700;
		std::ostringstream osstream;
		osstream << zoneNumber << "M";
		coord_param=osstream.str();
		found=true;
	}
	if(!found && (epsg==32761 || epsg==32661)) {
		coord_sys.assign("UPS");
		coord_param = (epsg==32761)? "S" : "N";
		found=true;
	}
	if(!found) {
		//anything else has to be processed by proj4
		coord_sys.assign("PROJ4");
		std::ostringstream osstream;
		osstream << epsg;
		coord_param=osstream.str();
	}
	setProj(coord_sys, coord_param);
}

/////////////////////////////////////////////////////private methods
/**
* @brief Method converting towards WGS84
* @param[in] i_easting easting of the coordinate to convert
* @param[in] i_northing northing of the coordinate to convert
* @param[out] o_latitude converted latitude
* @param[out] o_longitude converted longitude
*/
void Coords::convert_to_WGS84(double i_easting, double i_northing, double& o_latitude, double& o_longitude) const
{
	if((i_easting!=IOUtils::nodata) && (i_northing!=IOUtils::nodata)) {
		if(coordsystem=="UTM") UTM_to_WGS84(i_easting, i_northing, o_latitude, o_longitude);
		else if(coordsystem=="UPS") UPS_to_WGS84(i_easting, i_northing, o_latitude, o_longitude);
		else if(coordsystem=="CH1903") CH1903_to_WGS84(i_easting, i_northing, o_latitude, o_longitude);
		else if(coordsystem=="CH1903+") CH1903_to_WGS84(i_easting-2e6, i_northing-1e6, o_latitude, o_longitude);
		else if(coordsystem=="LOCAL") local_to_WGS84(i_easting, i_northing, o_latitude, o_longitude);
		else if(coordsystem=="PROJ4") PROJ4_to_WGS84(i_easting, i_northing, o_latitude, o_longitude);
		else if(coordsystem=="NULL") NULL_to_WGS84(i_easting, i_northing, o_latitude, o_longitude);
		else throw UnknownValueException("Unknown coordinate system \""+coordsystem+"\"", AT);
	} else {
		o_latitude = IOUtils::nodata;
		o_longitude = IOUtils::nodata;
	}
}

/**
* @brief Method converting towards WGS84
* @param[in] i_latitude latitude of the coordinate to convert
* @param[in] i_longitude longitude of the coordinate to convert
* @param[out] o_easting converted easting
* @param[out] o_northing converted northing
*/
void Coords::convert_from_WGS84(double i_latitude, double i_longitude, double& o_easting, double& o_northing) const
{
	if((i_latitude!=IOUtils::nodata) && (i_longitude!=IOUtils::nodata)) {
		if(coordsystem=="UTM") WGS84_to_UTM(i_latitude, i_longitude, o_easting, o_northing);
		else if(coordsystem=="UPS") WGS84_to_UPS(i_latitude, i_longitude, o_easting, o_northing);
		else if(coordsystem=="CH1903") WGS84_to_CH1903(i_latitude, i_longitude, o_easting, o_northing);
		else if(coordsystem=="CH1903+") {
			WGS84_to_CH1903(i_latitude, i_longitude, o_easting, o_northing);
			o_easting += 2e6;
			o_northing += 1e6;
		}
		else if(coordsystem=="LOCAL") WGS84_to_local(i_latitude, i_longitude, o_easting, o_northing);
		else if(coordsystem=="PROJ4") WGS84_to_PROJ4(i_latitude, i_longitude, o_easting, o_northing);
		else if(coordsystem=="NULL") WGS84_to_NULL(i_latitude, i_longitude, o_easting, o_northing);
		else throw UnknownValueException("Unknown coordinate system \""+coordsystem+"\"", AT);
	} else {
		o_easting = IOUtils::nodata;
		o_northing = IOUtils::nodata;
	}
}

/**
* @brief Parse a latitude or longitude
* It can be formatted as any of the following examples:
* - 46&deg; 48' 03" (with or without spaces, decimal or integer numbers)
* - 46d 48' 03" (with or without spaces, decimal or integer numbers)
* - 46 48' 03" (with spaces, decimal or integer numbers)
* - 48&deg; 48.02'(with or without spaces, decimal or integer numbers)
* - 46d 48.02'(with or without spaces, decimal or integer numbers)
* - 46 48.02'(with spaces, decimal or integer numbers)
* - 46.802&deg;
* - 46.802d
* - 46.802
* @param[in] dms string containing the coordinate
* @return coordinate in decimal
*/
double Coords::dms_to_decimal(const std::string& dms) {
	double d=IOUtils::nodata, m=IOUtils::nodata, s=IOUtils::nodata, decimal=IOUtils::nodata;

	if 	((sscanf(dms.c_str(), "%lf°%lf'%lf\"", &d, &m ,&s) < 3) &&
		(sscanf(dms.c_str(), "%lf° %lf' %lf\"", &d, &m ,&s) < 3) &&
		(sscanf(dms.c_str(), "%lfd%lf'%lf\"", &d, &m ,&s) < 3) &&
		(sscanf(dms.c_str(), "%lfd %lf' %lf\"", &d, &m ,&s) < 3) &&
		(sscanf(dms.c_str(), "%lf %lf' %lf\"", &d, &m ,&s) < 3) &&
		(sscanf(dms.c_str(), "%lf° %lf'", &d, &m) < 2) &&
		(sscanf(dms.c_str(), "%lf°%lf'", &d, &m) < 2) &&
		(sscanf(dms.c_str(), "%lfd %lf'", &d, &m) < 2) &&
		(sscanf(dms.c_str(), "%lfd%lf'", &d, &m) < 2) &&
		(sscanf(dms.c_str(), "%lf %lf'", &d, &m) < 2) &&
		(sscanf(dms.c_str(), "%lf°", &d) < 1) &&
		(sscanf(dms.c_str(), "%lfd", &d) < 1) &&
		(sscanf(dms.c_str(), "%lf", &d) < 1)) {
			throw InvalidFormatException("Can not parse given latitude or longitude: "+dms,AT);
	}

	decimal = fabs(d);
	if(m!=IOUtils::nodata) {
		decimal += m/60.;
	}
	if(s!=IOUtils::nodata) {
		decimal += s/3600.;
	}

	if(d<0.) return (-decimal);
	else return decimal;
}

/**
* @brief Parse a latitude-longitude pair
* It can be formatted as any of the following examples:
* - lat lon (without any spaces in the latitude or longitude string)
* - lat/lon
* - (lat;lon)
* - (lat,lon)
* @param[in] coordinates string containing the coordinates
* @param[out] lat parsed latitude
* @param[out] lon parsed longitude
*/
void Coords::parseLatLon(const std::string& coordinates, double&
lat, double& lon)
{
	const size_t len=64;
	char lat_str[len]=""; //each string must be able to accomodate the whole length to avoid buffer overflow
	char lon_str[len]="";

	if(coordinates.size()>=len) {
		throw InvalidFormatException("Given lat/lon string is too long! ",AT);
	}

	if     ((sscanf(coordinates.c_str(), "%[0-9.,°d'\"-] %[0-9.,°d'\"-]", lat_str, lon_str) < 2) &&
		(sscanf(coordinates.c_str(), "%[0-9.,°d'\"- ]/%[0-9.,°d'\"- ]", lat_str, lon_str) < 2) &&
		(sscanf(coordinates.c_str(), "(%[0-9.,°d'\"- ];%[0-9.,°d'\"- ])", lat_str, lon_str) < 2) &&
		(sscanf(coordinates.c_str(), "(%[0-9.°d'\"- ],%[0-9.°d'\"- ])", lat_str, lon_str) < 2)) {
			throw InvalidFormatException("Can not parse given lat/lon: "+coordinates,AT);
	}

	lat = dms_to_decimal(std::string(lat_str));
	lon = dms_to_decimal(std::string(lon_str));
}

/**
* @brief Converts a decimal latitude or longitude to degrees, minutes, seconds
* It formats its arguments as in the following example: 46&deg;48'03"
* @param[in] decimal decimal coordinate to convert
* @return string containing the formatted coordinate
*/
std::string Coords::decimal_to_dms(const double& decimal) {
	std::ostringstream dms;
	const double abs_dec = fabs(decimal);
	const int d = (int)floor(abs_dec);
	const int m = (int)floor( (abs_dec - (double)d)*60. );
	const double s = 3600.*(abs_dec - (double)d) - 60.*(double)m;

	if(decimal<0.)
		dms << "-";
	dms << d << "°" << m << "'" << fixed << setprecision(6) << s << "\"";
	return dms.str();
}

/**
* @brief Lenght of one degree of latitude
* This returns the lenght in meters of one degree of latitude around the given latitude
* (ie: latitude-.5, latitude+.5). See https://en.wikipedia.org/wiki/Latitude#The_length_of_a_degree_of_latitude
* @param[in] latitude latitude where to perform the computation
* @return lenght of one degree of latitude
*/
double Coords::lat_degree_lenght(const double& latitude) {
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double e2 = (a*a-b*b) / (a*a);	//ellispoid eccentricity, squared

	const double degree_length = (Cst::PI*a*(1.-e2)) / ( 180.*pow(1.-e2*Optim::pow2(sin(latitude*Cst::to_rad)), 1.5) );
	return fabs( degree_length );
}

/**
* @brief Lenght of one degree of longitude
* This returns the lenght in meters of one degree of longitude around the given latitude
* (ie: latitude-.5, latitude+.5). See https://en.wikipedia.org/wiki/Latitude#The_length_of_a_degree_of_latitude
* @param[in] latitude latitude where to perform the computation
* @return lenght of one degree of longitude
*/
double Coords::lon_degree_lenght(const double& latitude) {
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double e2 = (a*a-b*b) / (a*a);	//ellispoid eccentricity, squared

	const double degree_length = (Cst::PI*a*cos(latitude*Cst::to_rad)) / ( 180.*sqrt(1.-e2*Optim::pow2(sin(latitude*Cst::to_rad))) );
	return fabs( degree_length );
}

/**
* @brief Convert rotated lat/lon into geographic lat/lon
* Rotated coordinates are created by moving the North pole by a given offset along latitude and longitude.
* The goal is to put the equator through the center of the domain of interest, so a lat/lon grid can easily be
* approximated by a tangential cartesian coordinate system.
* (see http://www.cosmo-model.org/content/model/documentation/core/default.htm, part I, chapter 3.3 for more)
* @note To convert South Pole coordinates to North Pole coordinates, multiply the latitude by -1 and add 180 to the longitude.
* @param[in] lat_N North pole latitude offset
* @param[in] lon_N North pole longitude offset
* @param[in] lat_rot rotated latitude
* @param[in] lon_rot rotated longitude
* @param[out] lat_true geographic latitude
* @param[out] lon_true geographic longitude
*/
void Coords::rotatedToTrueLatLon(const double& lat_N, const double& lon_N, const double& lat_rot, const double& lon_rot, double &lat_true, double &lon_true)
{
	if(lat_N==IOUtils::nodata || lon_N==IOUtils::nodata || lat_rot==IOUtils::nodata || lon_rot==IOUtils::nodata) {
		lat_true = IOUtils::nodata;
		lon_true = IOUtils::nodata;
		return;
	}

	const double lat_pole_rad = lat_N*Cst::to_rad;
	const double lat_rot_rad = lat_rot*Cst::to_rad;
	const double lon_rot_rad = (fmod(lon_rot+180., 360.)-180.)*Cst::to_rad; //putting lon_rot_rad in [-PI; PI]

	lat_true = asin( sin(lat_rot_rad)*sin(lat_pole_rad) + cos(lat_rot_rad)*cos(lon_rot_rad)*cos(lat_pole_rad) ) * Cst::to_deg;
	lon_true = atan2( cos(lat_rot_rad)*sin(lon_rot_rad) , (sin(lat_pole_rad)*cos(lat_rot_rad)*cos(lon_rot_rad) - sin(lat_rot_rad)*cos(lat_pole_rad)) )*Cst::to_deg + lon_N;
	lon_true -= 180.; //HACK
	lon_true = fmod(lon_true+180., 360.)-180.; //putting lon_rot_rad in [-180; 180]
}

/**
* @brief Convert geographic lat/lon into rotated lat/lon
* Rotated coordinates are created by moving the North pole by a given offset along latitude and longitude.
* The goal is to put the equator through the center of the domain of interest, so a lat/lon grid can easily be
* approximated by a tangential cartesian coordinate system.
* (see http://www.cosmo-model.org/content/model/documentation/core/default.htm, part I, chapter 3.3 for more)
* @note To convert South Pole coordinates to North Pole coordinates, multiply the latitude by -1 and add 180 to the longitude.
* @param[in] lat_N North pole latitude offset
* @param[in] lon_N North pole longitude offset
* @param[in] lat_true geographic latitude
* @param[in] lon_true geographic longitude
* @param[out] lat_rot rotated latitude
* @param[out] lon_rot rotated longitude
*/
void Coords::trueLatLonToRotated(const double& lat_N, const double& lon_N, const double& lat_true, const double& lon_true, double &lat_rot, double &lon_rot)
{
	if(lat_N==IOUtils::nodata || lon_N==IOUtils::nodata || lat_true==IOUtils::nodata || lon_true==IOUtils::nodata) {
		lat_rot = IOUtils::nodata;
		lon_rot = IOUtils::nodata;
		return;
	}

	const double lon_norm = (fmod(lon_true+180., 360.)-180.)*Cst::to_rad; //putting lon_true in [-PI; PI]
	const double lat_true_rad = lat_true*Cst::to_rad;
	const double lat_pole_rad = lat_N*Cst::to_rad;
	const double lon_pole_rad = lon_N*Cst::to_rad;
	const double delta_lon_rad = lon_norm-lon_pole_rad;

	lat_rot = asin( sin(lat_true_rad)*sin(lat_pole_rad) + cos(lat_true_rad)*cos(lat_pole_rad)*cos(delta_lon_rad) ) * Cst::to_deg;
	lon_rot = atan2( cos(lat_true_rad)*sin(delta_lon_rad) , (cos(lat_true_rad)*sin(lat_pole_rad)*cos(delta_lon_rad) - sin(lat_true_rad)*cos(lat_pole_rad)) ) * Cst::to_deg;
	lon_rot += 180.; //HACK
	lon_rot = fmod(lon_rot+180., 360.)-180.; //putting lon_rot_rad in [-180; 180]
}

/**
* @brief Coordinate conversion: from WGS84 Lat/Long to Swiss grid
* See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf for more.
* @param[in] lat_in Decimal Latitude
* @param[in] long_in Decimal Longitude
* @param[out] east_out easting coordinate (Swiss system)
* @param[out] north_out northing coordinate (Swiss system)
*/
void Coords::WGS84_to_CH1903(double lat_in, double long_in, double& east_out, double& north_out) const
{
	//converts WGS84 coordinates (lat,long) to the Swiss coordinates. See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf
	//The elevation is supposed to be above sea level, so it does not require any conversion
	//lat and long must be decimal (and they will be converted to seconds)
	const double phi_p = (lat_in*3600. - 169028.66) / 10000.;
	const double lambda_p = (long_in*3600. - 26782.5) / 10000.;

	east_out = 600072.37
		+ 211455.93	* lambda_p
		- 10938.51	* lambda_p * phi_p
		- 0.36		* lambda_p * (phi_p*phi_p)
		- 44.54		* (lambda_p*lambda_p*lambda_p);

	north_out = 200147.07
		+ 308807.95	* phi_p
		+ 3745.25	* (lambda_p*lambda_p)
		+ 76.63		* (phi_p*phi_p)
		- 194.56	* (lambda_p*lambda_p) * phi_p
		+ 119.79	* (phi_p*phi_p*phi_p);

	/*// if necessary for the elevation, uncomment this block
	h_out = h_in - 49.55
		+ 2.73		* lambda_p
		+ 6.94		* phi_p;
	*/
}

/**
* @brief Coordinate conversion: from Swiss grid to WGS84 Lat/Long
* See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf for more.
* @param[in] east_in easting coordinate (Swiss system)
* @param[in] north_in northing coordinate (Swiss system)
* @param[out] lat_out Decimal Latitude
* @param[out] long_out Decimal Longitude
*/
void Coords::CH1903_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
{
	//converts Swiss coordinates to WGS84 coordinates (lat,long). See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf
	//The elevation is supposed to be above sea level, so it does not require any conversion
	//lat and long are decimal
	const double y_p = (east_in - 600000.) / 1000000.;
	const double x_p = (north_in - 200000.) / 1000000.;

	const double lambda_p = 2.6779094
		+ 4.728982	* y_p
		+ 0.791484	* y_p * x_p
		+ 0.1306	* y_p * (x_p*x_p)
		- 0.0436	* (y_p*y_p*y_p);

	const double phi_p = 16.9023892
		+ 3.238272	* x_p
		- 0.270978	* (y_p*y_p)
		- 0.002528	* (x_p*x_p)
		- 0.0447	* (y_p*y_p) * x_p
		- 0.0140	* (x_p*x_p*x_p);

	lat_out = phi_p * 100./36.;
	long_out = lambda_p * 100./36.;

	/*// if necessary for the elevation, uncomment this block
	h_out = h_in + 49.55
		- 12.60		* y_p
		- 22.64		* x_p;
	*/
}

int Coords::getUTMZone(const double i_latitude, const double i_longitude, std::string& zone_out) const
{//This routine determines the correct UTM letter designator for the given latitude
//UTM limits its coverage to [80S , 84N], outside of this, returns Y/Z/A/B for the zone

	//computing zone number, assuming longitude in [-180. ; 180[
	int ZoneNumber = int((i_longitude + 180.)/6.) + 1;

	// Special zones for Scandinavia
	if( i_latitude >= 72.0 && i_latitude < 84.0 ) {
		if(      i_longitude >= 0.0  && i_longitude <  9.0 ) ZoneNumber = 31;
		else if( i_longitude >= 9.0  && i_longitude < 21.0 ) ZoneNumber = 33;
		else if( i_longitude >= 21.0 && i_longitude < 33.0 ) ZoneNumber = 35;
		else if( i_longitude >= 33.0 && i_longitude < 42.0 ) ZoneNumber = 37;
	 }
	if( latitude >= 56.0 && i_latitude < 64.0 && i_longitude >= 3.0 && i_longitude < 12.0 ) {
		ZoneNumber = 32;
	}

	//getting zone letter
	char zoneLetter='Z';
	if     ((0 >= i_longitude) && (i_latitude >  84)) zoneLetter = 'Y';
	else if((0 <  i_longitude) && (i_latitude >  84)) zoneLetter = 'Z';
	else if((84 >= i_latitude) && (i_latitude >= 72)) zoneLetter = 'X';
	else if((72 > i_latitude) && (i_latitude >= 64)) zoneLetter = 'W';
	else if((64 > i_latitude) && (i_latitude >= 56)) zoneLetter = 'V';
	else if((56 > i_latitude) && (i_latitude >= 48)) zoneLetter = 'U';
	else if((48 > i_latitude) && (i_latitude >= 40)) zoneLetter = 'T';
	else if((40 > i_latitude) && (i_latitude >= 32)) zoneLetter = 'S';
	else if((32 > i_latitude) && (i_latitude >= 24)) zoneLetter = 'R';
	else if((24 > i_latitude) && (i_latitude >= 16)) zoneLetter = 'Q';
	else if((16 > i_latitude) && (i_latitude >= 8)) zoneLetter = 'P';
	else if(( 8 > i_latitude) && (i_latitude >= 0)) zoneLetter = 'N'; //first zone of Northern hemisphere
	else if(( 0 > i_latitude) && (i_latitude >= -8)) zoneLetter = 'M'; //first zone of Southern hemisphere
	else if((-8 > i_latitude) && (i_latitude >= -16)) zoneLetter = 'L';
	else if((-16 > i_latitude) && (i_latitude >= -24)) zoneLetter = 'K';
	else if((-24 > i_latitude) && (i_latitude >= -32)) zoneLetter = 'J';
	else if((-32 > i_latitude) && (i_latitude >= -40)) zoneLetter = 'H';
	else if((-40 > i_latitude) && (i_latitude >= -48)) zoneLetter = 'G';
	else if((-48 > i_latitude) && (i_latitude >= -56)) zoneLetter = 'F';
	else if((-56 > i_latitude) && (i_latitude >= -64)) zoneLetter = 'E';
	else if((-64 > i_latitude) && (i_latitude >= -72)) zoneLetter = 'D';
	else if((-72 > i_latitude) && (i_latitude >= -80)) zoneLetter = 'C';
	else if((0 >=  i_longitude) && (i_latitude <= -80)) zoneLetter = 'A';
	else if((0 <   i_longitude) && (i_latitude <= -80)) zoneLetter = 'B';

	std::ostringstream zone;
	zone << ZoneNumber << zoneLetter;
	zone_out = zone.str();

	return ZoneNumber;
}

/**
* @brief Coordinate conversion: from WGS84 Lat/Long to UTM grid
* See http://www.oc.nps.edu/oc2902w/maps/utmups.pdf for more.
* @param[in] lat_in Decimal Latitude
* @param[in] long_in Decimal Longitude
* @param[out] east_out easting coordinate (Swiss system)
* @param[out] north_out northing coordinate (Swiss system)
*/
void Coords::WGS84_to_UTM(double lat_in, double long_in, double& east_out, double& north_out) const
{//converts WGS84 coordinates (lat,long) to UTM coordinates.
//See USGS Bulletin 1532 or http://earth-info.nga.mil/GandG/publications/tm8358.2/TM8358_2.pdf
//also http://www.uwgb.edu/dutchs/usefuldata/UTMFormulas.HTM
//also http://www.oc.nps.edu/oc2902w/maps/utmups.pdf or Chuck Gantz (http://www.gpsy.com/gpsinfo/geotoutm/)
	//Geometric constants
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double e2 = (a*a-b*b) / (a*a);	//ellispoid eccentricity, squared
	const double eP2 = e2 / (1.-e2);	//second ellispoid eccentricity, squared (=(a²-b²)/b²)
	const double k0 = 0.9996;	//scale factor for the projection

	//getting posistion parameters
	std::string zone;
	long_in = fmod(long_in+360.+180., 360.) - 180.; //normalized to [-180. ; 180.[
	const double Long = long_in * Cst::to_rad;
	const double Lat = lat_in * Cst::to_rad;
	int zoneNumber = getUTMZone(lat_in, long_in, zone);
	short int in_zoneNumber;
	char in_zoneLetter;
	parseUTMZone(coordparam, in_zoneLetter, in_zoneNumber);
	if(in_zoneNumber!=zoneNumber) {
		std::cerr << "[W] requested UTM zone is not appropriate for the given coordinates. Normally, It should be zone ";
		std::cerr << zoneNumber << "\n";
		zoneNumber = in_zoneNumber;
	}
	const double long0 = (double)((zoneNumber - 1)*6 - 180 + 3) * Cst::to_rad; //+3 puts origin in middle of zone

	//Geometrical parameters
	const double nu = a / sqrt(1.-e2*Optim::pow2(sin(Lat))); //radius of curvature of the earth perpendicular to the meridian plane
	const double p = (Long-long0);

	//calculating first the coefficients of the series, then the Meridional Arc M itself
	const double n = (a-b)/(a+b);
	const double n2=n*n, n3=n*n*n, n4=n*n*n*n, n5=n*n*n*n*n, n6=n*n*n*n*n*n;
	const double A = a           * (1. - n + 5./4.*(n2 - n3) + 81./64.*(n4 - n5));
	const double B = (3./2.*a)   * (n - n2 + 7./8.*(n3 - n4) + 55./64.*(n5 - n6));
	const double C = (15./16.*a) * (n2 - n3 + 3./4.*(n4 - n5));
	const double D = (35./48.*a) * (n3 - n4 + 11./16.*(n5 - n6));
	const double E = (315./512.*a) * (n4 - n5); //correction of ~0.03mm
	const double M = A*Lat - B*sin(2.*Lat) + C*sin(4.*Lat) - D*sin(6.*Lat) + E*sin(8.*Lat);

	//calculating the coefficients for the series
	const double K1 = M*k0;
	const double K2 = 1./4.*k0*nu*sin(2.*Lat);
	const double K3 = (k0*nu*sin(Lat)*Optim::pow3(cos(Lat))*1./24.) * (5. - Optim::pow2(tan(Lat)) + 9.*eP2*Optim::pow2(cos(Lat)) + 4.*eP2*eP2*Optim::pow4(cos(Lat)));
	const double K4 = k0*nu*cos(Lat);
	const double K5 = (k0*nu*Optim::pow3(cos(Lat))*1./6.) * (1. - Optim::pow2(tan(Lat)) + eP2*Optim::pow2(cos(Lat)));

	north_out = K1 + K2*p*p + K3*p*p*p*p;
	east_out = K4*p + K5*p*p*p + 500000.0;

	if(Lat < 0)
		north_out += 10000000.0; //offset for southern hemisphere
}

/**
* @brief Coordinate conversion: from UTM grid to WGS84 Lat/Long
* See http://www.oc.nps.edu/oc2902w/maps/utmups.pdf for more.
* @param[in] east_in easting coordinate (UTM)
* @param[in] north_in northing coordinate (UTM)
* @param[out] lat_out Decimal Latitude
* @param[out] long_out Decimal Longitude
*/
void Coords::UTM_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
{//converts UTM coordinates to WGS84 coordinates (lat,long).
//See USGS Bulletin 1532 or http://earth-info.nga.mil/GandG/publications/tm8358.2/TM8358_2.pdf
//also http://www.uwgb.edu/dutchs/usefuldata/UTMFormulas.HTM
//also http://www.oc.nps.edu/oc2902w/maps/utmups.pdf or Chuck Gantz (http://www.gpsy.com/gpsinfo/geotoutm/)
	//Geometric constants
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double e2 = (a*a-b*b) / (a*a);	//ellispoid eccentricity, squared
	const double eP2 = e2 / (1.-e2);	//second ellispoid eccentricity, squared (=(a²-b²)/b²)
	const double k0 = 0.9996;		//scale factor for the projection

	//UTM Zone information
	short int zoneNumber;
	char zoneLetter;
	parseUTMZone(coordparam, zoneLetter, zoneNumber);

	//set reference parameters: central meridian of the zone, true northing and easting
	//please note that the special zones still use the reference meridian as given by their zone number (ie: even if it might not be central anymore)
	const int long0 = ((int)zoneNumber - 1)*6 - 180 + 3;  //+3 puts origin in "middle" of zone as required for the projection meridian (might not be the middle for special zones)
	if(zoneLetter<='M') {
		north_in -= 10000000.0; //offset used for southern hemisphere
	}
	east_in -= 500000.0; //longitude offset: x coordinate is relative to central meridian

	//calculating footprint latitude fp (it should be done using a few iterations)
	const double arc = north_in/k0; //Meridional arc
	const double mu = arc / (a*(1.-e2/4.-3.*e2*e2/64.-5.*e2*e2*e2/256.));
	const double e1 = (1.-b/a) / (1.+b/a); //simplification of [1 - (1 - e2)1/2]/[1 + (1 - e2)1/2]
	const double J1 = (3./2.*e1 - 27./32.*e1*e1*e1);
	const double J2 = (21./16.*e1*e1 - 55./32.*e1*e1*e1*e1);
	const double J3 = (151./96.*e1*e1*e1);
	const double J4 = (1097./512.*e1*e1*e1*e1);
	const double fp = mu + J1*sin(2.*mu) + J2*sin(4.*mu) + J3*sin(6.*mu) + J4*sin(8.*mu);

	//calculating the parameters
	const double C1 = eP2 * Optim::pow2(cos(fp));
	const double T1 = Optim::pow2( tan(fp) );
	const double R1 = a*(1.-e2) / pow((1.-e2*Optim::pow2(sin(fp))), 1.5);
	const double N1 = a / sqrt(1.-e2*Optim::pow2(sin(fp)));
	const double D = east_in / (N1*k0);

	//calculating the coefficients of the series for latitude and longitude
	const double Q1 = N1*tan(fp)/R1;
	const double Q2 = 0.5*D*D;
	const double Q3 = (5. + 3.*T1 + 10.*C1 - 4.*C1*C1 - 9.*eP2) * 1./24.*D*D*D*D;
	const double Q4 = (61. + 90.*T1 + 298.*C1 + 45.*T1*T1 - 3.*C1*C1 - 252.*eP2) * 1./720.*D*D*D*D*D*D;
	//const double Q4extra = (1385. + 3633.*T1 + 4095.*T1*T1 + 1575.*T1*T1*T1) * 1./40320.*D*D*D*D*D*D*D*D;

	const double Q5 = D;
	const double Q6 = (1. + 2.*T1 + C1) * 1./6.*D*D*D;
	const double Q7 = (5. - 2.*C1 + 28.*T1 - 3.*C1*C1 + 8.*eP2 + 24.*T1*T1) * 1./120.*D*D*D*D*D;
	//const double Q7extra = (61. + 662.*T1 + 1320.*T1*T1 +720.*T1*T1*T1) * 1./5040.*D*D*D*D*D*D*D;

	lat_out = (fp - Q1 * (Q2 - Q3 + Q4 /*+Q4extra*/))*Cst::to_deg;
	long_out = (double)long0 + ((Q5 - Q6 + Q7 /*-Q7extra*/)/cos(fp))*Cst::to_deg;
}

/**
* @brief Coordinate conversion: from WGS84 Lat/Long to Universal Polar Stereographic grid
* see J. Hager, J. Behensky, B. Drew, <i>THE UNIVERSAL GRIDS: Universal Transverse Mercator (UTM) and Universal Polar Stereographic (UPS)</i>, 1989, Defense Mapping Agency, DMATM 8358.2.
* This is valid above latitudes 84N or above 80S.
* @param[in] lat_in Decimal Latitude
* @param[in] long_in Decimal Longitude
* @param[out] east_out easting coordinate (Swiss system)
* @param[out] north_out northing coordinate (Swiss system)
*/
void Coords::WGS84_to_UPS(double lat_in, double long_in, double& east_out, double& north_out) const
{
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double e2 = (a*a-b*b) / (a*a);	//ellispoid eccentricity, squared
	const double k0 = 0.994;		//scale factor for the projection
	const double FN = 2000000.;		//false northing
	const double FE = 2000000.;		//false easting

	const double e = sqrt(e2);
	const double C0 = 2.*a / sqrt(1.-e2) * pow( (1.-e)/(1.+e) , e/2.);
	const double tan_Zz = pow( (1.+e*sin(lat_in*Cst::to_rad))/(1.-e*sin(lat_in*Cst::to_rad)) , e/2. ) * tan( Cst::PI/4. - lat_in*Cst::to_rad/2.);
	const double R = k0*C0*tan_Zz;

	if(coordparam=="N") {
		north_out = FN - R*cos(long_in*Cst::to_rad);
	} else if(coordparam=="S") {
		north_out = FN + R*cos(long_in*Cst::to_rad);
	} else {
		throw InvalidFormatException("Invalid UPS zone: "+coordparam+". It should be either N or S.",AT);
	}
	east_out = FE + R*sin(long_in*Cst::to_rad);
}

/**
* @brief Coordinate conversion: from Universal Polar Stereographic grid to WGS84 Lat/Long
* see J. Hager, J. Behensky, B. Drew, <i>THE UNIVERSAL GRIDS: Universal Transverse Mercator (UTM) and Universal Polar Stereographic (UPS)</i>, 1989, Defense Mapping Agency, DMATM 8358.2.
* This is valid above latitudes 84N or above 80S.
* @param[in] east_in easting coordinate (UTM)
* @param[in] north_in northing coordinate (UTM)
* @param[out] lat_out Decimal Latitude
* @param[out] long_out Decimal Longitude
*/
void Coords::UPS_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
{
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double e2 = (a*a-b*b) / (a*a);	//ellispoid eccentricity, squared
	const double k0 = 0.994;		//scale factor for the projection
	const double FN = 2000000.;		//false northing
	const double FE = 2000000.;		//false easting

	const double Delta_N = north_in - FN;
	const double Delta_E = east_in - FE;

	//which hemisphere are we in?
	bool northern_hemisphere;
	if(coordparam=="N") {
		northern_hemisphere = true;
	} else if(coordparam=="S") {
		northern_hemisphere = false;
	} else {
		throw InvalidFormatException("Invalid UPS zone: "+coordparam+". It should be either N or S.",AT);
	}

	//computing longitude
	if(Delta_N==0.) {
		if(Delta_E==0.) {
			long_out = 0.;
			if(northern_hemisphere) {
				lat_out = 90.;
			} else {
				lat_out = -90.;
			}
			return;
		} else {
			if(Delta_E>0.) long_out = 90.;
			else long_out = -90.;
		}
	} else {
		if(northern_hemisphere) {
			long_out = atan2( Delta_E , -Delta_N) * Cst::to_deg;
		} else {
			long_out = atan2( Delta_E , Delta_N) * Cst::to_deg;
		}
	}

	//computing latitude
	double R;
	if(Delta_N==0.) {
		R = fabs(Delta_E);
	} else if(Delta_E==0.) {
		R = fabs(Delta_N);
	} else {
		if(Delta_N>Delta_E) R = fabs( Delta_N / cos(long_out*Cst::to_rad));
		else R = fabs( Delta_E / sin(long_out*Cst::to_rad) );
	}
	const double e = sqrt(e2);
	const double C0 = 2.*a / sqrt(1.-e2) * pow( (1.-e)/(1.+e) , e/2.);
	const double tan_Zz = R / (k0*C0); //isometric colatitude
	const double chi = Cst::PI/2. - atan( tan_Zz )*2.;

	const double A = e2/2. + 5.*e2*e2/24. + e2*e2*e2/12. + 13.*e2*e2*e2*e2/360.;
	const double B = 7.*e2*e2/48. + 29.*e2*e2*e2/240. + 811.*e2*e2*e2*e2/11520.;
	const double C = 7.*e2*e2*e2/120. + 81.*e2*e2*e2*e2/1120.;
	const double D = 4279.*e2*e2*e2*e2/161280.;
	lat_out = chi + A*sin(2.*chi) + B*sin(4.*chi) + C*sin(6.*chi) + D*sin(8.*chi);
	lat_out *= Cst::to_deg;

	if(!northern_hemisphere) lat_out *=-1.;
}

void Coords::parseUTMZone(const std::string& zone_info, char& zoneLetter, short int& zoneNumber) const
{ //helper method: parse a UTM zone specification string into letter and number
	if ((sscanf(zone_info.c_str(), "%hd%c", &zoneNumber, &zoneLetter) < 2) &&
		(sscanf(zone_info.c_str(), "%hd %c)", &zoneNumber, &zoneLetter) < 2)) {
			throw InvalidFormatException("Can not parse given UTM zone: "+zone_info,AT);
	}
	if(zoneNumber<1 || zoneNumber>60) {
		throw InvalidFormatException("Invalid UTM zone: "+zone_info+" (zone should be between 1 and 60, zone letter in [C-M, N-X])",AT);
	}
	zoneLetter = (char)toupper(zoneLetter); //just in case... (sorry for the pun!)
	if(zoneLetter=='Y' || zoneLetter=='Z' || zoneLetter=='A' || zoneLetter=='B') {
			//Special zones for the poles: we should NOT use UTM in these regions!
			throw InvalidFormatException("Invalid UTM zone: "+zone_info+" (trying to use UTM in polar regions)",AT);
	}
}

/**
* @brief Coordinate conversion: from WGS84 Lat/Long to proj4 parameters
* @param lat_in Decimal Latitude
* @param long_in Decimal Longitude
* @param east_out easting coordinate (target system)
* @param north_out northing coordinate (target system)
*/
void Coords::WGS84_to_PROJ4(double lat_in, double long_in, double& east_out, double& north_out) const
{
#ifdef PROJ4
	const std::string src_param="+proj=latlong +datum=WGS84 +ellps=WGS84";
	const std::string dest_param="+init=epsg:"+coordparam;
	projPJ pj_latlong, pj_dest;
	double x=long_in*Cst::to_rad, y=lat_in*Cst::to_rad;

	if ( !(pj_dest = pj_init_plus(dest_param.c_str())) ) {
		pj_free(pj_dest);
		throw InvalidArgumentException("Failed to initalize Proj4 with given arguments: "+dest_param, AT);
	}
	if ( !(pj_latlong = pj_init_plus(src_param.c_str())) ) {
		pj_free(pj_latlong);
		pj_free(pj_dest);
		throw InvalidArgumentException("Failed to initalize Proj4 with given arguments: "+src_param, AT);
	}

	const int p = pj_transform(pj_latlong, pj_dest, 1, 1, &x, &y, NULL );
	if(p!=0) {
		pj_free(pj_latlong);
		pj_free(pj_dest);
		throw ConversionFailedException("PROJ4 conversion failed: "+p, AT);
	}
	east_out = x;
	north_out = y;
	pj_free(pj_latlong);
	pj_free(pj_dest);
#else
	(void)lat_in;
	(void)long_in;
	(void)east_out;
	(void)north_out;
	throw IOException("Not compiled with PROJ4 support", AT);
#endif
}

/**
* @brief Coordinate conversion: from proj4 parameters to WGS84 Lat/Long
* @param east_in easting coordinate (Swiss system)
* @param north_in northing coordinate (Swiss system)
* @param lat_out Decimal Latitude
* @param long_out Decimal Longitude
*/
void Coords::PROJ4_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
{
#ifdef PROJ4
	const std::string src_param="+init=epsg:"+coordparam;
	const std::string dest_param="+proj=latlong +datum=WGS84 +ellps=WGS84";
	projPJ pj_latlong, pj_src;
	double x=east_in, y=north_in;

	if ( !(pj_src = pj_init_plus(src_param.c_str())) ) {
		pj_free(pj_src);
		throw InvalidArgumentException("Failed to initalize Proj4 with given arguments: "+src_param, AT);
	}
	if ( !(pj_latlong = pj_init_plus(dest_param.c_str())) ) {
		pj_free(pj_latlong);
		pj_free(pj_src);
		throw InvalidArgumentException("Failed to initalize Proj4 with given arguments: "+dest_param, AT);
	}

	const int p = pj_transform(pj_src, pj_latlong, 1, 1, &x, &y, NULL );
	if(p!=0) {
		pj_free(pj_latlong);
		pj_free(pj_src);
		throw ConversionFailedException("PROJ4 conversion failed: "+p, AT);
	}
	long_out = x*RAD_TO_DEG;
	lat_out = y*RAD_TO_DEG;
	pj_free(pj_latlong);
	pj_free(pj_src);
#else
	(void)east_in;
	(void)north_in;
	(void)lat_out;
	(void)long_out;
	throw IOException("Not compiled with PROJ4 support", AT);
#endif
}

void Coords::distance(const Coords& destination, double& o_distance, double& o_bearing) const {
//HACK: this is the 2D distance, it does not work in 3D!!
	if(isSameProj(destination)) {
		//we can use simple cartesian grid arithmetic
		o_distance = sqrt( Optim::pow2(easting - destination.getEasting()) + Optim::pow2(northing - destination.getNorthing()) );
		o_bearing = atan2( northing - destination.getNorthing() , easting - destination.getEasting() );
		o_bearing = fmod( o_bearing*Cst::to_deg+360. , 360. );
	} else {
		switch(distance_algo) {
			case GEO_COSINE:
				o_distance = cosineDistance(latitude, longitude, destination.getLat(), destination.getLon(), o_bearing);
				break;
			case GEO_VINCENTY:
				o_distance = VincentyDistance(latitude, longitude, destination.getLat(), destination.getLon(), o_bearing);
				break;
		}
	}
}

/**
* @brief Coordinate conversion: from WGS84 Lat/Long to local grid as given in coordparam
* @param lat_in Decimal Latitude
* @param long_in Decimal Longitude
* @param east_out easting coordinate (target system)
* @param north_out northing coordinate (target system)
*/
void Coords::WGS84_to_local(double lat_in, double long_in, double& east_out, double& north_out) const
{
	double alpha;
	double dist;

	if((ref_latitude==IOUtils::nodata) || (ref_longitude==IOUtils::nodata)) {
		east_out = IOUtils::nodata;
		north_out = IOUtils::nodata;
		//throw InvalidArgumentException("No reference coordinate provided for LOCAL projection", AT);
	} else {
		switch(distance_algo) {
			case GEO_COSINE:
				dist = cosineDistance(ref_latitude, ref_longitude, lat_in, long_in, alpha);
				break;
			case GEO_VINCENTY:
				dist = VincentyDistance(ref_latitude, ref_longitude, lat_in, long_in, alpha);
				break;
		}

		east_out = dist*sin(alpha*Cst::to_rad);
		north_out = dist*cos(alpha*Cst::to_rad);
	}
}


/**
* @brief Coordinate conversion: from local grid as given in coordparam to WGS84 Lat/Long
* @param east_in easting coordinate (in local system)
* @param north_in northing coordinate (in local system)
* @param lat_out Decimal Latitude
* @param long_out Decimal Longitude
*/
void Coords::local_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
{
	const double dist = sqrt( Optim::pow2(east_in) + Optim::pow2(north_in) );
	const double bearing = fmod( atan2(east_in, north_in)*Cst::to_deg+360. , 360.);

	if((ref_latitude==IOUtils::nodata) || (ref_longitude==IOUtils::nodata)) {
		lat_out = IOUtils::nodata;
		long_out = IOUtils::nodata;
		//throw InvalidArgumentException("No reference coordinate provided for LOCAL projection", AT);
	} else {
		switch(distance_algo) {
			case GEO_COSINE:
				cosineInverse(ref_latitude, ref_longitude, dist, bearing, lat_out, long_out);
				break;
			case GEO_VINCENTY:
				VincentyInverse(ref_latitude, ref_longitude, dist, bearing, lat_out, long_out);
				break;
		}
	}
}

void Coords::NULL_to_WGS84(double /*east_in*/, double /*north_in*/, double& /*lat_out*/, double& /*long_out*/) const
{
	throw InvalidArgumentException("The projection has not been initialized!", AT);
}

void Coords::WGS84_to_NULL(double /*lat_in*/, double /*long_in*/, double& /*east_out*/, double& /*north_out*/) const
{
	throw InvalidArgumentException("The projection has not been initialized!", AT);
}

/**
* @brief Spherical law of cosine Distance calculation between points in WGS84 (decimal Lat/Long)
* See http://www.movable-type.co.uk/scripts/latlong.html for more
* @param[in] lat_ref Decimal Latitude (const double&)
* @param[in] lon_ref Decimal Longitude (const double&)
* @param[in] distance Distance in meters (const double&)
* @param[in] bearing bearing in degrees, 0 being north (const double&)
* @param[out] lat Decimal latitude of target point (double&)
* @param[out] lon Decimal longitude of target point (double&)
*/
void Coords::cosineInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon)
{
	const double lat_ref_rad = lat_ref*Cst::to_rad;
	const double bearing_rad = bearing*Cst::to_rad;

	if(IOUtils::checkEpsilonEquality(distance, 0., .01)) {
		//distance is too small, it could create numerical problems
		lat = lat_ref;
		lon = lon_ref;
		return;
	}

	lat = asin( sin(lat_ref_rad)*cos(distance/Cst::earth_R0) +
	            cos(lat_ref_rad)*sin(distance/Cst::earth_R0)*cos(bearing_rad) );
	lon = lon_ref*Cst::to_rad + atan2( sin(bearing_rad)*sin(distance/Cst::earth_R0)*cos(lat_ref_rad) ,
	                              cos(distance/Cst::earth_R0) - sin(lat_ref_rad)*sin(lat) );
	lon = fmod(lon+Cst::PI, 2.*Cst::PI) - Cst::PI;

	lat *= Cst::to_deg;
	lon *= Cst::to_deg;
}

/**
* @brief Spherical law of cosine Distance calculation between points in WGS84 (decimal Lat/Long)
* See http://www.movable-type.co.uk/scripts/latlong.html for more
* @param[in] lat1 Decimal Latitude (const double&)
* @param[in] lon1 Decimal Longitude (const double&)
* @param[in] lat2 Decimal Latitude (const double&)
* @param[in] lon2 Decimal Longitude (const double&)
* @param[out] alpha average bearing
* @return distance (double)
*/
double Coords::cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha)
{
	if(lat1==lat2 && lon1==lon2) {
		//distance is zero, it creates numerical problems -> skip calculation
		alpha = 0.;
		return 0.;
	}

	const double d = acos(
		sin(lat1*Cst::to_rad) * sin(lat2*Cst::to_rad)
		+ cos(lat1*Cst::to_rad) * cos(lat2*Cst::to_rad) * cos((lon2-lon1)*Cst::to_rad)
		) * Cst::earth_R0;

	alpha = atan2( sin((lon2-lon1)*Cst::to_rad)*cos(lat2*Cst::to_rad) ,
			cos(lat1*Cst::to_rad)*sin(lat2*Cst::to_rad) - sin(lat1*Cst::to_rad)*cos(lat2*Cst::to_rad)*cos((lon2-lon1)*Cst::to_rad)
 		) / Cst::to_rad;
	alpha = fmod((alpha+360.), 360.);

	return d;
}

/**
* @brief Vincenty Distance calculation between points in WGS84 (decimal Lat/Long)
* See T. Vincenty, "Closed formulas for the direct and reverse geodetic problems",
* Journal of Geodesy, 51, 3, 1977, DOI:10.1007/BF02521599,
* see http://www.springerlink.com/content/y7108u6862473583 for more
* @param[in] lat1 Decimal Latitude (const double&)
* @param[in] lon1 Decimal Longitude (const double&)
* @param[in] lat2 Decimal Latitude (const double&)
* @param[in] lon2 Decimal Longitude (const double&)
* @param[in] alpha average bearing (double&)
* @return distance (double)
*/
double Coords::VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha)
{
	const double thresh = 1.e-12;	//convergence absolute threshold
	const int n_max = 100;		//maximum number of iterations
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double f = (a - b) / a;	//ellispoid flattening

	const double L = (lon1 - lon2)*Cst::to_rad;
	const double U1 = atan( (1.-f)*tan(lat1*Cst::to_rad) );
	const double U2 = atan( (1.-f)*tan(lat2*Cst::to_rad) );

	double lambda = L, lambda_p=0.;
	double sin_sigma, cos_sigma, sigma, cos_alpha2, cos_2sigma_m;
	int n=0;
	do {
		sin_sigma = sqrt( Optim::pow2(cos(U2)*sin(lambda)) + Optim::pow2(cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda)) );
		if(sin_sigma==0.) return 0.; //co-incident points
		cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda);
		sigma = atan2(sin_sigma,cos_sigma);
		const double sin_alpha = cos(U1)*cos(U2)*sin(lambda) / sin_sigma;
		cos_alpha2 = 1. - Optim::pow2(sin_alpha);
		if(lat1==0. && lat2==0.) {
			cos_2sigma_m = 0.;
		} else {
			cos_2sigma_m = cos_sigma - 2.*sin(U1)*sin(U2)/cos_alpha2;
		}
		const double C = f/16. * cos_alpha2*(4.+f*(4.-3.*cos_alpha2));
		lambda_p = lambda;
		lambda = L + (1.-C)*f*sin_alpha*(
			sigma + C*sin_sigma*( cos_2sigma_m + C * cos_sigma * (-1.+2.*Optim::pow2(cos_2sigma_m)) )
			);
		n++;
	} while ( (n<n_max) && (fabs(lambda - lambda_p) > thresh) );

	if(n>n_max) {
		throw IOException("Distance calculation not converging", AT);
	}

	const double u2 = cos_alpha2 * (a*a - b*b) / (b*b);
	const double A = 1. + u2/16384. * ( 4096.+u2*(-768.+u2*(320.-175.*u2)) );
	const double B = u2/1024. * ( 256.+u2*(-128.+u2*(74.-47.*u2)) );
	const double delta_sigma = B*sin_sigma*( cos_2sigma_m+B/4.*( cos_sigma*(-1.+2.*Optim::pow2(cos_2sigma_m)) - B/6.*(cos_2sigma_m*(-3.+4.*Optim::pow2(sin_sigma))*(-3.+4.*Optim::pow2(cos_2sigma_m))) ) );

	const double s = b*A*(sigma - delta_sigma);	//distance between the two points

	//computation of the average forward bearing
	double alpha1 = atan2(cos(U2)*sin(lambda), cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(lambda)) / Cst::to_rad; //forward azimuth
	//double alpha2 = atan2(cos(U1)*sin(lambda), sin(U1)*cos(U2)-cos(U1)*sin(U2)*cos(lambda)) / Cst::to_rad; //reverse azimuth

	//trying to get a normal compass bearing... TODO: make sure it works and understand why
	alpha1 = fmod(-alpha1+360., 360.);
	//alpha2 = fmod(alpha2+180., 360.);

	//we only keep the forward bearing, otherwise the reverse projection will not produce the initial point
	alpha = alpha1;
	return s;
}

/**
* @brief Vincenty Inverse calculation giving WGS84 (decimal Lat/Long) position
* given a start location (lat,lon) a distance and a bearing
* See T. Vincenty, "Closed formulas for the direct and reverse geodetic problems",
* Journal of Geodesy, 51, 3, 1977, DOI:10.1007/BF02521599,
* see http://www.springerlink.com/content/y7108u6862473583 for more
* @param[in] lat_ref Decimal Latitude (const double&)
* @param[in] lon_ref Decimal Longitude (const double&)
* @param[in] distance Distance in meters (const double&)
* @param[in] bearing bearing in degrees, 0 being north (const double&)
* @param[out] lat Decimal latitude of target point (double&)
* @param[out] lon Decimal longitude of target point (double&)
*/
void Coords::VincentyInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon)
{//well, actually this is the DIRECT Vincenty formula
	const double thresh = 1.e-12;	//convergence absolute threshold
	const double a = ellipsoids[E_WGS84].a;	//major ellipsoid semi-axis, value for wgs84
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis, value for wgs84
	const double f = (a - b) / a;	//ellispoid flattening

	const double alpha1 = bearing*Cst::to_rad;
	const double tanU1 = (1.-f)*tan(lat_ref*Cst::to_rad);
	const double cosU1 = 1./sqrt(1.+tanU1*tanU1);
	const double sinU1 = tanU1*cosU1;
	const double sigma1 = atan2(tanU1,cos(alpha1));
	const double sinAlpha = cosU1*sin(alpha1);
	const double cos2alpha = 1. - sinAlpha*sinAlpha;
	const double u2 = cos2alpha * (a*a - b*b) / (b*b);
	const double A = 1. + u2/16384. * (4096. + u2*(-768.+u2*(320.-175.*u2)) );
	const double B = u2/1024. * (256. + u2*(-128.+u2*(74.-47.*u2)));

	double sigma = distance / (b*A);
	double sigma_p = 2.*Cst::PI;
	double cos2sigma_m = cos( 2.*sigma1 + sigma ); //required to avoid uninitialized value

	while (fabs(sigma - sigma_p) > thresh) {
		cos2sigma_m = cos( 2.*sigma1 + sigma );
		double delta_sigma = B*sin(sigma) * ( cos2sigma_m + B/4. * (
			cos(sigma)*(-1.+2.*cos2sigma_m*cos2sigma_m)
			-B/6. * cos2sigma_m * (-3.+4.*Optim::pow2(sin(sigma))) * (-3.+4.*cos2sigma_m*cos2sigma_m)
			) );
		sigma_p = sigma;
		sigma = distance / (b*A) + delta_sigma;
	}

	lat = atan2( sinU1*cos(sigma) + cosU1*sin(sigma)*cos(alpha1),
		     (1.-f) * sqrt( sinAlpha*sinAlpha + Optim::pow2(sinU1*sin(sigma) - cosU1*cos(sigma)*cos(alpha1)) )
		   );
	const double lambda = atan2( sin(sigma)*sin(alpha1), cosU1*cos(sigma) - sinU1*sin(sigma)*cos(alpha1) );
	const double C = f/16. * cos2alpha * (4.+f*(4.-3.*cos2alpha));
	const double L = lambda - (1.-C) * f * sinAlpha * (
				sigma + C * sin(sigma) * ( cos2sigma_m+C*cos(sigma) * (-1.+2.*cos2sigma_m*cos2sigma_m) )
				);

	lat = lat * Cst::to_deg;
	lon = lon_ref + (L*Cst::to_deg);
	//const double alpha2 = atan2( sinAlpha, -(sinU1*sin(sigma)-cosU1*cos(sigma)*cos(alpha1)) ); //reverse azimuth
}

void Coords::clearCoordinates() {
//sets safe defaults for all internal variables (except function pointers and maps)
	latitude = longitude = IOUtils::nodata;
	altitude = IOUtils::nodata;
	easting = northing = IOUtils::nodata;
	grid_i = grid_j = grid_k = IOUtils::inodata;
}

void Coords::setDefaultValues() {
//sets safe defaults for all internal variables (except function pointers and maps)
	clearCoordinates();
	ref_latitude = ref_longitude = IOUtils::nodata;
	distance_algo = GEO_COSINE;
}

} //end namespace
