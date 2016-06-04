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
#ifndef __COORDS_H__
#define __COORDS_H__

#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#include <string>
#include <iostream>

namespace mio {
	
/**
 * @page dev_coords Coordinate systems support developement guide
 * Geographic coordinates are transparently supported (ie converted to/from lat/lon as well as between cartesian coordinate systems) through
 * the Coords class. Therefore, coordinate system conversions are quite easy, as illustrated by the example below.
 * @code
 * Coords point1("CH1903","");
 * point1.setXY(785425. , 191124., 1400.);
 * std::cout << "Lat=" << point1.getLat() << " lon=" << point1.getLon() << "\n"
 * @endcode
 * Adding support for another cartesian coordinate system is simple, once you have the official conversion wgs84 to/from this system.
 * 
* @section coords_implementation Conversion implementation
* National cartesian coordinate systems often come with an official conversion formula to and from WGS84 lat/lon. This could be found by
* the national geographic authority. Then implement two <b>private</b> methods into the Coords class (replacing the 
* "XXX" by an abbreviation for the coordinate system you are implementing):
* @code
 * void WGS84_to_XXX(double lat_in, double long_in, double& east_out, double& north_out) const;
 * void XXX_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;
 * @endcode
 * It is highly recommended to add the link to the official conversion specification in the doxygen comments of these two methods.
 * 
 * @section coords_registering Registering the coordinate system
 * Once these two methods have been implemented, they need to be registered for the user to be able to use them. This is done by following
 * these steps:
 *      -# retrieve the appropriate EPSG code at http://www.epsg-registry.org/ for the coordinate system. It is a good idea to use the associated abbreviation
 *         as keyword/abbreviation in the source code as well as for the end user.
 *      -# map the coordinate system abbreviation to its EPSG code in Coords::setEPSG()
 *      -# map the coordinate system abbreviation to its EPSG code in Coords::getEPSG()
 *      -# link the coordinate system abbreviation to the conversion methods you have implemented in Coords::convert_to_WGS84() and Coords::convert_from_WGS84()
 * 
 * @section coords_documentation Documenting the coordinate system
 * Please update and expand the doxygen documentation at the begining of the Coords.cc file in order to specify the coordinate system keywords 
 * that has been used (ie the the coordinate system abbreviation). Feel free to add links to official documentation about this coordinate system. 
 * Please also properly document the conversion (with links to the official specification) in the conversion implementations. 
 * 
 * @section coord_testing Testing the implementation
 * Try to set a Coords object constructed with the chosen keywords to a set of Easting/Northing and then retrieve the lat/lon as well as the 
 * opposite. Make sure that for various points, including points close to the boundaries of the coordinate system definition, the conversion 
 * remains correct. Often the official specifications come with a set of test coordinates that you can use. Then try to set the coordinate
 * system by EPSG code and make sure Easting/Northing to and from WGS84 conversions work properly.
 * 
 */

/**
 * @class Coords
 * @brief A class to handle geographic coordinate systems.
 * This class offers an easy way to transparently convert between various coordinate systems. For any
 * given object, as soon as a latitude/longitude can be computed/read, it will be used as a reference.
 * This means that every subsequent change of projection system or read will be the conversion of this
 * reference lat/lon position (until a new "set" is called). See Coords::setProj for the supported coordinate systems.
 *
 * @ingroup data_str
 * @author Mathias Bavay
 * @date   2009-01-20
*/

class Coords {
public:
	///Keywords for selecting the algorithm for computing geodesic distances
	typedef enum GEO_DISTANCES {
		GEO_COSINE, ///< Spherical law of cosine (See http://www.movable-type.co.uk/scripts/latlong.html)
		GEO_VINCENTY ///< Vincenty ellispoid formula (See T. Vincenty, "Closed formulas for the direct and reverse geodetic problems", Journal of Geodesy, 51, 3, 1977, DOI:10.1007/BF02521599, or http://www.springerlink.com/content/y7108u6862473583 for more)
	} geo_distances;

	//Constructors
	Coords();
	Coords(const std::string& in_coordinatesystem, const std::string& in_parameters="");
	Coords(const std::string& in_coordinatesystem, const std::string& in_parameters, const std::string& coord_spec);
	Coords(const double& in_lat_ref, const double& in_long_ref);
	Coords(const Coords& c);

	//Operators
	Coords& operator=(const Coords&); ///<Assignement operator
	bool operator==(const Coords&) const; ///<Operator that tests for equality
	bool operator!=(const Coords&) const; ///<Operator that tests for inequality
	bool isNodata() const;
	void moveByXY(const double& x_displacement, const double& y_displacement);
	void moveByBearing(const double& i_bearing, const double& i_distance);

	static Coords merge(const Coords& coord1, const Coords& coord2);
	void merge(const Coords& coord2);

	//Getter methods
	double getEasting() const;
	double getNorthing() const;
	double getLat() const;
	double getLon() const;
	double getAltitude() const;
	int getGridI() const;
	int getGridJ() const;
	int getGridK() const;
	void getProj(std::string& proj_type, std::string& proj_args) const;
	std::string printLatLon() const;
	short int getEPSG() const;

	const std::string toString() const;
	friend std::iostream& operator<<(std::iostream& os, const Coords& coord);
	friend std::iostream& operator>>(std::iostream& is, Coords& coord);

	//Setter methods
	void setLatLon(const double in_latitude, const double in_longitude, const double in_altitude, const bool in_update=true);
	void setLatLon(const std::string& in_coordinates, const double in_altitude, const bool in_update=true);
	void setXY(const double in_easting, const double in_northing, const double in_altitude, const bool in_update=true);
	void setGridIndex(const int in_grid_i, const int in_grid_j, const int in_grid_k, const bool in_invalidate=true);
	void setAltitude(const double in_altitude, const bool in_update=true);
	void setProj(const std::string& in_coordinatesystem, const std::string& in_parameters="");
	void setLocalRef(const double in_ref_latitude, const double in_ref_longitude);
	void setLocalRef(const std::string in_coordparam);
	void setDistances(const geo_distances in_algo);
	void setEPSG(const int epsg);

	void check();
	double distance(const Coords& destination) const;
	bool isSameProj(const Coords& target) const;
	void copyProj(const Coords& source, const bool i_update=true);

	//Static helper methods
	static double dms_to_decimal(const std::string& dms);
	static std::string decimal_to_dms(const double& decimal);
	static void parseLatLon(const std::string& coordinates, double& lat, double& lon);

	static double lat_degree_lenght(const double& latitude);
	static double lon_degree_lenght(const double& latitude);

	static void rotatedToTrueLatLon(const double& lat_N, const double& lon_N, const double& lat_rot, const double& lon_rot, double &lat_true, double &lon_true);
	static void trueLatLonToRotated(const double& lat_N, const double& lon_N, const double& lat_true, const double& lon_true, double &lat_rot, double &lon_rot);

	static double cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha);
	static void cosineInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon);
	static double VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha);
	static void VincentyInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon);

 private:
	//Coordinates conversions
	void convert_to_WGS84(double i_easting, double i_northing, double& o_latitude, double& o_longitude) const;
	void convert_from_WGS84(double i_latitude, double i_longitude, double& o_easting, double& o_northing) const;

	void WGS84_to_CH1903(double lat_in, double long_in, double& east_out, double& north_out) const;
	void CH1903_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;
	void WGS84_to_UTM(double lat_in, double long_in, double& east_out, double& north_out) const;
	void UTM_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;
	void WGS84_to_UPS(double lat_in, double long_in, double& east_out, double& north_out) const;
	void UPS_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;
	void WGS84_to_PROJ4(double lat_in, double long_in, double& east_out, double& north_out) const;
	void PROJ4_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;
	void WGS84_to_local(double lat_in, double long_in, double& east_out, double& north_out) const;
	void local_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;
	void WGS84_to_NULL(double lat_in, double long_in, double& east_out, double& north_out) const;
	void NULL_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;

	void parseUTMZone(const std::string& zone_info, char& zoneLetter, short int& zoneNumber) const;

	//Distances calculations
	void distance(const Coords& destination, double& o_distance, double& o_bearing) const;

 private:
	void clearCoordinates();
	void setDefaultValues();
	int getUTMZone(const double latitude, const double longitude, std::string& zone_out) const;

 private:
	double ref_latitude, ref_longitude;
	double altitude; ///<altitude of the point (the altitude is currently NOT dependant on the projection)
	double latitude; ///<latitude of the point
	double longitude; ///<longitude of the point
	double easting; ///<east coordinate of the point in a cartesian grid
	double northing; ///<north coordinate of the point in a cartesian grid
	int grid_i; ///<grid index i (please notice that this index is NOT automatically regenerated NOR checked)
	int grid_j; ///<grid index j (please notice that this index is NOT automatically regenerated NOR checked)
	int grid_k; ///<grid index k (please notice that this index is NOT automatically regenerated NOR checked)

	std::string coordsystem;
	std::string coordparam;
	geo_distances distance_algo;

	///Keywords for selecting an ellipsoid to use
	enum ELLIPSOIDS_NAMES {
		E_WGS84, ///<Globally useable WGS84 ellipsoid
		E_GRS80, ///<GRS80 ellispoid, equivalent to WGS84 but used by America and Australia
		E_AIRY, ///<Airy ellispoid, good fit for the UK
		E_INTL1924, ///<International 1924 ellispoid, good for most of Europe
		E_CLARKE1880, ///<Clarke 1880, good for Africa
		E_GRS67 ///<GRS67 ellispoid, good for South America
	};
	struct ELLIPSOID {
		double a;
		double b;
	};
	static const struct ELLIPSOID ellipsoids[6];
};
} //end namespace

#endif
