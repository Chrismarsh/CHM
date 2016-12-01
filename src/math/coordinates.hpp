#pragma once
#include <cmath>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/squared_distance_2.h>
#include <boost/function.hpp>
namespace math
{
    namespace gis
    {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef K::Point_3 Point_3;
        typedef K::Point_2 Point_2;
        typedef K::Vector_2 Vector_2;

        /**
         * Computes a new point based on bearing (rad, CW from North) and distance for Lat long
         * http://www.movable-type.co.uk/scripts/latlong.html
         * @param src Source point in decimal degress
         * @param bearing Degrees, CW from north
         * @param distance Meters
         * @return
         */
        Point_2 point_from_bearing_latlong(Point_3 src, double bearing, double distance);

        /**
         * Computes a new point based on bearing (rad, CW from North) and distance for UTM
         * http://gis.stackexchange.com/a/76175
         * @param src Source point in UTM meters
         * @param bearing Degrees, CW from north
         * @param distance Meters
         * @return
         */
        Point_2 point_from_bearing_UTM(Point_3 src, double bearing, double distance);

        /**
         * Computes the distance between two lat long pairs using haversine formula
         * @param pt1
         * @param pt2
         * @return
         */
        double distance_latlong(Point_3 pt1, Point_3 pt2);
        double distance_UTM(Point_3 pt1, Point_3 pt2);

        extern boost::function<double(Point_3 pt1, Point_3 pt2)> distance;
        extern boost::function<Point_2(Point_3 src, double bearing, double distance)> point_from_bearing;

        /**
         * Converts a North-based bearing to polar coodinates.
         * @param bearing In degrees, N=0
         * @return Returns theta in radians
         */
        double bearing_to_polar(double bearing);

        /**
         * Convertes a polar unit vector to cartesian components
         * @param theta In radians
         * @return
         */
        Vector_2 polar_to_cartesian(double theta);

        /**
         * Converts a North-based bearing to cartesian vector components
         * @param bearing In degrees, N=0
         * @return x,y components
         */
        Vector_2 bearing_to_cartesian(double bearing);

        /**
         * Converts a cartesian x,y vector into compar bearing, N=0.
         * @param cart Degrees
         * @return
         */
        double cartesian_to_bearing(Vector_2 cart);

    }
}