#pragma once



#include <string>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>

namespace pt = boost::property_tree;


#include <tbb/concurrent_vector.h>

#include "interpolation.h"
#include "station.hpp"
#include "math/coordinates.hpp"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>


#include <CGAL/Euclidean_distance.h>

/**
 * Basin wide parameters such as transmissivity, solar elevation, solar aspect, etc.
 * 
 **/
class global
{

    struct Distance
    {

    };

    typedef CGAL::Simple_cartesian<double> Kernel;
    typedef boost::tuple<Kernel::Point_2, boost::shared_ptr<station> >  Point_and_station;
    typedef CGAL::Search_traits_2<Kernel> Traits_base;
    typedef CGAL::Search_traits_adapter<Point_and_station,
            CGAL::Nth_of_tuple_property_map<0, Point_and_station>,
            Traits_base>                                              Traits;

    typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_circle;
    typedef CGAL::Kd_tree<Traits> Tree;

    typedef CGAL::Orthogonal_k_neighbor_search<Traits> Neighbor_search;

     //want to let core modify date time, etc without showing a public interface. 
    //This is because global gets passed to all modules and a rogue module could do something dumb
    //const doesn't save us as we actually do want to modify things
    friend class core;
    
private:
    boost::posix_time::ptime _current_date;
    var _variables;
    Tree _dD_tree; //spatial query tree
    //each station where observations are
    std::vector< boost::shared_ptr<station> > _stations;
    int _dt; //seconds
    bool _is_geographic;


    //if radius selection for stations is chosen this holds that
    double radius;
    double N; // meters, radius for station search

public:

    boost::function< std::vector< boost::shared_ptr<station> > ( double, double) > get_stations;

    //approximate UTC offset
    int _utc_offset;
    bool is_geographic();
    global();
    int year();
    int day();
    interp_alg interp_algorithm;
    /*
     * Montboosth on [1,12]
     */
    int month();
    int hour();
    int min();
    int sec();
    int dt();
    boost::posix_time::ptime posix_time();
    uint64_t posix_time_int();


    double station_search_radius;
    bool first_time_step;
    std::string get_variable(std::string variable);

    /**
     * Inserts a new station
     * @param s
     */
    void insert_station(boost::shared_ptr<station> s);

    /**
 * Inserts new stations for a vector of stations
 * @param s
 */
    void insert_stations(tbb::concurrent_vector< boost::shared_ptr<station> >& stations);

    /**
     * Returns a set of stations within the search radius (meters) centered on the point x,y
     * @param x
     * @param y
     * @param radius
     * @return List stations that satisfy search criterion
     */
    std::vector< boost::shared_ptr<station> > get_stations_in_radius(double x, double y,double radius);

    /**
     * Returns the nearest station to x,y. Ignores elevation
     * @param x
     * @param y
     * @param N Number neighbours to find
     * @return
     */
    std::vector< boost::shared_ptr<station> > nearest_station(double x, double y,unsigned int N=1);


    std::vector< boost::shared_ptr<station> > stations();
//    template<class T>
//    T get_parameter_value(std::string key);

    pt::ptree parameters;
    
};
//
//template<class T>
//T global::get_parameter_value(std::string key)
//{
//    return parameters.get<T>(key);
//}
