#pragma once


#include <string>

#include <boost/utility.hpp>

#include "timeseries.hpp"

/**
* \class station
*
* \brief Concept of a met station.
*
* Allows the station to represent a timeseries that has a location (x,y), an elevation, and a station ID.
* As well, it wraps the iterators of the timeserires instance allowing for easy stepping and access.
*/
//namespace CGAL
//{
//    class Exact_predicates_inexact_constructions_kernel;
//}



class station : boost::noncopyable
{
public:
        
    /**
    *    Default constructor
    */
    station();

    /**
        Creates a new station with the specified attributes.

        \param ID Station name
        \param x UTM coord
        \param y UTM coord
        \param elevation station elevation
        */
    station(std::string ID, double x, double y, double elevation);

    /**
    * Default destructor
    */
    ~station();
    
    /**
    * Opens a metfile
    * \param file Fully qualified path to the file to open
    */
    void open(std::string file);


    /**
    * Returns the timestep iterator for the current timestep
    * \return A timstep iterator
    */
     timestep& now() ;

    /**
    * Increments the internal timestep iterator to the next timestep
    * \return False if next timestep is one past end of timeseries
    */
    bool next();


    /**
    * Returns the x UTM coordinate of the station
    * \return UTM coordinate
    */
    double x();

    /**
    * Returns the y UTM coordinate of the station
    * \return UTM coordinate
    */
    double y();


    /**
    * Sets the X coordinate
    * \param x UTM coordinate
    */
    void x(double x);


    /**
    * Sets the Y coordinate
    * \param y UTM coordinate
    */
    void y(double y);


    /**
    * Returns the station elevation (m)
    * \return Station elevation
    */
    double z();


    /**
    * Sets the station elevation (m).
    * \param elevation Station elevation
    */
    void z(double elevation);

    /**
    * Sets the ID of the station (name)
    * \param ID
    */
    void ID(std::string ID);

    /**
    * Returns the station ID
    * \return Station ID
    */
    std::string ID();

    /**
    * Resets the internal iterators to point to the begingin of the timeseries
    */
    void reset_itrs();

    /**
    * List all (including module provided) variables
    * \return Vector containing a list of variable names
    */
    std::vector<std::string> list_variables();

    /**
    * Returns the vector of all the dates in the timeseries
    * \return
    */
    timeseries::date_vec date_timeseries();

    /**
     * Adds a new variable to the underlying timeseries. Invalidates all internal iterators and resets them to the beging of the timeseries
     * @param var Variable name to add
     */
    void add_variable(std::string var);

    /**
    * Returns the underlying timeseries
    * \return
    */
    timeseries* raw_timeseries();

    /**
    * Returns the length of the timeseries. This is the number of records in the timeseries
    */
    size_t timeseries_length();

    /**
    * Returns the specificed variable of the station at this current timeseries\
    * \param variable variable name
    * \return value of the requested variable
    */
    double get(std::string variable);

    /**
     * Save the internal timeseres to a file
     * \param filename
     */
    void tofile(std::string file);

    /**
     * Obtain the triangle id that is closest to this station
     * @return
     */
    size_t closest_face();

    /**
     * Set the face id that the station is closest to
     * @param face
     */
    void set_closest_face(size_t face);

friend std::ostream& operator<<(std::ostream &strm, const station &s) ;

private:
        std::string _ID;
        timeseries * _obs;
        timeseries::iterator _itr;
        double _x;
        double _y;
        double _z;


    // There is a circular reference that is created when station knows about the triangulation
    // I could not get it to properly compile the mesh_elem reference
    // So this stores the id, which we can easily look up in domain modules
    // Perhaps rethink this if it becomes a problem later
        size_t _face; // holds the face that the station is closet to
   
};


