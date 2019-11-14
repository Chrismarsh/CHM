//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#pragma once


#include <string>

#include <boost/utility.hpp>

#include "variablestorage.hpp"
#include "timeseries.hpp"

/**
* \class station
*
* \brief Concept of a met station.
*
* Allows the station to represent a timeseries that has a location (x,y), an elevation, and a station ID.
* As well, it wraps the iterators of the timeserires instance allowing for easy stepping and access.
*/
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
    * Resets the internal iterators to point to the begining of the timeseries
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

    double& operator[](const uint64_t& hash);
    double& operator[](const std::string& variable);

friend std::ostream& operator<<(std::ostream &strm, const station &s) ;

private:
        std::string _ID;
        timeseries * _obs;
        timeseries::iterator _itr;
        double _x;
        double _y;
        double _z;

        variablestorage _timestep_data;

};


