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
* Allows the station to represent a timestep that has a location (x,y), an elevation, and a station ID.
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
    station(std::string ID, double x, double y, double elevation, std::set<std::string> variables = {});

    /**
    * Default destructor
    */
    ~station();
    
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
    * List all (including module provided) variables
    * \return Vector containing a list of variable names
    */
    std::vector<std::string> list_variables();

    /**
     * Initializes the station to store the specified variables. This will destroy all data
     * @param variables
     */
    void init(std::set<std::string> variables);



    /**
    * Returns the current hour, 24-hour format
    */
    int hour();

    /**
    * Returns the current minute
    */
    int min();

    /**
    * Returns the current second
    */
    int sec();

    /**
    *  Current month number, starting with Jan = 1.
    */
    int month();


    /**
    *  Current day number starting at 1
    */
    int day();


    /**
    * Current year
    */
    int year();

    /**
    * Returns the boost::gregoiran structure for time calculations.
    */
    boost::gregorian::date get_gregorian();

    /**
    * Returns the boost::posix_time structure for time calculations
    */
    boost::posix_time::ptime get_posix();

    /**
    * Returns the boost::posix_time structure for time calculations
    */
    void set_posix(boost::posix_time::ptime ts );

    /**
 * Returns true if the specified variable is available in the timeseries.
 */
    bool has(const std::string &variable);

    double& operator[](const uint64_t& hash);
    double& operator[](const std::string& variable);



friend std::ostream& operator<<(std::ostream &strm, const station &s) ;

private:
    std::string _ID;
    double _x;
    double _y;
    double _z;

    variablestorage _timestep_data;

    boost::posix_time::ptime _current_ts;

};


