#pragma once


#include <string>

#include <boost/utility.hpp>

#include "forcing_data.hpp"
/*
Class: Station
Concept of a met station.
Allows the station to have a location (x,y) and a station ID. 
As well, it wraps the iterators of the forcing_data instance.

Example:
>Station s("ExampleStation1","obs.txt",5,400);
>s.now().get<int>("RH");
>s.next();
*/
class station : boost::noncopyable
{
public:
        
    /*
    Function: Station
    Default constructor

    Parameters: 
    None

    Throws:
    Never

    Returns:   
    - 
     */
    station();

    /*
    Function: Station
    Creates a station from 

    Parameters: 
    std::string ID - Station ID
    std::string file - Met file to open
    unsigned int x - X coord
    unsigned int y - Y coord
    float elevation - Station elevation

    Throws:
    std::runtime_error - on error

    Returns:   
    - 
     */
    station(std::string ID, std::string file, unsigned int x, unsigned int y, float elevation);

    /*
    Function: open
    Opens a given met file

    Parameters: 
    std::string file - path to met file.

    Throws:
    std::runtime_error - on error

    Returns:   
    void 
     */
    void open(std::string file);


    /*
    Function: now
    Returns the data for the current time step

    Parameters: 
    None

    Throws:
    Never

    Returns:   
    forcing_data::timestep - The current timestep. See <timestep>
     */
    forcing_data::timestep now();


    /*
    Function: next
    Gets the next timestep	

    Parameters: 
    None

    Throws:
    Never

    Returns:   
    bool - Returns true if this is the last timestep
     */
    bool next();

    /*
    Function: get_x
    Gets the X coordinate. Begins a 0.

    Parameters: 
    None

    Throws:
    Never

    Returns:   
    unsigned int - X corrdinate of the station
     */
    unsigned int get_x();

    /*
    Function: get_y
    Gets the Y coordinate. Begins at 0

    Parameters: 
    None

    Throws:
    Never

    Returns:   
    unsigned int - Y corrdinate of the station
     */
    unsigned int get_y();


    /*
    Function: set_x
    Sets the X coordinate. Begins at 0

    Parameters: 
    unsigned int x - X coordinate

    Throws:
    Never

    Returns:   
    void - 
     */
    void set_x(unsigned int x);


    /*
    Function: set_y
    Sets the Y coordinate. Begins at 0

    Parameters: 
    unsigned int y - Y coordinate

    Throws:
    Never

    Returns:   
    void - 
     */
    void set_y(unsigned int y);


    /*
    Function: get_elevation
             Returns the station's elevation

    Parameters: 

    Throws:
            Never

    Returns:   
            float - Station elevation
     */
    float get_elevation();


    /*
    Function: set_z
            Sets the station elevation	 

    Parameters: 
            float elevation - Station elevation

    Throws:
            Never

    Returns:   
            void 
     */
    void set_z(float elevation);
        

    void set_ID(std::string ID);


    std::string get_ID();
    
friend std::ostream& operator<<(std::ostream &strm, const station &s) ;

private:
        std::string _ID;
        forcing_data* _obs;
        forcing_data::const_iterator _itr;
        unsigned int _x;
        unsigned int _y;
        float _z;
};


