#pragma once

#include <boost/date_time/posix_time/posix_time.hpp> // for boost::posix


#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_vector.h>
#include <vector>
#include "exception.hpp"
#include "variable_map.hpp"
#include "crc_hash_compare.hpp"
#include "logger.h"

/*
        Class: timestep
        A given time step. 

 */
class timestep 
{
private:
    friend class time_series;
    
    //two different. boost::variant solves this, but is very slow 
    // needs to be either a boost::posix_time or double, and is almost always a double
    typedef tbb::concurrent_vector< double > variable_vec; 
    typedef tbb::concurrent_vector<  boost::posix_time::ptime > date_variable; 

    typedef tbb::concurrent_hash_map<std::string, variable_vec::iterator, crc_hash_compare> itr_map;

    //holds the iterators for the current timestep. 
    //these are iterators into each vector in the variable hashmap
    itr_map _itrs; 
    date_variable::iterator _date_itr;

    
public:

    /*
    Function: timestep
    Default empty constructor

    Parameters: 

    Throws:
    Never

    Returns:   
    - 
     */
    timestep();


    /*
    Function: timestep
    Copy constructor

    Parameters: 
    const timestep & src - 

    Throws:
    Never

    Returns:   
    - 
     */
    timestep(const timestep& src);

    ~timestep();

    /*
    Function: to_string
    Converts the current time step to a std::string. There is no guarantee made about variable output order. 

    Parameters: 

    Throws:
    Never

    Returns:   
    std::string - Variables. Dates will be in the form 2008-Feb-23 23:59:59
     */
    std::string to_string();




    /*
    Function: month
             Gets the current month. At worst, will do O(n) where n=number of variables headers.
             If this function is called many times, or more calculations are needed it might be best either save the value or use "GetGreorian" or "get_posix" methods. 

    Parameters: 
            None

    Throws:
            Never

    Returns:   
            std::string - Current month number, starting with Jan = 1. If the date can't be determined (perhaps there is no date value) returns -1
     */
    int month();


    /*
    Function: day
             Gets the current day number. At worst, will do O(n) where n=number of variables headers.
             If this function is called many times, or more calculations are needed it might be best either save the value or use "GetGreorian" or "get_posix" methods. 

    Parameters: 
            None

    Throws:
            Never

    Returns:   
            int - Current day number starting at 1. -1 on error
     */
    int day();


    /*
    Function: year
             Gets the current year. At worst, will do O(n) where n=number of variables headers.
             If this function is called many times, or more calculations are needed it might be best either save the value or use "GetGreorian" or "get_posix" methods. 

    Parameters: 
            None

    Throws:
            Never

    Returns:   
            int - Current year. -1 on error
     */
    int year();

    /*
    Function: get_gregorian
             Returns the boost::gregoiran structure for time calculations. Performs at worst O(n) where n=number of variable headers
             See <http://www.boost.org/doc/libs/1_38_0/doc/html/date_time/gregorian.html#date_time.gregorian.date_class>

    Parameters: 
            None

    Throws:
            std::runtime_error is there is no date variable in memory

    Returns:   
            boost::gregorian::date - boost gregorian date structure.
     */
    boost::gregorian::date get_gregorian();

    /*
    Function: get_posix
              Returns the boost::posix_time structure for time calculations. Performs at worst O(n) where n=number of variable headers
              See <http://www.boost.org/doc/libs/1_38_0/doc/html/date_time/posix_time.html#date_time.posix_time.ptime_class>

    Parameters: 
            None

    Throws:
            std::runtime_error is there is no date variable in memory

    Returns:   
            boost::posix_time::ptime - boost posix time structure
     */
    boost::posix_time::ptime get_posix();

    /*
    Function: get
    Gets the value associated with a variable. The routine will do its best to cast to type T using boost::lexical_cast<T>()

    Parameters: 
    std::string varName - Variable name. Must be the same name that's in the headers, read form the file.

    Throws:
    std::runtime_error - if it cannot be cast to type T

    Returns:   
    boost::T - Value of the variable at the current timestep
     */
    double get(std::string varName) ;
    
    void set(std::string varName, double value);


};
