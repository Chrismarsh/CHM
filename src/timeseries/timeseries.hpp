#pragma once


#include <string>
#include <fstream>
#include <vector>
#include <ctime>
#include <algorithm>

#include <boost/date_time/posix_time/posix_time.hpp> // for boost::posix

#include <boost/variant.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/utility.hpp>
#include <boost/tuple/tuple.hpp>


#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_vector.h>

#include "regex_tokenizer.hpp"
#include "exception.hpp"
#include "logger.hpp"
#include "timestep.hpp"
#include "crc_hash_compare.hpp"


    
/**
\class timeseries
\brief Holds the meterological data.

This class holds meterological data. Each variable name is a string which acts as a key into a map, with the map entry containing a tbb::concurrent_vector allowing for parallel usage.

    "var1"        |     "var2"     |     "var3"   |
    ------------------------------------------------
    [...]         |      [...]     |      [...]   |
    vector 1      |     vector 2   |     vector 3 |
    [...]         |      [...]     |      [...]   |

 */
class timeseries : boost::noncopyable
{

public:
    
    //mesh elements need to see these
    //two different types: boost::variant solves this, but is very slow 
    // needs to be either a boost::posix_time or double, and is almost always a double
    typedef tbb::concurrent_vector< double > variable_vec;
    typedef tbb::concurrent_vector< boost::posix_time::ptime  > date_vec; 
    
    class iterator;
    

    timeseries();
 
    ~timeseries();


    /**
    * Return the timeseries for the given variable.
    * \param variable Variable name
    * \return A vector of this variable for the entire duration
    */
    variable_vec get_time_series(std::string variable);

    /**
    * Returns the datetime series
    * \return A boost::ptime vector
    */
    date_vec get_date_timeseries();

    /**
    * Returns a list of all the variable names in this timeseries
    * \return A vector of variable names
    */
    std::vector<std::string> list_variables();

    /**
    * Returns the length (number of elements) of the timeseries.
    */
    int get_timeseries_length();

    /**
    * Initializes an empty timeseries with the given variables and the given date timeseries
    * \param variables Set of variable names
    * \param datetime A complete datetime vector
    */
    void init(std::set<std::string> variables, date_vec datetime);
    

    /**
    * Returns an iterator at the start of the observations. This iterator may then be used to access a given variable.
    * \return An iterator
    */
    iterator begin();


    /**
    * Returns an iterator of one past the end of the observations
    * \return An iterator
    */
    iterator end();

    /**
    *  Opens an observation file. An observation file is organized in a tab, ",", or space delimited  columns, with
    each column representing an independent observation, and each row is a timesteps measurement.

    For example:

        Date			    Rh	Tair	Precip
        20080220T000000		50	-12		2
        20080221T000015		40	-10		0


    Some restrictions:
        - No more than 2147483647 steps. At 1s intervals, this equates to roughly 68 years.
        - Consistent units. You mustn't have mm on one line, then meters on the next, for the same observation
        - Has to be on a constant time step. The first interval is taken as the interval for the rest of the file
        - Missing values are not currently allowed - that is, each row must be complete with n entries where n is number of variables.
        - Whitespace, tab or comma delimited. Allows for mixed usage. ex 1234, 4543 890 is legal
        - Values must be numeric

    Integer styles:

         +1234
         -1234
         1234567890

    Float ing point:

         12.34
         12.
         .34
         12.345
         1234.45
         +12.34
         -12.34
         +1234.567e-89
         -1234.567e89

    Time:
        - Must be in one column in the following ISO 8601 date time form:
        - YYYYMMDDThhmmss   e.g., 20080131T235959
    \param path Fully qualified path
    */
    void open(std::string path);

    /**
    *  Writes the timeseries to file. Order of variable output not deterministic.
    *  \param file Full qualified path
    */
    void to_file(std::string file);

    /**
    * Returns a pair of iterators point to the start and end of the specified range. If the range is not found, a forcing_error exception is thrown.
    * \param start_time Start of range
    * \param end_time End of range
    * \return Pair of iterators that point to the start and end of the time series, inclusive.
    */
    boost::tuple<timeseries::iterator, timeseries::iterator> range(boost::posix_time::ptime start_time,boost::posix_time::ptime end_time);

    /**
    * Returns an iterator to the requested time.
    * \param time Time
    * \return iterator accessing the timeseries at the requested time
    */
    iterator find(boost::posix_time::ptime time);

    /**
    * Determines if the last opened file was successfully opened.
    * \return True on success
    */
    bool is_open();

    /**
    * Returns the last opened file name
    */
    std::string get_opened_file();

    /**
    * Calculates the miniumum value in a range [start,end] for the given variable
    * \param start Start iterator
    * \param end End iterator
    * \return min value
    */
    double range_min(timeseries::iterator& start, timeseries::iterator& end, std::string variable);

    /**
    * Calculates the maximum value in a range [start,end] for the given variable
    * \param start Start iterator
    * \param end End iterator
    * \return max value
    */
    double range_max(timeseries::iterator& start, timeseries::iterator& end, std::string variable);
    
private:
    typedef tbb::concurrent_hash_map<std::string, variable_vec, crc_hash_compare> ts_hashmap;


    // This is a hashmap interface, vector back end
    // "var1"        |     "var2"     |     "var3"   |      
    // ------------------------------------------------      
    //      [...]    |      [...]     |      [...]   |
    //    vector 1   |     vector 2   |     vector 3 |
    //      [...]    |      [...]     |      [...]   |
    ts_hashmap _variables;
    date_vec _date_vec;
    
    size_t _cols;
    size_t _rows;
    bool _isOpen;
    std::string _file;
    size_t _timeseries_length;

    
    //pushes variables back, only useful for reading from a file
    void push_back(double data, std::string variable);

    

};


/* Class: iterator
 Used to iterate over the timeseries instance.
 Thread safe.
 Dereference returns a timestep object.
  */
 class timeseries::iterator : public boost::iterator_facade<
                         timeseries::iterator,
                         timestep,
                         boost::bidirectional_traversal_tag> 
 {
 public:
    iterator();
    iterator(const iterator& src);
    ~iterator();
    iterator& operator=(const iterator& rhs);
 private:
     //the following satisfies the reqs for a boost::facade bidirectional iterator
     friend class boost::iterator_core_access;
     friend class timeseries;

     timestep& dereference() const;
     bool equal(iterator const& other) const;
     void increment();
     void decrement();
     void advance(timeseries::iterator::difference_type N);
     std::ptrdiff_t distance_to(iterator const& other) const;

     //iterators for the current step
     timestep* _currentStep;

 };


    