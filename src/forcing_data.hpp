#pragma once


#include <string>
#include <fstream>
#include <vector>


#include <boost/date_time/posix_time/posix_time.hpp> // for boost::posix
#include <boost/crc.hpp>      // for boost::crc_basic, boost::crc_optimal
#include <boost/cstdint.hpp>  // for boost::uint16_t
#include <boost/variant.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/utility.hpp>


#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_vector.h>

#include "regex_tokenizer.hpp"
#include "exception.hpp"

/*
        Class: forcing_data
                Holds the meterological data.
                forcing_data classes must currently not be copied or assigned from other instances.
*/
class time_series : boost::noncopyable
{
private:

        // Used to compare the hashes of the main data structure
        //not case sensitive
        class HashCompare
        {
        public:
                static size_t hash( const std::string& x);
                static bool equal(const std::string& s1, const std::string& s2);
        };

        //typedefs must go here after hashcompare decl

        typedef boost::variant< float,int,std::string, boost::posix_time::ptime > Datatypes;
        typedef tbb::concurrent_vector< Datatypes > Variable;
        typedef tbb::concurrent_hash_map<std::string, Variable ,HashCompare> ObsTable;


        //This is a hashmap interface, vector back end
        // "var1"        |     "var2"     |     "var3"   |      hashmap< std::string, tbb::concurrent_vector >
        // ------------------------------------------------      
        //      [...]    |      [...]     |      [...]   |
        //    vector 1   |     vector 2   |     vector 3 |
        //      [...]    |      [...]     |      [...]   |
        ObsTable m_variables;

        int m_cols;
        int m_rows;
        bool m_isOpen;
        std::string _file;

public:
        /*
        Class: timestep
        A given time step. Returned from dereferencing a forcing_data::const_iterator. Not to be directly instantiated. 
        Cannot directly modify the observation data

        Example:
        >forcing_data::const_iterator itr;
        >[...]
        >(*itr).get<std::string>("Date")
        */
        class timestep
        {
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
                template<typename T>
                T get(std::string varName)
                {
                        ConstItrMap::const_accessor a;
                        T out;

                        if(!m_itrs.find(a,varName))
                                throw std::runtime_error("Variable " + varName + " does not exist.");

                        try
                        {
                                out = boost::get<T>((a->second)[0]);
                        }
                        catch (boost::bad_get e)
                        {
                                //failed, *really* try to cast
                                //normally needed to cast out the date to a string
                                try
                                {
                                        out = boost::lexical_cast<T>((a->second)[0]);
                                }
                                catch (boost::bad_lexical_cast e)
                                {
                                        throw std::runtime_error("Variable " + varName + " could not be cast correctly");
                                }
                        }

                        return out;
                }


        private:
                friend class time_series;
                typedef tbb::concurrent_hash_map<std::string, Variable::const_iterator ,HashCompare> ConstItrMap;

                ConstItrMap m_itrs;//holds the iterators
        };



        class const_iterator;
        time_series();
        ~time_series();

        /*
        Function: begin
                 Returns an iterator at the start of the observations

        Parameters: 
                None

        Throws:
                Never

        Returns:   
                forcing_data::const_iterator - Iterator positioned at the begining of the observation
        */
        const_iterator begin();


        /*
        Function: end
                 Returns an iterator of one past the end of the observations

        Parameters: 
                None

        Throws:
                Never

        Returns:   
                forcing_data::const_iterator - Iterator positioned one past the end of the observations
        */
        const_iterator end();

        /*
        Function: open
                 Opens an observation file. An observation file is organized in a tab, ",", or space delimited  columns, with
                each column representing an independent observation, and each row is a timesteps measurement. For example:
                >Date					Rh	Tair	Precip	Notes
                >20080220T000000		50	-12		2		Station_1
                >20080221T000015		40	-10		0		Station_1
                >		[...]
                Some restrictions:
                        - No more than 2147483647 steps. At 1s intervals, this equates to roughly 68 years.
                        - Consistent units. You mustn't have mm on one line, then meters on the next, for the same observation
                        - Has to be on a constant time step. The first interval is taken as the interval for the rest of the file
                        - Missing values are not currently allowed - that is, each row must be complete with n entries where n is number of variables.
                        - Whitespace, tab or comma delimited. Allows for mixed usage. ex 1234, 4543 890 is legal
                        - Values can be one of
                                - String (special characters are allowed, no commas)
                                        Will match any thing with:
                                        >~`!@#$%^&*(){[}]|\:;"'<>.?/
                                        >OR has upper case / lower case letters
                                - Integer
                                        Will match the following styles:
                                        >+1234
                                        >-1234
                                        >1234567890
                                - Floating point
                                        Will match the following styles:
                                        >12.34
                                        >12.
                                        >.34
                                        >12.345
                                        >1234.45
                                        >+12.34
                                        >-12.34
                                        >+1234.567e-89
                                        >-1234.567e89
                                - Time
                                        - Must be in one column in the following ISO 8601 date time form:
                                                - 20080131T235959
                                                  YYYYMMDDThhmmss


        Parameters: 
                std::string path - Fully quantified path to the file

        Throws:
                std::runtime_error - On error

        Returns:   
                 void 
        */
        void open(std::string path);

        /*
        Function: to_file
                 Writes the internal observation structure to file. Order of columns is not preserved due to how the data is stored.

        Parameters: 
                std::string file - File to save to

        Throws:
                std::runtime_error On error

        Returns:   
                void 
        */
        void to_file(std::string file);



        /*
        Function: is_open
                 Determines if a file was successfully opened

        Parameters: 

        Throws:
                Never

        Returns:   
                bool - True if opened
        */
        bool is_open();
        
        std::string get_opened_file();


        /*
        Class: const_iterator
        Used to iterate over a forcing_data instance.
        Thread safe.
        Dereference returns a timestep object.
        Example:
        >forcing_data obs;
        >obs.open("Observations.txt");
        >
        >forcing_data::const_iterator itr;
        >
        >for(itr=obs.begin();itr != obs.end(); itr++)
        >{	
        >std::cout << std::endl << (*itr).to_string() << std::endl;		
        >}
        See also:
                <timestep>
        */
        class const_iterator : public boost::iterator_facade<
                const_iterator,
                Variable,
                boost::random_access_traversal_tag,
                timestep>
        {
        public:
                const_iterator();
                const_iterator(const const_iterator& src);
                ~const_iterator();
                const_iterator& operator=(const const_iterator& rhs);
        private:
                //the following satisfies the reqs for a boost::facade bidirectional iterator
                friend class boost::iterator_core_access;
                friend class time_series;

                const time_series::timestep& dereference() const;
                bool equal(const_iterator const& other) const;
                void increment();
                void decrement();

                //iterators for the current step
                time_series::timestep m_currentStep;

        };




};
