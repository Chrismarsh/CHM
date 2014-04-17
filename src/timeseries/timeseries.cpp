

#include "timeseries.hpp"

void time_series::open(std::string path)
{
    std::fstream file(path.c_str());
    std::string line = "";

    int cols = 0;

    //tokenizer
    regex_tokenizer token;
    //contains the column headers
    std::vector<std::string> header;

    if (!file.is_open())
        BOOST_THROW_EXCEPTION(file_read_error()
            << boost::errinfo_errno(errno)
            << boost::errinfo_file_name(path));


    bool done = false;
    token.set_regex("[^,\\r\\n\\s]+"); //anything but whitespace or ,
    //read in the file, skip any blank lines at the top of the file
    while (!done)
    {
        getline(file, line);
        header = token.tokenize<std::string>(line);

        //skip any line that is a blank, or that has text in it
        if (header.size() != 0)
            done = true;

    }

    //take that the number of headers is how many columns there should be
    _cols = header.size();

    for (std::vector<std::string>::const_iterator itr = header.begin();
            itr != header.end();
            itr++)
    {
        ts_hashmap::const_accessor a;
        if (!_variables.insert(a, *itr))
        {
            BOOST_THROW_EXCEPTION(forcing_insertion_error()
                    << errstr_info(std::string("Failed to insert ") + *itr)
                    << boost::errinfo_file_name(path)
                    );

        }
    }



    //a string defined as anything with a letter in it or that has a special character in this list
    //~`!@#$%^&*(){[}]|\:;"'<>.?/
    //regex_tokenizer string("[\\x21\\x22\\x23\\x24\\x25\\x26\\x27\\x28\\x29\\x2A\\x2F\\x3A\\x3B\\x3C\\x3D\\x3E\\x3F\\x40\\x5B\\x5C\\x5D\\x5E\\x5F\\x60\\x7B\\x7C\\x7D\\x88\\x98]|.*[A-Za-z]+.*"); 
    regex_tokenizer string(".*"); //anything 
    regex_tokenizer integer("^[-+]?\\d+$");
    regex_tokenizer floating("^[-+]?(?:[0-9]+\\.(?:[0-9]*)?|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?$");
    regex_tokenizer dateTime("[0-9]{8}T[0-9]{6}"); //iso standard time

    std::vector<std::string> values; //values on a line


    token.set_regex("[^,\\r\\n\\s]+"); //anything but whitespace or ,

    int lines = 0;

    //this is used to save the name of the date header so we can find it to check the timestep later
    std::string dateHeader = "";

    while (getline(file, line))
    {
        lines++;
        //how many cols, make sure that equals the number of headers read in.
        unsigned int cols_so_far = 0;
        values = token.tokenize<std::string>(line);

        //make sure it isn't a blank line
        if (values.size() != 0)
        {
            //get the col name
            std::vector<std::string>::const_iterator headerItr = header.begin();

            int cols_so_far = 0;
            //for each column
            for (std::vector<std::string>::const_iterator itr = values.begin();
                    itr != values.end();
                    itr++)
            {
//                std::vector<std::string> strings;
//                std::vector<std::string> ints;
                std::vector<std::string> doubles;
                std::vector<std::string> dates;

                //check to see if it's an integer first
//                if ((ints = integer.tokenize<std::string>(*itr)).size() == 1)
//                {
//                    ts_hashmap::accessor a;
//                    if (!_variables.find(a, *headerItr))
//                        BOOST_THROW_EXCEPTION(forcing_lookup_error()
//                            << errstr_info(std::string("Failed to find ") + *headerItr)
//                            << boost::errinfo_file_name(path)
//                            );
//
//                    try
//                    {
//                        a->second.push_back(boost::lexical_cast<int>(ints[0]));
//                    } catch (...)
//                    {
//                        BOOST_THROW_EXCEPTION(forcing_badcast()
//                                << errstr_info("Failed to cast " + ints[0] + " to an int.")
//                                << boost::errinfo_file_name(path)
//                                );
//                    }
//
//
//                } else 
                if ((doubles = floating.tokenize<std::string>(*itr)).size() == 1)
                {
                    ts_hashmap::accessor a;
                    if (!_variables.find(a, *headerItr))
                        BOOST_THROW_EXCEPTION(forcing_lookup_error()
                            << errstr_info(std::string("Failed to find ") + *headerItr)
                            << boost::errinfo_file_name(path)
                            );
                    try
                    {
                        a->second.push_back(boost::lexical_cast<double>(doubles[0]));
                    } catch (...)
                    {
                        BOOST_THROW_EXCEPTION(forcing_badcast()
                                << errstr_info("Failed to cast " + doubles[0] + " to a double.")
                                << boost::errinfo_file_name(path)
                                );
                    }


                } else if ((dates = dateTime.tokenize<std::string>(*itr)).size() == 1)
                {
//                    ts_hashmap::accessor a;
//
//                    if (dateHeader == "")
//                        dateHeader = *headerItr;
//
//                    if (!_variables.find(a, *headerItr))
//                        BOOST_THROW_EXCEPTION(forcing_lookup_error()
//                            << errstr_info(std::string("Failed to find ") + *headerItr)
//                            << boost::errinfo_file_name(path)
//                            );
                    _date_vec.push_back(boost::posix_time::from_iso_string(dates[0]));
//                    a->second.push_back(boost::posix_time::from_iso_string(dates[0]));



                }                    //if we get here, it's none of the above, so assume it's a string
//                else if ((strings = string.tokenize<std::string>(*itr)).size() == 1)
//                {
//
//                    ts_hashmap::accessor a;
//                    if (!_variables.find(a, *headerItr))
//                        BOOST_THROW_EXCEPTION(forcing_lookup_error()
//                            << errstr_info(std::string("Failed to find ") + *headerItr)
//                            << boost::errinfo_file_name(path)
//                            );
//
//                    a->second.push_back(strings[0]);
//                }
                else
                {
                    //something has gone horribly wrong
                    BOOST_THROW_EXCEPTION(forcing_no_regexmatch()
                            << errstr_info("Unable to match any regex for " + *itr)
                            << boost::errinfo_file_name(path)
                            );
                }

                //next header
                headerItr++;
                cols_so_far++;

            }
            if (cols_so_far != _cols)
            {
                BOOST_THROW_EXCEPTION(forcing_badcast()
                        << errstr_info("Expected " + boost::lexical_cast<std::string>(_cols) + "lines")
                        << boost::errinfo_file_name(path)
                        );
            }
            _rows++;
        }

    } //end of file read

    _isOpen = true;
    _file = path;

    //check to make sure what we have read in makes sense
    //Check for:
    //	- Each col has the same number of rows
    //	- Time steps are equal

    //get iters for each variables
    std::string* headerItems = new std::string[_variables.size()];

    int i = 0;
    //build a list of all the headers
    //unknown order
    for (ts_hashmap::iterator itr = _variables.begin(); itr != _variables.end(); itr++)
    {
        headerItems[i++] = itr->first;
    }

    //get and save each accessor
    bool bfirst = true;
    ts_hashmap::const_accessor first;


    for (unsigned int l = 0; l < _variables.size(); l++)
    {
        if (bfirst)
        {
            _variables.find(first, headerItems[0]);
            bfirst = false;
        }

        ts_hashmap::const_accessor a;
        if (!_variables.find(a, headerItems[l]))
            BOOST_THROW_EXCEPTION(forcing_lookup_error()
                << errstr_info(std::string("Failed to find ") + headerItems[l])
                << boost::errinfo_file_name(path)
                );

        //check all cols are the same sizes
        if (first->second.size() != a->second.size())
            BOOST_THROW_EXCEPTION(forcing_lookup_error()
                << errstr_info("Row " + a->first + " is a different size.")
                << boost::errinfo_file_name(path));
    }

    first.release();

    delete[] headerItems;

}

std::string time_series::get_opened_file()
{
    return _file;
}

time_series::time_series()
{
    _cols = 0;
    _rows = 0;
    _isOpen = false;
}

time_series::~time_series()
{

}

void time_series::to_file(std::string file)
{
    std::ofstream out;
    out.open(file.c_str());

    if (!out.is_open())
        BOOST_THROW_EXCEPTION(file_read_error()
            << boost::errinfo_errno(errno)
            << boost::errinfo_file_name(file));

    std::string* headerItems = new std::string[_variables.size()];

    //build a list of all the headers
    //unknown order
    int i = 0;
    variable::const_iterator *tItr = new variable::const_iterator[_variables.size()];
    for (ts_hashmap::iterator itr = _variables.begin(); itr != _variables.end(); itr++)
    {
        headerItems[i] = itr->first;
        out << "\t" << itr->first;

        //save vector iterators
        tItr[i] = itr->second.begin();
        i++;
    }
    out << std::endl;


    for (int k = 0; k < _rows; k++)
    {
        for (int j = 0; j < _cols; j++)
        {
            out << "\t" << *(tItr[j]);
            tItr[j]++;
        }
        out << std::endl;
    }

    delete[] tItr;
    delete[] headerItems;
}

bool time_series::is_open()
{
    return _isOpen;

}

time_series::const_iterator time_series::begin()
{
    const_iterator step;
    //build a list of all the headers
    //unknown order

    for (ts_hashmap::iterator itr = _variables.begin(); itr != _variables.end(); itr++)
    {
        timestep::const_itr_map::accessor a;
        if (!step._currentStep._itrs.insert(a, itr->first))
        {
            BOOST_THROW_EXCEPTION(forcing_insertion_error()
                    << errstr_info("Failed to insert " + itr->first)
                    );
        }
        a->second = itr->second.begin();

    }

    return step;

}

time_series::const_iterator time_series::end()
{
    const_iterator step;
    //build a list of all the headers
    //unknown order

    for (ts_hashmap::iterator itr = _variables.begin(); itr != _variables.end(); itr++)
    {
        timestep::const_itr_map::accessor a;

        if (!step._currentStep._itrs.insert(a, itr->first))
        {
            BOOST_THROW_EXCEPTION(forcing_insertion_error()
                    << errstr_info("Failed to insert " + itr->first)
                    );
        }
        a->second = itr->second.end();

    }

    return step;
}


const timestep& time_series::const_iterator::dereference() const
{
    return _currentStep;
}

bool time_series::const_iterator::equal(const_iterator const& other) const
{
    bool isEqual = false;

    //different sizes
    if (_currentStep._itrs.size() != other._currentStep._itrs.size())
    {
        return false;
    }

    unsigned int i = 0;

    std::string *thisHeaders = new std::string[_currentStep._itrs.size()];
    std::string *otherHeaders = new std::string[_currentStep._itrs.size()];

    variable::const_iterator *thisIters = new variable::const_iterator[_currentStep._itrs.size()];
    variable::const_iterator *otherIters = new variable::const_iterator[_currentStep._itrs.size()];

    //this' headers and iters
    for (timestep::const_itr_map::const_iterator itr = _currentStep._itrs.begin(); itr != _currentStep._itrs.end(); itr++)
    {
        thisHeaders[i] = itr->first;
        thisIters[i] = itr->second; //iterator
        i++;
    }

    i = 0;
    //other's headers and iters
    for (timestep::const_itr_map::const_iterator itr = other._currentStep._itrs.begin(); itr != other._currentStep._itrs.end(); itr++)
    {
        otherHeaders[i] = itr->first; //key
        otherIters[i] = itr->second; //iterator
        i++;
    }

    for (i = 0; i < _currentStep._itrs.size(); i++)
    {
        for (unsigned int j = 0; j < _currentStep._itrs.size(); j++)
        {
            if (thisHeaders[i] == otherHeaders[j])
            {
                if (thisIters[i] == otherIters[j])
                    isEqual = true;
            }
        }
    }

    delete[] thisHeaders;
    delete[] otherHeaders;
    delete[] thisIters;
    delete[] otherIters;

    return isEqual;

}

void time_series::const_iterator::increment()
{
    //walks the map locking each node so that the increment can happen
    //walk order is not guaranteed
    unsigned int size = _currentStep._itrs.size();
    std::string *headers = new std::string[size];
    timestep::const_itr_map::accessor *accesors = new timestep::const_itr_map::accessor[size];
    int i = 0;

    for (timestep::const_itr_map::iterator itr = _currentStep._itrs.begin(); itr != _currentStep._itrs.end(); itr++)
    {
        _currentStep._itrs.find(accesors[i], itr->first);
        (accesors[i]->second)++;
        i++;
        //invalid cached values
        _currentStep._month_cache = -1;
    
    }

    delete[] headers;
    delete[] accesors;
    
   

}

void time_series::const_iterator::decrement()
{
    //walks the map locking each node so that the increment can happen
    //walk order is not guaranteed
    unsigned int size = _currentStep._itrs.size();
    std::string *headers = new std::string[size];
    timestep::const_itr_map::accessor *accesors = new timestep::const_itr_map::accessor[size];
    int i = 0;

    for (timestep::const_itr_map::iterator itr = _currentStep._itrs.begin(); itr != _currentStep._itrs.end(); itr++)
    {
        _currentStep._itrs.find(accesors[i], itr->first);
        (accesors[i]->second)--;
        i++;
        //invalid cached values
        _currentStep._month_cache = -1;
    }

    delete[] headers;
    delete[] accesors;

}

time_series::const_iterator::const_iterator()
{

}

time_series::const_iterator::const_iterator(const const_iterator& src)
{
    _currentStep = timestep(src._currentStep);
}

time_series::const_iterator::~const_iterator()
{

}

time_series::const_iterator& time_series::const_iterator::operator=(const time_series::const_iterator& rhs)
{
    if (this == &rhs)
        return (*this);
    _currentStep = timestep(rhs._currentStep);
    return *this;
}

