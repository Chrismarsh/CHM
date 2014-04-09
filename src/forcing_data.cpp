

#include "forcing_data.hpp"

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
    m_cols = header.size();

    for (std::vector<std::string>::const_iterator itr = header.begin();
            itr != header.end();
            itr++)
    {
        ObsTable::const_accessor a;
        if (!m_variables.insert(a, *itr))
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
                std::vector<std::string> strings;
                std::vector<std::string> ints;
                std::vector<std::string> floats;
                std::vector<std::string> dates;

                //check to see if it's an integer first
                if ((ints = integer.tokenize<std::string>(*itr)).size() == 1)
                {
                    ObsTable::accessor a;
                    if (!m_variables.find(a, *headerItr))
                        BOOST_THROW_EXCEPTION(forcing_lookup_error()
                            << errstr_info(std::string("Failed to find ") + *headerItr)
                            << boost::errinfo_file_name(path)
                            );

                    try
                    {
                        a->second.push_back(boost::lexical_cast<int>(ints[0]));
                    } catch (...)
                    {
                        BOOST_THROW_EXCEPTION(forcing_badcast()
                                << errstr_info("Failed to cast " + ints[0] + " to an int.")
                                << boost::errinfo_file_name(path)
                                );
                    }


                } else if ((floats = floating.tokenize<std::string>(*itr)).size() == 1)
                {
                    ObsTable::accessor a;
                    if (!m_variables.find(a, *headerItr))
                        BOOST_THROW_EXCEPTION(forcing_lookup_error()
                            << errstr_info(std::string("Failed to find ") + *headerItr)
                            << boost::errinfo_file_name(path)
                            );
                    try
                    {
                        a->second.push_back(boost::lexical_cast<float>(floats[0]));
                    } catch (...)
                    {
                        BOOST_THROW_EXCEPTION(forcing_badcast()
                                << errstr_info("Failed to cast " + floats[0] + " to an int.")
                                << boost::errinfo_file_name(path)
                                );
                    }


                } else if ((dates = dateTime.tokenize<std::string>(*itr)).size() == 1)
                {
                    ObsTable::accessor a;

                    if (dateHeader == "")
                        dateHeader = *headerItr;

                    if (!m_variables.find(a, *headerItr))
                        BOOST_THROW_EXCEPTION(forcing_lookup_error()
                            << errstr_info(std::string("Failed to find ") + *headerItr)
                            << boost::errinfo_file_name(path)
                            );

                    a->second.push_back(boost::posix_time::from_iso_string(dates[0]));



                }                    //if we get here, it's none of the above, so assume it's a string
                else if ((strings = string.tokenize<std::string>(*itr)).size() == 1)
                {

                    ObsTable::accessor a;
                    if (!m_variables.find(a, *headerItr))
                        BOOST_THROW_EXCEPTION(forcing_lookup_error()
                            << errstr_info(std::string("Failed to find ") + *headerItr)
                            << boost::errinfo_file_name(path)
                            );

                    a->second.push_back(strings[0]);
                }
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
            if (cols_so_far != m_cols)
            {
                BOOST_THROW_EXCEPTION(forcing_badcast()
                        << errstr_info("Expected " + boost::lexical_cast<std::string>(m_cols) + "lines")
                        << boost::errinfo_file_name(path)
                        );
            }
            m_rows++;
        }

    } //end of file read

    m_isOpen = true;
    _file = path;

    //check to make sure what we have read in makes sense
    //Check for:
    //	- Each col has the same number of rows
    //	- Time steps are equal

    //get iters for each variables
    std::string* headerItems = new std::string[m_variables.size()];

    int i = 0;
    //build a list of all the headers
    //unknown order
    for (ObsTable::iterator itr = m_variables.begin(); itr != m_variables.end(); itr++)
    {
        headerItems[i++] = itr->first;
    }

    //get and save each accessor
    bool bfirst = true;
    ObsTable::const_accessor first;


    for (unsigned int l = 0; l < m_variables.size(); l++)
    {
        if (bfirst)
        {
            m_variables.find(first, headerItems[0]);
            bfirst = false;
        }

        ObsTable::const_accessor a;
        if (!m_variables.find(a, headerItems[l]))
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
    m_cols = 0;
    m_rows = 0;
    m_isOpen = false;
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

    std::string* headerItems = new std::string[m_variables.size()];

    //build a list of all the headers
    //unknown order
    int i = 0;
    Variable::const_iterator *tItr = new Variable::const_iterator[m_variables.size()];
    for (ObsTable::iterator itr = m_variables.begin(); itr != m_variables.end(); itr++)
    {
        headerItems[i] = itr->first;
        out << "\t" << itr->first;

        //save vector iterators
        tItr[i] = itr->second.begin();
        i++;
    }
    out << std::endl;


    for (int k = 0; k < m_rows; k++)
    {
        for (int j = 0; j < m_cols; j++)
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
    return m_isOpen;

}

time_series::const_iterator time_series::begin()
{
    const_iterator step;
    //build a list of all the headers
    //unknown order

    for (ObsTable::iterator itr = m_variables.begin(); itr != m_variables.end(); itr++)
    {
        timestep::ConstItrMap::accessor a;
        if (!step.m_currentStep.m_itrs.insert(a, itr->first))
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

    for (ObsTable::iterator itr = m_variables.begin(); itr != m_variables.end(); itr++)
    {
        timestep::ConstItrMap::accessor a;

        if (!step.m_currentStep.m_itrs.insert(a, itr->first))
        {
            BOOST_THROW_EXCEPTION(forcing_insertion_error()
                    << errstr_info("Failed to insert " + itr->first)
                    );
        }
        a->second = itr->second.end();

    }

    return step;
}

size_t time_series::HashCompare::hash(const std::string& x)
{
    boost::crc_32_type crc32;
    std::string xlower = boost::algorithm::to_lower_copy<std::string>(x);
    crc32.process_bytes(xlower.c_str(), xlower.length());

    return crc32.checksum();
}

bool time_series::HashCompare::equal(const std::string& s1, const std::string& s2)
{
    std::string ss1 = boost::algorithm::to_lower_copy<std::string>(s1);
    std::string ss2 = boost::algorithm::to_lower_copy<std::string>(s2);

    return ss1 == ss2;
}

const time_series::timestep& time_series::const_iterator::dereference() const
{
    return m_currentStep;
}

bool time_series::const_iterator::equal(const_iterator const& other) const
{
    bool isEqual = false;

    //different sizes
    if (m_currentStep.m_itrs.size() != other.m_currentStep.m_itrs.size())
    {
        return false;
    }

    unsigned int i = 0;

    std::string *thisHeaders = new std::string[m_currentStep.m_itrs.size()];
    std::string *otherHeaders = new std::string[m_currentStep.m_itrs.size()];

    Variable::const_iterator *thisIters = new Variable::const_iterator[m_currentStep.m_itrs.size()];
    Variable::const_iterator *otherIters = new Variable::const_iterator[m_currentStep.m_itrs.size()];

    //this' headers and iters
    for (timestep::ConstItrMap::const_iterator itr = m_currentStep.m_itrs.begin(); itr != m_currentStep.m_itrs.end(); itr++)
    {
        thisHeaders[i] = itr->first;
        thisIters[i] = itr->second; //iterator
        i++;
    }

    i = 0;
    //other's headers and iters
    for (timestep::ConstItrMap::const_iterator itr = other.m_currentStep.m_itrs.begin(); itr != other.m_currentStep.m_itrs.end(); itr++)
    {
        otherHeaders[i] = itr->first; //key
        otherIters[i] = itr->second; //iterator
        i++;
    }

    for (i = 0; i < m_currentStep.m_itrs.size(); i++)
    {
        for (unsigned int j = 0; j < m_currentStep.m_itrs.size(); j++)
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
    unsigned int size = m_currentStep.m_itrs.size();
    std::string *headers = new std::string[size];
    timestep::ConstItrMap::accessor *accesors = new timestep::ConstItrMap::accessor[size];
    int i = 0;

    for (timestep::ConstItrMap::iterator itr = m_currentStep.m_itrs.begin(); itr != m_currentStep.m_itrs.end(); itr++)
    {
        m_currentStep.m_itrs.find(accesors[i], itr->first);
        (accesors[i]->second)++;
        i++;
    }

    delete[] headers;
    delete[] accesors;

}

void time_series::const_iterator::decrement()
{
    //walks the map locking each node so that the increment can happen
    //walk order is not guaranteed
    unsigned int size = m_currentStep.m_itrs.size();
    std::string *headers = new std::string[size];
    timestep::ConstItrMap::accessor *accesors = new timestep::ConstItrMap::accessor[size];
    int i = 0;

    for (timestep::ConstItrMap::iterator itr = m_currentStep.m_itrs.begin(); itr != m_currentStep.m_itrs.end(); itr++)
    {
        m_currentStep.m_itrs.find(accesors[i], itr->first);
        (accesors[i]->second)--;
        i++;
    }

    delete[] headers;
    delete[] accesors;
}

time_series::const_iterator::const_iterator()
{

}

time_series::const_iterator::const_iterator(const const_iterator& src)
{
    m_currentStep = timestep(src.m_currentStep);
}

time_series::const_iterator::~const_iterator()
{

}

time_series::const_iterator& time_series::const_iterator::operator=(const time_series::const_iterator& rhs)
{
    if (this == &rhs)
        return (*this);
    m_currentStep = timestep(rhs.m_currentStep);
    return *this;
}

time_series::timestep::timestep(const timestep& src)
{
    m_itrs = ConstItrMap(src.m_itrs);
}

time_series::timestep::timestep()
{

}

time_series::timestep::~timestep()
{

}

std::string time_series::timestep::to_string()
{
    std::string s;

    for (time_series::timestep::ConstItrMap::const_iterator itr = m_itrs.begin(); itr != m_itrs.end(); itr++)
    {
        s += boost::lexical_cast<std::string>(*(itr->second)) + std::string("\t");
    }

    return s;
}

int time_series::timestep::month()
{
    int d = -1;
    boost::posix_time::ptime time;
    boost::gregorian::date date;


    for (time_series::timestep::ConstItrMap::const_iterator itr = m_itrs.begin(); itr != m_itrs.end(); itr++)
    {
        try
        {
            time = boost::get<boost::posix_time::ptime>(*(itr->second));
            date = boost::gregorian::from_string(boost::lexical_cast<std::string>(time.date()));
            return date.month();
        } catch (boost::bad_get e)
        {
            //keep going
        }

    }

    return d;
}

int time_series::timestep::day()
{
    int d = -1;
    boost::posix_time::ptime time;
    boost::gregorian::date date;


    for (time_series::timestep::ConstItrMap::const_iterator itr = m_itrs.begin(); itr != m_itrs.end(); itr++)
    {
        try
        {
            time = boost::get<boost::posix_time::ptime>(*(itr->second));
            date = boost::gregorian::from_string(boost::lexical_cast<std::string>(time.date()));
            return date.day();
        } catch (boost::bad_get e)
        {
            //keep going
        }

    }

    return d;
}

int time_series::timestep::year()
{
    int d = -1;
    boost::posix_time::ptime time;
    boost::gregorian::date date;


    for (time_series::timestep::ConstItrMap::const_iterator itr = m_itrs.begin(); itr != m_itrs.end(); itr++)
    {
        try
        {
            time = boost::get<boost::posix_time::ptime>(*(itr->second));
            date = boost::gregorian::from_string(boost::lexical_cast<std::string>(time.date()));
            return date.year();
        } catch (boost::bad_get e)
        {
            //keep going
        }

    }

    return d;
}

boost::gregorian::date time_series::timestep::get_gregorian()
{
    boost::posix_time::ptime time;
    boost::gregorian::date date;

    for (time_series::timestep::ConstItrMap::const_iterator itr = m_itrs.begin(); itr != m_itrs.end(); itr++)
    {
        try
        {
            time = boost::get<boost::posix_time::ptime>(*(itr->second));
            date = boost::gregorian::from_string(boost::lexical_cast<std::string>(time.date()));
            return date;
        } catch (boost::bad_get e)
        {
            //keep going
        }

    }

    //if we get this far, didn't find it, so bail
    BOOST_THROW_EXCEPTION(forcing_error()
            << errstr_info("No date variable found")
            );

}

boost::posix_time::ptime time_series::timestep::get_posix()
{
    boost::posix_time::ptime time;


    for (time_series::timestep::ConstItrMap::const_iterator itr = m_itrs.begin(); itr != m_itrs.end(); itr++)
    {
        try
        {
            time = boost::get<boost::posix_time::ptime>(*(itr->second));
            return time;
        } catch (boost::bad_get e)
        {
            //keep going
        }

    }

    //if we get this far, didn't find it, so bail
    BOOST_THROW_EXCEPTION(forcing_error()
            << errstr_info("No date variable found")
            );

}
