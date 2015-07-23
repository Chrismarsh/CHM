

#include "timeseries.hpp"

void timeseries::push_back(double data, std::string variable)
{
    ts_hashmap::accessor a;
    _variables.insert(a, variable);
    a->second.push_back(data);
}
void timeseries::init(std::set<std::string> variables, date_vec datetime)
{
   for (auto& v: variables)
   {
       ts_hashmap::accessor a;
       _variables.insert(a,v);
   }

   //preallocate all the memory required
   for (auto& itr : _variables)
   {
        size_t size = datetime.size();
        itr.second.resize(size);
//       for(int i = 0; i < size; i++)
//            itr.second.push_back(3.14159);
   }

   //setup date vector
   _date_vec = datetime;
}

 timeseries::date_vec timeseries::get_date_timeseries()
 {
     return _date_vec;
 }
std::vector<std::string> timeseries::list_variables()
{
    std::vector<std::string> vars;
    for(auto& itr : _variables)
    {
        vars.push_back(itr.first);
    }
    
    return vars;
}
timeseries::variable_vec timeseries::get_time_series(std::string variable)
{
    ts_hashmap::const_accessor a;
    if(!_variables.find(a, variable))
    {
        BOOST_THROW_EXCEPTION(forcing_error()
                                << errstr_info("Unable to find " + variable));
    }   
    return a->second;
}

boost::tuple<timeseries::iterator, timeseries::iterator> timeseries::range(boost::posix_time::ptime start_time,boost::posix_time::ptime end_time)
{
    //look for our requested timestep
    auto itr_find = std::find(_date_vec.begin(),_date_vec.end(),start_time);
    if ( itr_find == _date_vec.end())
    {
        BOOST_THROW_EXCEPTION(forcing_timestep_notfound()
                                << errstr_info("Timestep not found"));
    }
    
    //Find the first one
    //get offset from iterator
    int dist_start = std::distance(_date_vec.begin(), itr_find);
    
    iterator start_step;

    //iterate over the map of vectors and build a list of all the variable names
    //unknown order
    for (auto& itr : _variables)
    {
        //itr_map is holding the iterators into each vector
        timestep::itr_map::accessor a;
        //create the keyname for this variable and store the iterator
        if (!start_step._currentStep->_itrs.insert(a, itr.first))
        {
            BOOST_THROW_EXCEPTION(forcing_error()
                    << errstr_info("Failed to insert " + itr.first)
                    );
        }
        
        //insert the iterator
        a->second = itr.second.begin()+dist_start;
    }

    //set the date vector to be the begining of the internal data vector
    start_step._currentStep->_date_itr = _date_vec.begin()+dist_start;
    
    
    
    //ok we can cheat and start from where we currently are instead of two straight calls to find
    itr_find = std::find(_date_vec.begin()+dist_start,_date_vec.end(),end_time);
    //get offset from iterator
    int dist_end = std::distance(_date_vec.begin(), itr_find);
    ++dist_end; //get 1 past where we are going
    iterator end_step;

    //iterate over the map of vectors and build a list of all the variable names
    //unknown order
    for (auto& itr : _variables)
    {
        //itr_map is holding the iterators into each vector
        timestep::itr_map::accessor a;
        //create the keyname for this variable and store the iterator
        if (!end_step._currentStep->_itrs.insert(a, itr.first))
        {
            BOOST_THROW_EXCEPTION(forcing_error()
                    << errstr_info("Failed to insert " + itr.first)
                    );
        }
        
        //insert the iterator
        a->second = itr.second.begin()+dist_end;
    }

    //set the date vector to be the begining of the internal data vector
    end_step._currentStep->_date_itr = _date_vec.begin()+dist_end;
    
    return boost::tuple<timeseries::iterator, timeseries::iterator>(start_step,end_step);
    
    
}
timeseries::iterator timeseries::find(boost::posix_time::ptime time)
{
    //look for our requested timestep
    auto itr = std::find(_date_vec.begin(),_date_vec.end(),time);
    if ( itr == _date_vec.end())
    {
        BOOST_THROW_EXCEPTION(forcing_timestep_notfound()
                                << errstr_info("Timestep not found"));
    }
    
    //get offset from iterator
    int dist = std::distance(_date_vec.begin(), itr);
    
    iterator step;

    //iterate over the map of vectors and build a list of all the variable names
    //unknown order
    for (auto& itr : _variables)
    {
        //itr_map is holding the iterators into each vector
        timestep::itr_map::accessor a;
        //create the keyname for this variable and store the iterator
        if (!step._currentStep->_itrs.insert(a, itr.first))
        {
            BOOST_THROW_EXCEPTION(forcing_insertion_error()
                    << errstr_info("Failed to insert " + itr.first)
                    );
        }
        
        //insert the iterator
        a->second = itr.second.begin()+dist;
    }

    //set the date vector to be the begining of the internal data vector
    step._currentStep->_date_itr = _date_vec.begin()+dist;
    
    return step;
}

void timeseries::open(std::string path)
{
    std::fstream file(path.c_str());
    std::string line = "";

    //tokenizer
    regex_tokenizer token;
    //contains the column headers
    std::vector<std::string> header;

    if (!file.is_open())
        BOOST_THROW_EXCEPTION(file_read_error()
            << boost::errinfo_errno(errno)
            << boost::errinfo_file_name(path));

    LOG_DEBUG << "Parsing file " + path;
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
    regex_tokenizer floating("^[-+]?(?:[0-9]+\\.?(?:[0-9]*)?|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?$");
    regex_tokenizer dateTime("[0-9]{8}T[0-9]{6}"); //iso standard time

    std::vector<std::string> values; //values on a line


    token.set_regex("[^,\\r\\n\\s]+"); //anything but whitespace or ,

    int lines = 0;

    //this is used to save the name of the date header so we can find it to check the timestep later
    std::string dateHeader = "";

    while (getline(file, line))
    {
        lines++;
        

        values = token.tokenize<std::string>(line);

        //make sure it isn't a blank line
        if (values.size() != 0)
        {
            //get the col name
            std::vector<std::string>::const_iterator headerItr = header.begin();

            //how many cols, make sure that equals the number of headers read in.
            size_t cols_so_far = 0;
            //for each column
            for (std::vector<std::string>::const_iterator itr = values.begin();
                    itr != values.end();
                    itr++)
            {
                std::vector<std::string> doubles;
                std::vector<std::string> dates;

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
                        LOG_VERBOSE << "Found " << *headerItr << ": " << doubles[0];
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

                    LOG_VERBOSE << "Found " << *headerItr << ": " << dates[0];
                    _date_vec.push_back(boost::posix_time::from_iso_string(dates[0])); //from_iso_string
                    
                    //now we know where the date colum is, we remove it from the hashmap if we haven't already
                    ts_hashmap::accessor a;
                    if (_variables.find(a, *headerItr))
                        _variables.erase(a);
                  
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
    _timeseries_length = lines;

    //check to make sure what we have read in makes sense
    //Check for:
    //	- Each col has the same number of rows
    //	- Time steps are equal

    //get iters for each variables
    LOG_VERBOSE << "Read in " << _variables.size() << " variables";
    std::string* headerItems = new std::string[_variables.size()];

    int i = 0;
    //build a list of all the headers
    //unknown order
    for (ts_hashmap::iterator itr = _variables.begin(); itr != _variables.end(); itr++)
    {
        LOG_VERBOSE << itr->first;
        headerItems[i++] = itr->first;
    }

    //get and save each accessor
    size_t d_length = _date_vec.size();


    for (unsigned int l = 0; l < _variables.size(); l++)
    {
        //compare all columns to date length

        ts_hashmap::const_accessor a;
        if (!_variables.find(a, headerItems[l]))
            BOOST_THROW_EXCEPTION(forcing_lookup_error()
                << errstr_info(std::string("Failed to find ") + headerItems[l])
                << boost::errinfo_file_name(path)
                );

        //check all cols are the same size as the first col
        LOG_VERBOSE << "Column " + headerItems[l] + " length=" + boost::lexical_cast<std::string>( a->second.size()), + "expected=" + boost::lexical_cast<std::string>(d_length);
        if (d_length != a->second.size())
        {
            LOG_ERROR << "Col " + headerItems[l] + " is a different size. Expected size="+boost::lexical_cast<std::string>(d_length);
            BOOST_THROW_EXCEPTION(forcing_lookup_error()
                << errstr_info("Col " + headerItems[l] + " is a different size. Expected size="+boost::lexical_cast<std::string>(d_length))
                << boost::errinfo_file_name(path));
        }
        
    }

   

    delete[] headerItems;

}

int timeseries::get_timeseries_length()
{
    return _timeseries_length;
}

std::string timeseries::get_opened_file()
{
    return _file;
}

timeseries::timeseries()
{
    _cols = 0;
    _rows = 0;
    _isOpen = false;
    _timeseries_length=0;
}


timeseries::~timeseries()
{

}

void timeseries::to_file(std::string file)
{
    std::ofstream out;
    out.open(file.c_str());
//    out << std::fixed << std::setprecision(8);
    if (!out.is_open())
        BOOST_THROW_EXCEPTION(file_read_error()
            << boost::errinfo_errno(errno)
            << boost::errinfo_file_name(file));

    
    std::string* headerItems = new std::string[_variables.size()];

    //build a list of all the headers
    //unknown order
    int i = 0;
    out << "Date\t";
    variable_vec::const_iterator *tItr = new variable_vec::const_iterator[_variables.size()];
    for (ts_hashmap::iterator itr = _variables.begin(); itr != _variables.end(); itr++)
    {
        headerItems[i] = itr->first;
        out << "\t" << itr->first;

        //save vector iterators
        tItr[i] = itr->second.begin();
        _rows = itr->second.size();
        i++;
    }
    out << std::endl;

    
    for (size_t k = 0; k < _rows; k++)
    {
        out << boost::posix_time::to_iso_string(_date_vec.at(k)) << "\t";
//        out << _date_vec.at(k) << "\t";
        for (size_t j = 0; j < _variables.size(); j++)
        {
            out << "\t" << *(tItr[j]);
            tItr[j]++;
        }
        out << std::endl;
    }

    delete[] tItr;
    delete[] headerItems;
}

bool timeseries::is_open()
{
    return _isOpen;

}

double timeseries::range_max(timeseries::iterator& start, timeseries::iterator& end, std::string variable )
{
    auto m = std::max_element(start->get_itr(variable),++end->get_itr(variable));  //because _element is [first,last)
    return *m;
}

double timeseries::range_min(timeseries::iterator& start, timeseries::iterator& end, std::string variable )
{
    auto m = std::min_element(start->get_itr(variable),++end->get_itr(variable)); //because _element is [first,last)
    return *m;
}


//iterator implementation
//------------------------
timeseries::iterator timeseries::begin()
{
    iterator step;

    //iterate over the map of vectors and build a list of all the variable names
    //unknown order
    for (auto& itr : _variables)
    {
        //itr_map is holding the iterators into each vector
        timestep::itr_map::accessor a;
        //create the keyname for this variable and store the iterator
        if (!step._currentStep->_itrs.insert(a, itr.first))
        {
            BOOST_THROW_EXCEPTION(forcing_insertion_error()
                    << errstr_info("Failed to insert " + itr.first)
                    );
        }
        
        //insert the iterator
        a->second = itr.second.begin();
    }

    //set the date vector to be the begining of the internal data vector
    step._currentStep->_date_itr = _date_vec.begin();

//    for (auto& itr : _variables)
//   {
//       LOG_DEBUG << itr.first << ":";
//       for(auto& jtr : itr.second)
//       {
//           LOG_DEBUG << boost::lexical_cast<std::string>(jtr);
//       }
//       
//   }
   
//    LOG_DEBUG << step->get("t");
    return step;

}


timeseries::iterator timeseries::end()
{
    iterator step;
    //loop over the variable map and save the iterator to the end
    //unknown order that'll get the variables in.
    for (auto& itr : _variables)
    {
        timestep::itr_map::accessor a;

        if (!step._currentStep->_itrs.insert(a, itr.first))
        {
            BOOST_THROW_EXCEPTION(forcing_insertion_error()
                    << errstr_info("Failed to insert " + itr.first)
                    );
        }
        a->second = itr.second.end();

    }
    step._currentStep->_date_itr = _date_vec.end();
    return step;
}


timestep& timeseries::iterator::dereference() const
{
    return *_currentStep;
}

bool timeseries::iterator::equal(iterator const& other) const
{
    bool isEqual = false;

    //different sizes? try to bail early
    if (_currentStep->_itrs.size() != other._currentStep->_itrs.size())
    {
        return false;
    }
    
    //no point checking headers as the order built is undefined
    //check each iterator
    for (auto& itr : _currentStep->_itrs)
    {
        for (auto& jtr : other._currentStep->_itrs)
        {
            if (itr.second == jtr.second)
                isEqual = true;
        }
    }

    if (isEqual && !(_currentStep->_date_itr == other._currentStep->_date_itr))
        isEqual = false; //negate if the date vectors don't match
    
    return isEqual;

}

void timeseries::iterator::increment()
{
    //walks the map locking each node so that the increment can happen
    //walk order is not guaranteed
    unsigned int size = _currentStep->_itrs.size();
    std::string *headers = new std::string[size];
    timestep::itr_map::accessor *accesors = new timestep::itr_map::accessor[size];
    int i = 0;

    for (auto& itr : _currentStep->_itrs)
    {
        _currentStep->_itrs.find(accesors[i], itr.first);
        (accesors[i]->second)++;     
        i++;
    }
    
    _currentStep->_date_itr++;

    delete[] headers;
    delete[] accesors;
}

void timeseries::iterator::decrement()
{
    //walks the map locking each node so that the increment can happen
    //walk order is not guaranteed
    unsigned int size = _currentStep->_itrs.size();
    std::string *headers = new std::string[size];
    timestep::itr_map::accessor *accesors = new timestep::itr_map::accessor[size];
    int i = 0;

    for (auto& itr : _currentStep->_itrs)
    {
        _currentStep->_itrs.find(accesors[i], itr.first);
        (accesors[i]->second)--;
        
        i++;
    }
    _currentStep->_date_itr--;
    delete[] headers;
    delete[] accesors;

}

timeseries::iterator::iterator()
{
    _currentStep = boost::make_shared<timestep>();
}

timeseries::iterator::iterator(const iterator& src)
{
    _currentStep = boost::make_shared<timestep>(src._currentStep);
}

timeseries::iterator::~iterator()
{
   // delete _currentStep;
}

timeseries::iterator& timeseries::iterator::operator=(const timeseries::iterator& rhs)
{
    if (this == &rhs)
        return (*this);
    _currentStep = boost::make_shared<timestep>(rhs._currentStep);
    return *this;
}

std::ptrdiff_t timeseries::iterator::distance_to(timeseries::iterator const& other) const
{
    return std::distance(this->_currentStep->_date_itr,other._currentStep->_date_itr);
}

void timeseries::iterator::advance(timeseries::iterator::difference_type N)
{
    for (int i = 0;i<N;i++)
        this->increment();
}