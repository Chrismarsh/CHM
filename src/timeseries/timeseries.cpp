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



#include "timeseries.hpp"

void timeseries::push_back(double data, std::string variable)
{
    _variables[variable].push_back(data);

}

void timeseries::init_new_variable(std::string variable)
{
    size_t size = _date_vec.size();
    if (size == 0)
    {
        CHM_THROW_EXCEPTION(forcing_error, "Adding variable to uninitialized timeseries");
    }

    _variables[variable].assign(size,-9999.0);
}

void timeseries::init(std::set<std::string> variables, boost::posix_time::ptime start_time, boost::posix_time::ptime end_time, boost::posix_time::time_duration dt)
{
    size_t size = 0;
    boost::posix_time::ptime ts = start_time;

    //figure out how many timesteps we need
    // < end_time as we will do +1 to equal end_time
    while(ts < end_time)
    {
        ts =  start_time + dt*size;
        ++size;
    }

    // if we have exactly 1 timestep, special case this
    if(size == 0 && start_time == end_time)
        size = 1;

    _timeseries_length = size;

    for (auto& v: variables)
    {
        _variables[v].assign(_timeseries_length,-9999.0);
    }


    _date_vec.resize(_timeseries_length);
    for(size_t i=0; i <_timeseries_length;i++)
    {
        _date_vec[i] = start_time + dt*i;
    }

}
void timeseries::init(std::set<std::string> variables, date_vec datetime)
{
    size_t size = datetime.size();

   for (auto& v: variables)
   {
       _variables[v].assign(size,-9999.0);
   }

   //setup date vector
   _date_vec = datetime;
}

 timeseries::date_vec timeseries::get_date_timeseries()
 {
     return _date_vec;
 }
std::set<std::string> timeseries::list_variables()
{
    std::set<std::string> vars;
    for(auto& itr : _variables)
    {
        vars.insert(itr.first);
    }
    
    return vars;
}
double& timeseries::at(std::string variable, size_t idx)
{
    auto res = _variables.find(variable);
    if(res == _variables.end())
    {
        CHM_THROW_EXCEPTION(forcing_error, "Unable to find " + variable);
    }
    return const_cast<double&>(res->second.at(idx));
}

timeseries::variable_vec timeseries::get_time_series(std::string variable)
{
    auto res = _variables.find(variable);
    if(res == _variables.end())
    {
        CHM_THROW_EXCEPTION(forcing_error, "Unable to find " + variable);
    }   
    return res->second;
}

void timeseries::subset(boost::posix_time::ptime start,boost::posix_time::ptime end)
{
    //look for our requested timestep
    auto itrstart = std::find(_date_vec.begin(),_date_vec.end(),start);
    if ( itrstart == _date_vec.end())
    {
        CHM_THROW_EXCEPTION(forcing_timestep_notfound, "Start timestep not found");
    }

    //Find the first one
    //get offset from iterator
    auto dist_start = std::distance(_date_vec.begin(), itrstart);
    auto itrend = std::find(_date_vec.begin()+dist_start,_date_vec.end(),end);

    if(itrend == _date_vec.end())
    {
        SPDLOG_WARN("Requested end date is past last date. Setting date end = time series end.");
        itrend = std::next(_date_vec.begin(),  _date_vec.size() - 1); //skip to last item
    }
    else{
        itrend++;//need to include the last item we asked for, so step once more.
    }

    auto dist_end = std::distance(_date_vec.begin(), itrend);

    //iterate over the map of vectors and build a list of all the variable names
    //unknown order
    for (auto& itr : _variables)
    {
       auto start_itr = itr.second.begin() + dist_start;
       auto end_itr = itr.second.begin() + dist_end;

        std::vector<double> temp(start_itr,end_itr);
        //insert the iterator
        itr.second = temp;
    }

    auto start_itr =_date_vec.begin() + dist_start;
    auto end_itr = _date_vec.begin() + dist_end;
    date_vec temp(start_itr,end_itr);
    _date_vec = temp;

}
boost::tuple<timeseries::iterator, timeseries::iterator> timeseries::range(boost::posix_time::ptime start_time,boost::posix_time::ptime end_time)
{
    //look for our requested timestep
    auto itr_find = std::find(_date_vec.begin(),_date_vec.end(),start_time);
    if ( itr_find == _date_vec.end())
    {
        CHM_THROW_EXCEPTION(forcing_timestep_notfound, "Timestep not found");
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
//        timestep::itr_map::accessor a;
        //create the keyname for this variable and store the iterator

//        auto res = start_step._currentStep->_itrs.insert(itr.first);
//        if (!start_step._currentStep->_itrs.insert(a, itr.first))
//        {
//            BOOST_THROW_EXCEPTION(forcing_error()
//                    << errstr_info("Failed to insert " + itr.first)
//                    );
//        }
//
        start_step._currentStep->_itrs[itr.first]= itr.second.begin()+dist_start;
        //insert the iterator
//        res->second =
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

//        //create the keyname for this variable and store the iterator
//        if (!end_step._currentStep->_itrs.insert(itr.first))
//        {
//            BOOST_THROW_EXCEPTION(forcing_error()
//                    << errstr_info("Failed to insert " + itr.first)
//                    );
//        }
//
        //insert the iterator
        end_step._currentStep->_itrs[itr.first] = itr.second.begin()+dist_end;
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
        CHM_THROW_EXCEPTION(forcing_timestep_notfound, "Timestep not found");
    }
    
    //get offset from iterator
    int dist = std::distance(_date_vec.begin(), itr);
    
    iterator step;

    //iterate over the map of vectors and build a list of all the variable names
    //unknown order
    for (auto& itr : _variables)
    {
//        //create the keyname for this variable and store the iterator
//        if (!step._currentStep->_itrs.insert(itr.first))
//        {
//            BOOST_THROW_EXCEPTION(forcing_insertion_error()
//                    << errstr_info("Failed to insert " + itr.first)
//                    );
//        }
        
        //insert the iterator
        step._currentStep->_itrs[itr.first] = itr.second.begin()+dist;
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
    {
//        std::stringstream ss;
//        ss << "" << boost::errinfo_errno(errno) << boost::errinfo_file_name(path);
        CHM_THROW_EXCEPTION(file_read_error, boost::to_string(boost::errinfo_errno(errno)));
    };

    SPDLOG_DEBUG("Parsing file {}", path);
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

//    for (std::vector<std::string>::const_iterator itr = header.begin();
//            itr != header.end();
//            itr++)
//    {
//        if (!_variables.insert(*itr))
//        {
//            BOOST_THROW_EXCEPTION(forcing_insertion_error()
//                    << errstr_info(std::string("Failed to insert ") + *itr)
//                    << boost::errinfo_file_name(path)
//                    );
//
//        }
//    }



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
//                    auto res = _variables.find(*headerItr);
//                    if (res == _variables.end())
//                        BOOST_THROW_EXCEPTION(forcing_lookup_error()
//                            << errstr_info(std::string("Failed to find ") + *headerItr)
//                            << boost::errinfo_file_name(path)
//                            );
                    try
                    {
                        //LOG_VERBOSE << "Found " << *headerItr << ": " << doubles[0];
                        _variables[*headerItr].push_back(boost::lexical_cast<double>(doubles[0]));
                    } catch (...)
                    {
                        CHM_THROW_EXCEPTION(forcing_badcast, "Failed to cast " + doubles[0] + " to a double. " + path);
                    }
                } else if ((dates = dateTime.tokenize<std::string>(*itr)).size() == 1)
                {

                    //LOG_VERBOSE << "Found " << *headerItr << ": " << dates[0];
                    _date_vec.push_back(boost::posix_time::from_iso_string(dates[0])); //from_iso_string
                    
                    //now we know where the date colum is, we remove it from the hashmap if we haven't already
                    _variables.erase(*headerItr);
                  
                }
                else
                {
                    //something has gone horribly wrong
                    CHM_THROW_EXCEPTION(forcing_no_regexmatch,"Unable to match any regex for " + *itr + ". Line: " + std::to_string(lines) + path);

                }

                //next header
                headerItr++;
                cols_so_far++;

            }
            if (cols_so_far != _cols)
            {
                CHM_THROW_EXCEPTION(forcing_badcast,"Expected " + std::to_string(_cols) + " columns on line " + std::to_string(_rows) + path);
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
    SPDLOG_DEBUG("Read in {} variables", _variables.size());
    std::string* headerItems = new std::string[_variables.size()];

    int i = 0;
    //build a list of all the headers
    //unknown order
    for (ts_hashmap::iterator itr = _variables.begin(); itr != _variables.end(); itr++)
    {
        //LOG_VERBOSE << itr->first;
        headerItems[i++] = itr->first;
    }

    //get and save each accessor
    size_t d_length = _date_vec.size();


    for (unsigned int l = 0; l < _variables.size(); l++)
    {
        //compare all columns to date length
        auto res = _variables.find( headerItems[l]);
        if (res == _variables.end())
        {
            CHM_THROW_EXCEPTION(forcing_lookup_error, std::string("Failed to find ") + headerItems[l] + path);
        }

        //check all cols are the same size as the first col
        if (d_length != res->second.size())
        {
            SPDLOG_ERROR("Col {} is a different size. Expected size={}", headerItems[l], boost::lexical_cast<std::string>(d_length));
            CHM_THROW_EXCEPTION(forcing_lookup_error, "Col " + headerItems[l] + " is a different size. Expected size="+boost::lexical_cast<std::string>(d_length)+path);

        }
        
    }

    delete[] headerItems;

    //we can only check date-time consistency if we have more than 1 datetime
    if (_date_vec.size() > 1)
    {
        auto dt = _date_vec.at(1) - _date_vec.at(0);

        for (size_t i = 1; i < _date_vec.size(); ++i)
        {
            //using our calculated timestep, check what we think out timestep should be
            auto pred_ts = _date_vec.at(i - 1) + dt;
            auto &actual_ts = _date_vec.at(i);
            if (pred_ts != actual_ts)
            {
                //streams will pretty-print the boost time nicely
                std::stringstream expected_ts;
                expected_ts << pred_ts;
                std::stringstream act_ts;
                act_ts << actual_ts;

                CHM_THROW_EXCEPTION(forcing_lookup_error,"On line " + std::to_string(i + 1) +
                                                             " the timestep is inconsistent with dt. Expected "
                                                             + expected_ts.str() + " got " + act_ts.str() + path);
            }
        }
    }
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
#ifdef USE_SPARSEHASH
    _variables.set_empty_key("");
#endif
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
    {
        CHM_THROW_EXCEPTION(file_read_error, boost::to_string(boost::errinfo_errno(errno)) + file);
    }

    
    std::string* headerItems = new std::string[_variables.size()];

    //build a list of all the headers
    //unknown order
    int i = 0;
    out << "datetime";
    variable_vec::const_iterator *tItr = new variable_vec::const_iterator[_variables.size()];
    for (ts_hashmap::iterator itr = _variables.begin(); itr != _variables.end(); itr++)
    {
        headerItems[i] = itr->first;
        out << "," << itr->first;

        //save vector iterators
        tItr[i] = itr->second.begin();
        _rows = itr->second.size();
        i++;
    }
    out << std::endl;

    
    for (size_t k = 0; k < _rows; k++)
    {
        out << boost::posix_time::to_iso_string(_date_vec.at(k));
        for (size_t j = 0; j < _variables.size(); j++)
        {
            out << "," << *(tItr[j]);
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
//        auto res = step._currentStep->_itrs.insert(itr.first);
        //create the keyname for this variable and store the iterator
//        if (res == )
//        {
//            BOOST_THROW_EXCEPTION(forcing_insertion_error()
//                    << errstr_info("Failed to insert " + itr.first)
//                    );
//        }
//
        //insert the iterator
        step._currentStep->_itrs[itr.first] = itr.second.begin();
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
//        auto res = step._currentStep->_itrs.insert(itr.first);
//        if (!)
//        {
//            BOOST_THROW_EXCEPTION(forcing_insertion_error()
//                    << errstr_info("Failed to insert " + itr.first)
//                    );
//        }
        step._currentStep->_itrs[itr.first] = itr.second.end();

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
//    unsigned int size = _currentStep->_itrs.size();
//    std::string *headers = new std::string[size];
//    timestep::itr_map::accessor *accesors = new timestep::itr_map::accessor[size];
    int i = 0;

    for (auto& itr : _currentStep->_itrs)
    {
        itr.second++;
//        _currentStep->_itrs.find(accesors[i], itr.first);
//        (accesors[i]->second)++;
//        i++;
    }
    
    _currentStep->_date_itr++;
//
//    delete[] headers;
//    delete[] accesors;
}

void timeseries::iterator::decrement()
{
    //walks the map locking each node so that the increment can happen
    //walk order is not guaranteed
//    unsigned int size = _currentStep->_itrs.size();
//    std::string *headers = new std::string[size];
//    timestep::itr_map::accessor *accesors = new timestep::itr_map::accessor[size];
//    int i = 0;

    for (auto& itr : _currentStep->_itrs)
    {
//        _currentStep->_itrs.find(accesors[i], itr.first);
//        (accesors[i]->second)--;
        itr.second --;
//        i++;
    }
    _currentStep->_date_itr--;
//    delete[] headers;
//    delete[] accesors;

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