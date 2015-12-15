
#include "timestep.hpp"


timestep::timestep(const boost::shared_ptr<timestep> src)
{
   
    _itrs = itr_map(src->_itrs);
    _date_itr = date_variable::iterator(src->_date_itr);
}

timestep::timestep()
{
   
}

timestep::~timestep()
{

}

std::string timestep::to_string()
{
    std::stringstream s;
    s << get_posix() << "\t";
    for (timestep::itr_map::const_iterator itr = _itrs.begin(); itr != _itrs.end(); itr++)
    {
        s << boost::lexical_cast<std::string>(*(itr->second)) << std::string("\t");
    }

    return s.str();
}

int timestep::month()
{
    return _date_itr->date().month();
}


int timestep::day()
{
    return _date_itr->date().day();

}

int timestep::year()
{
    return _date_itr->date().year();
}

int timestep::hour()
{
    return boost::posix_time::to_tm(*_date_itr).tm_hour;
}
int timestep::min()
{
    return boost::posix_time::to_tm(*_date_itr).tm_min;
}
int timestep::sec()
{
    return boost::posix_time::to_tm(*_date_itr).tm_sec;
}


boost::gregorian::date timestep::get_gregorian()
{
    boost::gregorian::date date;

    date = boost::gregorian::from_string(boost::lexical_cast<std::string>(_date_itr->date()));
    return date;

}

boost::posix_time::ptime timestep::get_posix()
{
    return *_date_itr;
}

double timestep::get(std::string variable)
{
//    itr_map::const_accessor a;
//    if (!_itrs.find(a, variable))
//        BOOST_THROW_EXCEPTION( forcing_lookup_error() << errstr_info("Variable " + variable + " does not exist."));
    auto res = _itrs.find(variable);
    if(res == _itrs.end())
        BOOST_THROW_EXCEPTION( forcing_lookup_error() << errstr_info("Variable " + variable + " does not exist."));

    double out = res->second[0];

    return out;
}

void timestep::set(std::string variable, double value)
{
    auto res = _itrs.find(variable);
    if(res == _itrs.end())
        BOOST_THROW_EXCEPTION( forcing_lookup_error() << errstr_info("Variable " + variable + " does not exist."));

//    itr_map::accessor a;
//    if (!_itrs.find(a, variable))
//    {
//        BOOST_THROW_EXCEPTION( forcing_lookup_error() << errstr_info("Variable " + variable + " does not exist."));
//    }

    res->second[0] = value;

}

timestep::variable_vec::iterator timestep::get_itr(std::string varName)
{
    auto res = _itrs.find(varName);
    if(res == _itrs.end())
        BOOST_THROW_EXCEPTION( forcing_lookup_error() << errstr_info("Variable " + varName + " does not exist."));

//    itr_map::const_accessor a;
//    if (!_itrs.find(a, varName))
//        BOOST_THROW_EXCEPTION( forcing_lookup_error() << errstr_info("Variable " + varName + " does not exist."));


    return res->second;

}

