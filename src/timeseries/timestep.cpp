
#include "timestep.hpp"


timestep::timestep(const timestep& src)
{
   
    _itrs = const_itr_map(src._itrs);
    _date_itr = date_variable::const_iterator(src._date_itr);
}

timestep::timestep()
{
   
}

timestep::~timestep()
{

}

std::string timestep::to_string()
{
    std::string s;

    for (timestep::const_itr_map::const_iterator itr = _itrs.begin(); itr != _itrs.end(); itr++)
    {
        s += boost::lexical_cast<std::string>(*(itr->second)) + std::string("\t");
    }

    return s;
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