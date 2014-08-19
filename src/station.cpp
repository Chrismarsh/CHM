#include "station.hpp"

station::~station()
{
    delete _obs;
}

station::station()
{
        _x = 0;
        _y = 0;
        _z = 0.0;
        _obs = NULL;

}

station::station( std::string ID, std::string file, unsigned int x, unsigned int y, float elevation )
{
        _ID = ID;
        _x = x;
        _y = y;
        _z = elevation;
        _obs = NULL; //initialized in openfile

        open(file);
}

void station::open( std::string file )
{
    try
    {
        _obs = new time_series();
        _obs->open(file);

        _itr = _obs->begin();
    }
    catch(exception_base& e)
    {
        //e << errstr_info( std::string("Station:") + _ID ); //hope like hell we've got the ID at this point
        throw;
    }
}

time_series::date_vec station::get_date_timeseries()
{
    return _obs->get_date_timeseries();
}

int station::get_timeseries_length()
{
    return _obs->get_timeseries_length();
}
 
std::vector<std::string> station::list_variables()
{
    return _obs->list_variables();
}


 timestep station::now() 
{
    return *_itr;
}
double station::get(std::string variable)
{
    return _itr->get(variable);
}
void station::reset_itrs()
{
    _itr = _obs->begin();
}

bool station::next()
{
        ++_itr;
        if(_itr == _obs->end())
                return false;
        else
                return true;
}

unsigned int station::get_x()
{
        return _x;
}

unsigned int station::get_y()
{
        return _y;
}

void station::set_x( unsigned int x )
{
        _x = x;
}

void station::set_y( unsigned int y )
{
        _y = y;
}

float station::get_z()
{
        return _z;
}

void station::set_z( float elevation )
{
        _z = elevation;
}

void station::set_ID(std::string ID)
{
    _ID = ID;
}

std::string station::get_ID()
{
    return _ID;
}

std::ostream& operator<<(std::ostream &strm, const station &s)
{
    if(s._obs)
        return strm << "ID=" << s._ID << " (x,y,z)=("<<s._x<<","<<s._y<<","<<s._z<<") ,"<<"forcing="<<s._obs->get_opened_file();
    else
        return strm << "ID=" << s._ID << "; (x,y,z)=("<<s._x<<","<<s._y<<","<<s._z<<"); "<<"forcing=Not opened";
}