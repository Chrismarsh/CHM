#include "global.hpp"

global::global()
{
    first_time_step = true;
    _utc_offset = 0;
    _is_point_mode = false;
}

bool global::is_geographic()
{
    return _is_geographic;
}
int global::year()
{
    return _current_date.date().year();
}
int global::day()
{
    return _current_date.date().day();
}
int global::month()
{
    return _current_date.date().month();
}
int global::hour()
{
    return boost::posix_time::to_tm(_current_date).tm_hour;
}
int global::min()
{
    return boost::posix_time::to_tm(_current_date).tm_min;
}
int global::sec()
{
    return boost::posix_time::to_tm(_current_date).tm_sec;
}
boost::posix_time::ptime global::posix_time()
{
    return _current_date;
}

uint64_t global::posix_time_int()
{

    const boost::posix_time::ptime epoch = boost::posix_time::from_time_t(0);
    boost::posix_time::time_duration duration = _current_date - epoch;
    return duration.total_seconds();
}

std::string global::get_variable(std::string variable)
{
    return _variables(variable);
    
}

int global::dt()
{
    return _dt;
}

void global::insert_station(boost::shared_ptr<station> s)
{
    _stations.push_back(s);
    _dD_tree.insert( boost::make_tuple(Kernel::Point_2(s->x(),s->y()),s) );
//    _dD_NN_tree.insert( boost::make_tuple(Kernel::Point_2(s->x(),s->y()),s) );
}

void global::insert_stations(tbb::concurrent_vector< boost::shared_ptr<station> >& stations)
{
    _stations.resize( stations.size() );
#pragma omp for
    for(size_t i = 0; i< stations.size(); ++i)
    {
        boost::shared_ptr<station> s (stations.at(i));
        _stations.at(i) = s;
        _dD_tree.insert( boost::make_tuple(Kernel::Point_2(s->x(),s->y()),s) );
    }

//    _dD_NN_tree.insert( boost::make_tuple(Kernel::Point_2(s->x(),s->y()),s) );
}

std::vector< boost::shared_ptr<station> > global::get_stations_in_radius(double x, double y, double radius )
{
    // define exact circular range query  (fuzziness=0)
    Kernel::Point_2 center(x, y);
    Fuzzy_circle exact_range(center, radius);

    std::vector<boost::tuple<Kernel::Point_2, boost::shared_ptr<station> > > result;
    _dD_tree.search(std::back_inserter(result), exact_range);

    std::vector< boost::shared_ptr<station> > stations;

    for (auto& itr : result)
    {
        stations.push_back( boost::get<1>(itr));
    }
    return stations;

}

std::vector< boost::shared_ptr<station> > global::nearest_station(double x, double y,unsigned int N)
{
    Kernel::Point_2 query(x,y);
    Neighbor_search search(_dD_tree, query, N);

    std::vector< boost::shared_ptr<station> > stations;

    for (auto itr : search)
    {
        stations.push_back( boost::get<1>(itr.first));
    }
    return stations;

}

std::vector< boost::shared_ptr<station> > global::stations()
{
    return _stations;
}

bool  global::is_point_mode()
{
    return _is_point_mode;
}