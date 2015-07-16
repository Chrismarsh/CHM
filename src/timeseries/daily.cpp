
#include "daily.hpp"

namespace daily
{

    double min(timeseries& ts,timeseries::iterator& now, std::string variable)
    {
        auto start = start_of_day(ts,now);
        auto end = end_of_day(ts, now);

        return ts.range_min(start, end, variable);
    }

    double max(timeseries& ts,timeseries::iterator& now, std::string variable)
    {
        auto start = start_of_day(ts,now);
        auto end = end_of_day(ts, now);

        return ts.range_max(start, end, variable);
    }

    timeseries::iterator start_of_day(timeseries& ts, timeseries::iterator& now)
    {

        //today but midnight
        boost::posix_time::ptime time ( boost::gregorian::date(now->year(),now->month(),now->day()),
                                        boost::posix_time::time_duration(0,0,0));

        return ts.find(time);

    }

    timeseries::iterator end_of_day(timeseries& ts, timeseries::iterator& now)
    {

        //today but almost midnight
        boost::posix_time::ptime time ( boost::gregorian::date(now->year(),now->month(),now->day()),
                boost::posix_time::time_duration(23,0,0));

        timeseries::iterator itr = ts.find(time);

        int day = now->day();
        while(itr != ts.end() && itr->day() != day ) //stop at the last element in this day
            ++itr;

        return itr;

    }
}
