#include "daily.hpp"

namespace daily
{
    double mean(time_series::iterator& start, time_series::iterator& end, std::string variable )
    {
        double m = 0.0;
        double count = 0.0;
        for(auto& itr = start; itr < end; itr++)
        {
            m += itr->get(variable);
            count++;
        }
        return m/count;
    }
    
    double max(time_series::iterator& start, time_series::iterator& end, std::string variable )
    {
        auto m = std::max_element(start->get_itr(variable),end->get_itr(variable));
        return *m;
    }
    
}
