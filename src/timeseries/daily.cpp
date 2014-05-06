#include "daily.hpp"

namespace daily
{
    double mean(time_series::const_iterator now, std::string variable )
    {
        bool found_start = false;
        int day = now->day();//our 'current' step
        while (!found_start)
        {
            --now;
            if (now->day() != day) //found the end of the previous day)
            {
                now++; //get back to the start of our day
                found_start = true;
            }
        }
        
        //first one
        double mean = now->get(variable);
        double counter = 1;
        
        ++now;
        while(now->day() != day)
        {
            mean += now->get(variable);
            counter++;
            ++now;
        }
        
        return mean/counter;
    }
    
    
}
