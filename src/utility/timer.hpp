#pragma once

#include <chrono>

#include <stack>

//Title: clock

class timer
{
public:
    typedef std::chrono::milliseconds ms;
    typedef std::chrono::nanoseconds ns;
    typedef std::chrono::seconds s;

    /*
    Function: clock
            Creates a clock.

    Parameters: 
    None

    Throws:
    Never

    Returns:   
		
     */
    timer()
    {
        _resolution = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::duration(1)).count();
    }

    /*
    Function: Tic
            Begins a timer

    Parameters: 
            None

    Throws:
            Never

    Returns:   
				
     */
    void tic()
    {
        _timers.push(std::chrono::high_resolution_clock::now());
    }

    /*
    Function: Toc
            Returns the time elapsed in micro-seconds since the last Tic() call

    Parameters: 
            None

    Throws:
            Never

    Returns:   
            double - Time elapsed, micro-seconds. -1.0 if Tic() hasn't been called yet.
     */
    template<class K>
    double toc()
    {
        auto end = std::chrono::high_resolution_clock::now();
        auto start = _timers.top();_timers.pop();
        auto dur = std::chrono::duration_cast<K>(end - start);
        return dur.count();
    }

private:
    std::stack<std::chrono::high_resolution_clock::time_point> _timers;
    double _resolution; //in ns
};

