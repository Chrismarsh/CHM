#pragma once

#include <chrono>

#include <stack>

//Title: clock
    typedef std::chrono::milliseconds ms;
    typedef std::chrono::nanoseconds ns;
    typedef std::chrono::seconds s;
class timer
{
public:


    /*
    Function: clock
            Creates a clock.

    Parameters: 
    None

    Throws:
    Never

    Returns:   
		
     */
    timer();

    /*
    Function: Tic
            Begins a timer

    Parameters: 
            None

    Throws:
            Never

    Returns:   
				
     */
    void tic();

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
template <typename T>
    double toc();

private:
    std::stack<std::chrono::high_resolution_clock::time_point> _timers;
    double _resolution; //in ns
};


template <typename T>
double timer::toc()
{
    auto end = std::chrono::high_resolution_clock::now();
    auto start = _timers.top();
    _timers.pop();
    auto dur = std::chrono::duration_cast< T>(end - start);
    return dur.count();
}