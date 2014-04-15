#pragma once

#ifdef _WIN32   
    #include <windows.h>
#elif __MACH__
    #include <sys/time.h>
    #include <mach/mach.h>
#else          
    #include <time.h>
#endif


#include <stack>

//Title: clock

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
    double toc();

private:


#ifdef _WIN32
    LARGE_INTEGER m_frequency;
    std::stack<LARGE_INTEGER> m_timers;
#else
    std::stack<timespec> m_timers;
#endif
};
