
#include "timer.hpp"

void timer::tic()
{

#ifdef _WIN32
    LARGE_INTEGER startCount;
    startCount.QuadPart = 0;
    QueryPerformanceCounter(&startCount);
    m_timers.push(startCount);
#elif __MACH__
  timespec startCount;
  struct timeval now;
  int rv = gettimeofday(&now, NULL);
 
  startCount.tv_sec = now.tv_sec;
  
  m_timers.push(startCount);
  
#else
    timespec startCount;
    clock_gettime(CLOCK_MONOTONIC, &startCount); //CLOCK_PROCESS_CPUTIME_ID
    m_timers.push(startCount);
#endif
}

timer::timer()
{

#ifdef _WIN32
    QueryPerformanceFrequency(&m_frequency);
#endif

}

double timer::toc()
{



#ifdef _WIN32
    LARGE_INTEGER endCount;
    endCount.QuadPart = 0;

    QueryPerformanceCounter(&endCount);

    if (m_timers.empty())
        return -1.0;

    m_start = (m_timers.top()).QuadPart * (1000000.0 / m_frequency.QuadPart);
    m_timers.pop();
    m_end = endCount.QuadPart * (1000000.0 / m_frequency.QuadPart);
    return m_end - m_start;
 #elif __MACH__


    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    
    timespec end;
    end.tv_sec = now.tv_sec;
//    end.tv_nsec = now.tv_nsec;
    
    timespec start = m_timers.top();
    m_timers.pop();

    timespec temp;
    temp.tv_sec = end.tv_sec - start.tv_sec;
//    temp.tv_nsec = end.tv_nsec - start.tv_nsec;

    return temp.tv_sec;
#else

    timespec end;
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    
    timespec temp;
    timespec start = m_timers.top();
    m_timers.pop();
    
    //following:
    //http://www.guyrutenberg.com/2007/09/22/profiling-code-using-clock_gettime/
    if ((end.tv_nsec - start.tv_nsec) < 0)
    {
        temp.tv_sec = end.tv_sec - start.tv_sec - 1;
        temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
    } else
    {
        temp.tv_sec = end.tv_sec - start.tv_sec;
        temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }

    return temp.tv_sec;
#endif

    
}
