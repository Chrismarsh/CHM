#include "timer.hpp"


timer::timer()
{
    _resolution = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::duration(1)).count();
}

void timer::tic()
{
    _timers.push(std::chrono::high_resolution_clock::now());
}


