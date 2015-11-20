#pragma once

#include <curses.h>
#include <signal.h>
#include <functional>


#include <sys/ioctl.h>
#include "exception.hpp"
#include "logger.hpp"
#include "ncstream.h"
//ncurses UI
class ui
{
public:
    ui();
    ~ui();
   void init();

    void write_timestep(std::string ts);
    void write_progress(int prog);
    void write_time_estimate(std::string time);
    void write_meantime(std::string s);
    void handle_sig_winch(int sig);
private:
    WINDOW* topbar;  // top bar
    WINDOW* top;     //runtime information
    WINDOW* bottom;   //logging information

    WINDOW* create_window(int height, int width, int startx, int starty);

    //Timestep location
    int ts_x;
    int ts_y;

    //Progress location
    int prog_x;
    int prog_y;

    //time estimate location
    int time_x;
    int time_y;

    //mean time location
    int meantime_x;
    int meantime_y;

    struct sigaction sa;
    SCREEN *scr;
    //swap out cout to bottom window
    nc_window_streambuf* buf;
};


