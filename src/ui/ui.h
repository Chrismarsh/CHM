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

    //write the current timestep
    void write_timestep(std::string ts);

    //current model progress percentage
    void write_progress(int prog);

    //how much longer the model will take
    void write_time_estimate(std::string time);

    //the average time per timestep
    void write_meantime(std::string s);

    //supplemental information
    void write_model_name(std::string s);

    //current working directory
    void write_cwd(std::string s);

    //writes the mesh details
    void write_mesh_details(int tri);

    //writes the mesh details
    void write_modules(std::string s);

    //handle window resize. currently broken :(
    // TODO: fix sig_winch handler
    void handle_sig_winch(int sig);

    void end();
private:
    WINDOW* topbar;  // top bar
    WINDOW* bottombar;  // top bar
    WINDOW* top_left;     //runtime information
    WINDOW* top_right;     //runtime information
    WINDOW* bottom;   //logging information

    WINDOW* create_window(int height, int width, int startx, int starty);

    bool is_init;
    /*
     * TOP LEFT PANE
     */
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

    //mean time location
    int topbar_x;
    int topbar_y;

    /*
     * TOP RIGHT PANE
     */
    //mesh information
    int mesh_name_x;
    int mesh_name_y;

    int module_x;
    int module_y;

    std::string topbar_title;

    struct sigaction sa;
    SCREEN *scr;
    //swap out cout to bottom window
    nc_window_streambuf* buf;
};


