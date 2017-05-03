//
// Created by chris on 19/11/15.
//

#include "ui.h"
std::function<void(int)> sig_callback;
void sig_winch_wrapper(int sig)
{
    sig_callback(sig);
}

ui::ui()
{
    buf = nullptr;
    scr = nullptr;
    topbar  = nullptr;
    top_left = nullptr;
    top_right     = nullptr;
    bottom  = nullptr;
    bottombar  = nullptr;

    //layout:
//    [ topbar ]
//    [  top | right  ]
//    [     bottom   ]
    /* [ CHM 0.1a  -- Model Name     ]
     * Current time step is: 2012-May-08 00:00:00
     * Progress 70%
     * Mean time per step 50ms
     * Estimated completion is: xxxxxx
     */

    //Timestep location
    ts_x = 2;
    ts_y = 1;

    //Progress location
    prog_x = 2;
    prog_y = ts_y+2;

    //time estimate location
    meantime_x = 2;
    meantime_y = prog_y+2;

    //time estimate location
    time_x = 2;
    time_y = meantime_y+2;

    topbar_x = 3;
    topbar_y = 0;

    mesh_name_x = 1;
    mesh_name_y = 1;

    module_y = 3;
    module_x = 1;

    topbar_title = " CHM v0.1a ";

    is_init = false;

};



ui::~ui()
{


}


void ui::end()
{
    if(is_init)
    {
        LOG_DEBUG << "Cleaning up ncurses";
        endwin();
        //gnome terminal isn't resetting properly if this is called from within a python script (e.g., xvalidation)
//        std::system("stty sane");
    }

}

void ui::write_modules(std::string s)
{
    if(is_init)
    {
        mvwprintw(top_right,module_y, module_x, "Modules: %s",s.c_str());
        wrefresh(top_right);
    }
}
void ui::write_mesh_details(int tri)
{
    if(is_init)
    {
        mvwprintw(top_right,mesh_name_y, mesh_name_x, "# Triangles: %d",tri);
        wrefresh(top_right);
    }

}

void ui::write_model_name(std::string s)
{
    if(is_init)
    {
        mvwprintw(topbar,topbar_y, topbar_x, "%s -- %s",topbar_title.c_str(),s.c_str());
        wrefresh(topbar);
    }

}

void ui::write_timestep(std::string ts)
{
    if(is_init)
    {
        mvwprintw(top_left, ts_y, ts_x, "Current timestep is: %s", ts.c_str());
        wrefresh(top_left);
    }

}
void ui::write_progress(int prog)
{
    if(is_init)
    {
        mvwprintw(top_left, prog_y, prog_x, "Progress %d%", prog);
        wrefresh(top_left);
    }

}
void ui::write_cwd(std::string s)
{
    if(is_init)
    {
        mvwprintw(bottombar, 0, 0, "%s", s.c_str());
        wrefresh(bottombar);
    }

}
void ui::write_time_estimate(std::string time)
{
    if(is_init)
    {
        mvwprintw(top_left, time_y, time_x, "Estimated completion %s", time.c_str());
        wrefresh(top_left);
    }
}

void ui::write_meantime(std::string s)
{
    if(is_init)
    {
        mvwprintw(top_left, meantime_y, meantime_x, "Mean timestep time %s", s.c_str());
        wrefresh(top_left);
    }

}
void ui::init()
{
    FILE *fd = fopen("/dev/tty", "r+");
    if(!fd)
        BOOST_THROW_EXCEPTION(model_init_error() << errstr_info("Unable to initialize ncurses!"));

    scr = newterm(NULL, fd, fd);
    noecho(); //don't echo input
    curs_set(0); //no cursor
    start_color(); //colors!!!!!
    refresh();

    int topbar_width = COLS;
    int topbar_height = 1;
    int topbar_startx = 0;
    int topbar_starty = 0;
    topbar = create_window(topbar_height, topbar_width, topbar_startx, topbar_starty);

    int bottombar_width = COLS;
    int bottombar_height = 1;
    int bottombar_startx = 0;
    int bottombar_starty = LINES-bottombar_height;
    bottombar = create_window(bottombar_height, bottombar_width, bottombar_startx, bottombar_starty);

    int top_width = COLS/2;
    int top_height = topbar_height+time_y+1+bottombar_height; // make sure we can fit everything
    int top_startx = 0;
    int top_starty = topbar_height;
    top_left = create_window(top_height, top_width, top_startx, top_starty);

    int top_rwidth = COLS/2;
    int top_rheight = topbar_height+time_y+1+bottombar_height; // make sure we can fit everything
    int top_rstartx = top_width;
    int top_rstarty = topbar_height;
    top_right = create_window(top_rheight, top_rwidth, top_rstartx, top_rstarty);


    int bottom_width = COLS;
    int bottom_height = LINES-topbar_height-top_height-bottombar_height;
    int bottom_startx = 0;
    int bottom_starty = top_starty+top_height;
    bottom = create_window(bottom_height, bottom_width, bottom_startx, bottom_starty);

    if (!topbar || !top_left || !bottom )
    {
        BOOST_THROW_EXCEPTION(model_init_error() << errstr_info("Unable to initialize ncurses!"));
    }

    is_init = true;

    buf = new nc_window_streambuf( bottom, std::cout );

    //color topbar
    init_pair(1,COLOR_BLACK, COLOR_WHITE);
    wbkgd(topbar, COLOR_PAIR(1));

    //color bottombar
    init_pair(1,COLOR_BLACK, COLOR_WHITE);
    wbkgd(bottombar, COLOR_PAIR(1));

    //allow bottom to scroll
    scrollok(bottom,TRUE);

    mvwaddstr(topbar,topbar_y, topbar_x, topbar_title.c_str());

//    box(top_left, 0, 'bs');
//    box(top_right, 0, 'bs');
    //(top, 0, 'bs');
   // wborder(bottom,0,0,' ',' ',0,0,0,0);

    sig_callback = std::bind(&ui::handle_sig_winch,this,std::placeholders::_1);
    memset(&sa, 0, sizeof(struct sigaction));
    sa.sa_handler = sig_winch_wrapper;
   // sigaction(SIGWINCH, &sa, NULL);

    wrefresh(topbar);
    wrefresh(bottombar);
    wrefresh(top_left);
    wrefresh(top_right);
    wrefresh(bottom);


}

WINDOW* ui::create_window(int height, int width, int startx, int starty)
{
    WINDOW* local_win =  newwin(height, width, starty, startx);

    wrefresh(local_win);

    return local_win;
}

void ui::handle_sig_winch(int sig)
{
//    endwin();
//    init();
//    write_status("resized");
//    std::ifstream myfile ("CHM.log");
//    std::string line;
//    if (myfile.is_open())
//    {
//        while ( getline (myfile,line) )
//        {
//            std::cout << line << '\n';
//        }
//        myfile.close();
//    }
//    return;
    std::system("resize &> /dev/null");

    struct winsize ws;
    ioctl(1, TIOCGWINSZ, &ws);
    int lines,cols;
    cols = COLS;//ws.ws_col;
    lines = LINES;// ws.ws_row;

    wresize(stdscr,lines,cols);

    int topbar_width = cols;
    int topbar_height = 1;
    wresize(topbar,topbar_height,topbar_width);

    mvwaddstr(topbar,topbar_y, topbar_x, " CHM v0.1a");


    int top_width = cols;
    int top_height = topbar_height+time_y+1; // make sure we can fit everything
    wresize(top_left, top_height, top_width);

    int bottom_width = cols;
    int bottom_height = lines-topbar_height-top_height;
    wresize(bottom,bottom_height,bottom_width);
    refresh();

    wrefresh(topbar);
    wrefresh(top_left);
    wrefresh(bottom);


}
