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
    topbar  = nullptr;
    top     = nullptr;
    bottom  = nullptr;


    //layout:

    /*
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

};



ui::~ui()
{
    endwin();
    //gnome terminal isn't resetting properly if this is called from within a python script (e.g., xvalidation)
    std::system("stty sane;");

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
    mvwaddstr(topbar,0,(cols-7)/2,"CHM version 0.1a");

    int top_width = cols;
    int top_height = topbar_height+time_y+1; // make sure we can fit everything
    wresize(top,top_height,top_width);

    int bottom_width = cols;
    int bottom_height = lines-topbar_height-top_height;
    wresize(bottom,bottom_height,bottom_width);
    refresh();

    wrefresh(topbar);
    wrefresh(top);
    wrefresh(bottom);


}

void ui::write_timestep(std::string ts)
{
    if(top)
    {
        mvwprintw(top,ts_y, ts_x, "Current timestep is: %s",ts.c_str());
        wrefresh(top);
    }

}
void ui::write_progress(int prog)
{
    if(top)
    {
        mvwprintw(top,prog_y, prog_x, "Progress %d%",prog);
        wrefresh(top);
    }

}
void ui::write_time_estimate(std::string time)
{
    if(top)
    {
        mvwprintw(top,time_y, time_x, "Estimated completion %s",time.c_str());
        wrefresh(top);
    }
}

void ui::write_meantime(std::string s)
{
    if(top)
    {
        mvwprintw(top,meantime_y, meantime_x, "Mean timestep time %s", s.c_str());
        wrefresh(top);
    }

}
void ui::init()
{
    FILE *fd = fopen("/dev/tty", "r+");
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

    int top_width = COLS;
    int top_height = topbar_height+time_y+1; // make sure we can fit everything
    int top_startx = 0;
    int top_starty = topbar_height;
    top = create_window(top_height, top_width, top_startx, top_starty);

    int bottom_width = COLS;
    int bottom_height = LINES-topbar_height-top_height;
    int bottom_startx = 0;
    int bottom_starty = top_starty+top_height;
    bottom = create_window(bottom_height, bottom_width, bottom_startx, bottom_starty);

    if (!topbar || !top || !bottom )
    {
        BOOST_THROW_EXCEPTION(model_init_error() << errstr_info("Unable to initialize ncurses!"));
    }

    buf = new nc_window_streambuf( bottom, std::cout );

    //color topbar
    init_pair(1,COLOR_BLACK, COLOR_WHITE);
    wbkgd(topbar, COLOR_PAIR(1));

    //allow bottom to scroll
    scrollok(bottom,TRUE);

    mvwaddstr(topbar,0,(COLS-7)/2,"CHM version 0.1a");

    box(top, 0, 'bs');
    //(top, 0, 'bs');
   // wborder(bottom,0,0,' ',' ',0,0,0,0);

    sig_callback = std::bind(&ui::handle_sig_winch,this,std::placeholders::_1);
    memset(&sa, 0, sizeof(struct sigaction));
    sa.sa_handler = sig_winch_wrapper;
   // sigaction(SIGWINCH, &sa, NULL);

    wrefresh(topbar);
    wrefresh(top);
    wrefresh(bottom);
    refresh();

}

WINDOW* ui::create_window(int height, int width, int startx, int starty)
{
    WINDOW* local_win =  newwin(height, width, starty, startx);

    wrefresh(local_win);

    return local_win;
}