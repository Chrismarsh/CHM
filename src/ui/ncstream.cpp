//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

//
// Created by chris on 19/11/15.
//

#include "ncstream.h"


nc_window_streambuf::nc_window_streambuf( WINDOW * p,
                                          unsigned long curses_attr )
        : m_pnl(p), m_flags(curses_attr), m_os(0),m_old(0)
{
    // Tell parent class that we want to call overflow() for each
    // input char:
    this->setp( 0, 0 );
    this->setg( 0, 0, 0 );
    scrollok(p,true);
    //mvwinch( p, p->_maxy, 0 );
    mvwinch( p, 0, 0 );
}

nc_window_streambuf::nc_window_streambuf( WINDOW * p,
                                          std::ostream & os,
                                          unsigned long curses_attr )
        : m_pnl(p), m_flags(curses_attr), m_os(&os),m_old(os.rdbuf())
{
    this->setp( 0, 0 );
    this->setg( 0, 0, 0 );
    os.rdbuf( this );
    scrollok(p,true);
    //mvwinch(p, p->_maxy, 0 );
    mvwinch( p, 0, 0 );
}

void nc_window_streambuf::copy( const nc_window_streambuf & rhs )
{
    if( this != &rhs )
    {
        this->m_pnl = rhs.m_pnl;
        this->m_flags = rhs.m_flags;
        this->m_os = rhs.m_os;
        this->m_old = rhs.m_old;
    }
}

nc_window_streambuf::nc_window_streambuf( const nc_window_streambuf & rhs )
{
    this->copy(rhs);
}

nc_window_streambuf & nc_window_streambuf::operator=( const nc_window_streambuf & rhs )
{
    this->copy(rhs);
    return *this;
}


nc_window_streambuf::~nc_window_streambuf()
{
    if( this->m_os )
    {
        this->m_os->rdbuf( this->m_old );
    }
}

int nc_window_streambuf::overflow( int c )
{
    int ret = c;
    if( c != EOF )
    {
        if( this->m_flags )
        {
            wattron(this->m_pnl, this->m_flags );
            if( ERR == waddch(this->m_pnl, (chtype)c ) ) ret = EOF;
            wattroff( this->m_pnl, this->m_flags );
        }
        else if( ERR == waddch( this->m_pnl, (chtype)c ) ) ret = EOF;
    }
    if( (EOF==c) || std::isspace(c) )
    {
        if( EOF == this->sync() ) ret = EOF;
    }
    return ret;
}

int nc_window_streambuf::sync()
{
    if( stdscr && this->m_pnl )
    {
        //return (ERR == wrefresh( this->m_pnl )) ? EOF : 0;
        // ^^^ doing wrefresh here can hose panels :(
//        if( 0 != panel_below(0) )
//        {
//            update_panels();
//            return doupdate();
//        }
//        else
//        {
            return (ERR == wrefresh( this->m_pnl ))
                   ? EOF
                   : 0;
//        }
    }
    return EOF;
}

