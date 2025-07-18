// 	Copyright (C) 2011  Chris Marsh
//
// 	This program is free software: you can redistribute it and/or modify
// 	it under the terms of the GNU General Public License as published by
// 	the Free Software Foundation, either version 3 of the License, or
// 	(at your option) any later version.
//
// 	This program is distributed in the hope that it will be useful,
// 	but WITHOUT ANY WARRANTY; without even the implied warranty of
// 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// 	GNU General Public License for more details.
//
// 	You should have received a copy of the GNU General Public License
// 	along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <netcdf.h>
#include <netcdf_par.h>
#include <boost/mpi.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <string>

#include <boost/mpi.hpp>

#include "triangulation.hpp"
#include "timer.hpp"
#include "global.hpp"

typedef boost::shared_ptr<triangulation> mesh;

class ugrid
{
public:
    ugrid(mesh m, boost::shared_ptr<global> g);
    ~ugrid();

    /**
     * Check a netcdf C call's return value and convert to exception if needed
     * @param status
     */
    void nc_chk_ret(int status);

    void close_ugrid();
    void write_ugrid(std::vector<std::string> output_variables );
    void open_ugrid(std::vector<std::string> output_variables, std::string fname);
    void init_ugrid(std::vector<std::string> output_variables, std::string fname);

private:
    //holds the file id for the ugrid output netcdf
    int _ugrid_fid;

    //maps the variable string to the netcdf id to write to file
    std::map<std::string, int> _ugrid_id_var;

    std::string _fname;

    mesh _mesh;
    boost::shared_ptr<global> _global;

    boost::mpi::environment _mpi_env;
    boost::mpi::communicator _comm_world;
};

