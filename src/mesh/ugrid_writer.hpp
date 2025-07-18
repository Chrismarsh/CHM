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


class ugrid_writer
{
public:
    ugrid_writer(mesh m, boost::shared_ptr<global> g,  bool write_parameters, std::string fname);
    ~ugrid_writer();

    /**
     * Check a netcdf C call's return value and convert to exception if needed
     * @param status
     */
    void nc_chk_ret(int status);

    void close_ugrid();
    void write_ugrid(const std::vector<std::string>& output_variables);
    void open_ugrid(const std::vector<std::string>& output_variables);
    void init_ugrid(const std::vector<std::string>& output_variables);

private:
    //holds the file id for the ugrid output netcdf
    int _ugrid_fid;

    //maps the variable string to the netcdf id to write to file
    std::map<std::string, int> _ugrid_id_var;

    std::string _fname;

    // track the number of outputs that have been done to correctly compute the offset in the nc
    size_t _time_index;

    mesh _mesh;
    boost::shared_ptr<global> _global;

    boost::mpi::environment _mpi_env;
    boost::mpi::communicator _comm_world;

    bool _write_parameters;
};

