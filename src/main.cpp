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

#include <iostream>
#include <string>
#include <boost/filesystem.hpp>

#define BOOST_SPIRIT_THREADSAFE

#include "core.hpp"

int main (int argc, char *argv[])
{
    core kernel;

    int ret = 0;
    try
    {
        kernel.init(argc, argv) ;

        kernel.run();

        kernel.end();
    }
    catch(chm_done& e)
    {
        kernel.end();
    }
    catch(const module_error& e)
    {
        ret = -1;
        SPDLOG_ERROR(boost::diagnostic_information(e));
    }
    catch(const exception_base& e)
    {
        ret = -1;
        SPDLOG_ERROR(boost::diagnostic_information(e));
    }
    catch(const std::exception& e)
    {
        ret = -1;
        SPDLOG_ERROR(boost::diagnostic_information(e));
    }
    catch(...)
    {
       SPDLOG_ERROR("Unknown exception");
       ret = -1;
    }

    // if we have an exception, ensure we tear down all of the MPI. In Non MPI mode this won't do anything
    if(ret == -1)
        kernel.end(true);
    else
        kernel.end();
    try
    {

#if BOOST_VERSION < 107400
        boost::filesystem::copy_file(kernel.log_file_path,kernel.o_path / "CHM.log", boost::filesystem::copy_option::overwrite_if_exists);
#else
        boost::filesystem::copy_file(kernel.log_file_path,kernel.o_path / "CHM.log", boost::filesystem::copy_options::overwrite_existing);
#endif
}

    catch(...)
    {

    }

    return ret;
}
