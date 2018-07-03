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

#include "core.h"




int main (int argc, char *argv[])
{
    core kernel;

    try
    {
        kernel.init(argc, argv) ;
        
        kernel.run();

        kernel.end();
    }
    catch( boost::exception& e)
    {
        kernel.end(); //make sure endwin() is called so ncurses cleans up

        LOG_ERROR << boost::diagnostic_information(e);

        return -1;
    }

    boost::filesystem::copy_file(kernel.log_file_path,kernel.o_path / "CHM.log", boost::filesystem::copy_option::overwrite_if_exists);

	return 0;
}
