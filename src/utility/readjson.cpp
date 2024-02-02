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

#include "readjson.hpp"

pt::ptree read_json(const std::string& path)
{

    // the imbue appears to fix some oddness in how the json parser parses the str.
    // http://stackoverflow.com/q/35275314/410074
    // putting it here fixed a segfault very deep in the json parse

    std::ifstream in(path);
    in.imbue(std::locale());

    std::stringstream json_file;
    if (in.is_open())
    {
        std::string line;
        while ( getline (in,line) )
        {
            json_file << line << '\n';
        }
        in.close();
    }
    else
    {
        CHM_THROW_EXCEPTION(config_error, "Unable to open " + path);

    }

    //strip out json comments for ptress
    bool whitespace = true;
    auto stripped = stripComments(json_file.str(), whitespace);
//    auto json_file_stripped = std::stringstream(stripped);
    std::stringstream json_file_stripped (stripped);
    json_file_stripped.imbue(std::locale());
    pt::ptree config;
    try
    {
        pt::read_json( json_file_stripped, config);
    }
    catch (pt::json_parser_error &e)
    {
        CHM_THROW_EXCEPTION(config_error, "Error reading file: " + path + " on line: " + std::to_string(e.line()) + " with error: " + e.message());
    }

    return config;
}