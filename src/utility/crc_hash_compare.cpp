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

#include "crc_hash_compare.hpp"

size_t crc_hash_compare::hash( const std::string& x )
{
        boost::crc_32_type crc32;
        //std::string xlower = boost::algorithm::to_lower_copy<std::string>(x);
        crc32.process_bytes(x.c_str(),x.length());

        return crc32.checksum();
}

bool  crc_hash_compare::equal( const std::string& s1, const std::string& s2 )
{
       // std::string ss1 = boost::algorithm::to_lower_copy<std::string>(s1);
       // std::string ss2 = boost::algorithm::to_lower_copy<std::string>(s2);

        return s1 == s2;
}

