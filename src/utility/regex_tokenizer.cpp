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

#include "regex_tokenizer.hpp"

void regex_tokenizer::set_regex( std::string exp, bool case_sensitive )
{
        try
        {
                if(!case_sensitive)
                        re.assign(exp,boost::regex_constants::icase);
                else
                        re.assign(exp);
        }
        catch ( ... )
        {
                CHM_THROW_EXCEPTION(module_not_found, std::string("Invalid regular expression ") + exp);
        }

}

std::string regex_tokenizer::get_regex()
{
    return re.expression();
}



regex_tokenizer::regex_tokenizer() :
        FLOAT_REGEX("-?[0-9]+\\.?[0-9]*")
{

}

regex_tokenizer::regex_tokenizer( std::string regex, bool case_senstive /*= true*/ ) :
        FLOAT_REGEX("-?[0-9]+\\.?[0-9]*")
{
        try
        {
                set_regex(regex,case_senstive);
        }
        catch ( ... )
        {
                //nothing. don't want a zombie
        }

}
regex_tokenizer::~regex_tokenizer()
{

}
