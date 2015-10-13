//
// Created by Chris Marsh on 2015-10-13.
//

#pragma once



#include <string>
#include <cstdio>


//http://stackoverflow.com/a/26197300/410074
template <typename... Ts>
std::string str_format (const std::string &fmt, Ts... vs)
{
    unsigned required = std::snprintf(nullptr, 0, fmt.c_str(), vs...) + 1;
        // See comments: the +1 is necessary, while the first parameter
        //               can also be set to nullptr

    char bytes[required];
    std::snprintf(bytes, required, fmt.c_str(), vs...);

    return std::string(bytes);
}