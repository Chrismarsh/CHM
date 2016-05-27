//strips comments from a json file. needed as since boost 1.59 ptree no longer parses comments
//from here
//https://github.com/andrew-d/json-strip

#pragma once

#include <string>

std::string stripComments(const std::string& str, bool whitespace = true);