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

#include "jsonstrip.hpp"

#include <cctype>   // for isspace
#include <iostream> // for debugging

#define COMMENT_TYPE_NONE   0
#define COMMENT_TYPE_SINGLE 1
#define COMMENT_TYPE_MULTI  2

#undef STATE_DEBUG


typedef void (*stripFunction)(const std::string& str, size_t start, size_t end, std::string& out);


void stripWithoutWhitespace(const std::string& str, size_t start, size_t end, std::string& out) {
    // Do nothing.
    ((void)str);
    ((void)start);
    ((void)end);
    ((void)out);
    return;
}

void stripWithWhitespace(const std::string& str, size_t start, size_t end, std::string& out) {
    for (size_t i = start; i < end; i++) {
        char ch = str[i];

        if (isspace(ch)) {
            out.push_back(ch);
        } else {
            out.push_back(' ');
        }
    }
}


template <typename StripFunc>
std::string stripCommentsImpl(const std::string& str, StripFunc strip) {
    std::string ret;
    ret.reserve(str.length());

    char currentChar, nextChar;
    bool insideString = false;
    int commentType = COMMENT_TYPE_NONE;

    size_t offset = 0;
    for (size_t i = 0; i < str.length(); i++) {
        currentChar = str[i];

        if (i < str.length() - 1) {
            nextChar = str[i + 1];
        } else {
            nextChar = '\0';
        }

        // If we're not in a comment, check for a quote.
        if (commentType == COMMENT_TYPE_NONE && currentChar == '"') {
            bool escaped = false;

            // If the previous character was a single slash, and the one before
            // that was not (i.e. the previous character is escaping this quote
            // and is not itself escaped), then the quote is escaped.
            if (i >= 2 && str[i - 1] == '\\' && str[i - 2] != '\\') {
                escaped = true;
            }

            if (!escaped) {
                insideString = !insideString;
            }
        }

        if (insideString) {
            continue;
        }

        if (commentType == COMMENT_TYPE_NONE && currentChar == '/' && nextChar == '/') {
#if STATE_DEBUG
            std::cerr << "[no comment] found single comment, adding slice from " << offset << " to " << i << std::endl;
            std::cerr << "str:\n" << ret << std::endl;
#endif

            ret.append(str, offset, i - offset);
            offset = i;
            commentType = COMMENT_TYPE_SINGLE;

            // Skip second '/'
            i++;
        } else if (commentType == COMMENT_TYPE_SINGLE && currentChar == '\r' && nextChar == '\n') {
#if STATE_DEBUG
            std::cerr << "[single comment] exiting comment (1), stripping slice from " << offset << " to " << i << std::endl;
            std::cerr << "str:\n" << ret << std::endl;
#endif

            // Skip '\r'
            i++;

            commentType = COMMENT_TYPE_NONE;
            strip(str, offset, i, ret);
            offset = i;

            continue;
        } else if (commentType == COMMENT_TYPE_SINGLE && currentChar == '\n') {
#if STATE_DEBUG
            std::cerr << "[single comment] exiting comment (2), stripping slice from " << offset << " to " << i << std::endl;
            std::cerr << "str:\n" << ret << std::endl;
#endif

            commentType = COMMENT_TYPE_NONE;
            strip(str, offset, i, ret);
            offset = i;
        } else if (commentType == COMMENT_TYPE_NONE && currentChar == '/' && nextChar == '*') {
#if STATE_DEBUG
            std::cerr << "[no comment] found multi comment, adding slice from " << offset << " to " << i << std::endl;
            std::cerr << "str:\n" << ret << std::endl;
#endif

            ret.append(str, offset, i - offset);
            offset = i;
            commentType = COMMENT_TYPE_MULTI;

            // Skip the '*'
            i++;
            continue;
        } else if (commentType == COMMENT_TYPE_MULTI && currentChar == '*' && nextChar == '/') {
#if STATE_DEBUG
            std::cerr << "[multi comment] exiting comment, stripping slice from " << offset << " to " << (i + 1) << std::endl;
            std::cerr << "str:\n" << ret << std::endl;
#endif

            // Skip '*'
            i++;

            commentType = COMMENT_TYPE_NONE;
            strip(str, offset, i + 1, ret);
            offset = i + 1;
            continue;
        }
    }

#if STATE_DEBUG
    std::cerr << "[done] adding substring from " << offset << " to end" << std::endl;
#endif

    ret.append(str, offset, str.length() - offset);
    ret.shrink_to_fit();
    return ret;
}


std::string stripComments(const std::string& str, bool whitespace) {
    if (whitespace) {
        return stripCommentsImpl(str, stripWithWhitespace);
    } else {
        return stripCommentsImpl(str, stripWithoutWhitespace);
    }
}