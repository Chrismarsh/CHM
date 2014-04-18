#pragma once
#include <boost/crc.hpp>      // for boost::crc_basic, boost::crc_optimal
#include <boost/cstdint.hpp>  // for boost::uint16_t
#include <string>
#include <algorithm>
//#include <boost/utility.hpp>
#include <boost/algorithm/string/case_conv.hpp>
   // Used to compare the hashes of the main data structure
   //not case sensitive
    class crc_hash_compare
    {
    public:
            static size_t hash( const std::string& x);
            static bool equal(const std::string& s1, const std::string& s2);
    };
    