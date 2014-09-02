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

