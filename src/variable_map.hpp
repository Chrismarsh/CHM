#pragma once

#include <string>

#include <tbb/concurrent_hash_map.h>
#include "utility/crc_hash_compare.hpp"

typedef tbb::concurrent_hash_map<std::string,std::string,crc_hash_compare> var_hashmap;

class var
{
public:
    var();
    ~var();
    std::string operator()(std::string variable);
    void init_from_file(std::string path);
    
   var_hashmap _varmap;
};



//this relates internal variable names to the names in met files/other input
//namespace variables
//{
//  const  std::string RH = "RH";
//  const  std::string Tair = "t";
//  const  std::string datetime = "datetime";
//};


//#define RH "rh"
//#define TAIR "t"
//#define DATE "datetime"