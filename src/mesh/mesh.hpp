#pragma once

#include "triangulation.h"
#include "logger.h"
#include "json_spirit.h"

#include "libmaw.h"

#include <boost/shared_ptr.hpp>
#include <boost/crc.hpp>      // for boost::crc_basic, boost::crc_optimal
#include <boost/cstdint.hpp>  // for boost::uint16_t

#include <tbb/concurrent_hash_map.h>

#include <boost/iterator.hpp>

class mesh
{
public:
    mesh(boost::shared_ptr<maw::matlab_engine> engine);
    ~mesh();
    void add_mesh(std::string file, std::string ID);
    void plot(std::string ID);
    
    size_t size();
    
    typedef triangle mesh_elem;
    
    mesh_elem& operator()(size_t t);
    
private:
    boost::shared_ptr<maw::matlab_engine> _engine;
    boost::shared_ptr<maw::graphics> _gfx;
    
    // Used to compare the hashes of the main data structure
        //not case sensitive
    class hash_compare
    {
    public:
            static size_t hash( const std::string& x);
            static bool equal(const std::string& s1, const std::string& s2);
    };
    
    
    typedef tbb::concurrent_hash_map<std::string, boost::shared_ptr< triangulation > , hash_compare> mesh_hash;
    boost::shared_ptr< triangulation > _mesh;
    
};
