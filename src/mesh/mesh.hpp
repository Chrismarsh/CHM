#pragma once

#include "triangulation.h"
#include "logger.h"
#include "json_spirit.h"

#include "libmaw.h"

#include "crc_hash_compare.hpp"

#include <boost/shared_ptr.hpp>


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
        
    typedef tbb::concurrent_hash_map<std::string, boost::shared_ptr< triangulation > , crc_hash_compare> mesh_hash;
    boost::shared_ptr< triangulation > _mesh;
    
};
