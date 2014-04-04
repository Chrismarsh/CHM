#pragma once

#include "triangulation.h"
#include "logger.h"
#include "json_spirit.h"

#include "libmaw.h"

#include <boost/shared_ptr.hpp>
class mesh
{
public:
    mesh(boost::shared_ptr<maw::matlab_engine> engine);
    ~mesh();
    void add_mesh(std::string file, std::string ID);
    
private:

    

    boost::shared_ptr<maw::matlab_engine> _engine;
    boost::shared_ptr<maw::graphics> _gfx;
    std::vector< boost::shared_ptr< triangulation > > _meshes;
    
};
