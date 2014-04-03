#pragma once

#include "triangulation.h"
#include "logger.h"
#include "json_spirit.h"

#include "libmaw.h"

#include <boost/shared_ptr.hpp>
class mesh
{
public:
    mesh();
    ~mesh();
    /** Read in a mesh.config file.  This is a JSON based file. 
	 * Within this file are a collection of meshes that are expected to have the same number of x,y
	 * points. This is done so that, for example, elevation, forest cover, sky-view factor, etc 
	 * may be added individually. Generation of the meshes should be done via the utilities for this.
	 * An example of mesh.config is:
	\code	{
		"meshes":
		{
			"DEM":
			{
				"file": "mesh.asc"
			}
			,
			"Veg":
			{
				"file": "veg.asc"
			},
			"svf":
			{
				"file": "svf.asc"
			}
		}	
		}
	\endcode
	 * @param file The location of the mesh.config file to read in
	**/
	void read_mesh(std::string file);
private:
    maw::matlab_engine* _engine;
    maw::graphics* _gfx;
    
    //for each mesh, holds information about that mesh
    struct _mesh_config
    {
        std::string name;
        std::string file;
    };

    std::vector<_mesh_config> _meshes;
    
};
