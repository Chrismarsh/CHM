#pragma once

//std
#include <string>
#include <functional>
#include <fstream>
#include <iostream>

//json library
#include <json_spirit.h>





/// Primary unstructured mesh class
/**
 * Defines an unstructured mesh. This mesh is defined by x,y,z and any ancillary information at each vertex
 **/
class mesh
{

public:
	//Constructor
	mesh();

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
	//for each mesh, holds information about that mesh
	struct _mesh_config
	{
	    std::string name;
	    std::string file;
	};

	std::vector<_mesh_config> _meshes;

};

