#pragma once


#include <string>
#include <fstream>
#include <vector>
#include <errno.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>



#include "logger.h"
#include "exception.hpp"


#include "json_spirit.h"

#include "module_factory.hpp"
#include "station.hpp"
#include "mesh.hpp"
#include "interp_t_air.hpp"
#include "interp_rh.hpp"
#include "timer.hpp"
/// The main model core
/**
 * The main model core, handles initialization of the model
**/

class core
{
public:
    
    
  /**
   * Reads the main JSON configuration file. It assumes the base of the JSON is an object. That is, the file
   * starts with { ... }.
    Within this file are a collection of meshes that are expected to have the same number of x,y
    * points. This is done so that, for example, elevation, forest cover, sky-view factor, etc 
    * may be added individually. Generation of the meshes should be done via the utilities for this.
    * An example of mesh.config is:
    * \code	
    *  {
    *    "meshes":
    *    {
    *            "DEM":
    *            {
    *                    "file": "mesh.asc"
    *            }
    *            ,
    *            "Veg":
    *            {
    *                    "file": "veg.asc"
    *            },
    *            "svf":
    *            {
    *                    "file": "svf.asc"
    *            }
    *    }	
    *   }
    *   \endcode
   * @param file The file to open
  **/
    void read_config_file(std::string file);

    core();
    ~core();

    bool is_debug();


    void run();

	     

private:
    bool _is_debug;

    log_level _log_level;
    boost::shared_ptr< text_sink > _log_sink;
    module_factory _mfactory;
    
    
    tbb::concurrent_vector< boost::shared_ptr<station> > _stations;
    
    boost::shared_ptr<mesh> _mesh;
    
    boost::shared_ptr<maw::matlab_engine> _engine;
    
    
    
    //holds all the modules that are to be run on each mesh element
    std::vector< module > _modules;
    
    
    
};


