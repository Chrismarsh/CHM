#pragma once


#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <errno.h>
#include <utility> // std::pair
#include <set>

//graph
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>



#include "logger.h"
#include "exception.hpp"

#include "triangulation.h"


#include "json_spirit.h"

#include "module_factory.hpp"
#include "station.hpp"

#include "interp_t_air.hpp"
#include "interp_rh.hpp"
#include "timer.hpp"
#include "global.hpp"


typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, 
      boost::property<boost::vertex_color_t, boost::default_color_type>
    > Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef std::pair<int, int> Edge;

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
   * Within this file are a collection of meshes that are expected to have the same number of x,y
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
    void config_debug(const json_spirit::Value& value);
    void config_modules(const json_spirit::Value& value);
    void config_forcing(const json_spirit::Value& value);
    void config_meshes(const json_spirit::Value& value);
    void config_matlab(const json_spirit::Value& value);
    void config_output(const json_spirit::Value& value);
    void config_global(const json_spirit::Value& value);
    
    /**
     * Initializes the logger and Matlab engine
     */
    core();
    ~core();

    void run();


	     

private:
    //current level of the logger. Defaults to debug, but can be changed via configuration settings
    log_level _log_level;
    //a text file log
    boost::shared_ptr< text_sink > _log_sink;
    
    //module factory for creating the specified modules
    module_factory _mfactory;
    
    //each station where observations are    
    tbb::concurrent_vector< boost::shared_ptr<station> > _stations;
    
    //main mesh object
    boost::shared_ptr< triangulation > _mesh;
    
    //matlab engine
    boost::shared_ptr<maw::matlab_engine> _engine;
    
       
    //holds all the modules that are to be run on each mesh element
    //pair as we also need to store the make order
    std::vector< std::pair<module,size_t> > _modules;
    
    boost::shared_ptr<global> _global;
    
    //calculates the order modules are to be run in
    void _determine_module_dep();
    
        
    //holds a unique list of all variables provided by all the modules;
    std::set<std::string> _module_provided_variable_list;
    
    
    class output_info
    {
    public:
        output_info()
        {
            plot = false;
            out_file = "";
            northing = 0;
            easting = 0;
            face = NULL;
        }
        enum output_type
        {
            timeseries,
            mesh
        };
        output_type type;
        bool plot;
        std::string out_file;
        double northing;
        double easting;
        std::vector<std::string> variables;
        mesh_elem face;
    };
    
    
    std::vector<output_info> _outputs;
};



