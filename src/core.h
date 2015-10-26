#pragma once


#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <errno.h>
#include <utility> // std::pair
#include <set>
#include <chrono>
#include <map>
#include <stdio.h>
#include <cstdlib>

//graph
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>
namespace pt = boost::property_tree;
namespace po = boost::program_options;

#include "logger.hpp"
#include "exception.hpp"

#include "triangulation.hpp"



#include "module_factory.hpp"
#include "station.hpp"

#include "timer.hpp"
#include "global.hpp"
#include "str_format.h"

struct vertex{
    std::string name;
};

struct edge{
    std::string variable;
};

//typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
//      boost::property<boost::vertex_color_t, boost::default_color_type>
//    > Graph;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,boost::property<boost::vertex_index_t,int,vertex>,edge> Graph;


//from http://stackoverflow.com/questions/11369115/how-to-print-a-graph-in-graphviz-with-multiple-properties-displayed
template <class VariableMap>
class edge_writer {
public:
    edge_writer(VariableMap v) : vm(v) {}
    template <class Edge>
    void operator()(ostream &out, const Edge& e) const {
        out << "[label=\"" << vm[e] << "\", edgetype=" << vm[e] << "]";
    }
private:
    VariableMap vm;
};

template <class VariableMap>
inline edge_writer<VariableMap>
make_edge_writer(VariableMap v) {
    return edge_writer<VariableMap>(v);
}

//typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
//typedef std::pair<int, int> Edge;


/**
 * The main model core, handles initialization of the model
**/

class core
{
    friend class CoreTest;
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
    void init(int argc, char **argv);
    void config_debug(const pt::ptree& value);
    void config_modules(const pt::ptree& value,const pt::ptree& config);
    void config_meshes(const pt::ptree& value);
    void config_forcing(const pt::ptree& value);

    void config_matlab(const pt::ptree& value);
    void config_output(const pt::ptree& value);
    void config_global(const pt::ptree& value);

    // .first = config file to use
    // .second = extra options, if any.
    typedef boost::tuple<
            std::string,
            std::vector<std::pair<std::string,std::string>>,
            std::vector<std::string>
    > cmdl_opt;
    cmdl_opt config_cmdl_options(int argc, char **argv);

    /**
     * Initializes the logger and Matlab engine
     */
    core();
    ~core();

    void run();
    pt::ptree _cfg;

protected:
    //current level of the logger. Defaults to debug, but can be changed via configuration settings
    log_level _log_level;

    //a text file log
    boost::shared_ptr< text_sink > _log_sink;
    
    //module factory for creating the specified modules
    module_factory _mfactory;


    //main mesh object
    boost::shared_ptr< triangulation > _mesh;
    
#ifdef MATLAB
    //matlab engine
    boost::shared_ptr<maw::matlab_engine> _engine;
#endif
       
    //holds all the modules that are to be run on each mesh element
    //pair as we also need to store the make order
    std::vector< std::pair<module,size_t> > _modules;
    std::vector< std::vector < module> > _chunked_modules;
    
    boost::shared_ptr<global> _global;
    
    //calculates the order modules are to be run in
    void _determine_module_dep();
    
        
    //holds a unique list of all variables provided by all the met files;
    std::set<std::string> _provided_var_met_files;
    //unique list of all variables provided by all the modules
    std::set<std::string> _provided_var_module;
    
    
    class output_info
    {
    public:
        output_info()
        {
            fname = "";
            northing = 0;
            easting = 0;
            face = NULL;
        }
        enum output_type
        {
            timeseries,
            mesh
        };
        enum mesh_outputs
        {
            vtp,
            vtu,
            ascii
        };
        output_type type;
        std::vector<mesh_outputs> mesh_output_formats;
        std::string fname;
        double northing;
        double easting;
        std::vector<std::string> variables;
        mesh_elem face;
    };
    
    
    std::vector<output_info> _outputs;
};



