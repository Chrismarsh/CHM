#pragma once

#define BOOST_SPIRIT_THREADSAFE

//vtk includes
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkVersion.h>
#include <vtkStringArray.h>


//std includes
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
#include <chrono>
#include <algorithm>

//boost includes
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <tbb/concurrent_vector.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/bind.hpp>
namespace pt = boost::property_tree;
namespace po = boost::program_options;

#include <ogr_spatialref.h>
#ifdef _OPENMP
#include <omp.h>
#endif
//includes from CHM
#include "logger.hpp"
#include "exception.hpp"
#include "triangulation.hpp"
#include "filter_factory.h"
#include "module_factory.hpp"
#include "station.hpp"
#include "timer.hpp"
#include "global.hpp"
#include "str_format.h"
#include "ui.h"
#include "interpolation.h"
#include "readjson.hpp"
#include "version.h"
#include "math/coordinates.hpp"

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
    void config_options(const pt::ptree &value);
    void config_modules(const pt::ptree& value,const pt::ptree& config,std::vector<std::string> remove,std::vector<std::string> add);
    void config_meshes(const pt::ptree& value);
    void config_forcing(pt::ptree& value);
    void config_module_overrides(const pt::ptree& value);
    void config_parameters(pt::ptree &value);
    void config_matlab(const pt::ptree& value);
    void config_output(const pt::ptree& value);
    void config_global(const pt::ptree& value);

    // .first = config file to use
    // .second = extra options, if any.
    typedef boost::tuple<
            std::string, //config file path to load. defaults to CHM.config
            std::vector<std::pair<std::string,std::string>>, //insert or overide config value
            std::vector<std::string>, //remove config value
            std::vector<std::string>, //remove module
            std::vector<std::string>  // add module

    > cmdl_opt;

    cmdl_opt config_cmdl_options(int argc, char **argv);

    /**
     * Initializes the logger and Matlab engine
     */
    core();
    ~core();

    void run();
    void end();
    pt::ptree _cfg;

protected:
    //current level of the logger. Defaults to debug, but can be changed via configuration settings
    log_level _log_level;

    // if the users passes in a config file path that isn't the currently directory
    // e.g. CHM -f /some/other/path/CHM.json
    // then we need to affix every file IO (excep the log ?) with this path.
    boost::filesystem::path cwd_dir;
    boost::filesystem::path o_path; //path to output folder

    bool _output_station_ptv; //should we ouput the station ptv file? if we have no output section, then don't do this.

    //this is called via system call when the model is done to notify the user
    std::string _notification_script;

    //a text file log
    boost::shared_ptr< text_sink > _log_sink;
    boost::shared_ptr< text_sink > _cout_log_sink;
    
    //module factory for creating the specified modules
    module_factory _mfactory;
    filter_factory _filtfactory;

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
    std::vector< std::pair<std::string,std::string> > _overrides;
    boost::shared_ptr<global> _global;
    
    //calculates the order modules are to be run in
    void _determine_module_dep();

    interp_alg _interpolation_method;
        
    //holds a unique list of all variables provided by all the met files;
    std::set<std::string> _provided_var_met_files;
    //unique list of all variables provided by all the modules
    std::set<std::string> _provided_var_module;

    //unique set of all the paramters provided by the meshes
    std::set<std::string> _provided_parameters;
    std::set<std::string> _provided_initial_conditions;

    boost::posix_time::ptime* _start_ts;
    boost::posix_time::ptime* _end_ts;

    //a full timeseires per triangle
    bool _per_triangle_timeseries;

    struct point_mode_info
    {
        bool enable;
        std::string output;
        std::string forcing;

    } point_mode;


    class output_info
    {
    public:
        output_info()
        {
            frequency=1;
            fname = "";
            latitude = 0;
            longitude = 0;
            face = nullptr;
            name = "";
        }
        enum output_type
        {
            time_series,
            mesh
        };
        enum mesh_outputs
        {
            vtp,
            vtu,
            ascii
        };

        output_type type;
        std::string name;
        std::vector<mesh_outputs> mesh_output_formats;
        std::string fname;
        double latitude;
        double longitude;
        std::set<std::string> variables;
        mesh_elem face;
        timeseries ts;
        size_t frequency;

    };

    bool _enable_ui;
    ui _ui;
    std::vector<output_info> _outputs;
};



