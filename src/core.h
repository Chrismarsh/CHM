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
   * @param file The file to open
  **/
    void read_config_file(std::string file);

    core();
    ~core();

    bool is_debug();


    void run();
// 	log_level get_log_level();
// 	void set_log_level(log_level level);
	     

private:
    bool _is_debug;
    std::string _module_config;
    log_level _log_level;
    boost::shared_ptr< text_sink > _log_sink;
    module_factory _mfactory;
    
    
    std::vector < boost::shared_ptr<station> > _stations;
    
    boost::shared_ptr<mesh> _mesh;
    
};


