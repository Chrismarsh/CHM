
#include "config.h"

void config::open(std::string file)
{
    try
    {
        pt::read_json(file, *(reinterpret_cast<pt::ptree*>(this)));
    }
    catch (pt::json_parser_error &e)
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info( "Error reading file: " + e.filename() + " on line: " + std::to_string(e.line()) + "with error: " + e.message()));
    }


}

config::config()
{

}

config::~config()
{

}
