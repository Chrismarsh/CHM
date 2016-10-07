#include "readjson.hpp"

pt::ptree read_json(const std::string& path)
{

    // the imbue appears to fix some oddness in how the json parser parses the str.
    // http://stackoverflow.com/q/35275314/410074
    // putting it here fixed a segfault very deep in the json parse

    std::ifstream in(path);
    in.imbue(std::locale());

    std::stringstream json_file;
    if (in.is_open())
    {
        std::string line;
        while ( getline (in,line) )
        {
            json_file << line << '\n';
        }
        in.close();
    }
    else
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info("Unable to open " + path));

    }

    //strip out json comments for ptress
    bool whitespace = true;
    auto stripped = stripComments(json_file.str(), whitespace);
//    auto json_file_stripped = std::stringstream(stripped);
    std::stringstream json_file_stripped (stripped);
    json_file_stripped.imbue(std::locale());
    pt::ptree config;
    try
    {
        pt::read_json( json_file_stripped, config);
    }
    catch (pt::json_parser_error &e)
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info(
                "Error reading file: " + path + " on line: " + std::to_string(e.line()) + " with error: " +
                e.message()));
    }

    return config;
}