#pragma once

#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>

#include <string>

struct exception_base: virtual std::exception, virtual boost::exception { };

//matlab
struct matlab_error: virtual exception_base {};
struct matlab_null_return : virtual matlab_error{};
struct matlab_engine_failure : virtual matlab_error{};

// IO errors
struct io_error: virtual exception_base { };
struct file_read_error: virtual io_error { };
struct file_write_error: virtual io_error { };


// Module errors
struct module_error : virtual exception_base { };
struct module_not_found : virtual module_error { };
struct no_modules_defined : virtual module_error {};
struct module_data_error : virtual module_error {};

//regex tokenizer errors
struct tokenizer_error : virtual exception_base{};
struct invalid_rexp : virtual tokenizer_error {} ;
struct bad_lexical_cast : virtual tokenizer_error {};

//forcing data errors
struct forcing_error : virtual exception_base{};
struct forcing_insertion_error : virtual forcing_error{};
struct forcing_lookup_error : virtual forcing_error{};
struct forcing_badcast : virtual forcing_error{};
struct forcing_no_regexmatch : virtual forcing_error{};
struct forcing_timestep_mismatch : virtual forcing_error{};
struct forcing_timestep_notfound : virtual forcing_error{};
struct forcing_no_stations : virtual forcing_error{};

//interpolation errors
struct interp_error : virtual exception_base{};
struct interp_unknown_type : virtual interp_error {};

//configuration errors
struct config_error : virtual exception_base{};
struct missing_value_error : virtual exception_base{};

typedef boost::error_info<struct errstr_info_,std::string> errstr_info;

//mesh errors
struct mesh_error : virtual exception_base{};
struct mesh_insertion_error : virtual mesh_error{};
struct mesh_lookup_error : virtual mesh_error{};

struct interpolation_error : virtual exception_base{};

//boost::errinfo_type_info_name("LOL")
