#pragma once

#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>

#include <string>

struct exception_base: virtual std::exception, virtual boost::exception { };

// IO errors
struct io_error: virtual exception_base { };
struct file_read_error: virtual io_error { };
struct file_write_error: virtual io_error { };


// Module errors
struct module_error : virtual exception_base { };
struct module_not_found : virtual module_error { };
struct no_modules_defined : virtual module_error {};



typedef boost::error_info<struct errstr_info_,std::string> errstr_info;






//boost::errinfo_type_info_name("LOL")
