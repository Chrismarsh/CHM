#pragma once

#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>

#include <string>

struct exception_base: virtual std::exception, virtual boost::exception { };
struct io_error: virtual exception_base { };
struct file_read_error: virtual io_error { };

