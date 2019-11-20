//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#pragma once

#include <iostream>
#include <fstream>

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/version.hpp>
#if (BOOST_VERSION / 100 % 1000) < 56
    #include <boost/log/utility/empty_deleter.hpp>
#else
    #include <boost/core/null_deleter.hpp>
#endif

#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/exception_handler.hpp>

#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/severity_feature.hpp>
#include <boost/log/sources/record_ostream.hpp>

#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>



#include <boost/date_time/posix_time/posix_time.hpp>

#include <boost/log/attributes/named_scope.hpp>
#include <boost/log/expressions/formatters/named_scope.hpp>

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace sinks = boost::log::sinks;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;

 

 enum log_level
 {
     verbose,
     debug,
     warning,
     info,
     error
 };

// The formatting logic for the severity level
template< typename CharT, typename TraitsT >
inline std::basic_ostream< CharT, TraitsT >& operator<< (std::basic_ostream< CharT, TraitsT >& strm, log_level lvl)
{
    static const char* const str[] =
    {
        "verbose",
	    "debug",
	    "warning",
        "info",
	    "error"
    };
    if (static_cast< std::size_t >(lvl) < (sizeof(str) / sizeof(*str)))
	strm << str[lvl];
    else
	strm << static_cast< int >(lvl);
    return strm;
}

typedef src::severity_logger_mt<log_level> sev_logger;
typedef sinks::synchronous_sink< sinks::text_ostream_backend> text_sink;

BOOST_LOG_INLINE_GLOBAL_LOGGER_DEFAULT(logger, sev_logger)
BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", log_level)
//BOOST_LOG_ATTRIBUTE_KEYWORD(scope, "Scope", attrs::named_scope::value_type)


#define	LOG_VERBOSE       BOOST_LOG_NAMED_SCOPE(__PRETTY_FUNCTION__) BOOST_LOG_SEV(logger::get(), verbose) 
#define	LOG_DEBUG       BOOST_LOG_NAMED_SCOPE(__PRETTY_FUNCTION__) BOOST_LOG_SEV(logger::get(), debug) 
#define	LOG_WARNING 	BOOST_LOG_NAMED_SCOPE(__PRETTY_FUNCTION__) BOOST_LOG_SEV(logger::get(), warning) 
#define	LOG_ERROR 	BOOST_LOG_NAMED_SCOPE(__PRETTY_FUNCTION__) BOOST_LOG_SEV(logger::get(), error) 
#define	LOG_INFO 	BOOST_LOG_NAMED_SCOPE(__PRETTY_FUNCTION__) BOOST_LOG_SEV(logger::get(), info) 