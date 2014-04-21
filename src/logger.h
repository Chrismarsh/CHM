#pragma once

#include <iostream>
#include <fstream>

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>


#include <boost/log/utility/empty_deleter.hpp>
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


#define	LOG_DEBUG       BOOST_LOG_NAMED_SCOPE(__PRETTY_FUNCTION__) BOOST_LOG_SEV(logger::get(), log_level::debug) 
#define	LOG_WARNING 	BOOST_LOG_NAMED_SCOPE(__PRETTY_FUNCTION__) BOOST_LOG_SEV(logger::get(), log_level::warning) 
#define	LOG_ERROR 	BOOST_LOG_NAMED_SCOPE(__PRETTY_FUNCTION__) BOOST_LOG_SEV(logger::get(), log_level::error) 
#define	LOG_INFO 	BOOST_LOG_NAMED_SCOPE(__PRETTY_FUNCTION__) BOOST_LOG_SEV(logger::get(), log_level::info) 