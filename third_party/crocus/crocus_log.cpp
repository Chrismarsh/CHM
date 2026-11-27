// SPDX-License-Identifier: GPL-3.0-or-later
#include "crocus_log.hpp"

#include <cstddef>
#include <utility>

namespace crocus
{
namespace
{
thread_local LogScope::Writer log_writer = nullptr;
thread_local const void* log_context = nullptr;
} // namespace

LogScope::LogScope(Writer writer, const void* context) noexcept
    : previous_writer_(std::exchange(log_writer, writer)), previous_context_(std::exchange(log_context, context))
{
}

LogScope::~LogScope()
{
    log_writer = previous_writer_;
    log_context = previous_context_;
}

// Returning zero asks Fortran to write the original diagnostic to its standard unit.
extern "C" int chm_crocus_write_log(int level, const char* text, std::size_t size) noexcept
{
    if (!log_writer)
        return 0;
    try
    {
        return log_writer(static_cast<LogLevel>(level), std::string_view(text, size), log_context) ? 1 : 0;
    }
    catch (...)
    {
        return 0;
    }
}
} // namespace crocus
