// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include "crocus_bridge.hpp"
#include "crocus_log.hpp"
#include "logger.hpp"
#include <array>
#include <cstddef>

namespace crocus
{
struct ColumnLogContext
{
    std::size_t face_id, step;
    const Forcing& met;
};

inline bool log_to_chm(LogLevel level, std::string_view message, const void* context)
{
    constexpr std::array levels{spdlog::level::debug, spdlog::level::info, spdlog::level::warn, spdlog::level::err,
                                spdlog::level::critical};
    const auto severity = levels.at(static_cast<std::size_t>(level));
    auto logger = spdlog::default_logger();
    if (!logger)
        return false;
    if (!logger->should_log(severity))
        return level != LogLevel::critical;
    if (context)
    {
        const auto& column = *static_cast<const ColumnLogContext*>(context);
        const auto& m = column.met;
        SPDLOG_LOGGER_CALL(logger, severity, "[Crocus] face={} step={} time={:04}-{:02}-{:02}T{:02}:{:02}:{:02} {}",
                           column.face_id, column.step, m.year, m.month, m.day, m.seconds / 3600, m.seconds / 60 % 60,
                           m.seconds % 60, message);
    }
    else
    {
        SPDLOG_LOGGER_CALL(logger, severity, "[Crocus] {}", message);
    }
    if (level == LogLevel::critical)
        logger->flush();
    return true;
}
} // namespace crocus
