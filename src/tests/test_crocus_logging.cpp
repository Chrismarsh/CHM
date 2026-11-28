// SPDX-License-Identifier: GPL-3.0-or-later
#include "crocus_log_sink.hpp"
#include <gtest/gtest.h>
#include <spdlog/sinks/ostream_sink.h>
#include <sstream>

TEST(Crocus, FlushesNativeFatalProfileToRankLogger)
{
    struct Buffer : std::stringbuf
    {
        int flushes = 0;
        int sync() override
        {
            ++flushes;
            return std::stringbuf::sync();
        }
    } buffer;
    std::ostream stream(&buffer);
    auto sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(stream);
    auto logger = std::make_shared<spdlog::logger>("rank 7", sink);
    logger->set_pattern("[%n] [%l] %v");
    auto previous = spdlog::default_logger();
    spdlog::set_default_logger(logger);
    crocus::Forcing met{2017, 10, 1, 3600};
    const crocus::ColumnLogContext column{1842, 1, met};
    crocus::log_to_chm(crocus::LogLevel::critical, "SNOWCRO: erreur tempe snow\nTEMP (K) 70.6", &column);
    spdlog::set_default_logger(previous);

    EXPECT_NE(buffer.str().find("[rank 7] [critical] [Crocus] face=1842 step=1 time=2017-10-01T01:00:00 "
                                "SNOWCRO: erreur tempe snow\nTEMP (K) 70.6"),
              std::string::npos);
    EXPECT_GT(buffer.flushes, 0);
}
