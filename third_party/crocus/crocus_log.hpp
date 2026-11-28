// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include <string_view>

namespace crocus
{
// Values shared with crocus_logging.f90.
enum class LogLevel : int
{
    debug = 0,
    info,
    warning,
    error,
    critical
};

// The sink is synchronous and applies only to native calls on the current thread.
class LogScope
{
  public:
    // Return false to request standard-output/error fallback.
    using Writer = bool (*)(LogLevel, std::string_view, const void*);
    explicit LogScope(Writer writer, const void* context = nullptr) noexcept;
    ~LogScope();
    LogScope(const LogScope&) = delete;
    LogScope& operator=(const LogScope&) = delete;

  private:
    Writer previous_writer_;
    const void* previous_context_;
};
} // namespace crocus
