/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file PGLog.hpp
/// @brief Structured logging facade and shared terminal progress display for PowerGrid.

#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/fmt/fmt.h>
#include <indicators/block_progress_bar.hpp>
#include <indicators/cursor_control.hpp>
#include <indicators/color.hpp>

#include <memory>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Logging
// ---------------------------------------------------------------------------

/// @brief Initialise the PowerGrid logger with two sinks.
inline void PG_LOG_INIT(const std::string& level_str = "info",
                        const std::string& log_path  = "powergrid.log") {
    auto level = spdlog::level::from_str(level_str);

    auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
    stderr_sink->set_level(level);

    auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
        log_path, 10UL * 1024UL * 1024UL, 3);
    file_sink->set_level(spdlog::level::debug);

    std::vector<spdlog::sink_ptr> sinks{stderr_sink, file_sink};
    auto logger = std::make_shared<spdlog::logger>("powergrid", sinks.begin(), sinks.end());
    logger->set_level(spdlog::level::debug);
    spdlog::set_default_logger(logger);
}

#define PG_DEBUG(...) SPDLOG_DEBUG(__VA_ARGS__)
#define PG_INFO(...)  SPDLOG_INFO(__VA_ARGS__)
#define PG_WARN(...)  SPDLOG_WARN(__VA_ARGS__)
#define PG_ERROR(...) SPDLOG_ERROR(__VA_ARGS__)

// ---------------------------------------------------------------------------
// Progress bars
// ---------------------------------------------------------------------------

using PGProgressBar = indicators::BlockProgressBar;

inline std::vector<std::shared_ptr<PGProgressBar>>& PG_BARS() {
    static std::vector<std::shared_ptr<PGProgressBar>> bars;
    return bars;
}

inline size_t PG_PROGRESS_ADD(std::shared_ptr<PGProgressBar> bar) {
    PG_BARS().push_back(bar);
    return PG_BARS().size() - 1;
}

inline void PG_PROGRESS_TICK(size_t idx, const std::string& postfix = "") {
    auto& bar = *PG_BARS()[idx];
    if (!postfix.empty())
        bar.set_option(indicators::option::PostfixText{postfix});
    bar.tick();
}

inline void PG_PROGRESS_DONE(size_t idx) {
    PG_BARS()[idx]->mark_as_completed();
}
