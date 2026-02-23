#pragma once

#include <ftxui/dom/elements.hpp>
#include <ftxui/screen/color.hpp>
#include <string>
#include <vector>
#include <deque>
#include <cmath>
#include <ctime>
#include <chrono>
#include "PGViewState.hpp"
#include "Sparkline.hpp"

using namespace ftxui;

/// @brief Format seconds as "Xm Ys" or "Xh Ym" for display.
inline std::string format_eta(double seconds) {
    if (seconds < 0 || std::isnan(seconds) || std::isinf(seconds))
        return "";
    int total = static_cast<int>(std::round(seconds));
    if (total < 60)
        return std::to_string(total) + "s";
    int m = total / 60;
    int s = total % 60;
    if (m < 60)
        return std::to_string(m) + "m" + (s > 0 ? std::to_string(s) + "s" : "");
    int h = m / 60;
    m = m % 60;
    return std::to_string(h) + "h" + std::to_string(m) + "m";
}

/// @brief Format wall clock completion time as "~HH:MM".
inline std::string format_wall_clock(double eta_sec) {
    if (eta_sec < 0 || std::isnan(eta_sec) || std::isinf(eta_sec))
        return "";
    auto now = std::chrono::system_clock::now();
    auto finish = now + std::chrono::seconds(static_cast<long>(std::round(eta_sec)));
    auto tt = std::chrono::system_clock::to_time_t(finish);
    struct tm local;
    localtime_r(&tt, &local);
    char buf[8];
    std::snprintf(buf, sizeof(buf), "~%02d:%02d", local.tm_hour, local.tm_min);
    return buf;
}

/// @brief Format a double in scientific notation, e.g. "1.23e-2".
inline std::string format_sci(double v) {
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%.2e", v);
    return buf;
}

/// @brief Render the progress bar section.
///
/// Returns an Element containing all active (non-completed) progress bars,
/// each with gauge, percentage, postfix, ETA, and sparklines for metrics.
inline Element RenderProgress(const PGViewState& state) {
    // Must be called with state.mu locked
    Elements bar_rows;

    for (size_t i = 0; i < state.bars.size(); ++i) {
        const auto& bar = state.bars[i];
        if (bar.completed) continue;

        // Label (padded to 6 chars)
        std::string label = bar.label;
        while (label.size() < 6) label += ' ';

        // Percentage
        int pct = static_cast<int>(bar.fraction() * 100.0f);
        std::string pct_str = std::to_string(pct) + "%";

        // ETA info
        double eta = bar.eta_seconds();
        std::string eta_str = (eta >= 0) ? "ETA " + format_eta(eta) : "";
        std::string wall_str = (eta >= 0) ? format_wall_clock(eta) : "";

        // Build the bar row
        bar_rows.push_back(
            hbox({
                text(label) | bold,
                text(" "),
                gauge(bar.fraction()) | flex | color(Color::Cyan),
                text(" " + pct_str),
                text("  " + bar.postfix),
                text("  " + eta_str) | dim,
                text(" " + wall_str) | dim,
            })
        );

        // Metrics sparkline row (only if we have metrics data)
        if (!bar.metrics.empty()) {
            // Extract error_norm and penalty series
            std::deque<double> err_series, pen_series;
            for (const auto& m : bar.metrics) {
                err_series.push_back(m.error_norm);
                pen_series.push_back(m.penalty);
            }

            std::string err_spark = render_sparkline(err_series);
            std::string pen_spark = render_sparkline(pen_series);
            std::string err_val = format_sci(bar.metrics.back().error_norm);
            std::string pen_val = format_sci(bar.metrics.back().penalty);

            bar_rows.push_back(
                hbox({
                    text("  err: ") | dim,
                    text(err_val) | color(Color::Yellow),
                    text("  " + err_spark),
                    text("    pen: ") | dim,
                    text(pen_val) | color(Color::Magenta),
                    text("  " + pen_spark),
                })
            );
        }
    }

    if (bar_rows.empty()) {
        return text(""); // No active bars
    }

    return vbox(std::move(bar_rows));
}

/// @brief Color for a log level.
inline Color log_level_color(const std::string& level) {
    if (level == "error" || level == "critical") return Color::Red;
    if (level == "warn")  return Color::Yellow;
    if (level == "debug") return Color::GrayDark;
    if (level == "info")  return Color::Green;
    return Color::White; // raw
}

/// @brief Format a log level tag like "[INFO ]".
inline std::string log_level_tag(const std::string& level) {
    if (level == "info")     return "[INFO ]";
    if (level == "debug")    return "[DEBUG]";
    if (level == "warn")     return "[WARN ]";
    if (level == "error")    return "[ERROR]";
    if (level == "critical") return "[CRIT ]";
    if (level == "trace")    return "[TRACE]";
    return "[     ]"; // raw
}

/// @brief Render the scrolling log section.
///
/// Returns an Element containing the last N log lines, auto-scrolled to bottom.
inline Element RenderLogs(const PGViewState& state, size_t max_visible = 200) {
    // Must be called with state.mu locked
    Elements log_lines;

    size_t start = 0;
    if (state.logs.size() > max_visible)
        start = state.logs.size() - max_visible;

    for (size_t i = start; i < state.logs.size(); ++i) {
        const auto& entry = state.logs[i];
        Color c = log_level_color(entry.level);
        std::string tag = log_level_tag(entry.level);

        log_lines.push_back(
            hbox({
                text(entry.timestamp) | dim,
                text(" " + tag + " ") | color(c),
                text(entry.message),
            })
        );
    }

    if (log_lines.empty()) {
        log_lines.push_back(text("Waiting for input...") | dim);
    }

    return vbox(std::move(log_lines)) | vscroll_indicator | frame | flex;
}
