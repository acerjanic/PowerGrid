#pragma once

#include <ftxui/dom/elements.hpp>
#include <ftxui/dom/canvas.hpp>
#include <ftxui/dom/node.hpp>
#include <ftxui/screen/color.hpp>
#include <string>
#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <chrono>
#include "PGViewState.hpp"
#include "Sparkline.hpp"
#include "GraphicsProtocol.hpp"

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

// ---------------------------------------------------------------------------
// Image preview rendering
// ---------------------------------------------------------------------------

// Forward declaration for the overlay-aware variant
class ImagePlaceholderNode;

/// @brief Custom FTXUI Node that reserves space for an inline terminal image.
///
/// During Render(), it fills the area with spaces and records the absolute
/// terminal position so the ImageOverlay thread can emit escape sequences
/// at the correct coordinates.
class ImagePlaceholderNode : public ftxui::Node {
public:
    ImagePlaceholderNode(int cols, int rows,
                         const ImagePlaneView& plane,
                         ImageOverlay* overlay)
        : cols_(cols), rows_(rows), plane_(plane), overlay_(overlay) {}

    void ComputeRequirement() override {
        requirement_.min_x = cols_;
        requirement_.min_y = rows_;
    }

    void SetBox(ftxui::Box box) override {
        Node::SetBox(box);
        // Record the position for the overlay thread.
        // Box uses 0-based coordinates; terminal escape sequences use 1-based.
        region_.term_row = box.y_min + 1;
        region_.term_col = box.x_min + 1;
        region_.cols = box.x_max - box.x_min + 1;
        region_.rows = box.y_max - box.y_min + 1;
        region_.pixels = plane_.pixels.data();
        region_.px_w = plane_.nx;
        region_.px_h = plane_.ny;
    }

    void Render(ftxui::Screen& screen) override {
        // Fill the area with spaces to claim it from FTXUI
        for (int y = box_.y_min; y <= box_.y_max; ++y) {
            for (int x = box_.x_min; x <= box_.x_max; ++x) {
                if (x >= 0 && x < screen.dimx() && y >= 0 && y < screen.dimy()) {
                    screen.PixelAt(x, y).character = " ";
                }
            }
        }
    }

    /// Get the recorded image region for the overlay thread.
    const ImageRegion& region() const { return region_; }

private:
    int cols_, rows_;
    const ImagePlaneView& plane_;
    ImageOverlay* overlay_;
    ImageRegion region_;
};

/// @brief Render a single image plane as a labeled FTXUI Canvas element.
///
/// When proto == None, uses DrawBlock (half-block characters) for 2:1 vertical
/// sub-pixel resolution. When a graphics protocol is available, creates a
/// placeholder node and registers the region with the overlay thread.
///
/// @param plane    The image plane data.
/// @param proto    The detected graphics protocol (None = Canvas fallback).
/// @param overlay  Pointer to the overlay thread (nullptr if proto == None).
/// @param[out] regions  Collects ImageRegion entries for the overlay thread.
inline Element RenderSinglePlane(const ImagePlaneView& plane,
                                  GraphicsProto proto = GraphicsProto::None,
                                  std::vector<std::shared_ptr<ImagePlaceholderNode>>* placeholders = nullptr) {
    if (plane.pixels.empty() || plane.nx == 0 || plane.ny == 0)
        return text("");

    if (proto != GraphicsProto::None && placeholders) {
        // Graphics protocol path: placeholder node
        int cols = static_cast<int>(plane.nx);
        int rows = static_cast<int>(plane.ny) / 2; // Half-block: 2 pixels per row
        if (rows < 1) rows = 1;

        auto node = std::make_shared<ImagePlaceholderNode>(cols, rows, plane, nullptr);
        placeholders->push_back(node);

        return vbox({
            text(plane.label) | bold | hcenter,
            node,
        }) | border;
    }

    // Canvas fallback (works on all terminals)
    int canvas_w = static_cast<int>(plane.nx);
    int canvas_h = static_cast<int>(plane.ny);

    Canvas c(canvas_w, canvas_h);
    for (int y = 0; y < canvas_h; y++) {
        for (int x = 0; x < canvas_w; x++) {
            float v = plane.pixels[static_cast<size_t>(y) * plane.nx + static_cast<size_t>(x)];
            uint8_t gray = static_cast<uint8_t>(std::clamp(v, 0.0f, 1.0f) * 255.0f);
            c.DrawBlock(x, y, true, Color(gray, gray, gray));
        }
    }

    return vbox({
        text(plane.label) | bold | hcenter,
        canvas(c),
    }) | border;
}

/// @brief Render the complete content area (images + logs).
///
/// Handles three cases:
/// - No image: full-width scrolling logs
/// - 2D (1 plane): side-by-side image + logs
/// - 3D (3 planes): 2×2 quad layout — Axial, Coronal, Sagittal, Logs
///
/// Must be called with state.mu locked.
inline Element RenderContentArea(const PGViewState& state,
                                  GraphicsProto proto = GraphicsProto::None,
                                  std::vector<std::shared_ptr<ImagePlaceholderNode>>* placeholders = nullptr) {
    // No image: full-width logs
    if (!state.has_image || state.latest_image.planes.empty()) {
        return RenderLogs(state);
    }

    const auto& img = state.latest_image;
    std::string iter_label = "iter " + std::to_string(img.iter);

    if (img.planes.size() == 1) {
        // 2D: side-by-side — image | logs
        return vbox({
            hbox({
                RenderSinglePlane(img.planes[0], proto, placeholders) | flex,
                separator(),
                RenderLogs(state) | flex,
            }),
            text(iter_label) | dim | hcenter,
        });
    }

    // 3D MPR: 2×2 quad layout
    // Top row:    Axial    | Coronal
    // Bottom row: Sagittal | Logs
    Element top_left  = (img.planes.size() > 0)
        ? RenderSinglePlane(img.planes[0], proto, placeholders)
        : text("");
    Element top_right = (img.planes.size() > 1)
        ? RenderSinglePlane(img.planes[1], proto, placeholders)
        : text("");
    Element bot_left  = (img.planes.size() > 2)
        ? RenderSinglePlane(img.planes[2], proto, placeholders)
        : text("");
    Element bot_right = RenderLogs(state);

    return vbox({
        hbox({
            top_left | flex,
            separator(),
            top_right | flex,
        }),
        separator(),
        hbox({
            bot_left | flex,
            separator(),
            bot_right | flex,
        }),
        text(iter_label) | dim | hcenter,
    });
}
