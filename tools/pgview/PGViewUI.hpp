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
#include <functional>
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

/// @brief Custom FTXUI Node that reserves space for an inline terminal image.
///
/// During Render(), it fills the area with spaces (so FTXUI doesn't draw
/// there) and pushes an ImageRegion with a COPY of the pixel data into the
/// shared region vector. The main loop then emits escape sequences at
/// those recorded positions via a Post() task.
class ImagePlaceholderNode : public ftxui::Node {
public:
    /// @param cols      Requested width in terminal columns.
    /// @param rows      Requested height in terminal rows.
    /// @param plane     Image plane (pixels are copied during Render).
    /// @param sink      Shared vector where regions are collected.
    ImagePlaceholderNode(int cols, int rows,
                         const ImagePlaneView& plane,
                         std::vector<ImageRegion>* sink)
        : cols_(cols), rows_(rows), plane_(plane), sink_(sink) {}

    void ComputeRequirement() override {
        requirement_.min_x = cols_;
        requirement_.min_y = rows_;
    }

    void SetBox(ftxui::Box box) override {
        Node::SetBox(box);
    }

    void Render(ftxui::Screen& screen) override {
        int box_w = box_.x_max - box_.x_min + 1;
        int box_h = box_.y_max - box_.y_min + 1;

        // Render half-block characters as a Canvas-like backdrop.
        // When a graphics protocol is active, the inline image will be
        // overlaid on top via escape sequences. The half-block rendering
        // eliminates the "flash" problem: when image data changes, FTXUI's
        // differential renderer writes new half-blocks (which already look
        // like the image), so the transition from old-inline-image →
        // new-half-blocks → new-inline-image is nearly invisible.
        // Without this, the transition was old-image → spaces → new-image,
        // creating an obvious flash.
        if (!plane_.pixels.empty() && plane_.nx > 0 && plane_.ny > 0) {
            double scale_x = static_cast<double>(plane_.nx) / box_w;
            double scale_y = static_cast<double>(plane_.ny) / (box_h * 2);
            int max_px = static_cast<int>(plane_.nx - 1);
            int max_py = static_cast<int>(plane_.ny - 1);

            for (int ty = 0; ty < box_h; ++ty) {
                for (int tx = 0; tx < box_w; ++tx) {
                    int sx = box_.x_min + tx;
                    int sy = box_.y_min + ty;
                    if (sx < 0 || sx >= screen.dimx() ||
                        sy < 0 || sy >= screen.dimy())
                        continue;

                    // Map terminal cell to image pixels (2 pixel rows per cell)
                    int px = std::min(static_cast<int>(tx * scale_x), max_px);
                    int py_top = std::min(static_cast<int>((ty * 2) * scale_y), max_py);
                    int py_bot = std::min(static_cast<int>((ty * 2 + 1) * scale_y), max_py);

                    float v_top = std::clamp(
                        plane_.pixels[static_cast<size_t>(py_top) * plane_.nx + px],
                        0.0f, 1.0f);
                    float v_bot = std::clamp(
                        plane_.pixels[static_cast<size_t>(py_bot) * plane_.nx + px],
                        0.0f, 1.0f);

                    uint8_t g_top = static_cast<uint8_t>(v_top * 255.0f);
                    uint8_t g_bot = static_cast<uint8_t>(v_bot * 255.0f);

                    auto& pixel = screen.PixelAt(sx, sy);
                    pixel.character = "\u2580"; // upper half block
                    pixel.foreground_color = Color(g_top, g_top, g_top);
                    pixel.background_color = Color(g_bot, g_bot, g_bot);
                }
            }
        } else {
            // No image data — fill with spaces
            for (int y = box_.y_min; y <= box_.y_max; ++y) {
                for (int x = box_.x_min; x <= box_.x_max; ++x) {
                    if (x >= 0 && x < screen.dimx() && y >= 0 && y < screen.dimy()) {
                        screen.PixelAt(x, y).character = " ";
                    }
                }
            }
        }

        // Push a region with OWNED pixel data into the collection vector.
        // At this point box_ is valid (SetBox was called before Render).
        if (sink_ && !plane_.pixels.empty()) {
            ImageRegion r;
            r.term_row = box_.y_min + 1; // FTXUI box is 0-based, ANSI is 1-based
            r.term_col = box_.x_min + 1;
            r.cols = box_w;
            r.rows = box_h;
            r.pixels = plane_.pixels; // copy
            r.px_w = plane_.nx;
            r.px_h = plane_.ny;
            sink_->push_back(std::move(r));
        }
    }

private:
    int cols_, rows_;
    const ImagePlaneView& plane_;
    std::vector<ImageRegion>* sink_;
};

/// @brief Render a single image plane as a labeled FTXUI element.
///
/// When proto == None, uses DrawBlock (half-block characters) for 2:1 vertical
/// sub-pixel resolution. When a graphics protocol is available, creates a
/// placeholder node that reserves space and registers an ImageRegion (with
/// copied pixel data) into the region_sink during FTXUI's Render pass.
///
/// @param plane        The image plane data.
/// @param proto        The detected graphics protocol (None = Canvas fallback).
/// @param region_sink  Collects ImageRegion entries for escape sequence emission.
inline Element RenderSinglePlane(const ImagePlaneView& plane,
                                  GraphicsProto proto = GraphicsProto::None,
                                  std::vector<ImageRegion>* region_sink = nullptr) {
    if (plane.pixels.empty() || plane.nx == 0 || plane.ny == 0)
        return text("");

    if (proto != GraphicsProto::None && region_sink) {
        // Graphics protocol path: placeholder node that reserves screen space.
        // The actual image will be drawn via escape sequences after FTXUI flushes.
        int cols = static_cast<int>(plane.nx);
        int rows = static_cast<int>(plane.ny) / 2; // ~2 pixels per terminal row
        if (rows < 1) rows = 1;

        auto node = std::make_shared<ImagePlaceholderNode>(cols, rows, plane, region_sink);

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
                                  std::vector<ImageRegion>* region_sink = nullptr) {
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
                RenderSinglePlane(img.planes[0], proto, region_sink) | flex,
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
        ? RenderSinglePlane(img.planes[0], proto, region_sink)
        : text("");
    Element top_right = (img.planes.size() > 1)
        ? RenderSinglePlane(img.planes[1], proto, region_sink)
        : text("");
    Element bot_left  = (img.planes.size() > 2)
        ? RenderSinglePlane(img.planes[2], proto, region_sink)
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

// ---------------------------------------------------------------------------
// Post-render callback wrapper
// ---------------------------------------------------------------------------

/// @brief FTXUI Node wrapper that calls a callback after its child renders.
///
/// Used to schedule image emission after FTXUI finishes its Render pass
/// (when all ImagePlaceholderNode positions are finalized) but before
/// FTXUI flushes the screen buffer. A Post() call from the callback
/// ensures the emit task runs on the NEXT iteration — after the flush.
class PostRenderNode : public ftxui::Node {
public:
    PostRenderNode(Element child, std::function<void()> callback)
        : callback_(std::move(callback)) {
        children_.push_back(std::move(child));
    }

    void ComputeRequirement() override {
        Node::ComputeRequirement();
        requirement_ = children_[0]->requirement();
    }

    void SetBox(ftxui::Box box) override {
        Node::SetBox(box);
        children_[0]->SetBox(box);
    }

    void Render(ftxui::Screen& screen) override {
        children_[0]->Render(screen);
        if (callback_) callback_();
    }

private:
    std::function<void()> callback_;
};

/// @brief Wrap an Element with a post-render callback.
///
/// The callback fires after the wrapped element (and all its children)
/// have completed their Render() calls.
inline Element with_post_render(Element inner, std::function<void()> callback) {
    return std::make_shared<PostRenderNode>(std::move(inner), std::move(callback));
}
