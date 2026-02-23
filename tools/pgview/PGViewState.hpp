#pragma once

#include <string>
#include <vector>
#include <deque>
#include <mutex>
#include <chrono>

/// @brief Log entry for display.
struct LogEntry {
    std::string timestamp;  // "HH:MM:SS" local time
    std::string level;      // "info", "debug", "warn", "error"
    std::string message;
};

/// @brief A single image plane for display in the TUI.
struct ImagePlaneView {
    std::string label;   // "Image", "Axial", "Coronal", "Sagittal"
    size_t nx = 0, ny = 0;
    std::vector<float> pixels; // Row-major, [0, 1] normalized magnitudes
};

/// @brief A complete image frame (one or more planes) for TUI preview.
struct ImageFrame {
    size_t bar_id = 0;
    size_t iter = 0;
    std::vector<ImagePlaneView> planes;  // 1 plane (2D) or 3 planes (3D MPR)
};

/// @brief Per-bar convergence metric sample.
struct MetricSample {
    double error_norm = 0.0;
    double penalty = 0.0;
};

/// @brief State for a single progress bar.
struct ProgressBarState {
    std::string label;
    size_t max = 0;
    size_t current = 0;
    std::string postfix;
    bool completed = false;
    std::chrono::steady_clock::time_point start_time;

    /// Ring buffer of convergence metrics (last N samples).
    std::deque<MetricSample> metrics;
    static constexpr size_t kMaxMetrics = 40;

    void push_metric(double err, double pen) {
        metrics.push_back({err, pen});
        if (metrics.size() > kMaxMetrics)
            metrics.pop_front();
    }

    /// Compute ETA in seconds, or -1 if not enough data.
    double eta_seconds() const {
        if (current == 0 || completed) return -1.0;
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - start_time).count();
        double remaining = elapsed * static_cast<double>(max - current)
                         / static_cast<double>(current);
        return remaining;
    }

    /// Progress fraction [0, 1].
    float fraction() const {
        if (max == 0) return 0.0f;
        return static_cast<float>(current) / static_cast<float>(max);
    }
};

/// @brief Thread-safe shared state between the reader thread and the UI thread.
struct PGViewState {
    mutable std::mutex mu;

    /// App info from the "start" lifecycle message.
    std::string app_name;
    std::string app_version;

    /// Progress bars indexed by their id.
    std::vector<ProgressBarState> bars;

    /// Log ring buffer.
    std::deque<LogEntry> logs;
    static constexpr size_t kMaxLogs = 2000;

    /// Set when the child process sends "exit" or stdin closes.
    bool finished = false;
    int exit_code = 0;

    /// Latest image preview for TUI display.
    ImageFrame latest_image;
    bool has_image = false;

    /// Session timing — set by set_start() and set_exit().
    std::chrono::steady_clock::time_point session_start_time;
    std::chrono::steady_clock::time_point session_end_time;
    bool has_session_start = false;
    bool has_session_end = false;

    // --- Mutated by the reader thread, read by the UI thread ---

    void set_start(const std::string& app, const std::string& ver) {
        std::lock_guard<std::mutex> lock(mu);
        app_name = app;
        app_version = ver;
        session_start_time = std::chrono::steady_clock::now();
        has_session_start = true;
    }

    void set_exit(int code) {
        std::lock_guard<std::mutex> lock(mu);
        finished = true;
        exit_code = code;
        session_end_time = std::chrono::steady_clock::now();
        has_session_end = true;
    }

    void add_log(const std::string& ts, const std::string& level,
                 const std::string& msg) {
        std::lock_guard<std::mutex> lock(mu);
        // Convert ISO timestamp to HH:MM:SS for display
        std::string display_ts = ts;
        if (ts.size() >= 19) {
            // "2026-02-23T14:30:01.123Z" → "14:30:01"
            display_ts = ts.substr(11, 8);
        }
        logs.push_back({display_ts, level, msg});
        if (logs.size() > kMaxLogs)
            logs.pop_front();
    }

    void add_raw_log(const std::string& line) {
        std::lock_guard<std::mutex> lock(mu);
        auto now = std::chrono::system_clock::now();
        auto tt = std::chrono::system_clock::to_time_t(now);
        struct tm local;
        localtime_r(&tt, &local);
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%02d:%02d:%02d",
                      local.tm_hour, local.tm_min, local.tm_sec);
        logs.push_back({buf, "raw", line});
        if (logs.size() > kMaxLogs)
            logs.pop_front();
    }

    void progress_add(size_t id, const std::string& label, size_t max) {
        std::lock_guard<std::mutex> lock(mu);
        if (id >= bars.size())
            bars.resize(id + 1);
        bars[id].label = label;
        bars[id].max = max;
        bars[id].current = 0;
        bars[id].completed = false;
        bars[id].start_time = std::chrono::steady_clock::now();
        bars[id].metrics.clear();
    }

    void progress_tick(size_t id, size_t current, const std::string& postfix) {
        std::lock_guard<std::mutex> lock(mu);
        if (id >= bars.size()) return;
        bars[id].current = current;
        if (!postfix.empty())
            bars[id].postfix = postfix;
    }

    void progress_done(size_t id) {
        std::lock_guard<std::mutex> lock(mu);
        if (id >= bars.size()) return;
        bars[id].completed = true;
    }

    void add_metrics(size_t bar_id, double error_norm, double penalty) {
        std::lock_guard<std::mutex> lock(mu);
        if (bar_id >= bars.size()) return;
        bars[bar_id].push_metric(error_norm, penalty);
    }

    void set_image(size_t bar_id, size_t iter, std::vector<ImagePlaneView> planes) {
        std::lock_guard<std::mutex> lock(mu);
        latest_image = {bar_id, iter, std::move(planes)};
        has_image = true;
    }

    // --- Summary helpers (call with mu locked) ---

    /// Total elapsed seconds from start to exit (or now if no exit yet).
    double total_elapsed_seconds() const {
        if (!has_session_start) return 0.0;
        auto end = has_session_end ? session_end_time
                                   : std::chrono::steady_clock::now();
        return std::chrono::duration<double>(end - session_start_time).count();
    }

    /// Count completed PCG solves (each = one reconstructed image).
    size_t num_completed_images() const {
        size_t count = 0;
        for (const auto& bar : bars) {
            if (bar.completed && bar.label == "PCG")
                ++count;
        }
        return count;
    }

    /// Total PCG iterations across all completed PCG bars.
    size_t total_pcg_iterations() const {
        size_t total = 0;
        for (const auto& bar : bars) {
            if (bar.label == "PCG")
                total += bar.current;
        }
        return total;
    }

    /// Get the final error norm from the last completed PCG bar, or -1 if none.
    double last_error_norm() const {
        for (auto it = bars.rbegin(); it != bars.rend(); ++it) {
            if (it->label == "PCG" && it->completed && !it->metrics.empty())
                return it->metrics.back().error_norm;
        }
        return -1.0;
    }
};
