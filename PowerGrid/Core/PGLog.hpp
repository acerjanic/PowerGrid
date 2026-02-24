/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file PGLog.hpp
/// @brief Structured logging facade, JSONL output for TUI, and progress bar display.
///
/// Default mode: TUI/JSONL — log messages and progress updates are emitted as
/// single-line JSON objects on stderr, suitable for piping to `pgview`.
/// Pass `no_tui=true` to PG_LOG_INIT to fall back to classic spdlog+indicators.

#pragma once

#include "Version.h"

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/fmt/fmt.h>
#include <indicators/block_progress_bar.hpp>
#include <indicators/dynamic_progress.hpp>
#include <indicators/cursor_control.hpp>
#include <indicators/color.hpp>

#include <memory>
#include <string>
#include <vector>
#include <mutex>
#include <cstdio>
#include <ctime>
#include <chrono>

// POSIX headers for auto-spawning pgview
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#include <fcntl.h>
#include <cstdlib>
#include <climits>

#ifdef __APPLE__
#include <mach-o/dyld.h> // _NSGetExecutablePath
#endif

// ---------------------------------------------------------------------------
// TUI mode flag
// ---------------------------------------------------------------------------

/// @brief Query or set the TUI/JSONL output mode.
///
/// Default is true (TUI/JSONL enabled). Call with `false` to disable before
/// PG_LOG_INIT, or let PG_LOG_INIT set it via the `no_tui` parameter.
inline bool& PG_TUI_MODE(bool set_value = true, bool do_set = false) {
    static bool mode = true; // default: TUI/JSONL enabled
    if (do_set)
        mode = set_value;
    return mode;
}

// ---------------------------------------------------------------------------
// JSONL stderr sink for spdlog
// ---------------------------------------------------------------------------

namespace pg_detail {

/// @brief Escape a string for JSON output (handles ", \, newlines, tabs).
inline std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:   out += c;      break;
        }
    }
    return out;
}

/// @brief Get the spdlog level name as a lowercase string.
inline const char* level_name(spdlog::level::level_enum lvl) {
    switch (lvl) {
        case spdlog::level::trace:    return "trace";
        case spdlog::level::debug:    return "debug";
        case spdlog::level::info:     return "info";
        case spdlog::level::warn:     return "warn";
        case spdlog::level::err:      return "error";
        case spdlog::level::critical: return "critical";
        default:                      return "unknown";
    }
}

/// @brief Format current time as ISO 8601 with milliseconds.
inline std::string iso8601_now() {
    auto now = std::chrono::system_clock::now();
    auto tt = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;
    struct tm utc;
    gmtime_r(&tt, &utc);
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02d.%03dZ",
                  utc.tm_year + 1900, utc.tm_mon + 1, utc.tm_mday,
                  utc.tm_hour, utc.tm_min, utc.tm_sec,
                  static_cast<int>(ms.count()));
    return buf;
}

/// @brief spdlog sink that writes JSONL log messages to stderr.
///
/// Each log message becomes a single-line JSON object:
/// {"type":"log","ts":"...","level":"info","file":"...","line":N,"msg":"..."}
template <typename Mutex>
class jsonl_stderr_sink : public spdlog::sinks::base_sink<Mutex> {
protected:
    void sink_it_(const spdlog::details::log_msg& msg) override {
        // Extract source file basename
        std::string file = "unknown";
        int line = 0;
        if (!msg.source.empty()) {
            file = msg.source.filename;
            line = msg.source.line;
            // Extract basename only
            auto pos = file.find_last_of("/\\");
            if (pos != std::string::npos)
                file = file.substr(pos + 1);
        }

        std::string payload = fmt::format(
            R"({{"type":"log","ts":"{}","level":"{}","file":"{}","line":{},"msg":"{}"}})",
            iso8601_now(),
            level_name(msg.level),
            json_escape(file),
            line,
            json_escape(std::string(msg.payload.data(), msg.payload.size()))
        );
        payload += '\n';
        std::fwrite(payload.data(), 1, payload.size(), stderr);
        std::fflush(stderr);
    }

    void flush_() override {
        std::fflush(stderr);
    }
};

using jsonl_stderr_sink_mt = jsonl_stderr_sink<std::mutex>;

// ---------------------------------------------------------------------------
// Auto-spawn pgview helper
// ---------------------------------------------------------------------------

/// @brief Get the directory containing the currently running executable.
inline std::string get_exe_dir() {
#ifdef __APPLE__
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);
    std::vector<char> buf(size);
    if (_NSGetExecutablePath(buf.data(), &size) == 0) {
        // Resolve symlinks
        char resolved[PATH_MAX];
        if (realpath(buf.data(), resolved)) {
            std::string path(resolved);
            auto pos = path.find_last_of('/');
            if (pos != std::string::npos)
                return path.substr(0, pos);
        }
    }
#else
    // Linux: /proc/self/exe
    char resolved[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", resolved, sizeof(resolved) - 1);
    if (len > 0) {
        resolved[len] = '\0';
        std::string path(resolved);
        auto pos = path.find_last_of('/');
        if (pos != std::string::npos)
            return path.substr(0, pos);
    }
#endif
    return "";
}

/// @brief Locate the pgview binary. Search order:
/// 1. Same directory as the current executable (installed layout: both in bin/)
/// 2. Relative to exe dir: ../tools/pgview/pgview (CMake build tree layout)
/// 3. PGVIEW_PATH environment variable (explicit override)
/// 4. PATH search
inline std::string find_pgview() {
    std::string exe_dir = get_exe_dir();

    if (!exe_dir.empty()) {
        // 1. Same directory (installed layout)
        std::string candidate = exe_dir + "/pgview";
        if (::access(candidate.c_str(), X_OK) == 0)
            return candidate;

        // 2. Build tree layout: exe is in build/apps/, pgview in build/tools/pgview/
        candidate = exe_dir + "/../tools/pgview/pgview";
        char resolved[PATH_MAX];
        if (realpath(candidate.c_str(), resolved) && ::access(resolved, X_OK) == 0)
            return resolved;
    }

    // 3. PGVIEW_PATH environment variable
    const char* pgview_env = std::getenv("PGVIEW_PATH");
    if (pgview_env && ::access(pgview_env, X_OK) == 0)
        return pgview_env;

    // 4. Search PATH
    const char* path_env = std::getenv("PATH");
    if (path_env) {
        std::string path_str(path_env);
        size_t start = 0;
        while (start < path_str.size()) {
            size_t end = path_str.find(':', start);
            if (end == std::string::npos) end = path_str.size();
            std::string dir = path_str.substr(start, end - start);
            std::string candidate = dir + "/pgview";
            if (::access(candidate.c_str(), X_OK) == 0)
                return candidate;
            start = end + 1;
        }
    }

    return ""; // Not found
}

/// @brief Manages the pgview child process lifecycle.
///
/// Spawns pgview via fork/exec, redirects parent's stderr to a pipe
/// feeding pgview's stdin. pgview renders TUI on stdout (the terminal).
struct TuiProcess {
    pid_t child_pid = -1;
    int pipe_write_fd = -1;
    int original_stderr_fd = -1;
    bool active = false;

    /// @brief Spawn pgview and redirect stderr to it.
    /// @returns true on success, false if pgview not found or fork fails.
    bool spawn(const std::string& pgview_path) {
        if (pgview_path.empty()) return false;

        // Check if stdout is a terminal (pgview needs a tty for its TUI)
        if (!isatty(STDOUT_FILENO)) return false;

        // Create a pipe: parent writes to pipe_fds[1], child reads from pipe_fds[0]
        int pipe_fds[2];
        if (pipe(pipe_fds) != 0) return false;

        // Save the original stderr fd before redirecting
        original_stderr_fd = dup(STDERR_FILENO);
        if (original_stderr_fd < 0) {
            close(pipe_fds[0]);
            close(pipe_fds[1]);
            return false;
        }

        // Ignore SIGPIPE so recon continues if pgview exits early
        signal(SIGPIPE, SIG_IGN);

        child_pid = fork();
        if (child_pid < 0) {
            // Fork failed — restore and clean up
            close(pipe_fds[0]);
            close(pipe_fds[1]);
            close(original_stderr_fd);
            original_stderr_fd = -1;
            return false;
        }

        if (child_pid == 0) {
            // === Child process (pgview) ===
            // stdin = pipe read end
            dup2(pipe_fds[0], STDIN_FILENO);
            close(pipe_fds[0]);
            close(pipe_fds[1]);

            // Close the saved original stderr in the child
            if (original_stderr_fd >= 0)
                close(original_stderr_fd);

            // stdout stays connected to the terminal (for TUI rendering)
            // stderr goes to /dev/null to avoid loops
            int devnull = open("/dev/null", O_WRONLY);
            if (devnull >= 0) {
                dup2(devnull, STDERR_FILENO);
                close(devnull);
            }

            execl(pgview_path.c_str(), "pgview", nullptr);
            // If exec fails, exit quietly
            _exit(127);
        }

        // === Parent process (reconstruction) ===
        close(pipe_fds[0]); // Close read end in parent

        // Redirect stderr to pipe write end
        dup2(pipe_fds[1], STDERR_FILENO);
        close(pipe_fds[1]); // Close the extra fd, STDERR_FILENO now points to the pipe

        pipe_write_fd = STDERR_FILENO; // stderr IS the pipe now
        active = true;
        return true;
    }

    /// @brief Non-blocking check whether pgview is still running.
    ///
    /// If pgview has exited, restores the original stderr fd and marks
    /// this process as inactive. Returns true if pgview is still alive.
    bool check_alive() {
        if (!active || child_pid <= 0) return false;

        int status = 0;
        pid_t result = waitpid(child_pid, &status, WNOHANG);
        if (result == child_pid) {
            // pgview has exited — restore original stderr
            clearerr(stderr);
            if (original_stderr_fd >= 0) {
                dup2(original_stderr_fd, STDERR_FILENO);
                close(original_stderr_fd);
                original_stderr_fd = -1;
            }
            child_pid = -1;
            active = false;
            return false;
        }
        return true; // Still running
    }

    /// @brief Close the pipe and wait for pgview to exit.
    void wait() {
        if (!active) return;

        // Flush stderr (the pipe)
        std::fflush(stderr);

        // Restore the original stderr
        if (original_stderr_fd >= 0) {
            dup2(original_stderr_fd, STDERR_FILENO);
            close(original_stderr_fd);
            original_stderr_fd = -1;
        }

        // Wait for pgview to finish (user presses 'q')
        if (child_pid > 0) {
            int status = 0;
            waitpid(child_pid, &status, 0);
            child_pid = -1;
        }

        active = false;
    }

    /// @brief Access the singleton TuiProcess.
    static TuiProcess& instance() {
        static TuiProcess proc;
        return proc;
    }
};

} // namespace pg_detail

// ---------------------------------------------------------------------------
// Logging initialisation
// ---------------------------------------------------------------------------

/// @brief Initialise the PowerGrid logger.
///
/// @param level_str  Log level string ("debug", "info", "warn", "error").
/// @param log_path   Path to the rotating log file (always human-readable).
/// @param no_tui     If true, disable JSONL and use classic spdlog+indicators.
inline void PG_LOG_INIT(const std::string& level_str = "info",
                        const std::string& log_path  = "powergrid.log",
                        bool no_tui = false) {
    // Set the TUI mode flag
    PG_TUI_MODE(!no_tui, true);

    auto level = spdlog::level::from_str(level_str);

    // File sink — always human-readable
    auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
        log_path, 10UL * 1024UL * 1024UL, 3);
    file_sink->set_level(spdlog::level::debug);

    std::vector<spdlog::sink_ptr> sinks;

    if (PG_TUI_MODE()) {
        // TUI mode: JSONL on stderr
        auto jsonl_sink = std::make_shared<pg_detail::jsonl_stderr_sink_mt>();
        jsonl_sink->set_level(level);
        sinks = {jsonl_sink, file_sink};
    } else {
        // Classic mode: colored stderr
        auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
        stderr_sink->set_level(level);
        sinks = {stderr_sink, file_sink};
    }

    auto logger = std::make_shared<spdlog::logger>("powergrid", sinks.begin(), sinks.end());
    logger->set_level(spdlog::level::debug);
    spdlog::set_default_logger(logger);
}

#define PG_DEBUG(...) SPDLOG_DEBUG(__VA_ARGS__)
#define PG_INFO(...)  SPDLOG_INFO(__VA_ARGS__)
#define PG_WARN(...)  SPDLOG_WARN(__VA_ARGS__)
#define PG_ERROR(...) SPDLOG_ERROR(__VA_ARGS__)

// ---------------------------------------------------------------------------
// Lifecycle helpers (TUI mode only)
// ---------------------------------------------------------------------------

/// @brief Reinitialise the spdlog logger for classic (non-TUI) mode.
///
/// Called when pgview fails to spawn, to swap the JSONL sink for a colored
/// stderr sink so the user sees human-readable output instead of raw JSON.
inline void PG_LOG_REINIT_CLASSIC(const std::string& log_path = "powergrid.log") {
    PG_TUI_MODE(false, true);

    auto level = spdlog::get("powergrid")
                     ? spdlog::get("powergrid")->level()
                     : spdlog::level::info;

    auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
        log_path, 10UL * 1024UL * 1024UL, 3);
    file_sink->set_level(spdlog::level::debug);

    auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
    stderr_sink->set_level(level);

    std::vector<spdlog::sink_ptr> sinks = {stderr_sink, file_sink};
    auto logger = std::make_shared<spdlog::logger>("powergrid", sinks.begin(), sinks.end());
    logger->set_level(spdlog::level::debug);
    spdlog::set_default_logger(logger);
}

/// @brief Emit a "start" lifecycle event and auto-spawn pgview in TUI mode.
///
/// When TUI mode is active, this function:
/// 1. Locates the pgview binary (next to this executable, or on PATH)
/// 2. Forks pgview as a child process
/// 3. Redirects this process's stderr to pgview's stdin via a pipe
/// 4. pgview renders its TUI on stdout (the terminal)
///
/// If pgview is not found or stdout isn't a terminal, automatically falls
/// back to classic spdlog+indicators mode (human-readable colored output).
inline void PG_TUI_START(const std::string& app_name,
                         const std::string& version = POWERGRID_VERSION_STRING) {
    if (!PG_TUI_MODE()) return;

    // Try to auto-spawn pgview
    std::string pgview_path = pg_detail::find_pgview();
    auto& tui = pg_detail::TuiProcess::instance();
    bool spawned = false;
    if (!pgview_path.empty()) {
        spawned = tui.spawn(pgview_path);
    }

    if (spawned) {
        // pgview is running — emit JSONL start message through the pipe
        auto msg = fmt::format(
            R"({{"type":"start","app":"{}","version":"{}"}})",
            pg_detail::json_escape(app_name),
            pg_detail::json_escape(version));
        msg += '\n';
        std::fwrite(msg.data(), 1, msg.size(), stderr);
        std::fflush(stderr);
    } else {
        // pgview not available — fall back to classic spdlog+indicators
        PG_LOG_REINIT_CLASSIC();
        PG_INFO("Starting {} v{}", app_name, version);
    }
}

/// @brief Emit an "exit" lifecycle event and wait for pgview to finish.
///
/// If pgview is still running, emits the exit JSONL message, closes the
/// pipe, and waits for the user to press 'q' in the TUI before returning.
/// If pgview already died (fell back to classic mode), just logs the exit.
inline void PG_TUI_EXIT(int code) {
    auto& tui = pg_detail::TuiProcess::instance();

    if (PG_TUI_MODE() && tui.active) {
        // pgview is running — emit JSONL exit and wait
        auto msg = fmt::format(R"({{"type":"exit","code":{}}})", code);
        msg += '\n';
        std::fwrite(msg.data(), 1, msg.size(), stderr);
        std::fflush(stderr);
        tui.wait();
    } else {
        // Classic mode or pgview already died
        PG_INFO("Reconstruction finished (exit code {})", code);
    }
}

// ---------------------------------------------------------------------------
// Progress bars — dual mode: JSONL (TUI) or indicators (classic)
// ---------------------------------------------------------------------------

using PGProgressBar = indicators::BlockProgressBar;

/// @brief Singleton manager for stacked progress bar display.
///
/// In TUI mode: emits JSONL progress messages, does not use indicators display.
/// In classic mode: delegates to indicators::DynamicProgress for ANSI rendering.
///
/// If pgview dies mid-run, automatically switches to classic mode and
/// re-registers all active bars with the indicators display.
struct PGProgressManager {
    /// Owns the bars so they outlive the DynamicProgress reference wrappers.
    std::vector<std::shared_ptr<PGProgressBar>> owned_bars;
    /// Whether each bar has been registered with indicators::DynamicProgress.
    std::vector<bool> registered_with_display;
    /// Stacked display container — references bars in owned_bars.
    /// Only used in classic (non-TUI) mode.
    indicators::DynamicProgress<PGProgressBar> display;
    /// Mutex for thread-safe output.
    std::mutex mu;
    /// Next progress bar ID for JSONL messages.
    size_t next_id = 0;

    PGProgressManager() {
        display.set_option(indicators::option::HideBarWhenComplete{true});
    }

    /// Check if pgview is still alive; if not, switch to classic mode.
    /// Must be called with mu held.
    void ensure_tui_or_fallback() {
        if (!PG_TUI_MODE()) return; // Already in classic mode

        auto& tui = pg_detail::TuiProcess::instance();
        if (!tui.active) return; // pgview was never spawned

        if (!tui.check_alive()) {
            // pgview has died — stderr is restored by check_alive()
            // Switch spdlog to classic mode
            PG_LOG_REINIT_CLASSIC();
            PG_INFO("pgview exited, switching to classic output");

            // Re-register all active (non-completed) bars with indicators
            for (size_t i = 0; i < owned_bars.size(); ++i) {
                if (!registered_with_display[i]) {
                    display.push_back(*owned_bars[i]);
                    registered_with_display[i] = true;
                }
            }
        }
    }

    /// Register a new bar. In TUI mode, emits JSONL "add" and returns an index.
    /// In classic mode, adds to indicators::DynamicProgress.
    size_t add(std::shared_ptr<PGProgressBar> bar,
               const std::string& label,
               size_t max_progress) {
        std::lock_guard<std::mutex> lock(mu);
        ensure_tui_or_fallback();

        size_t id = next_id++;
        owned_bars.push_back(bar);

        if (PG_TUI_MODE()) {
            registered_with_display.push_back(false);
            // Emit JSONL add message
            auto msg = fmt::format(
                R"({{"type":"progress","action":"add","id":{},"label":"{}","max":{}}})",
                id,
                pg_detail::json_escape(label),
                max_progress);
            msg += '\n';
            std::fwrite(msg.data(), 1, msg.size(), stderr);
            std::fflush(stderr);
        } else {
            // Classic mode: register with indicators display
            display.push_back(*bar);
            registered_with_display.push_back(true);
        }
        return id;
    }

    /// Tick a bar. In TUI mode, emits JSONL "tick".
    void tick(size_t idx, size_t current, const std::string& postfix = "") {
        std::lock_guard<std::mutex> lock(mu);
        ensure_tui_or_fallback();

        if (PG_TUI_MODE()) {
            auto msg = fmt::format(
                R"({{"type":"progress","action":"tick","id":{},"current":{},"postfix":"{}"}})",
                idx, current,
                pg_detail::json_escape(postfix));
            msg += '\n';
            std::fwrite(msg.data(), 1, msg.size(), stderr);
            std::fflush(stderr);
        } else {
            if (idx < owned_bars.size()) {
                auto& bar = display[idx];
                if (!postfix.empty())
                    bar.set_option(indicators::option::PostfixText{postfix});
                bar.tick();
            }
        }
    }

    /// Mark a bar as completed.
    void done(size_t idx) {
        std::lock_guard<std::mutex> lock(mu);
        ensure_tui_or_fallback();

        if (PG_TUI_MODE()) {
            auto msg = fmt::format(
                R"({{"type":"progress","action":"done","id":{}}})", idx);
            msg += '\n';
            std::fwrite(msg.data(), 1, msg.size(), stderr);
            std::fflush(stderr);
        } else {
            if (idx < owned_bars.size())
                display[idx].mark_as_completed();
        }
    }

    /// Emit metrics for convergence monitoring.
    /// In TUI mode: emits JSONL. In classic mode: logs a readable summary.
    void metrics(size_t bar_id, size_t iter, double error_norm, double penalty) {
        std::lock_guard<std::mutex> lock(mu);
        ensure_tui_or_fallback();

        if (PG_TUI_MODE()) {
            auto msg = fmt::format(
                R"({{"type":"metrics","bar_id":{},"iter":{},"error_norm":{:.6e},"penalty":{:.6e}}})",
                bar_id, iter, error_norm, penalty);
            msg += '\n';
            std::fwrite(msg.data(), 1, msg.size(), stderr);
            std::fflush(stderr);
        } else {
            // Classic mode: log metrics as readable text
            SPDLOG_INFO("iter {:3d}  err={:.4e}  penalty={:.4e}", iter, error_norm, penalty);
        }
    }
};

/// Access the singleton progress manager.
inline PGProgressManager& PG_PROGRESS() {
    static PGProgressManager mgr;
    return mgr;
}

/// @brief Register a new progress bar with label and max count.
///
/// In TUI mode, emits a JSONL "add" message. In classic mode, adds to indicators.
/// Returns an index for use with PG_PROGRESS_TICK / PG_PROGRESS_DONE.
inline size_t PG_PROGRESS_ADD(std::shared_ptr<PGProgressBar> bar,
                              const std::string& label = "",
                              size_t max_progress = 0) {
    return PG_PROGRESS().add(bar, label, max_progress);
}

/// @brief Advance bar at @p idx, reporting current absolute count.
///
/// @param idx      Bar index from PG_PROGRESS_ADD.
/// @param current  Absolute progress count (1-based, for ETA computation).
/// @param postfix  Optional postfix text.
inline void PG_PROGRESS_TICK(size_t idx, size_t current,
                             const std::string& postfix = "") {
    PG_PROGRESS().tick(idx, current, postfix);
}

/// @brief Mark bar at @p idx as completed (hidden in classic mode, "done" in TUI).
inline void PG_PROGRESS_DONE(size_t idx) {
    PG_PROGRESS().done(idx);
}

/// @brief Emit convergence metrics for a progress bar (TUI mode only).
///
/// @param bar_id      Bar index from PG_PROGRESS_ADD.
/// @param iter        Current iteration number (1-based).
/// @param error_norm  Data-fidelity error norm (||yi - Ax||).
/// @param penalty     Roughness penalty value (R.Penalty(x)).
inline void PG_METRICS(size_t bar_id, size_t iter,
                       double error_norm, double penalty) {
    PG_PROGRESS().metrics(bar_id, iter, error_norm, penalty);
}

// ---------------------------------------------------------------------------
// Image preview for TUI — live reconstruction magnitude display
// ---------------------------------------------------------------------------

namespace pg_detail {

/// @brief Base64 encode a byte buffer (no line wrapping).
inline std::string base64_encode(const uint8_t* data, size_t len) {
    static const char table[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    out.reserve(((len + 2) / 3) * 4);
    for (size_t i = 0; i < len; i += 3) {
        uint32_t n = static_cast<uint32_t>(data[i]) << 16;
        if (i + 1 < len) n |= static_cast<uint32_t>(data[i + 1]) << 8;
        if (i + 2 < len) n |= static_cast<uint32_t>(data[i + 2]);
        out += table[(n >> 18) & 0x3F];
        out += table[(n >> 12) & 0x3F];
        out += (i + 1 < len) ? table[(n >> 6) & 0x3F] : '=';
        out += (i + 2 < len) ? table[n & 0x3F] : '=';
    }
    return out;
}

/// @brief A single image plane for preview (sender side).
struct ImagePlane {
    std::string label;  // "Image", "Axial", "Coronal", "Sagittal"
    size_t nx, ny;
    std::vector<float> pixels;  // Row-major, [0, 1] normalized
};

/// @brief Box-filter downsample (or nearest-neighbor upscale) a 2D float array.
///
/// Supports non-square targets to preserve aspect ratio.
///
/// @param src       Source pixels in row-major order.
/// @param src_w     Source width.
/// @param src_h     Source height.
/// @param target_w  Target width.
/// @param target_h  Target height.
/// @returns         Downsampled array of target_w × target_h floats.
inline std::vector<float> box_downsample(const std::vector<float>& src,
                                          size_t src_w, size_t src_h,
                                          size_t target_w, size_t target_h) {
    std::vector<float> dst(target_w * target_h, 0.0f);
    if (src_w == 0 || src_h == 0) return dst;

    // If source is smaller than target in both dims, nearest-neighbor upscale
    if (src_w <= target_w && src_h <= target_h) {
        for (size_t ty = 0; ty < target_h; ty++) {
            size_t sy = ty * src_h / target_h;
            for (size_t tx = 0; tx < target_w; tx++) {
                size_t sx = tx * src_w / target_w;
                dst[ty * target_w + tx] = src[sy * src_w + sx];
            }
        }
        return dst;
    }

    // Box filter: average all source pixels that map to each target pixel
    for (size_t ty = 0; ty < target_h; ty++) {
        size_t sy0 = ty * src_h / target_h;
        size_t sy1 = (ty + 1) * src_h / target_h;
        if (sy1 == sy0) sy1 = sy0 + 1;
        for (size_t tx = 0; tx < target_w; tx++) {
            size_t sx0 = tx * src_w / target_w;
            size_t sx1 = (tx + 1) * src_w / target_w;
            if (sx1 == sx0) sx1 = sx0 + 1;
            float sum = 0.0f;
            size_t count = 0;
            for (size_t sy = sy0; sy < sy1 && sy < src_h; sy++) {
                for (size_t sx = sx0; sx < sx1 && sx < src_w; sx++) {
                    sum += src[sy * src_w + sx];
                    count++;
                }
            }
            dst[ty * target_w + tx] = (count > 0) ? sum / static_cast<float>(count) : 0.0f;
        }
    }
    return dst;
}

/// @brief Compute aspect-ratio-preserving target dimensions.
///
/// Fits src_w × src_h into max_dim × max_dim while preserving aspect ratio.
/// The longer dimension gets max_dim; the shorter is scaled proportionally.
inline std::pair<size_t, size_t> fit_aspect(size_t src_w, size_t src_h,
                                             size_t max_dim) {
    if (src_w == 0 || src_h == 0) return {max_dim, max_dim};
    if (src_w >= src_h) {
        size_t tw = max_dim;
        size_t th = std::max(size_t(1), src_h * max_dim / src_w);
        return {tw, th};
    } else {
        size_t th = max_dim;
        size_t tw = std::max(size_t(1), src_w * max_dim / src_h);
        return {tw, th};
    }
}

/// @brief Extract preview planes from a complex image vector.
///
/// For 2D (Nz==1): returns 1 plane ("Image") — the full Nx×Ny magnitude.
/// For 3D (Nz>1):  returns 3 planes for MPR display:
///   - "Axial"    — z = Nz/2, showing Nx × Ny
///   - "Coronal"  — y = Ny/2, showing Nx × Nz
///   - "Sagittal" — x = Nx/2, showing Ny × Nz
///
/// Each plane is box-filter downsampled to fit within preview_max_dim while
/// preserving aspect ratio, then normalized to [0, 1] using a global max.
///
/// @tparam T1           Floating-point precision (float or double).
/// @param x_data        Pointer to complex image data (column-major, Nx×Ny×Nz).
/// @param n_elem        Number of complex elements.
/// @param Nx, Ny, Nz    Image dimensions.
/// @param preview_max_dim  Maximum preview dimension (default 128).
/// @returns             Vector of ImagePlane structs.
template <typename T1>
inline std::vector<ImagePlane> extract_preview_planes(
        const std::complex<T1>* x_data, size_t n_elem,
        size_t Nx, size_t Ny, size_t Nz,
        size_t preview_max_dim = 128) {

    std::vector<ImagePlane> planes;
    if (Nx == 0 || Ny == 0) return planes;

    // Armadillo uses column-major storage:
    // For a 3D volume stored as Nx*Ny*Nz column vector,
    // element (ix, iy, iz) = x_data[ix + iy*Nx + iz*Nx*Ny]

    auto mag = [](std::complex<T1> c) -> float {
        return static_cast<float>(std::abs(c));
    };

    if (Nz <= 1) {
        // 2D: single plane
        std::vector<float> slice(Nx * Ny);
        for (size_t iy = 0; iy < Ny; iy++)
            for (size_t ix = 0; ix < Nx; ix++)
                slice[iy * Nx + ix] = mag(x_data[ix + iy * Nx]);

        auto [tw, th] = fit_aspect(Nx, Ny, preview_max_dim);
        auto ds = box_downsample(slice, Nx, Ny, tw, th);
        planes.push_back({"Image", tw, th, std::move(ds)});
    } else {
        // 3D: extract three orthogonal slices

        // Axial — z = Nz/2, shows Nx × Ny
        size_t zc = Nz / 2;
        std::vector<float> axial(Nx * Ny);
        for (size_t iy = 0; iy < Ny; iy++)
            for (size_t ix = 0; ix < Nx; ix++)
                axial[iy * Nx + ix] = mag(x_data[ix + iy * Nx + zc * Nx * Ny]);

        // Coronal — y = Ny/2, shows Nx × Nz (rows=Nz, cols=Nx)
        size_t yc = Ny / 2;
        std::vector<float> coronal(Nx * Nz);
        for (size_t iz = 0; iz < Nz; iz++)
            for (size_t ix = 0; ix < Nx; ix++)
                coronal[iz * Nx + ix] = mag(x_data[ix + yc * Nx + iz * Nx * Ny]);

        // Sagittal — x = Nx/2, shows Ny × Nz (rows=Nz, cols=Ny)
        size_t xc = Nx / 2;
        std::vector<float> sagittal(Ny * Nz);
        for (size_t iz = 0; iz < Nz; iz++)
            for (size_t iy = 0; iy < Ny; iy++)
                sagittal[iz * Ny + iy] = mag(x_data[xc + iy * Nx + iz * Nx * Ny]);

        auto [tw_ax, th_ax] = fit_aspect(Nx, Ny, preview_max_dim);
        auto [tw_co, th_co] = fit_aspect(Nx, Nz, preview_max_dim);
        auto [tw_sa, th_sa] = fit_aspect(Ny, Nz, preview_max_dim);

        auto ds_ax = box_downsample(axial, Nx, Ny, tw_ax, th_ax);
        auto ds_co = box_downsample(coronal, Nx, Nz, tw_co, th_co);
        auto ds_sa = box_downsample(sagittal, Ny, Nz, tw_sa, th_sa);

        planes.push_back({"Axial", tw_ax, th_ax, std::move(ds_ax)});
        planes.push_back({"Coronal", tw_co, th_co, std::move(ds_co)});
        planes.push_back({"Sagittal", tw_sa, th_sa, std::move(ds_sa)});
    }

    // Normalize all planes to [0, 1] using global max for consistent windowing
    float global_max = 0.0f;
    for (const auto& p : planes)
        for (float v : p.pixels)
            if (v > global_max) global_max = v;

    if (global_max > 0.0f) {
        float inv_max = 1.0f / global_max;
        for (auto& p : planes)
            for (float& v : p.pixels)
                v *= inv_max;
    }

    return planes;
}

} // namespace pg_detail (continued)

// ---------------------------------------------------------------------------
// Image preview context — set by app, read by solver
// ---------------------------------------------------------------------------

/// @brief Singleton storing image dimensions for preview extraction.
///
/// The solver template doesn't know Nx/Ny/Nz directly. The app sets these
/// before calling reconSolve, and the solver reads them for preview extraction.
struct PGImagePreviewContext {
    size_t Nx = 0, Ny = 0, Nz = 1;
    static PGImagePreviewContext& instance() {
        static PGImagePreviewContext ctx;
        return ctx;
    }
};

/// @brief Set the image dimensions for preview extraction.
///
/// Call this in the app before entering the solver loop.
inline void PG_SET_IMAGE_DIMS(size_t nx, size_t ny, size_t nz = 1) {
    auto& ctx = PGImagePreviewContext::instance();
    ctx.Nx = nx;
    ctx.Ny = ny;
    ctx.Nz = nz;
}

/// @brief Emit an image preview for display in pgview.
///
/// In TUI mode: base64-encodes each plane and emits JSONL with a "planes" array.
/// In classic mode: no-op (images aren't displayable on a scrolling terminal).
///
/// @param bar_id  Bar index from PG_PROGRESS_ADD.
/// @param iter    Current iteration number (1-based).
/// @param planes  Vector of image planes to send.
inline void PG_IMAGE_PREVIEW(size_t bar_id, size_t iter,
                              const std::vector<pg_detail::ImagePlane>& planes) {
    if (!PG_TUI_MODE() || planes.empty()) return;

    auto& mgr = PG_PROGRESS();
    std::lock_guard<std::mutex> lock(mgr.mu);
    mgr.ensure_tui_or_fallback();
    if (!PG_TUI_MODE()) return;

    // Build JSONL with planes array
    std::string msg = fmt::format(
        R"({{"type":"image_preview","bar_id":{},"iter":{},"planes":[)", bar_id, iter);

    for (size_t i = 0; i < planes.size(); i++) {
        const auto& p = planes[i];
        std::string b64 = pg_detail::base64_encode(
            reinterpret_cast<const uint8_t*>(p.pixels.data()),
            p.pixels.size() * sizeof(float));
        if (i > 0) msg += ',';
        msg += fmt::format(
            R"({{"label":"{}","nx":{},"ny":{},"data":"{}"}})",
            pg_detail::json_escape(p.label), p.nx, p.ny, b64);
    }
    msg += "]}\n";

    std::fwrite(msg.data(), 1, msg.size(), stderr);
    std::fflush(stderr);
}
