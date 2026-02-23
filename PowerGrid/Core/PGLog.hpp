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

/// @brief Emit a "start" lifecycle event and auto-spawn pgview in TUI mode.
///
/// When TUI mode is active, this function:
/// 1. Locates the pgview binary (next to this executable, or on PATH)
/// 2. Forks pgview as a child process
/// 3. Redirects this process's stderr to pgview's stdin via a pipe
/// 4. pgview renders its TUI on stdout (the terminal)
///
/// If pgview is not found or stdout isn't a terminal, falls back to raw
/// JSONL output on stderr (still parseable by an external pgview).
inline void PG_TUI_START(const std::string& app_name,
                         const std::string& version = "1.1.0") {
    if (!PG_TUI_MODE()) return;

    // Try to auto-spawn pgview
    std::string pgview_path = pg_detail::find_pgview();
    auto& tui = pg_detail::TuiProcess::instance();
    if (!pgview_path.empty()) {
        tui.spawn(pgview_path);
        // If spawn fails, we continue — JSONL still goes to stderr
    }

    // Emit the start message (goes to pipe if pgview spawned, else raw stderr)
    auto msg = fmt::format(
        R"({{"type":"start","app":"{}","version":"{}"}})",
        pg_detail::json_escape(app_name),
        pg_detail::json_escape(version));
    msg += '\n';
    std::fwrite(msg.data(), 1, msg.size(), stderr);
    std::fflush(stderr);
}

/// @brief Emit an "exit" lifecycle event and wait for pgview to finish.
///
/// After emitting the exit JSONL message, this closes the pipe to pgview
/// and waits for the user to press 'q' in the TUI before returning.
/// If pgview was not spawned, this is a no-op beyond the JSONL emission.
inline void PG_TUI_EXIT(int code) {
    if (!PG_TUI_MODE()) return;

    // Emit the exit message
    auto msg = fmt::format(R"({{"type":"exit","code":{}}})", code);
    msg += '\n';
    std::fwrite(msg.data(), 1, msg.size(), stderr);
    std::fflush(stderr);

    // Wait for pgview to finish (user presses 'q')
    auto& tui = pg_detail::TuiProcess::instance();
    tui.wait();
}

// ---------------------------------------------------------------------------
// Progress bars — dual mode: JSONL (TUI) or indicators (classic)
// ---------------------------------------------------------------------------

using PGProgressBar = indicators::BlockProgressBar;

/// @brief Singleton manager for stacked progress bar display.
///
/// In TUI mode: emits JSONL progress messages, does not use indicators display.
/// In classic mode: delegates to indicators::DynamicProgress for ANSI rendering.
struct PGProgressManager {
    /// Owns the bars so they outlive the DynamicProgress reference wrappers.
    std::vector<std::shared_ptr<PGProgressBar>> owned_bars;
    /// Stacked display container — references bars in owned_bars.
    /// Only used in classic (non-TUI) mode.
    indicators::DynamicProgress<PGProgressBar> display;
    /// Mutex for thread-safe JSONL output.
    std::mutex jsonl_mutex;
    /// Next progress bar ID for JSONL messages.
    size_t next_id = 0;

    PGProgressManager() {
        display.set_option(indicators::option::HideBarWhenComplete{true});
    }

    /// Register a new bar. In TUI mode, emits JSONL "add" and returns an index.
    /// In classic mode, adds to indicators::DynamicProgress.
    size_t add(std::shared_ptr<PGProgressBar> bar,
               const std::string& label,
               size_t max_progress) {
        size_t id = next_id++;
        owned_bars.push_back(bar);

        if (PG_TUI_MODE()) {
            // Emit JSONL add message
            std::lock_guard<std::mutex> lock(jsonl_mutex);
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
        }
        return id;
    }

    /// Tick a bar. In TUI mode, emits JSONL "tick".
    void tick(size_t idx, size_t current, const std::string& postfix = "") {
        if (PG_TUI_MODE()) {
            std::lock_guard<std::mutex> lock(jsonl_mutex);
            auto msg = fmt::format(
                R"({{"type":"progress","action":"tick","id":{},"current":{},"postfix":"{}"}})",
                idx, current,
                pg_detail::json_escape(postfix));
            msg += '\n';
            std::fwrite(msg.data(), 1, msg.size(), stderr);
            std::fflush(stderr);
        } else {
            auto& bar = display[idx];
            if (!postfix.empty())
                bar.set_option(indicators::option::PostfixText{postfix});
            bar.tick();
        }
    }

    /// Mark a bar as completed.
    void done(size_t idx) {
        if (PG_TUI_MODE()) {
            std::lock_guard<std::mutex> lock(jsonl_mutex);
            auto msg = fmt::format(
                R"({{"type":"progress","action":"done","id":{}}})", idx);
            msg += '\n';
            std::fwrite(msg.data(), 1, msg.size(), stderr);
            std::fflush(stderr);
        } else {
            display[idx].mark_as_completed();
        }
    }

    /// Emit metrics (TUI mode only).
    void metrics(size_t bar_id, size_t iter, double error_norm, double penalty) {
        if (!PG_TUI_MODE()) return;
        std::lock_guard<std::mutex> lock(jsonl_mutex);
        auto msg = fmt::format(
            R"({{"type":"metrics","bar_id":{},"iter":{},"error_norm":{:.6e},"penalty":{:.6e}}})",
            bar_id, iter, error_norm, penalty);
        msg += '\n';
        std::fwrite(msg.data(), 1, msg.size(), stderr);
        std::fflush(stderr);
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
