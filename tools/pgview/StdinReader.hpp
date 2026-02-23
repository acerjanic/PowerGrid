#pragma once

#include <iostream>
#include <string>
#include <thread>
#include <functional>
#include <cstdio>
#include <unistd.h>
#include <cstring>
#include <nlohmann/json.hpp>
#include "PGViewState.hpp"
#include "Base64.hpp"

using json = nlohmann::json;

/// @brief Background thread that reads JSONL and updates PGViewState.
///
/// When data_fd >= 0 (auto-spawn mode), reads from that file descriptor
/// (the pipe from the reconstruction process). When data_fd < 0 (manual
/// pipe mode), reads from std::cin (which is stdin).
///
/// Non-JSON lines are treated as raw log messages (graceful fallback).
/// Calls `on_update()` after each state mutation to trigger UI redraw.
class StdinReader {
public:
    /// @param data_fd  File descriptor to read JSONL from. If < 0, uses std::cin.
    StdinReader(PGViewState& state, std::function<void()> on_update, int data_fd = -1)
        : state_(state), on_update_(std::move(on_update)), data_fd_(data_fd) {}

    /// Start the reader thread. Non-blocking.
    void start() {
        thread_ = std::thread([this]() { run(); });
    }

    /// Wait for the reader thread to finish (input closed).
    void join() {
        if (thread_.joinable())
            thread_.join();
    }

private:
    void run() {
        if (data_fd_ >= 0) {
            // Read from the pipe fd (auto-spawn mode)
            run_from_fd();
        } else {
            // Read from std::cin (manual pipe mode)
            run_from_cin();
        }

        // Input closed — mark finished if not already done via "exit" message
        {
            std::lock_guard<std::mutex> lock(state_.mu);
            if (!state_.finished) {
                state_.finished = true;
            }
        }
        on_update_();
    }

    /// Read JSONL lines from std::cin (stdin).
    void run_from_cin() {
        std::string line;
        while (std::getline(std::cin, line)) {
            process_line(line);
        }
    }

    /// Read JSONL lines from a raw file descriptor.
    void run_from_fd() {
        // Wrap fd in a FILE* for buffered line reading
        FILE* fp = fdopen(data_fd_, "r");
        if (!fp) return;

        // Use POSIX getline() for dynamic buffering — image_preview
        // messages can be 20KB+ (base64-encoded image data).
        char* linebuf = nullptr;
        size_t linecap = 0;
        ssize_t len;
        while ((len = ::getline(&linebuf, &linecap, fp)) > 0) {
            std::string line(linebuf, static_cast<size_t>(len));
            // Strip trailing newline
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
                line.pop_back();
            process_line(line);
        }
        free(linebuf);
        // Don't fclose — main.cpp owns the fd and will close it
    }

    /// Process a single line of input (JSON or raw text).
    void process_line(const std::string& line) {
        if (line.empty()) return;

        // Try to parse as JSON
        if (line.front() == '{') {
            try {
                auto j = json::parse(line);
                process_json(j);
                on_update_();
                return;
            } catch (...) {
                // Fall through to raw log
            }
        }

        // Not valid JSON — treat as raw log line
        state_.add_raw_log(line);
        on_update_();
    }

    void process_json(const json& j) {
        auto type = j.value("type", "");

        if (type == "log") {
            state_.add_log(
                j.value("ts", ""),
                j.value("level", "info"),
                j.value("msg", "")
            );
        } else if (type == "progress") {
            auto action = j.value("action", "");
            size_t id = j.value("id", 0u);

            if (action == "add") {
                state_.progress_add(
                    id,
                    j.value("label", ""),
                    j.value("max", 0u)
                );
            } else if (action == "tick") {
                state_.progress_tick(
                    id,
                    j.value("current", 0u),
                    j.value("postfix", "")
                );
            } else if (action == "done") {
                state_.progress_done(id);
            }
        } else if (type == "metrics") {
            state_.add_metrics(
                j.value("bar_id", 0u),
                j.value("error_norm", 0.0),
                j.value("penalty", 0.0)
            );
        } else if (type == "image_preview") {
            size_t bar_id = j.value("bar_id", 0u);
            size_t iter = j.value("iter", 0u);
            std::vector<ImagePlaneView> planes;
            if (j.contains("planes") && j["planes"].is_array()) {
                for (const auto& p : j["planes"]) {
                    ImagePlaneView plane;
                    plane.label = p.value("label", "");
                    plane.nx = p.value("nx", 0u);
                    plane.ny = p.value("ny", 0u);
                    std::string b64 = p.value("data", "");
                    auto raw = base64_decode(b64);
                    size_t expected = plane.nx * plane.ny * sizeof(float);
                    if (raw.size() == expected && plane.nx > 0 && plane.ny > 0) {
                        plane.pixels.resize(plane.nx * plane.ny);
                        std::memcpy(plane.pixels.data(), raw.data(), expected);
                        planes.push_back(std::move(plane));
                    }
                }
            }
            if (!planes.empty())
                state_.set_image(bar_id, iter, std::move(planes));
        } else if (type == "start") {
            state_.set_start(
                j.value("app", "PowerGrid"),
                j.value("version", "")
            );
        } else if (type == "exit") {
            state_.set_exit(j.value("code", 0));
        }
    }

    PGViewState& state_;
    std::function<void()> on_update_;
    int data_fd_;
    std::thread thread_;
};
