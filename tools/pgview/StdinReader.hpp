#pragma once

#include <iostream>
#include <string>
#include <thread>
#include <functional>
#include <nlohmann/json.hpp>
#include "PGViewState.hpp"

using json = nlohmann::json;

/// @brief Background thread that reads JSONL from stdin and updates PGViewState.
///
/// Non-JSON lines are treated as raw log messages (graceful fallback).
/// Calls `on_update()` after each state mutation to trigger UI redraw.
class StdinReader {
public:
    StdinReader(PGViewState& state, std::function<void()> on_update)
        : state_(state), on_update_(std::move(on_update)) {}

    /// Start the reader thread. Non-blocking.
    void start() {
        thread_ = std::thread([this]() { run(); });
    }

    /// Wait for the reader thread to finish (stdin closed).
    void join() {
        if (thread_.joinable())
            thread_.join();
    }

private:
    void run() {
        std::string line;
        while (std::getline(std::cin, line)) {
            if (line.empty()) continue;

            // Try to parse as JSON
            if (line.front() == '{') {
                try {
                    auto j = json::parse(line);
                    process_json(j);
                    on_update_();
                    continue;
                } catch (...) {
                    // Fall through to raw log
                }
            }

            // Not valid JSON — treat as raw log line
            state_.add_raw_log(line);
            on_update_();
        }

        // stdin closed — mark finished if not already done via "exit" message
        {
            std::lock_guard<std::mutex> lock(state_.mu);
            if (!state_.finished) {
                state_.finished = true;
            }
        }
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
    std::thread thread_;
};
