/// @file main.cpp
/// @brief pgview — TUI viewer for PowerGrid reconstruction JSONL output.
///
/// Usage: PowerGridIsmrmrd [args] 2>&1 | pgview
///     or: auto-spawned by PG_TUI_START() with stdin as pipe
///
/// Reads JSONL from stdin (or a saved pipe fd), renders scrolling logs and
/// progress bars with ETA, sparklines for convergence metrics, and wall
/// clock completion time. Press 'q' or Escape to quit.

#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/component/event.hpp>
#include <ftxui/dom/elements.hpp>
#include <atomic>
#include <thread>
#include <unistd.h>
#include <fcntl.h>

#include "PGViewState.hpp"
#include "StdinReader.hpp"
#include "PGViewUI.hpp"

using namespace ftxui;

int main() {
    // When auto-spawned, stdin is a pipe carrying JSONL data.
    // FTXUI needs stdin to be the terminal for keyboard input.
    // Solution: save the pipe fd, then reopen /dev/tty as stdin.
    int data_fd = -1;
    if (!isatty(STDIN_FILENO)) {
        // stdin is a pipe — save it and replace with /dev/tty
        data_fd = dup(STDIN_FILENO);
        int tty_fd = open("/dev/tty", O_RDONLY);
        if (tty_fd >= 0) {
            dup2(tty_fd, STDIN_FILENO);
            close(tty_fd);
        }
    }
    // If stdin was already a terminal (manual pipe: ... | pgview),
    // data_fd stays -1 and StdinReader will use stdin directly.

    PGViewState state;
    auto screen = ScreenInteractive::Fullscreen();

    // Start the stdin reader thread, posting custom events on updates
    // Pass data_fd so it reads JSONL from the pipe (not from stdin/terminal)
    StdinReader reader(state, [&screen]() {
        screen.Post(Event::Custom);
    }, data_fd);
    reader.start();

    // Build the UI renderer
    auto renderer = Renderer([&]() {
        std::lock_guard<std::mutex> lock(state.mu);

        // Title bar
        std::string title = "pgview";
        if (!state.app_name.empty()) {
            title = state.app_name;
            if (!state.app_version.empty())
                title += " v" + state.app_version;
        }
        if (state.finished) {
            title += " (exit " + std::to_string(state.exit_code) + ")";
        }

        auto title_bar = hbox({
            text(" " + title + " ") | bold | inverted,
            filler(),
            text(" q=quit ") | dim,
        });

        // Progress section
        Element progress_section = RenderProgress(state);

        // Log section
        Element log_section = RenderLogs(state);

        // Compose layout
        Elements layout;
        layout.push_back(title_bar);
        layout.push_back(separator());

        // Only show progress section if there are active bars
        bool has_active_bars = false;
        for (const auto& bar : state.bars) {
            if (!bar.completed) {
                has_active_bars = true;
                break;
            }
        }
        if (has_active_bars) {
            layout.push_back(progress_section);
            layout.push_back(separator());
        }

        layout.push_back(log_section);

        return vbox(std::move(layout)) | border;
    });

    // Handle keyboard events
    auto component = CatchEvent(renderer, [&](Event event) {
        if (event == Event::Character('q') || event == Event::Escape) {
            screen.Exit();
            return true;
        }
        return false;
    });

    screen.Loop(component);

    // Wait for reader thread to finish
    reader.join();

    // Clean up saved pipe fd
    if (data_fd >= 0)
        close(data_fd);

    return 0;
}
