/// @file main.cpp
/// @brief pgview — TUI viewer for PowerGrid reconstruction JSONL output.
///
/// Usage: PowerGridIsmrmrd [args] 2>&1 | pgview
///     or: auto-spawned by PG_TUI_START() with stdin as pipe
///
/// Reads JSONL from stdin (or a saved pipe fd), renders scrolling logs and
/// progress bars with ETA, sparklines for convergence metrics, and wall
/// clock completion time. Supports inline terminal images via iTerm2 and
/// Kitty protocols, with Canvas half-block fallback for other terminals.
/// Auto-exits when the reconstruction finishes and prints a summary.

#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/component/event.hpp>
#include <ftxui/dom/elements.hpp>
#include <atomic>
#include <thread>
#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <fcntl.h>

#include "PGViewState.hpp"
#include "StdinReader.hpp"
#include "PGViewUI.hpp"

using namespace ftxui;

// ---------------------------------------------------------------------------
// Summary printing (after TUI exits)
// ---------------------------------------------------------------------------

/// @brief Format a duration in seconds as human-readable string.
static std::string format_duration(double seconds) {
    if (seconds < 0 || std::isnan(seconds)) return "N/A";
    if (seconds < 60.0) {
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%.1fs", seconds);
        return buf;
    }
    int total = static_cast<int>(std::round(seconds));
    int h = total / 3600;
    int m = (total % 3600) / 60;
    int s = total % 60;
    if (h > 0)
        return std::to_string(h) + "h " + std::to_string(m) + "m " + std::to_string(s) + "s";
    return std::to_string(m) + "m " + std::to_string(s) + "s";
}

/// @brief Print a reconstruction summary to the terminal after the TUI exits.
static void print_summary(PGViewState& state) {
    std::lock_guard<std::mutex> lock(state.mu);

    // App name
    std::string app = state.app_name.empty() ? "PowerGrid" : state.app_name;
    if (!state.app_version.empty())
        app += " v" + state.app_version;

    // Status text + ANSI color
    const char* color_start = "";
    const char* color_end = "\033[0m";
    std::string status;
    if (state.has_session_end && state.exit_code == 0) {
        color_start = "\033[32m"; // green
        status = "completed successfully";
    } else if (state.has_session_end) {
        color_start = "\033[31m"; // red
        status = "failed (exit code " + std::to_string(state.exit_code) + ")";
    } else {
        color_start = "\033[33m"; // yellow
        status = "input closed unexpectedly";
    }

    // Total elapsed time
    double elapsed = state.total_elapsed_seconds();
    size_t images = state.num_completed_images();
    size_t iters = state.total_pcg_iterations();

    // Print
    std::string line(50, '-');

    std::printf("\n%s\n", line.c_str());
    std::printf("  %s %s%s%s\n", app.c_str(), color_start, status.c_str(), color_end);
    std::printf("\n");
    std::printf("  Total time:       %s\n", format_duration(elapsed).c_str());

    if (images > 0) {
        std::printf("  Images:           %zu\n", images);
        std::printf("  Time per image:   %s\n", format_duration(elapsed / static_cast<double>(images)).c_str());
        std::printf("  PCG iterations:   %zu\n", iters);
    }

    double last_err = state.last_error_norm();
    if (last_err >= 0) {
        std::printf("  Final error norm: %.4e\n", last_err);
    }

    std::printf("%s\n\n", line.c_str());
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

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

    // Detect terminal graphics protocol BEFORE entering fullscreen.
    GraphicsProto gfx_proto = detect_graphics_protocol();

    PGViewState state;
    auto screen = ScreenInteractive::Fullscreen();

    // Get an exit closure to call from the reader thread when reconstruction ends
    auto exit_tui = screen.ExitLoopClosure();

    // Start the stdin reader thread, posting custom events on updates.
    // When the reconstruction finishes (exit message or pipe close),
    // schedule the TUI to exit automatically.
    StdinReader reader(state, [&screen, &state, exit_tui]() {
        screen.Post(Event::Custom);  // trigger redraw

        // Auto-exit when reconstruction finishes
        std::lock_guard<std::mutex> lock(state.mu);
        if (state.finished) {
            screen.Post(exit_tui);
        }
    }, data_fd);
    reader.start();

    // Region vector populated during FTXUI's Render pass by ImagePlaceholderNode.
    // Consumed by the Post() task to emit escape sequences.
    // Both accesses happen on the main thread (no mutex needed).
    std::vector<ImageRegion> rendered_regions;

    // Build the UI renderer
    auto renderer = Renderer([&state, gfx_proto, &rendered_regions]() {
        std::lock_guard<std::mutex> lock(state.mu);

        // Clear regions from the previous frame. They were already consumed
        // by the Post() task, but clear anyway for safety.
        rendered_regions.clear();

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

        // Content area: images (quad/side-by-side) + logs, or full-width logs.
        // When a graphics protocol is available, RenderContentArea creates
        // ImagePlaceholderNode elements that push ImageRegion entries (with
        // copied pixel data) into rendered_regions during FTXUI's Render pass.
        layout.push_back(
            RenderContentArea(state, gfx_proto,
                              (gfx_proto != GraphicsProto::None)
                                  ? &rendered_regions : nullptr)
        );

        return vbox(std::move(layout)) | border;
    });

    // Handle keyboard events (manual quit) and graphics protocol emission
    auto component = CatchEvent(renderer, [&](Event event) {
        if (event == Event::Character('q') || event == Event::Escape) {
            screen.Exit();
            return true;
        }

        // After each state update, schedule image emission for the next
        // event loop iteration. The Post() task runs AFTER the current
        // frame is flushed to the terminal. FTXUI's differential rendering
        // won't touch the placeholder cells (spaces that haven't changed),
        // so the inline images persist until the next full redraw.
        if (event == Event::Custom && gfx_proto != GraphicsProto::None) {
            screen.Post([gfx_proto, &rendered_regions]() {
                if (!rendered_regions.empty()) {
                    emit_inline_images(gfx_proto, rendered_regions);
                }
            });
        }

        return false;
    });

    screen.Loop(component);

    // Close pipe fd to unblock reader thread if still running
    // (e.g., if user pressed 'q' before the reconstruction finished)
    if (data_fd >= 0) {
        close(data_fd);
        data_fd = -1;
    }

    // Wait for reader thread to finish
    reader.join();

    // Print summary to the terminal
    print_summary(state);

    return 0;
}
