#pragma once

#include <deque>
#include <string>
#include <algorithm>
#include <cmath>

/// @brief Render a sparkline from a sequence of doubles.
///
/// Maps values to Unicode block characters ▁▂▃▄▅▆▇█ (8 levels)
/// scaled between the min and max of the input buffer.
///
/// @param values  Recent metric samples.
/// @param width   Maximum number of characters in the sparkline.
/// @returns       UTF-8 string of sparkline characters.
inline std::string render_sparkline(const std::deque<double>& values, size_t width = 40) {
    if (values.empty()) return "";

    // Unicode block elements: 8 levels
    static const char* blocks[] = {
        "\u2581", // ▁
        "\u2582", // ▂
        "\u2583", // ▃
        "\u2584", // ▄
        "\u2585", // ▅
        "\u2586", // ▆
        "\u2587", // ▇
        "\u2588", // █
    };

    // Determine how many values to use (last `width` values)
    size_t n = std::min(values.size(), width);
    size_t start = values.size() - n;

    // Find min/max for scaling
    double vmin = values[start];
    double vmax = values[start];
    for (size_t i = start; i < values.size(); ++i) {
        vmin = std::min(vmin, values[i]);
        vmax = std::max(vmax, values[i]);
    }

    std::string result;
    result.reserve(n * 4); // UTF-8 block chars are ~3 bytes each

    double range = vmax - vmin;
    for (size_t i = start; i < values.size(); ++i) {
        int idx = 0;
        if (range > 0) {
            double norm = (values[i] - vmin) / range;
            idx = static_cast<int>(norm * 7.0);
            idx = std::clamp(idx, 0, 7);
        } else {
            idx = 4; // flat line in the middle
        }
        result += blocks[idx];
    }
    return result;
}
