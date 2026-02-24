#pragma once

/// @file GraphicsProtocol.hpp
/// @brief Terminal inline image support (iTerm2 + Kitty protocols) with Canvas fallback.
///
/// Provides:
/// - GraphicsProto enum + detect_graphics_protocol()
/// - Minimal uncompressed PNG encoder (for iTerm2)
/// - Raw RGBA packer (for Kitty)
/// - emit_inline_images() — synchronous function to write escape sequences

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <unistd.h>

// ---------------------------------------------------------------------------
// Graphics protocol detection
// ---------------------------------------------------------------------------

enum class GraphicsProto { None, ITerm2, Kitty };

/// @brief Detect the terminal's graphics protocol capability.
///
/// Call this BEFORE entering FTXUI fullscreen. Checks TERM_PROGRAM and TERM
/// environment variables.
inline GraphicsProto detect_graphics_protocol() {
    const char* tp = std::getenv("TERM_PROGRAM");
    const char* term = std::getenv("TERM");

    if (tp) {
        std::string s(tp);
        if (s == "iTerm.app")  return GraphicsProto::ITerm2;
        if (s == "WezTerm")    return GraphicsProto::Kitty;
        if (s == "ghostty")    return GraphicsProto::Kitty;
    }
    if (term) {
        std::string t(term);
        if (t == "xterm-kitty")   return GraphicsProto::Kitty;
        if (t == "xterm-ghostty") return GraphicsProto::Kitty;
    }
    return GraphicsProto::None;
}

// ---------------------------------------------------------------------------
// Base64 encoder (for emitting image data to the terminal)
// ---------------------------------------------------------------------------

namespace detail {

inline std::string base64_encode(const uint8_t* data, size_t len) {
    static const char table[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    out.reserve(((len + 2) / 3) * 4);
    for (size_t i = 0; i < len; i += 3) {
        uint32_t a = data[i];
        uint32_t b = (i + 1 < len) ? data[i + 1] : 0;
        uint32_t c = (i + 2 < len) ? data[i + 2] : 0;
        uint32_t triple = (a << 16) | (b << 8) | c;
        out += table[(triple >> 18) & 0x3F];
        out += table[(triple >> 12) & 0x3F];
        out += (i + 1 < len) ? table[(triple >> 6) & 0x3F] : '=';
        out += (i + 2 < len) ? table[triple & 0x3F] : '=';
    }
    return out;
}

} // namespace detail

// ---------------------------------------------------------------------------
// CRC-32 and Adler-32 (for PNG)
// ---------------------------------------------------------------------------

namespace detail {

inline uint32_t crc32(const uint8_t* data, size_t len) {
    uint32_t crc = 0xFFFFFFFF;
    for (size_t i = 0; i < len; ++i) {
        crc ^= data[i];
        for (int j = 0; j < 8; ++j) {
            if (crc & 1)
                crc = (crc >> 1) ^ 0xEDB88320;
            else
                crc >>= 1;
        }
    }
    return crc ^ 0xFFFFFFFF;
}

inline uint32_t adler32(const uint8_t* data, size_t len) {
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < len; ++i) {
        a = (a + data[i]) % 65521;
        b = (b + a) % 65521;
    }
    return (b << 16) | a;
}

} // namespace detail

// ---------------------------------------------------------------------------
// Minimal uncompressed PNG encoder (grayscale 8-bit) — for iTerm2
// ---------------------------------------------------------------------------

/// @brief Encode a grayscale 8-bit image as an uncompressed PNG.
///
/// Uses zlib "stored" blocks (no compression), so no zlib/libpng dependency.
inline std::vector<uint8_t> encode_grayscale_png(const uint8_t* gray,
                                                  uint32_t width,
                                                  uint32_t height) {
    size_t raw_row = 1 + width;
    size_t raw_size = raw_row * height;

    std::vector<uint8_t> raw(raw_size);
    for (uint32_t y = 0; y < height; ++y) {
        raw[y * raw_row] = 0x00; // filter: None
        std::memcpy(&raw[y * raw_row + 1], &gray[y * width], width);
    }

    // zlib stored blocks
    std::vector<uint8_t> zlib_data;
    zlib_data.push_back(0x78);
    zlib_data.push_back(0x01);

    size_t offset = 0;
    while (offset < raw_size) {
        size_t chunk = std::min(raw_size - offset, size_t(65535));
        bool last = (offset + chunk >= raw_size);
        zlib_data.push_back(last ? 0x01 : 0x00);
        uint16_t len = static_cast<uint16_t>(chunk);
        uint16_t nlen = ~len;
        zlib_data.push_back(len & 0xFF);
        zlib_data.push_back((len >> 8) & 0xFF);
        zlib_data.push_back(nlen & 0xFF);
        zlib_data.push_back((nlen >> 8) & 0xFF);
        zlib_data.insert(zlib_data.end(), raw.begin() + offset,
                         raw.begin() + offset + chunk);
        offset += chunk;
    }

    uint32_t adler = detail::adler32(raw.data(), raw.size());
    zlib_data.push_back((adler >> 24) & 0xFF);
    zlib_data.push_back((adler >> 16) & 0xFF);
    zlib_data.push_back((adler >> 8) & 0xFF);
    zlib_data.push_back(adler & 0xFF);

    // Assemble PNG
    std::vector<uint8_t> png;
    png.reserve(8 + 25 + zlib_data.size() + 12 + 12);

    auto write_be32 = [&](uint32_t v) {
        png.push_back((v >> 24) & 0xFF);
        png.push_back((v >> 16) & 0xFF);
        png.push_back((v >> 8) & 0xFF);
        png.push_back(v & 0xFF);
    };

    const uint8_t sig[] = {137, 80, 78, 71, 13, 10, 26, 10};
    png.insert(png.end(), sig, sig + 8);

    // IHDR
    {
        uint8_t ihdr[13];
        ihdr[0] = (width >> 24) & 0xFF;  ihdr[1] = (width >> 16) & 0xFF;
        ihdr[2] = (width >> 8) & 0xFF;   ihdr[3] = width & 0xFF;
        ihdr[4] = (height >> 24) & 0xFF;  ihdr[5] = (height >> 16) & 0xFF;
        ihdr[6] = (height >> 8) & 0xFF;   ihdr[7] = height & 0xFF;
        ihdr[8] = 8; ihdr[9] = 0; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;

        write_be32(13);
        size_t ts = png.size();
        png.push_back('I'); png.push_back('H'); png.push_back('D'); png.push_back('R');
        png.insert(png.end(), ihdr, ihdr + 13);
        write_be32(detail::crc32(&png[ts], 4 + 13));
    }

    // IDAT
    {
        uint32_t idat_len = static_cast<uint32_t>(zlib_data.size());
        write_be32(idat_len);
        size_t ts = png.size();
        png.push_back('I'); png.push_back('D'); png.push_back('A'); png.push_back('T');
        png.insert(png.end(), zlib_data.begin(), zlib_data.end());
        write_be32(detail::crc32(&png[ts], 4 + idat_len));
    }

    // IEND
    {
        write_be32(0);
        size_t ts = png.size();
        png.push_back('I'); png.push_back('E'); png.push_back('N'); png.push_back('D');
        write_be32(detail::crc32(&png[ts], 4));
    }

    return png;
}

// ---------------------------------------------------------------------------
// ImageRegion — owns its pixel data for safe cross-frame usage
// ---------------------------------------------------------------------------

/// @brief Describes one image region to render on the terminal.
///
/// Owns a copy of the pixel data so it remains valid after the state mutex
/// is released and FTXUI's element tree is rebuilt.
struct ImageRegion {
    int term_row = 0;   ///< 1-based terminal row
    int term_col = 0;   ///< 1-based terminal column
    int cols = 0;       ///< Width in terminal columns
    int rows = 0;       ///< Height in terminal rows
    std::vector<float> pixels; ///< Owned copy of pixel data
    size_t px_w = 0;    ///< Image width in pixels
    size_t px_h = 0;    ///< Image height in pixels
};

// ---------------------------------------------------------------------------
// Synchronous image emission — no thread needed
// ---------------------------------------------------------------------------

namespace detail {

inline void emit_iterm2(std::string& out, const ImageRegion& r) {
    std::vector<uint8_t> gray(r.px_w * r.px_h);
    for (size_t i = 0; i < r.px_w * r.px_h; ++i) {
        gray[i] = static_cast<uint8_t>(
            std::min(std::max(r.pixels[i], 0.0f), 1.0f) * 255.0f);
    }

    auto png = encode_grayscale_png(gray.data(),
                                     static_cast<uint32_t>(r.px_w),
                                     static_cast<uint32_t>(r.px_h));

    std::string b64 = base64_encode(png.data(), png.size());

    out += "\033]1337;File=inline=1"
           ";width=" + std::to_string(r.cols) +
           ";height=" + std::to_string(r.rows) +
           ";preserveAspectRatio=0"
           ":" + b64 + "\a";
}

inline void emit_kitty(std::string& out, const ImageRegion& r) {
    // Pack as RGBA
    std::vector<uint8_t> rgba(r.px_w * r.px_h * 4);
    for (size_t i = 0; i < r.px_w * r.px_h; ++i) {
        uint8_t g = static_cast<uint8_t>(
            std::min(std::max(r.pixels[i], 0.0f), 1.0f) * 255.0f);
        rgba[i * 4 + 0] = g;
        rgba[i * 4 + 1] = g;
        rgba[i * 4 + 2] = g;
        rgba[i * 4 + 3] = 255;
    }

    std::string b64 = base64_encode(rgba.data(), rgba.size());

    const size_t chunk_size = 4096;
    size_t offset = 0;
    bool first = true;

    while (offset < b64.size()) {
        size_t remaining = b64.size() - offset;
        size_t chunk = std::min(remaining, chunk_size);
        bool last = (offset + chunk >= b64.size());

        out += "\033_G";
        if (first) {
            out += "a=T,f=32"
                   ",s=" + std::to_string(r.px_w) +
                   ",v=" + std::to_string(r.px_h) +
                   ",C=1";
            first = false;
        }
        out += last ? ",m=0;" : ",m=1;";
        out += b64.substr(offset, chunk);
        out += "\033\\";

        offset += chunk;
    }
}

} // namespace detail

/// @brief Emit inline images to the terminal using escape sequences.
///
/// Call this from the main thread (e.g., a Post() task) AFTER FTXUI has
/// flushed its frame. FTXUI's differential rendering won't overwrite
/// placeholder cells (spaces) that haven't changed between frames, so
/// the images persist until the next full redraw.
///
/// For Kitty protocol, old images are deleted first to avoid stacking.
inline void emit_inline_images(GraphicsProto proto,
                                const std::vector<ImageRegion>& regions) {
    if (proto == GraphicsProto::None || regions.empty()) return;

    std::string output;

    // Kitty: delete all placements before re-emitting (prevents stacking)
    if (proto == GraphicsProto::Kitty) {
        output += "\033_Ga=d;\033\\";
    }

    // Save cursor
    output += "\0337";

    for (const auto& r : regions) {
        if (r.pixels.empty() || r.px_w == 0 || r.px_h == 0) continue;

        // Position cursor
        output += "\033[" + std::to_string(r.term_row) + ";" +
                  std::to_string(r.term_col) + "H";

        if (proto == GraphicsProto::ITerm2) {
            detail::emit_iterm2(output, r);
        } else if (proto == GraphicsProto::Kitty) {
            detail::emit_kitty(output, r);
        }
    }

    // Restore cursor
    output += "\0338";

    // Write all at once to minimize flicker
    ::write(STDOUT_FILENO, output.data(), output.size());
}
