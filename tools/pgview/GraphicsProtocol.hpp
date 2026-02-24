#pragma once

/// @file GraphicsProtocol.hpp
/// @brief Terminal inline image support (iTerm2 + Kitty protocols) with Canvas fallback.
///
/// Provides:
/// - GraphicsProto enum + detect_graphics_protocol()
/// - Minimal uncompressed PNG encoder (for iTerm2)
/// - Raw RGBA packer (for Kitty)
/// - ImageOverlay — thread that emits terminal escape sequences after FTXUI renders

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>
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
// CRC-32 (for PNG)
// ---------------------------------------------------------------------------

namespace detail {

inline uint32_t crc32(const uint8_t* data, size_t len) {
    // Standard CRC-32 (ISO 3309 / PNG)
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
/// The resulting PNG is valid but larger than a compressed one.
inline std::vector<uint8_t> encode_grayscale_png(const uint8_t* gray,
                                                  uint32_t width,
                                                  uint32_t height) {
    // Each PNG row: filter_byte(0x00) + width gray bytes
    size_t raw_row = 1 + width;
    size_t raw_size = raw_row * height;

    // Build raw (unfiltered) image data with filter byte = 0 per row
    std::vector<uint8_t> raw(raw_size);
    for (uint32_t y = 0; y < height; ++y) {
        raw[y * raw_row] = 0x00; // filter: None
        std::memcpy(&raw[y * raw_row + 1], &gray[y * width], width);
    }

    // Build zlib stream with stored blocks
    // zlib header: 0x78 0x01 (deflate, no compression)
    // For stored blocks: each block is max 65535 bytes
    // Block header: BFINAL(1 byte) LEN(2LE) NLEN(2LE) DATA
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

    // Adler-32 checksum of the raw data
    uint32_t adler = detail::adler32(raw.data(), raw.size());
    zlib_data.push_back((adler >> 24) & 0xFF);
    zlib_data.push_back((adler >> 16) & 0xFF);
    zlib_data.push_back((adler >> 8) & 0xFF);
    zlib_data.push_back(adler & 0xFF);

    // --- Assemble PNG ---
    std::vector<uint8_t> png;
    png.reserve(8 + 25 + zlib_data.size() + 12 + 12);

    // Helper: write big-endian 32-bit
    auto write_be32 = [&](uint32_t v) {
        png.push_back((v >> 24) & 0xFF);
        png.push_back((v >> 16) & 0xFF);
        png.push_back((v >> 8) & 0xFF);
        png.push_back(v & 0xFF);
    };

    // PNG signature
    const uint8_t sig[] = {137, 80, 78, 71, 13, 10, 26, 10};
    png.insert(png.end(), sig, sig + 8);

    // IHDR chunk (13 bytes data)
    {
        uint8_t ihdr[13];
        ihdr[0] = (width >> 24) & 0xFF;  ihdr[1] = (width >> 16) & 0xFF;
        ihdr[2] = (width >> 8) & 0xFF;   ihdr[3] = width & 0xFF;
        ihdr[4] = (height >> 24) & 0xFF;  ihdr[5] = (height >> 16) & 0xFF;
        ihdr[6] = (height >> 8) & 0xFF;   ihdr[7] = height & 0xFF;
        ihdr[8] = 8;  // bit depth
        ihdr[9] = 0;  // color type: grayscale
        ihdr[10] = 0; // compression
        ihdr[11] = 0; // filter
        ihdr[12] = 0; // interlace

        write_be32(13); // length
        size_t type_start = png.size();
        png.push_back('I'); png.push_back('H'); png.push_back('D'); png.push_back('R');
        png.insert(png.end(), ihdr, ihdr + 13);
        write_be32(detail::crc32(&png[type_start], 4 + 13));
    }

    // IDAT chunk
    {
        uint32_t idat_len = static_cast<uint32_t>(zlib_data.size());
        write_be32(idat_len);
        size_t type_start = png.size();
        png.push_back('I'); png.push_back('D'); png.push_back('A'); png.push_back('T');
        png.insert(png.end(), zlib_data.begin(), zlib_data.end());
        write_be32(detail::crc32(&png[type_start], 4 + idat_len));
    }

    // IEND chunk
    {
        write_be32(0);
        size_t type_start = png.size();
        png.push_back('I'); png.push_back('E'); png.push_back('N'); png.push_back('D');
        write_be32(detail::crc32(&png[type_start], 4));
    }

    return png;
}

// ---------------------------------------------------------------------------
// Raw RGBA packer — for Kitty
// ---------------------------------------------------------------------------

/// @brief Pack grayscale float pixels into RGBA bytes (gray, gray, gray, 255).
inline std::vector<uint8_t> pack_grayscale_rgba(const float* pixels,
                                                 size_t width, size_t height) {
    std::vector<uint8_t> rgba(width * height * 4);
    for (size_t i = 0; i < width * height; ++i) {
        uint8_t g = static_cast<uint8_t>(std::min(std::max(pixels[i], 0.0f), 1.0f) * 255.0f);
        rgba[i * 4 + 0] = g;
        rgba[i * 4 + 1] = g;
        rgba[i * 4 + 2] = g;
        rgba[i * 4 + 3] = 255;
    }
    return rgba;
}

// ---------------------------------------------------------------------------
// ImageOverlay — thread for emitting inline images after FTXUI renders
// ---------------------------------------------------------------------------

/// @brief Describes one image region to render on the terminal.
struct ImageRegion {
    int term_row = 0;   ///< 1-based terminal row
    int term_col = 0;   ///< 1-based terminal column
    int cols = 0;       ///< Width in terminal columns
    int rows = 0;       ///< Height in terminal rows
    const float* pixels = nullptr; ///< Pointer to pixel data
    size_t px_w = 0;    ///< Image width in pixels
    size_t px_h = 0;    ///< Image height in pixels
};

/// @brief Thread that emits inline terminal images after FTXUI finishes rendering.
///
/// Usage:
///   1. Create ImageOverlay with the detected protocol
///   2. During FTXUI rendering, record image regions via set_regions()
///   3. After FTXUI flushes its frame, call notify() to trigger emission
///   4. The overlay thread emits escape sequences to draw images at recorded positions
class ImageOverlay {
public:
    explicit ImageOverlay(GraphicsProto proto)
        : proto_(proto) {
        if (proto_ != GraphicsProto::None) {
            thread_ = std::thread([this]() { run(); });
        }
    }

    ~ImageOverlay() {
        stop();
    }

    /// @brief Set the image regions to render on the next frame.
    void set_regions(std::vector<ImageRegion> regions) {
        std::lock_guard<std::mutex> lock(mu_);
        pending_regions_ = std::move(regions);
    }

    /// @brief Notify the overlay thread that FTXUI has finished rendering a frame.
    void notify() {
        {
            std::lock_guard<std::mutex> lock(mu_);
            frame_ready_ = true;
        }
        cv_.notify_one();
    }

    /// @brief Stop the overlay thread.
    void stop() {
        {
            std::lock_guard<std::mutex> lock(mu_);
            stop_ = true;
            frame_ready_ = true; // unblock the wait
        }
        cv_.notify_one();
        if (thread_.joinable())
            thread_.join();
    }

    GraphicsProto protocol() const { return proto_; }

private:
    void run() {
        while (true) {
            std::vector<ImageRegion> regions;
            {
                std::unique_lock<std::mutex> lock(mu_);
                cv_.wait(lock, [this]() { return frame_ready_; });
                if (stop_) return;
                frame_ready_ = false;
                regions = std::move(pending_regions_);
                pending_regions_.clear();
            }

            if (regions.empty()) continue;

            // Small delay to let FTXUI flush its frame to the terminal
            std::this_thread::sleep_for(std::chrono::milliseconds(5));

            // Build the complete escape sequence for all regions
            std::string output;

            // Save cursor position
            output += "\0337"; // ESC 7 — save cursor

            for (const auto& r : regions) {
                if (!r.pixels || r.px_w == 0 || r.px_h == 0) continue;

                // Position cursor at the image region
                output += "\033[" + std::to_string(r.term_row) + ";" +
                          std::to_string(r.term_col) + "H";

                if (proto_ == GraphicsProto::ITerm2) {
                    emit_iterm2(output, r);
                } else if (proto_ == GraphicsProto::Kitty) {
                    emit_kitty(output, r);
                }
            }

            // Restore cursor position
            output += "\0338"; // ESC 8 — restore cursor

            // Write all at once to minimize flicker
            if (!output.empty()) {
                ::write(STDOUT_FILENO, output.data(), output.size());
            }
        }
    }

    void emit_iterm2(std::string& out, const ImageRegion& r) {
        // Convert float pixels to grayscale uint8
        std::vector<uint8_t> gray(r.px_w * r.px_h);
        for (size_t i = 0; i < r.px_w * r.px_h; ++i) {
            gray[i] = static_cast<uint8_t>(
                std::min(std::max(r.pixels[i], 0.0f), 1.0f) * 255.0f);
        }

        // Encode as PNG
        auto png = encode_grayscale_png(gray.data(),
                                         static_cast<uint32_t>(r.px_w),
                                         static_cast<uint32_t>(r.px_h));

        // Base64 encode
        std::string b64 = detail::base64_encode(png.data(), png.size());

        // iTerm2 inline image protocol
        // ESC ] 1337 ; File=inline=1;width=Ncols;height=Nrows;preserveAspectRatio=0 : BASE64 BEL
        out += "\033]1337;File=inline=1"
               ";width=" + std::to_string(r.cols) +
               ";height=" + std::to_string(r.rows) +
               ";preserveAspectRatio=0"
               ":" + b64 + "\a";
    }

    void emit_kitty(std::string& out, const ImageRegion& r) {
        // Pack as RGBA
        auto rgba = pack_grayscale_rgba(r.pixels, r.px_w, r.px_h);

        // Base64 encode
        std::string b64 = detail::base64_encode(rgba.data(), rgba.size());

        // Kitty graphics protocol: chunked transmission
        // f=32 = raw RGBA, s=width, v=height, C=1 = don't move cursor
        // Chunks of 4096 base64 chars: m=1 (more), m=0 (last)
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

    GraphicsProto proto_;
    std::thread thread_;
    std::mutex mu_;
    std::condition_variable cv_;
    bool frame_ready_ = false;
    bool stop_ = false;
    std::vector<ImageRegion> pending_regions_;
};
