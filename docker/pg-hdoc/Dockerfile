FROM ghcr.io/acerjanic/powergrid-hpcsdk:HPCSDK_25.1

# Install hdoc build dependencies.
# hdoc's nix flake targets LLVM 14; use the same version here for compatibility.
RUN apt-get update && apt-get install -y \
        llvm-14-dev \
        libclang-14-dev \
        clang-14 \
        meson \
        ninja-build \
        xxd \
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/*

# Make llvm-config-14 the default llvm-config so meson's LLVM detection finds it.
RUN update-alternatives --install /usr/bin/llvm-config llvm-config \
        /usr/bin/llvm-config-14 100

# Build hdoc from source (tag 1.4.1).
# All non-LLVM dependencies are vendored in subprojects/.
RUN git clone https://github.com/hdoc/hdoc.git /tmp/hdoc \
    && cd /tmp/hdoc \
    && git checkout 1.4.1 \
    && git submodule update --init --recursive \
    && meson setup build \
    && ninja -C build \
    && cp /tmp/hdoc/build/hdoc /usr/local/bin/hdoc \
    && rm -rf /tmp/hdoc

WORKDIR /root/PowerGrid
