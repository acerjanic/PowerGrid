#!/bin/bash
# ===========================================================================
# benchmark_recon.sh
# Benchmark PowerGrid reconstruction at different matrix sizes.
# Reports median wall-clock time across multiple runs.
#
# Usage: ./scripts/benchmark_recon.sh [BUILD_DIR] [NUM_RUNS]
# ===========================================================================
set -e

BUILD=${1:-build_native}
NRUNS=${2:-3}
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Resolve build directory to absolute path
BUILD=$(cd "$BUILD" && pwd)

echo "============================================"
echo " PowerGrid Reconstruction Benchmark"
echo "============================================"
echo "Build dir: $BUILD"
echo "Runs per config: $NRUNS"
echo ""
printf "%-10s %-6s %-7s %-5s %-6s %-12s %-20s\n" \
       "Size" "Coils" "Spokes" "Nro" "Iters" "Time (s)" "App"
printf "%-10s %-6s %-7s %-5s %-6s %-12s %-20s\n" \
       "----" "-----" "------" "---" "-----" "--------" "---"

bench() {
    local label=$1 app=$2 size=$3 nc=$4 nk=$5 nro=$6 niter=$7
    shift 7
    local gen_extra=("$@")

    local file="$TMPDIR/${label}_${size}.h5"
    local outdir="$TMPDIR/${label}_${size}_out"

    # Generate synthetic data (only once per config)
    "$BUILD/GenSyntheticISMRMRD" -o "$file" -x "$size" -y "$size" \
        -c "$nc" -k "$nk" -r "$nro" "${gen_extra[@]}" 2>/dev/null

    local times=()
    for ((r=1; r<=NRUNS; r++)); do
        rm -rf "$outdir"
        mkdir -p "$outdir"

        local t
        t=$( { /usr/bin/time -p "$BUILD/$app" -i "$file" -o "$outdir/" \
               -n "$niter" -B 0.001 "${APP_ARGS[@]}" ; } 2>&1 | grep ^real | awk '{print $2}')
        times+=("$t")
    done

    # Sort and pick median
    IFS=$'\n' sorted=($(printf '%s\n' "${times[@]}" | sort -g)); unset IFS
    local median=${sorted[$((NRUNS/2))]}

    printf "%-10s %-6s %-7s %-5s %-6s %-12s %-20s\n" \
           "${size}x${size}" "$nc" "$nk" "$nro" "$niter" "$median" "$app"
}

# -----------------------------------------------------------------------
# SENSE + NUFFT benchmarks (no field map, L=1)
# -----------------------------------------------------------------------
echo ""
echo "--- SENSE + NUFFT (L=1, no field map) ---"
APP_ARGS=(-F NUFFT -t 1)
bench sense PowerGridIsmrmrd 32  4  64   64  10
bench sense PowerGridIsmrmrd 64  4  128  128 10
bench sense PowerGridIsmrmrd 128 8  256  256 10

# -----------------------------------------------------------------------
# SENSE + NUFFT + TimeSegmentation benchmarks (L=4, +/- 50 Hz)
# -----------------------------------------------------------------------
echo ""
echo "--- SENSE + NUFFT + TimeSeg (L=4, +/-50 Hz) ---"
APP_ARGS=(-F NUFFT -t 4)
bench sense_fm PowerGridIsmrmrd 32  4  64   64  10 -f 50
bench sense_fm PowerGridIsmrmrd 64  4  128  128 10 -f 50
bench sense_fm PowerGridIsmrmrd 128 8  256  256 10 -f 50

# -----------------------------------------------------------------------
# pcSENSE benchmarks (2 shots)
# -----------------------------------------------------------------------
echo ""
echo "--- pcSENSE (2 shots) ---"
APP_ARGS=(-s 2)
bench pcsense PowerGridPcSense 32  4  32   32  10 -s 2
bench pcsense PowerGridPcSense 64  4  64   64  10 -s 2

echo ""
echo "============================================"
echo " Benchmark complete."
echo "============================================"
