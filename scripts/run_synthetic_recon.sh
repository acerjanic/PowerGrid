#!/bin/bash
# ===========================================================================
# run_synthetic_recon.sh
# End-to-end test: generate synthetic ISMRMRD data, reconstruct with CLI apps,
# validate that output files are produced.
#
# Usage: ./scripts/run_synthetic_recon.sh [BUILD_DIR]
# ===========================================================================
set -e

BUILD=${1:-build_native}
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

PASS=0
FAIL=0

# Resolve build directory to absolute path
BUILD=$(cd "$BUILD" && pwd)

echo "============================================"
echo " PowerGrid Synthetic Reconstruction Tests"
echo "============================================"
echo "Build dir: $BUILD"
echo "Temp dir:  $TMPDIR"
echo ""

run_test() {
    local name="$1"; shift
    echo "--- $name ---"
    if "$@"; then
        echo "  PASS: $name"
        ((PASS++)) || true
    else
        echo "  FAIL: $name (exit code $?)"
        ((FAIL++)) || true
    fi
    echo ""
}

check_output() {
    local outdir="$1"
    local pattern="$2"
    local count
    count=$(find "$outdir" -name "$pattern" 2>/dev/null | wc -l)
    if [ "$count" -gt 0 ]; then
        echo "  Output files found: $count matching $pattern"
        return 0
    else
        echo "  WARNING: no output files matching $pattern in $outdir"
        return 1
    fi
}

# -----------------------------------------------------------------------
# Test 1: SENSE + NUFFT (no field map, L=1)
# -----------------------------------------------------------------------
echo "Generating SENSE test data..."
"$BUILD/GenSyntheticISMRMRD" -o "$TMPDIR/sense.h5" -x 64 -y 64 -c 4 -k 128 -r 128

mkdir -p "$TMPDIR/recon1"
run_test "PowerGridIsmrmrd (NUFFT, no field map, L=1)" \
    /usr/bin/time -p "$BUILD/PowerGridIsmrmrd" \
    -i "$TMPDIR/sense.h5" -o "$TMPDIR/recon1/" -F NUFFT -t 1 -n 10 -B 0.001

check_output "$TMPDIR/recon1" "*.nii"

# -----------------------------------------------------------------------
# Test 2: SENSE + NUFFT + field map (L=4)
# -----------------------------------------------------------------------
echo "Generating SENSE + field map test data..."
"$BUILD/GenSyntheticISMRMRD" -o "$TMPDIR/sense_fm.h5" -x 64 -y 64 -c 4 -k 128 -r 128 -f 50

mkdir -p "$TMPDIR/recon2"
run_test "PowerGridIsmrmrd (NUFFT, field map, L=4)" \
    /usr/bin/time -p "$BUILD/PowerGridIsmrmrd" \
    -i "$TMPDIR/sense_fm.h5" -o "$TMPDIR/recon2/" -F NUFFT -t 4 -n 10 -B 0.001

check_output "$TMPDIR/recon2" "*.nii"

# -----------------------------------------------------------------------
# Test 3: pcSENSE (2 shots)
# -----------------------------------------------------------------------
echo "Generating pcSENSE test data (2 shots)..."
"$BUILD/GenSyntheticISMRMRD" -o "$TMPDIR/pcsense.h5" -x 64 -y 64 -c 4 -k 128 -r 64 -s 2

mkdir -p "$TMPDIR/recon3"
run_test "PowerGridPcSense (2 shots)" \
    /usr/bin/time -p "$BUILD/PowerGridPcSense" \
    -i "$TMPDIR/pcsense.h5" -o "$TMPDIR/recon3/" -s 2 -n 10 -B 0.001

check_output "$TMPDIR/recon3" "*.nii"

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
echo "============================================"
echo " Results: $PASS passed, $FAIL failed"
echo "============================================"

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
