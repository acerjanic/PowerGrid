#!/usr/bin/env bash
# Generate PowerGrid API documentation via hdoc.
#
# Intended to run inside the powergrid-hdoc Docker image, which provides
# both hdoc and cmake.  See docker/pg-hdoc/Dockerfile.
#
# Usage (from the repo root on the host):
#   docker run --rm -v "$(pwd)":/root/PowerGrid powergrid-hdoc \
#       bash scripts/generate-docs.sh
#
# Output:
#   docs/api/index.html
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${ROOT}"

echo "==> Generating compile_commands.json..."
cmake -B build \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DOPENACC_GPU=OFF \
    -DMPISupport=OFF

echo "==> Running hdoc..."
hdoc --verbose

echo ""
echo "Docs written to: ${ROOT}/docs/api/index.html"
