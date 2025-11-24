#!/bin/bash
# Build script for FS3
# Usage: ./scripts/build.sh [Release|Debug] [additional cmake args]

set -e  # Exit on error

# Determine build type (default: Release)
BUILD_TYPE="${1:-Release}"
shift || true  # Remove first argument if it exists

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Build and install directories
BUILD_DIR="$PROJECT_ROOT/build"
INSTALL_DIR="$PROJECT_ROOT/install"

echo "=========================================="
echo "FS3 Build Script"
echo "=========================================="
echo "Build Type: $BUILD_TYPE"
echo "Project Root: $PROJECT_ROOT"
echo "Build Directory: $BUILD_DIR"
echo "Install Directory: $INSTALL_DIR"
echo "Additional CMake Args: $@"
echo "=========================================="

# Remove old build directory
if [ -d "$BUILD_DIR" ]; then
    echo "Removing old build directory..."
    rm -rf "$BUILD_DIR"
fi

# Create build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure with CMake
echo ""
echo "Configuring with CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    "$@"

# Build
echo ""
echo "Building..."
cmake --build . --parallel $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# Install
echo ""
echo "Installing to $INSTALL_DIR..."
cmake --install .

echo ""
echo "=========================================="
echo "Build completed successfully!"
echo "=========================================="
echo "Libraries installed to: $INSTALL_DIR/lib"
echo "Headers installed to: $INSTALL_DIR/include"
echo "CMake config installed to: $INSTALL_DIR/lib/cmake/FS3"
echo "=========================================="
