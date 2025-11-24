#!/bin/bash
# Install script for FS3 (from existing build)
# Usage: ./scripts/install.sh [custom-install-path]

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"

# Use custom install path if provided, otherwise use default
INSTALL_DIR="${1:-$PROJECT_ROOT/install}"

echo "=========================================="
echo "FS3 Install Script"
echo "=========================================="
echo "Install Directory: $INSTALL_DIR"
echo "=========================================="

if [ ! -d "$BUILD_DIR" ]; then
    echo "Error: Build directory does not exist!"
    echo "Please run ./scripts/build.sh first."
    exit 1
fi

cd "$BUILD_DIR"

# Install with custom prefix if specified
if [ -n "$1" ]; then
    echo "Installing to custom location: $INSTALL_DIR"
    cmake --install . --prefix "$INSTALL_DIR"
else
    echo "Installing to configured location: $INSTALL_DIR"
    cmake --install .
fi

echo ""
echo "=========================================="
echo "Installation completed successfully!"
echo "=========================================="
echo "To use FS3 in your project, add to CMakeLists.txt:"
echo ""
echo "  list(APPEND CMAKE_PREFIX_PATH \"$INSTALL_DIR\")"
echo "  find_package(FS3 REQUIRED)"
echo "  target_link_libraries(your_target PRIVATE FS3::FS3)"
echo "=========================================="
