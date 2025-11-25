#!/bin/bash
# Documentation script for FS3
# Usage: ./scripts/docs.sh [generate|open|clean|rebuild|pdf]

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"
DOCS_OUTPUT="$BUILD_DIR/docs/html"

ACTION="${1:-generate}"

case "$ACTION" in
    generate)
        echo "Generating documentation..."
        if [ ! -d "$BUILD_DIR" ]; then
            echo "Error: Build directory does not exist!"
            echo "Please run ./scripts/build.sh first."
            exit 1
        fi
        cd "$BUILD_DIR"
        cmake --build . --target docs
        echo "Documentation generated at: $DOCS_OUTPUT/index.html"
        ;;
    
    open)
        if [ ! -f "$DOCS_OUTPUT/index.html" ]; then
            echo "Documentation not found. Generating first..."
            "$0" generate
        fi
        echo "Opening documentation..."
        if command -v xdg-open &> /dev/null; then
            xdg-open "$DOCS_OUTPUT/index.html"
        elif command -v open &> /dev/null; then
            open "$DOCS_OUTPUT/index.html"
        elif command -v start &> /dev/null; then
            start "$DOCS_OUTPUT/index.html"
        else
            echo "Could not detect browser opener. Please open manually:"
            echo "$DOCS_OUTPUT/index.html"
        fi
        ;;
    
    clean)
        echo "Cleaning documentation..."
        if [ -d "$BUILD_DIR/docs" ]; then
            rm -rf "$BUILD_DIR/docs"
            echo "Documentation cleaned."
        else
            echo "No documentation to clean."
        fi
        ;;
    
    pdf)
        echo "Generating PDF documentation..."
        if [ ! -d "$BUILD_DIR" ]; then
            echo "Error: Build directory does not exist!"
            echo "Please run ./scripts/build.sh first."
            exit 1
        fi
        cd "$BUILD_DIR"
        if [ -d "docs/latex" ]; then
            cd docs/latex
            make
            echo "PDF generated at: $(pwd)/refman.pdf"
        else
            echo "LaTeX documentation not configured. Please enable in Doxyfile."
            exit 1
        fi
        ;;
    
    *)
        echo "Usage: $0 [generate|open|clean|rebuild|pdf]"
        echo ""
        echo "Commands:"
        echo "  generate  - Generate Doxygen documentation"
        echo "  open      - Open documentation in default browser"
        echo "  clean     - Remove generated documentation"
        echo "  pdf       - Generate PDF documentation (requires LaTeX)"
        exit 1
        ;;
esac
