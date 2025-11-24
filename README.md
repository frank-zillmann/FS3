# FS³ - Fast and Flexible Framework for Simple Simulations of Separation-Processes
**Created by Frank Zillmann during his Bachelor's Thesis: High-Performance Process Simulations for Magnetic Separation of Biomolecules**

See a [**Simulation of a hIgG Purification Process using the Rotor Stator High Gradient Magnetic Separation**](https://github.com/frank-zillmann/HGMS-hIgG-separation) as an example application of this framework.

## Setup and Build Instructions

### Prerequisites

#### Required Dependencies
- **C++ Compiler**: GCC 10+ or Clang 12+ with C++20 support
- **CMake**: Version 3.16 or higher
- **SUNDIALS**: Version 6.0+ with components: cvode, nvecserial, sunlinsolband, sunlinsolklu, sunmatrixsparse
- **ZLIB**: Development libraries

#### Optional Dependencies
- **Doxygen**: For generating documentation (with graphviz for diagrams)
- **LaTeX**: For PDF documentation generation

#### Platform-Specific Installation

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install build-essential cmake libsundials-dev zlib1g-dev
# Optional for documentation:
sudo apt install doxygen graphviz texlive-latex-base
```

**macOS (using Homebrew):**
```bash
brew install cmake sundials zlib
# Optional for documentation:
brew install doxygen graphviz
# For PDF: brew install --cask mactex
```

**Windows (using vcpkg):**
```powershell
# Install vcpkg first if you haven't: https://vcpkg.io/
vcpkg install sundials zlib

# When building, set CMAKE_TOOLCHAIN_FILE:
# -DCMAKE_TOOLCHAIN_FILE=[vcpkg root]/scripts/buildsystems/vcpkg.cmake
```

**Note:** CMake is typically installed separately as a system tool. On Windows, download from [cmake.org](https://cmake.org/download/) or use `winget install Kitware.CMake`.

### Clone and Setup
```bash
git clone --recursive https://github.com/frank-zillmann/FS3.git
cd FS3
```

**Note:** The `--recursive` flag ensures git submodules (Eigen and cnpy) are cloned.

### Build with Scripts (Recommended)

All scripts are located in the `scripts/` directory and handle the complete build process:

```bash
# Build in Release mode (optimized, with logging and benchmarking)
./scripts/build.sh Release

# Build in Release mode (benchmarking only, no logging)
./scripts/build.sh Release -DLOG_ENABLED=OFF

# Build in Release mode (no logging, no benchmarking - fastest)
./scripts/build.sh Release -DLOG_ENABLED=OFF -DBENCHMARK_ENABLED=OFF

# Build in Debug mode (with logging, no benchmarking - for debugging)
./scripts/build.sh Debug -DBENCHMARK_ENABLED=OFF

# Quick recompile after code changes (no CMake reconfiguration)
./scripts/recompile.sh

# Install to custom location
./scripts/install.sh /path/to/install

# Generate documentation
./scripts/docs.sh generate

# Open documentation in browser
./scripts/docs.sh open

# Other documentation commands
./scripts/docs.sh clean      # Remove generated docs
./scripts/docs.sh rebuild    # Clean and regenerate
./scripts/docs.sh pdf        # Generate PDF (requires LaTeX)
```

### VS Code Tasks

If you're using VS Code, press **`Ctrl + Shift + B`** (or **`Cmd + Shift + B`** on macOS) and select:

- **CMake: Release Build (Logging + Benchmarking)** - Optimized with logging and benchmarking
- **CMake: Release Build (Benchmarking only)** - Optimized with benchmarking, no logging
- **CMake: Release Build (no Logging, no Benchmarking)** - Optimized, fastest execution
- **CMake: Debug Build (Logging only)** - Debug build with AddressSanitizer, for debugging
- **CMake: Recompile** - Quick rebuild without CMake reconfiguration
- **Generate Documentation** - Generate Doxygen documentation
- **Open Documentation** - Open documentation in browser
- **Install FS3 Locally** - Install to `./install` directory

### CMake Build Options

All options can be passed to the build script:
```bash
./scripts/build.sh Release -DOPTION=VALUE
```

**Build Configuration:**
- `-DCMAKE_BUILD_TYPE=Release|Debug` - Build type (default: Release)
- `-DCMAKE_INSTALL_PREFIX=/path` - Installation directory (default: ./install)
- `-DCMAKE_PREFIX_PATH=/path/to/sundials` - Custom SUNDIALS location

**Feature Toggles:**
- `-DLOG_ENABLED=ON|OFF` - Enable logging functionality (default: ON)
- `-DBENCHMARK_ENABLED=ON|OFF` - Enable benchmarking functionality (default: ON)
- `-DENABLE_ASAN=ON|OFF` - Enable AddressSanitizer in Debug builds (default: ON)
- `-DBUILD_DOCS=ON|OFF` - Generate Doxygen documentation (default: ON if Doxygen found)

**Logging Configuration:**
- `-DLOG_FIRST_N_CALLS=1000` - Log first N function calls (default: 1000)
- `-DLOG_EVERY_N_CALLS=1000` - Log every Nth call after initial logging (default: 1000)
- `-DOUTPUT_DIR=/path/to/output` - Custom output directory for logs and observations

**Output Structure:**
When using FS3 in applications, the output directory (default: `./run_yyyy-mm-dd_hh-mm-ss`) contains:
- `logs/` - General logging output
- `obs/` - Observation/simulation data
- `bench/` - Benchmarking output (if enabled)

### Using FS3 in Your Project

After installation, add to your `CMakeLists.txt`:

```cmake
# Find FS3
list(APPEND CMAKE_PREFIX_PATH "/path/to/FS3/install")
find_package(FS3 REQUIRED)

# Link against FS3
add_executable(my_simulation main.cpp)
target_link_libraries(my_simulation PRIVATE FS3::FS3)
```

FS3 automatically provides:
- All required headers (FS3, Eigen, cnpy)
- Required libraries (FS3, cnpy, SUNDIALS, ZLIB)
- Proper include directories and compile definitions

### Development Notes

**For better debugging with Eigen types:**
Install [gdb-eigen-printers](https://github.com/gilleswaeber/gdb-eigen-printers) for pretty-printing Eigen matrices in GDB.

### Logging and Benchmarking Macros

**Logging Macros:**
- `LOG(file, msg)` - Writes to `logs/<file>` in the output directory. Only active when `LOG_ENABLED=ON`.
  ```cpp
  LOG("simulation.log", "Starting simulation at t=" << current_time);
  ```

**Benchmarking Macros:**
- `LOG_BENCHMARK(file, msg)` - Writes to `bench/<file>` in the output directory. Only active when `BENCHMARK_ENABLED=ON`.
  ```cpp
  LOG_BENCHMARK("performance.log", "Iteration " << i << " took " << duration << " ms");
  ```

- `BENCHMARK(var, { code })` - Measures wall time (ms) for code block. The duration is stored in `var`. When `BENCHMARK_ENABLED=OFF`, code still executes but without timing overhead.
  ```cpp
  double elapsed_ms;
  BENCHMARK(elapsed_ms, {
      // Code to benchmark
      solver.solve(system);
  });
  LOG_BENCHMARK("timings.log", "Solver took " << elapsed_ms << " ms");
  ```

**Note:** When disabled via CMake options, these macros compile to no-ops or minimal code, ensuring zero runtime overhead in production builds.

