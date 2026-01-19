# FS³ - Fast and Flexible Framework for Simple Simulations of Separation-Processes
**Created by Frank Zillmann during his Bachelor's Thesis: High-Performance Process Simulations for Magnetic Separation of Biomolecules**

See a [**Simulation of a hIgG Purification Process using the Rotor Stator High Gradient Magnetic Separation**](https://github.com/frank-zillmann/HGMS-hIgG-separation) as an example application of this framework.

## Prerequisites

### Required Dependencies
- **C++ Compiler**: C++20 support
- **CMake**: Version 3.16 or higher
- **SUNDIALS**: with components: cvode, nvecserial, sunlinsolband
- **ZLIB**: Required by cnpy for reading/writing compressed .npz files (NumPy archives)

### Optional Dependencies
- **Doxygen**: For generating documentation (with graphviz for diagrams)
- **LaTeX**: For PDF documentation generation

### Platform-Specific Installation

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install build-essential cmake libsundials-dev zlib1g-dev
```

**macOS (using Homebrew):**
```bash
brew install cmake sundials zlib
```

**Windows (using vcpkg):**
- Install Cmake from https://cmake.org/download/
- Install vcpkg from https://vcpkg.io/
```powershell
vcpkg install sundials zlib

# When building, set CMAKE_TOOLCHAIN_FILE:
# -DCMAKE_TOOLCHAIN_FILE=[vcpkg root]/scripts/buildsystems/vcpkg.cmake
```

**Note:** On macOS/Windows, the C++ compiler is installed separately from package managers (Xcode Command Line Tools for macOS, Visual Studio for Windows).

### Clone and Setup
```bash
git clone https://github.com/frank-zillmann/FS3.git --recursive
cd FS3
```

**Note:** The `--recursive` flag ensures git submodules (Eigen and cnpy) are cloned.

## Build with Scripts or VS Code Tasks

- Use the scripts under `scripts/` directory or the provided VS Code tasks (press **`Ctrl + Shift + B`**)
- `scripts/build.sh` script handles CMake configuration, building, and installation 
- Cmake options can be passed as arguments
- `./scripts/recompile.sh` for quick recompilation after code changes (no CMake reconfiguration)
- Example usages (for these settings VS Code tasks are available):

```bash
# Build in Release mode, no logging, no bechmarking - fastest and default
./scripts/build.sh

# Build in Release mode with logging enabled
./scripts/build.sh -DLOG_ENABLED=ON

# Build in Release mode with benchmarking enabled
./scripts/build.sh -DBENCHMARK_ENABLED=ON

# Build in Debug mode (includes AdressSanitizer) with logging enabled
./scripts/build.sh -DCMAKE_BUILD_TYPE=Debug -DLOG_ENABLED=ON

# Custom install location
./scripts/build.sh -DCMAKE_INSTALL_PREFIX=/custom/install/path

# Quick recompile after code changes (no CMake reconfiguration)
./scripts/recompile.sh
```

## CMake Build Options

All these options can be passed to the build script:
```bash
./scripts/build.sh -DOPTION=VALUE
```

**Build Configuration:**
- `-DCMAKE_BUILD_TYPE=Release|Debug` - Build type (default: Release)
- `-DCMAKE_INSTALL_PREFIX=/path` - Installation directory (default: ./install)
- `-DCMAKE_PREFIX_PATH=/path/to/sundials` - Custom SUNDIALS location

**Feature Toggles:**
- `-DLOG_ENABLED=ON|OFF` - Enable logging functionality (default: OFF)
- `-DBENCHMARK_ENABLED=ON|OFF` - Enable benchmarking functionality (default: OFF)
- `-DBUILD_DOCS=ON|OFF` - Generate Doxygen documentation (default: ON if Doxygen found)

**Logging Configuration:**
- `-DLOG_FIRST_N_CALLS=1000` - Log first N function calls (default: 1000)
- `-DLOG_EVERY_N_CALLS=1000` - Log every Nth call after initial logging (default: 1000)
- `-DOUTPUT_DIR=/path/to/output` - Custom output directory for logs and observations

When using FS3 in applications, the output directory (default: `./run_yyyy-mm-dd_hh-mm-ss`) contains:
- `logs/` - General logging output
- `obs/` - Observation/simulation data
- `bench/` - Benchmarking output (if enabled)

## Build Documentation

```bash
# Generate Doxygen documentation
./scripts/docs.sh generate

# Open documentation in browser
./scripts/docs.sh open

# Remove generated docs
./scripts/docs.sh clean

# Generate PDF (requires LaTeX)
./scripts/docs.sh pdf
```

## Using FS3 in C++

After installation, add to your `CMakeLists.txt`:

```cmake
# Find FS3
list(APPEND CMAKE_PREFIX_PATH "/path/to/FS3/install")
find_package(FS3 REQUIRED)

# Link against FS3
add_executable(my_simulation main.cpp)
target_link_libraries(my_simulation PRIVATE FS3::FS3)
```

**Logging Macros:**
- `LOG(file, msg)` - Writes to `logs/<file>` in the output directory. Only active when `LOG_ENABLED=ON`.
  ```cpp
  LOG("example.log", "Log message");
  ```

**Benchmarking Macros:**
- `LOG_BENCHMARK(file, msg)` - Writes to `bench/<file>` in the output directory. Only active when `BENCHMARK_ENABLED=ON`.
  ```cpp
  LOG_BENCHMARK("benchmark.log", "Benchmark message");
  ```

- `BENCHMARK(var, { code })` - Measures time in seconds for code block. The duration is stored in `var`. When `BENCHMARK_ENABLED=OFF`, code still executes but without timing overhead.
  ```cpp
  double t;default_components=python 
  BENCHMARK(t, {
      // Code to benchmark
  });
  LOG_BENCHMARK("benchmark.log", "Code took " << t << " s");
  ```

Note: For better debugging with Eigen, I recommend to install [gdb-eigen-printers](https://github.com/gilleswaeber/gdb-eigen-printers).

## Using FS3 in Python

**Prerequisites:**
- Python 3.8+ with development headers
- numpy (installed automatically as dependency)

**Option 1: Build wheel (recommended for distribution)**
```bash
pip install scikit-build-core nanobind
pip wheel . --no-deps -w install/dist/
pip install install/dist/fs3-*.whl
```

**Option 2: Direct install to conda/virtualenv**
```bash
# Builds and installs Python module to site-packages
./scripts/build.sh -DBUILD_PYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX # for conda
./scripts/build.sh -DBUILD_PYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV # for virtualenv
```

**Verify installation:**
```python
import fs3
print(fs3.Component("H+", charge=1, molar_mass_kg_per_mol=0.001))
```

**Note:** When `BUILD_PYTHON_BINDINGS=ON`, the C++ targets are not installed by default to avoid pollution of site-packages with C++ files.
Python users get only:
```
site-packages/fs3/
  fs3.cpython-<version>-<platform>.so  # Compiled extension
  fs3.pyi                               # Type hints for IDEs
```

