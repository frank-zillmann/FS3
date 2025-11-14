# FS³ - Fast and Flexible Framework for Simple Simulations of Separation-Processes
**Created by Frank Zillmann during his Bachelor's Thesis: High-Performance Process Simulations for Magnetic Separation of Biomolecules**

See a [**Simulation of a hIgG Purification Process using the Rotor Stator High Gradient Magnetic Separation**](https://github.com/frank-zillmann/HGMS-hIgG-separation) as an example application of this framework.

## Setup and Build Instructions

### Prerequisites
TODO: Coming soon

### Clone
Clone the repository using `git clone <repo-url>`

Then open the folder in VSCode (recommended) or use `cd <repo-name>` and continue with the scripts in command line (alternative).

### VS Code Tasks (recommended)
Press `Ctrl + Shift + B` and select:
- `CMake: Release Build (with Logging)` - Clean Release build with logging enabled (fast)
- `CMake: Release Build (no Logging)` - Clean Release build with logging disabled (fastest)
- `CMake: Complete Debug Build` - Clean build for debugging (slow)
- `CMake: Recompile` - Recompile without new CMake configuration (enough, if you made small code changes, but did not change CMake Settings or added/removed files)
- `Run: HGMS_hIgG_separation` - Run the simulation of the HGMS hIgG separation process
- `Generate Documentation` - Generate Doxygen documentation
- `Open Documentation` - Open documentation in browser

### Scripts (alternative e.g. not using VS Code)
```bash
# Build commands
./scripts/build.sh Release
./scripts/build.sh Release -DLOG_ENABLED=OFF
./scripts/build.sh Release -DBENCHMARK_ENABLED=OFF
./scripts/build.sh Debug
./scripts/recompile.sh

# Run simulation
./scripts/run-hgms.sh

# Documentation
./scripts/docs.sh generate
./scripts/docs.sh open
./scripts/docs.sh clean
./scripts/docs.sh rebuild
./scripts/docs.sh pdf
```

### CMake Build Options

- **Build Type**: `-DCMAKE_BUILD_TYPE=Release|Debug` (default: Release)
- **Custom Sundials**: `-DCMAKE_PREFIX_PATH=/path/to/sundials`
- **Logging**: `-DLOG_ENABLED=ON|OFF` (default: ON) - Enable/disable logging functionality
- **Benchmarking**: `-DBENCHMARK_ENABLED=ON|OFF` (default: ON) - Enable/disable benchmarking functionality (independent of logging)
- **Output Directory**: `-DOUTPUT_DIR=/path/to/output` (default: ./run_yyyy-yy-dd_hh-mm-ss) - Set custom output directory for outputs (logs + observations)
        - The run directory contains three subfolders: `logs/`, `obs/`, and `bench/` (for benchmarking output)
- **Log First N Calls**: `-DLOG_FIRST_N_CALLS=1000` (default: 1000) - Log first N function calls
- **Log Every N Calls**: `-DLOG_EVERY_N_CALLS=1000` (default: 1000) - Log every Nth call after initial logging
- **Build Documentation**: `-DBUILD_DOCS=ON|OFF` (default: ON) - Generate Doxygen documentation
- **AddressSanitizer**: `-DENABLE_ASAN=ON|OFF` (default: ON for Debug builds) - Enable memory error detection

### Notes

- For better debugging with the Eigen Types I can advise to install [gdb-eigen-printers](https://github.com/gilleswaeber/gdb-eigen-printers)

### Logging and Benchmarking Macros

- `LOG(file, msg)`: Writes to the `logs/` subfolder inside the timestamped run directory. Controlled by `LOG_ENABLED`.
- `LOG_BENCHMARK(file, msg)`: Same usage as `LOG`, but writes to the `bench/` subfolder. Controlled by `BENCHMARK_ENABLED`.
- `BENCHMARK(var, { code })`: Measures wall time (ms) for the enclosed code block, storing the duration in `var`. Controlled by `BENCHMARK_ENABLED`. When disabled, the code still executes, but no timing overhead is added.

