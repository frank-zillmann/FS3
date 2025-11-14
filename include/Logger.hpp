#ifndef LOGGER_HPP
#define LOGGER_HPP

// Requires C++17 for <filesystem>
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>

// ––––––––––––––––––
// Configuration Flags
// ––––––––––––––––––
// Enable/disable logging & benchmarking (default: enabled)
#ifndef LOG_ENABLED
#define LOG_ENABLED 1
#endif

// Enable/disable benchmarking (default: enabled)
#ifndef BENCHMARK_ENABLED
#define BENCHMARK_ENABLED 1
#endif

// Optional: compile-time top-level output directory.
// Example: g++ -DOUTPUT_DIR="\"/home/myuser\"" ...
#ifdef OUTPUT_DIR
// nothing else required here; used in init()
#endif

#ifndef LOG_FIRST_N_CALLS
#define LOG_FIRST_N_CALLS 1000  // log first 1000 calls
#endif

#ifndef LOG_EVERY_N_CALLS
#define LOG_EVERY_N_CALLS 1000  // log every 1000th call
#endif

namespace logger {

// Forward declarations
inline void register_cleanup_handlers();

// resolved log folder (OUTPUT_DIR/run_TIMESTAMP/logs or ./run_TIMESTAMP/logs)
inline std::string log_folder;
inline std::string obs_folder;
inline std::string bench_folder;
inline std::string run_folder;  // top-level timestamped run folder (parent of logs/ obs/ bench/)

// protect one-time init & run_folder changes
inline std::mutex init_mutex;

// protect log_streams map and writes
inline std::mutex log_streams_mutex;
inline std::unordered_map<std::string, std::shared_ptr<std::ofstream>> log_streams;

// single write mutex used to make each LOG line atomic (coarse but simple)
inline std::mutex write_mutex;

// create timestamped folder name: run_DD-MM-YYYY_HH-MM
inline std::string make_timestamped_folder_name() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream name;
    name << "run_" << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    return name.str();
}

// Resolve output directory. If OUTPUT_DIR is defined at compile-time, the
// log folder will be OUTPUT_DIR/<timestamp>/logs. Otherwise it's ./<timestamp>/logs.
inline void init_folders() {
    std::lock_guard<std::mutex> lock(init_mutex);
    if (!log_folder.empty() && !obs_folder.empty() && !bench_folder.empty() && !run_folder.empty())
        return;  // already inited

#ifdef OUTPUT_DIR
    std::filesystem::path top = OUTPUT_DIR;
#else
    std::filesystem::path top = std::filesystem::current_path();
#endif

    std::string time_stamped_folder_name = make_timestamped_folder_name();

    // set run folder and create subfolders
    std::filesystem::path run = top / time_stamped_folder_name;
    run_folder = run.string();

    // create timestamped log folder inside run
    std::filesystem::path path_logs = run / "logs";
    std::filesystem::create_directories(path_logs);
    log_folder = path_logs.string();

    // create timestamped obs folder inside run
    std::filesystem::path path_obs = run / "obs";
    std::filesystem::create_directories(path_obs);
    obs_folder = path_obs.string();

    // create timestamped benchmark folder inside run
    std::filesystem::path path_bench = run / "bench";
    std::filesystem::create_directories(path_bench);
    bench_folder = path_bench.string();
}

// Return the current log directory (creates it if needed)
inline std::string log_dir() {
    init_folders();
    return log_folder;
}

// Return the current obs directory (creates it if needed)
inline std::string obs_dir() {
    init_folders();
    return obs_folder;
}

// Return the current benchmark directory (creates it if needed)
inline std::string bench_dir() {
    init_folders();
    return bench_folder;
}

// Ensure a stream (topic filename relative to run_folder) is open and cached.
// topic examples: "log.txt", "reactions.log", "physics/solver.log"
inline std::shared_ptr<std::ofstream> ensure_log_stream(const std::string& relpath) {
    std::lock_guard<std::mutex> lock(log_streams_mutex);
    // key is the complete relative path under run_folder (e.g. "logs/file.log", "bench/file.log")
    const std::string& key = relpath;
    auto it = log_streams.find(key);
    if (it != log_streams.end()) return it->second;

    // ensure run folder exists and cleanup handlers are registered
    init_folders();
    register_cleanup_handlers();

    // Place all logs relative to run_folder
    std::filesystem::path p = std::filesystem::path(run_folder) / key;
    if (p.has_parent_path()) {
        std::filesystem::create_directories(p.parent_path());
    }

    auto ofs = std::make_shared<std::ofstream>(p.string(), std::ios::app | std::ios::binary);
    if (!ofs->is_open()) {
        // fallback to run_folder/log.txt
        std::filesystem::path fallback = std::filesystem::path(run_folder) / "unnamed.log";
        ofs->open(fallback.string(), std::ios::app | std::ios::binary);
    }

    log_streams.emplace(key, ofs);
    return ofs;
}

// convenience accessor returning reference (creates stream if missing)
inline std::ofstream& log_file(const std::string& relpath) { return *ensure_log_stream(relpath); }

// No separate ensure function for benchmark; both logs and benchmarks share the same cache.

// flush all open log_streams to disk
inline void flush_all_logs() {
    std::lock_guard<std::mutex> lock(log_streams_mutex);
    for (auto& kv : log_streams) {
        if (kv.second && kv.second->is_open()) kv.second->flush();
    }
}

// close all log_streams (useful on shutdown)
inline void close_all_logs() {
    std::lock_guard<std::mutex> lock(log_streams_mutex);
    for (auto& kv : log_streams) {
        if (kv.second && kv.second->is_open()) {
            kv.second->close();
        }
    }
    log_streams.clear();
}

// flush all open benchmark streams to disk
// No separate flushing/closing for benchmark streams; they share the same cache.

// Register cleanup handlers to flush/close logs on exit
inline void register_cleanup_handlers() {
    static bool registered = false;
    if (!registered) {
        // Register normal exit handler
        std::atexit([]() {
            flush_all_logs();
            close_all_logs();
        });

        // Register signal handlers for abnormal termination
        // SIGABRT is triggered by assert() failures
        auto signal_handler = [](int signal) {
            flush_all_logs();
            close_all_logs();
            // Re-raise the signal to allow default handling
            std::signal(signal, SIG_DFL);
            std::raise(signal);
        };

        std::signal(SIGABRT, signal_handler);  // assert() failures
        std::signal(SIGTERM, signal_handler);  // termination request
        std::signal(SIGINT, signal_handler);   // Ctrl+C
        std::signal(SIGSEGV, signal_handler);  // segmentation fault (may not always work)

        registered = true;
    }
    // Also ensure folders are initialized when this is called
    init_folders();
}

}  // namespace logger

// ––––––––––––––––––
// Macros
// ––––––––––––––––––

// LOG(file, msg)
// file  : relative filename inside the run folder (can include subdirs)
// msg   : <<-style stream expression (e.g. "Started id=" << id)
#if LOG_ENABLED
#define LOG(file, msg)                                                 \
    do {                                                               \
        std::ostringstream _oss;                                       \
        _oss << msg;                                                   \
        std::lock_guard<std::mutex> _lg(logger::write_mutex);          \
        logger::log_file(std::string("logs/") + (file)) << _oss.str(); \
    } while (0)
#else
#define LOG(file, msg) \
    do {               \
    } while (0)
#endif

// LOG_BENCHMARK(file, msg)
// Identical to LOG except for activation macro and expected path prefix (e.g. "bench/...")
#if BENCHMARK_ENABLED
#define LOG_BENCHMARK(file, msg)                                        \
    do {                                                                \
        std::ostringstream _oss;                                        \
        _oss << msg;                                                    \
        std::lock_guard<std::mutex> _lg(logger::write_mutex);           \
        logger::log_file(std::string("bench/") + (file)) << _oss.str(); \
    } while (0)
#else
#define LOG_BENCHMARK(file, msg) \
    do {                         \
    } while (0)
#endif

// BENCHMARK(var, { code })
// - Declares `realtype var` and runs `code`
// - If BENCHMARK_ENABLED == 1: var is set to elapsed time in seconds (s) as a realtype
// - If BENCHMARK_ENABLED == 0: code executes but no timing machinery is used
// Usage:
//   BENCHMARK(dt_s, { /* work to measure */ });
//   LOG_BENCHMARK("timings.log", "rhs time (s): " << dt_s);
#if BENCHMARK_ENABLED
#define BENCHMARK(var, code)                                                      \
    do {                                                                          \
        auto _bench_start = std::chrono::high_resolution_clock::now();            \
        code;                                                                     \
        auto _bench_end = std::chrono::high_resolution_clock::now();              \
        var = std::chrono::duration<realtype>(_bench_end - _bench_start).count(); \
    } while (0)
#else
#define BENCHMARK(var, code) \
    do {                     \
        code;                \
    } while (0)
#endif

#endif  // LOGGER_HPP
