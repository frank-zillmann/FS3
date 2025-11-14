
cmake_minimum_required(VERSION 3.16)
set(PROJECT_NAME "FS3")
project(${PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message("Building with CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# Option to set LOG_ENABLED preprocessor flag
option(LOG_ENABLED "Enable logging functionality" ON)
message("LOG_ENABLED: ${LOG_ENABLED}")

# Option to set BENCHMARK_ENABLED preprocessor flag (independent of logging)
option(BENCHMARK_ENABLED "Enable benchmarking functionality" ON)
message("BENCHMARK_ENABLED: ${BENCHMARK_ENABLED}")

set(LOG_FIRST_N_CALLS "1000" CACHE STRING "Log first N calls")
message("LOG_FIRST_N_CALLS: ${LOG_FIRST_N_CALLS}")

set(LOG_EVERY_N_CALLS "1000" CACHE STRING "Log every Nth call")
message("LOG_EVERY_N_CALLS: ${LOG_EVERY_N_CALLS}")

# Option to set OUTPUT_DIR preprocessor flag
set(OUTPUT_DIR "" CACHE STRING "Set output directory for logs")
if(OUTPUT_DIR)
    message("OUTPUT_DIR: ${OUTPUT_DIR}")
endif()

# Add -march=native to Release only (keep portability in check!)
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS_RELEASE
      "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  message("CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

# Option to enable AddressSanitizer
option(ENABLE_ASAN "Enable AddressSanitizer" ON)

if(ENABLE_ASAN AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Enable comprehensive AddressSanitizer flags
    set(ASAN_FLAGS "-fsanitize=address -fsanitize=leak -fsanitize=undefined -fno-omit-frame-pointer -fno-optimize-sibling-calls")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${ASAN_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${ASAN_FLAGS}")

    message("CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
    message("CMAKE_EXE_LINKER_FLAGS_DEBUG: ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
endif()

# SUNDIALS
find_package(SUNDIALS REQUIRED COMPONENTS cvode nvecserial sunlinsolband sunlinsolklu sunmatrixsparse)
if (NOT SUNDIALS_FOUND)
    message(FATAL_ERROR "SUNDIALS not found. Please install it or set the SUNDIALS_DIR environment variable.")
else()
    message(STATUS "SUNDIALS found: ${SUNDIALS_LIBRARIES}")
endif()

# EIGEN
add_subdirectory(external/eigen)
link_libraries(Eigen3::Eigen)

# CNPY
add_subdirectory(external/cnpy)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/cnpy)
link_libraries(cnpy)

# Configure preprocessor definitions
if(LOG_ENABLED)
    add_compile_definitions(LOG_ENABLED=1)
else()
    add_compile_definitions(LOG_ENABLED=0)
endif()

if(BENCHMARK_ENABLED)
    add_compile_definitions(BENCHMARK_ENABLED=1)
else()
    add_compile_definitions(BENCHMARK_ENABLED=0)
endif()

add_compile_definitions(LOG_FIRST_N_CALLS=${LOG_FIRST_N_CALLS})
add_compile_definitions(LOG_EVERY_N_CALLS=${LOG_EVERY_N_CALLS})

if(OUTPUT_DIR)
    add_compile_definitions(OUTPUT_DIR="${OUTPUT_DIR}")
endif()

# Add subdirectories
add_subdirectory(src)

# Doxygen documentation
option(BUILD_DOCS "Build documentation" ON)
if(BUILD_DOCS)
    find_package(Doxygen REQUIRED)
    if(DOXYGEN_FOUND)
        # Set project version for Doxygen
        set(PROJECT_VERSION "1.0.0")
        
        # Configure Doxyfile
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in 
                       ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
        
        # Add custom target for documentation
        add_custom_target(docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM
        )
        
        # Add custom target for opening documentation
        add_custom_target(docs-open
            COMMAND xdg-open ${CMAKE_CURRENT_BINARY_DIR}/docs/html/index.html
            DEPENDS docs
            COMMENT "Opening documentation in browser"
        )
        
        message(STATUS "Doxygen found. Use 'make docs' to generate documentation.")
    else()
        message(STATUS "Doxygen not found. Documentation will not be generated.")
    endif()
endif()
