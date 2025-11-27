
cmake_minimum_required(VERSION 3.16)

set(PROJECT_NAME "FS3")
project(${PROJECT_NAME} LANGUAGES CXX)
set(PROJECT_URL "https://github.com/frank-zillmann/FS3")
set(
  PROJECT_DESCRIPTION
  "Fast and Flexible Framework for Simple Simulations of Separation-Processes"
)

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
find_package(SUNDIALS REQUIRED COMPONENTS cvode nvecserial sunlinsolband)
if (NOT SUNDIALS_FOUND)
    message(FATAL_ERROR "SUNDIALS not found. Please install it or set the SUNDIALS_DIR environment variable.")
else()
    message(STATUS "SUNDIALS found: ${SUNDIALS_LIBRARIES}")
endif()

# ZLIB (required by cnpy)
find_package(ZLIB REQUIRED)
if (NOT ZLIB_FOUND)
    message(FATAL_ERROR "ZLIB not found. Please install it.")
endif()

# EIGEN
add_subdirectory(external/eigen)

# CNPY - EXCLUDE_FROM_ALL prevents installing cnpy's binaries/headers
add_subdirectory(external/cnpy EXCLUDE_FROM_ALL)

# Build FS3
file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)
# Exclude example, experiment and test files
list(FILTER SRC_FILES EXCLUDE REGEX ".*/(Examples|Experiments|Tests)/.*")

add_library(${PROJECT_NAME} STATIC ${SRC_FILES})

target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/cnpy>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
        # These dependencies are part of FS3's public API
        SUNDIALS::arkode
        SUNDIALS::cvode
        SUNDIALS::nvecserial
        SUNDIALS::sunlinsolband
        Eigen3::Eigen
        cnpy-static
        ZLIB::ZLIB
)

# Configure preprocessor definitions
if(LOG_ENABLED)
    target_compile_definitions(${PROJECT_NAME} PUBLIC LOG_ENABLED=1)
else()
    target_compile_definitions(${PROJECT_NAME} PUBLIC LOG_ENABLED=0)
endif()

if(BENCHMARK_ENABLED)
    target_compile_definitions(${PROJECT_NAME} PUBLIC BENCHMARK_ENABLED=1)
else()
    target_compile_definitions(${PROJECT_NAME} PUBLIC BENCHMARK_ENABLED=0)
endif()

target_compile_definitions(${PROJECT_NAME} PUBLIC LOG_FIRST_N_CALLS=${LOG_FIRST_N_CALLS})
target_compile_definitions(${PROJECT_NAME} PUBLIC LOG_EVERY_N_CALLS=${LOG_EVERY_N_CALLS})

if(OUTPUT_DIR)
    target_compile_definitions(${PROJECT_NAME} PUBLIC OUTPUT_DIR="${OUTPUT_DIR}")
endif()

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

# Install / package support
include(GNUInstallDirs)

# Install all targets to satisfy CMake export requirements
# cnpy-static is included in export because FS3 links to it, but since it's PRIVATE
# and statically linked, consumers don't directly use it
install(
    TARGETS ${PROJECT_NAME} eigen cnpy-static
    EXPORT ${PROJECT_NAME}Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install public headers
install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install cnpy headers
install(
    FILES ${CMAKE_CURRENT_SOURCE_DIR}/external/cnpy/cnpy.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Export target set for find_package
install(
    EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

# Generate config + version files
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)