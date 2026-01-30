cmake_minimum_required(VERSION 3.16)

set(PROJECT_NAME "FS3")
set(PROJECT_VERSION "0.0.1")
set(PROJECT_URL "https://github.com/frank-zillmann/FS3")
set(
  PROJECT_DESCRIPTION
  "Fast and Flexible Framework for Simple Simulations of Separation-Processes"
)

project(${PROJECT_NAME}
    VERSION ${PROJECT_VERSION}
    DESCRIPTION ${PROJECT_DESCRIPTION}
    HOMEPAGE_URL ${PROJECT_URL}
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Building with CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# Option to set LOG_ENABLED preprocessor flag
option(LOG_ENABLED "Enable logging functionality" OFF)
message(STATUS "LOG_ENABLED: ${LOG_ENABLED}")

# Option to set BENCHMARK_ENABLED preprocessor flag (independent of logging)
option(BENCHMARK_ENABLED "Enable benchmarking functionality" OFF)
message(STATUS "BENCHMARK_ENABLED: ${BENCHMARK_ENABLED}")

set(LOG_FIRST_N_CALLS "1000" CACHE STRING "Log first N calls")
message(STATUS "LOG_FIRST_N_CALLS: ${LOG_FIRST_N_CALLS}")

set(LOG_EVERY_N_CALLS "1000" CACHE STRING "Log every Nth call")
message(STATUS "LOG_EVERY_N_CALLS: ${LOG_EVERY_N_CALLS}")

# Option to set OUTPUT_DIR preprocessor flag
set(OUTPUT_DIR "" CACHE STRING "Set output directory for logs")
if(OUTPUT_DIR)
    message(STATUS "OUTPUT_DIR: ${OUTPUT_DIR}")
endif()

# Add -march=native to Release only (keep portability in check!)
# Note: Disable for portable wheels with -DUSE_NATIVE_ARCH=OFF
option(USE_NATIVE_ARCH "Use -march=native for Release builds" ON)
if(CMAKE_BUILD_TYPE STREQUAL "Release" AND USE_NATIVE_ARCH)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

# ==================== Dependencies ====================
include(FetchContent)

# SUNDIALS - System dependency (too complex to build inline)
find_package(SUNDIALS REQUIRED COMPONENTS cvode nvecserial sunlinsolband)
message(STATUS "SUNDIALS found: ${SUNDIALS_VERSION}")

# ZLIB - System dependency (required by cnpy)
find_package(ZLIB REQUIRED)
message(STATUS "ZLIB found: ${ZLIB_VERSION_STRING}")

# EIGEN - Header-only library, always fetched via FetchContent
message(STATUS "Fetching Eigen via FetchContent...")
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 5.0.1
    GIT_SHALLOW TRUE
)
# Just populate, don't add_subdirectory (avoids export/install issues)
FetchContent_Populate(eigen)
set(EIGEN_INCLUDE_DIR "${eigen_SOURCE_DIR}")
message(STATUS "Eigen fetched to: ${EIGEN_INCLUDE_DIR}")

# CNPY - Small library for NumPy file I/O (compiled directly into FS3)
message(STATUS "Fetching cnpy via FetchContent...")
FetchContent_Declare(
    cnpy
    GIT_REPOSITORY https://github.com/rogersce/cnpy.git
    GIT_TAG master
    GIT_SHALLOW TRUE
)
# Just populate, don't add_subdirectory, compile cnpy.cpp directly into FS3
FetchContent_Populate(cnpy)

# ==================== Build FS3 Library ====================
file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)
# Exclude example, experiment and test files (currently not needed)
# list(FILTER SRC_FILES EXCLUDE REGEX ".*/(Examples|Experiments|Tests)/.*")

# Add cnpy source directly to FS3 (avoids export issues)
list(APPEND SRC_FILES "${cnpy_SOURCE_DIR}/cnpy.cpp")

add_library(${PROJECT_NAME} STATIC ${SRC_FILES})

# Sanitizer flags: apply after target creation, use PUBLIC for propagation to applications
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(ASAN_FLAGS
        -fsanitize=address
        -fsanitize=leak
        -fsanitize=undefined
        -fno-omit-frame-pointer
        -fno-optimize-sibling-calls
    )
    target_compile_options(${PROJECT_NAME} PUBLIC $<$<CONFIG:Debug>:${ASAN_FLAGS}>)
    target_link_options(${PROJECT_NAME} PUBLIC $<$<CONFIG:Debug>:${ASAN_FLAGS}>)
    message(STATUS "Exporting sanitizer flags for Debug builds (PUBLIC): ${ASAN_FLAGS}")
endif()

target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        # cnpy headers only needed internally
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/cnpy>
        $<BUILD_INTERFACE:${cnpy_SOURCE_DIR}>
)

# Add Eigen include path (header-only, fetched via FetchContent)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC
    $<BUILD_INTERFACE:${EIGEN_INCLUDE_DIR}>
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
        # These dependencies are part of FS3's public API
        SUNDIALS::arkode
        SUNDIALS::cvode
        SUNDIALS::nvecserial
        SUNDIALS::sunlinsolband
    PRIVATE
        # ZLIB is required by cnpy (compiled directly into FS3)
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
    target_compile_definitions(${PROJECT_NAME} PUBLIC OUTPUT_DIR=${OUTPUT_DIR})
endif()

# ==================== Doxygen Documentation (Optional) ====================
option(BUILD_DOCS "Build documentation" OFF)
if(BUILD_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)        
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in 
                       ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
        
        add_custom_target(docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM
        )
        
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

# ==================== Python Bindings (Optional) ====================
option(BUILD_PYTHON_BINDINGS "Build Python bindings using nanobind" OFF)
if(BUILD_PYTHON_BINDINGS)
    message(STATUS "Python bindings enabled")

    # Enable PIC (needed when linking into Python module)
    set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    
    find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
    
    # Fetch nanobind
    FetchContent_Declare(
        nanobind
        GIT_REPOSITORY https://github.com/wjakob/nanobind.git
        GIT_TAG v2.10.2
        GIT_SHALLOW TRUE
    )
    message(STATUS "Fetching nanobind (this may take a moment on first build)...")
    FetchContent_MakeAvailable(nanobind)

    # Create the Python module
    nanobind_add_module(fs3 python/bindings.cpp)
    target_link_libraries(fs3 PRIVATE ${PROJECT_NAME})
    
    # Install the Python module
    if(DEFINED SKBUILD)
        # Building a wheel with scikit-build-core: install to root so 'import fs3' works
        install(TARGETS fs3 LIBRARY DESTINATION . COMPONENT python)
        install(FILES python/fs3.pyi DESTINATION . COMPONENT python)
    else()
        # Direct CMake install: detect Python version and install to lib/pythonX.Y/site-packages
        execute_process(
            COMMAND ${Python_EXECUTABLE} -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
            OUTPUT_VARIABLE PYTHON_VERSION_STRING
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set(PYTHON_SITE_PACKAGES "lib/python${PYTHON_VERSION_STRING}/site-packages")
        install(TARGETS fs3 LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT python)
        install(FILES python/fs3.pyi DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT python)
        message(STATUS "Python module will be installed to: ${CMAKE_INSTALL_PREFIX}/${PYTHON_SITE_PACKAGES}")
    endif()
    
    message(STATUS "Python module 'fs3' will be built")

# ==================== C++ Installation ====================
else()
    include(GNUInstallDirs)

    # Install FS3 library only (cnpy is statically linked and not exported)
    install(
        TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    # Install public headers
    install(
        DIRECTORY ${PROJECT_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT Development
    )

    # Export target set for find_package
    install(
        EXPORT ${PROJECT_NAME}Targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
        COMPONENT Development
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
        COMPONENT Development
    )
endif()
