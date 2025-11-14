file(GLOB_RECURSE SRC_FILES
    "*.cpp"
)
# Exclude example, experiment and test files
list(FILTER SRC_FILES EXCLUDE REGEX ".*/Examples/.*")
list(FILTER SRC_FILES EXCLUDE REGEX ".*/Experiments/.*")
list(FILTER SRC_FILES EXCLUDE REGEX ".*/Tests/.*")

add_library(${PROJECT_NAME} STATIC ${SRC_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(${PROJECT_NAME} PUBLIC
    # ${LAPACK_LIBRARIES}
    SUNDIALS::arkode
    SUNDIALS::cvode
    SUNDIALS::nvecserial
    # SUNDIALS::sunlinsollapackband
    SUNDIALS::sunlinsolband
    # SUNDIALS::sunlinsolklu
    # SUNDIALS::sunmatrixsparse
    Eigen3::Eigen
)

