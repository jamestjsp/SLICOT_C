# CMakeLists.txt for SLICOT C wrapper tests

# GoogleTest is already configured in the main CMakeLists.txt
# No need to redeclare it here

# Find all test source files
file(GLOB TEST_SOURCES "*.cpp")

if(NOT TEST_SOURCES)
    message(STATUS "No test source files found in tests directory.")
    return()
endif()

message(STATUS "Found test source files: ${TEST_SOURCES}")

# Add the test executable
add_executable(slicot_c_wrapper_tests ${TEST_SOURCES})

# Include directories for the test executable
target_include_directories(slicot_c_wrapper_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/include  # For SLICOT headers
)

# Define macros for export/import or static build based on SLICOT_BUILD_SHARED_LIBS
if(NOT SLICOT_BUILD_SHARED_LIBS)
    # Define SLICOT_STATIC when linking against a static library
    target_compile_definitions(slicot_c_wrapper_tests PRIVATE SLICOT_STATIC)
    message(STATUS "Test build: Using static linking (SLICOT_STATIC defined)")
else()
    message(STATUS "Test build: Using dynamic linking")
endif()

# Link against the main SLICOT library (which includes the C wrappers),
# GoogleTest libraries, and BLAS/LAPACK
target_link_libraries(slicot_c_wrapper_tests PRIVATE
    slicot          # Main library including Fortran code and C wrappers
    GTest::gtest    # GoogleTest main library
    GTest::gmock    # GoogleMock (often used with gtest)
    GTest::gtest_main # Provides main() function
    ${SLICOT_LINK_TARGETS} # Link BLAS/LAPACK found in parent scope
)

# Add test_utils.cpp to the test sources
target_sources(slicot_c_wrapper_tests
    PRIVATE
        test_utils.cpp
)

# Handle test data in a more elegant way
# Define a preprocessor macro that points to the test data directory
# This allows tests to find data files regardless of working directory
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/test_config.h
)

# Add the binary directory to include paths to find the generated config
target_include_directories(slicot_c_wrapper_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/include  # For SLICOT headers
    ${CMAKE_CURRENT_BINARY_DIR}  # For generated test_config.h
)

# No need to copy data if using absolute paths from test_config.h
# But keep it for backward compatibility or just in case
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Include CTest and GoogleTest modules for testing
include(GoogleTest)
enable_testing()

# Make sure the test working directory is set correctly
gtest_discover_tests(slicot_c_wrapper_tests
    PROPERTIES
        TIMEOUT 120  # Set a timeout for each test
    DISCOVERY_TIMEOUT 120  # Set a timeout for test discovery
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}  # Set the working directory for tests
)
