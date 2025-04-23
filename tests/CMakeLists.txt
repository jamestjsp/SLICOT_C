# CMakeLists.txt for SLICOT C wrapper tests

# First, download and configure GoogleTest if not already done
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

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

# Copy test data files to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Include CTest and GoogleTest modules for testing
include(GoogleTest)
enable_testing()

# Discover tests in the executable with more explicit options
gtest_discover_tests(slicot_c_wrapper_tests
    PROPERTIES
        TIMEOUT 120  # Set a timeout for each test
    DISCOVERY_TIMEOUT 120  # Set a timeout for test discovery
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}  # Set the working directory for tests
)
