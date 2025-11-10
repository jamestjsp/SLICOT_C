cmake_minimum_required(VERSION 3.15)

# ---
# NOTE: The SLICOT C wrapper generator always places Fortran CHARACTER string length arguments at the end of the argument list.
# This is portable for all compilers except Intel Fortran (ifort/ifx), which by default passes string lengths immediately after each string.
# The CMake configuration below sets the Intel Fortran flag -assume nomixed_str_len_arg to ensure string lengths are passed at the end,
# making the wrappers portable across compilers. Do not change this unless you know what you are doing.
# ---

# Project definition
project(SLICOT VERSION 5.9.0 LANGUAGES Fortran C CXX) # Add CXX language
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Options
option(SLICOT_BUILD_SHARED_LIBS "Build shared libraries" ON) # Default to shared for DLLs
option(SLICOT_BUILD_EXAMPLES "Build the SLICOT example programs and tests" OFF)
option(SLICOT_BUILD_TESTS "Build the SLICOT C wrapper tests using GoogleTest" ON) # Default to ON
option(SLICOT_LINK_AUX_WRAPPERS "Link auxiliary LAPACK wrappers (for incomplete LAPACK)" ON) # Add option for aux wrappers
# For vcpkg integration, use CMakePresets.json with CMAKE_TOOLCHAIN_FILE set to $env{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake
# See: https://learn.microsoft.com/en-us/vcpkg/get_started/get-started
option(SLICOT_USE_VCPKG "Use vcpkg for BLAS and LAPACK dependencies" OFF) # Option to use vcpkg
option(SLICOT_BUILD_C_WRAPPERS "Build C wrapper library for SLICOT" ON) # Option to build C wrappers

# Setup build properties
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set library type based on option
if(SLICOT_BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
else()
    set(BUILD_SHARED_LIBS OFF)
endif()

# Include standard modules
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)
include(FetchContent) # Include FetchContent module

# Include SLICOT helper functions
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SLICOTHelpers.cmake)

# Enable testing globally
enable_testing()

# --- GoogleTest ---
if(SLICOT_BUILD_TESTS)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip # Or specify a commit hash
        # GIT_REPOSITORY https://github.com/google/googletest.git
        # GIT_TAG release-1.14.0 # Or specify a commit hash
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

# Configure BLAS/LAPACK detection using helper function
slicot_configure_blas_lapack()
message(STATUS "Searching for BLAS/LAPACK vendor: ${BLA_VENDOR}")

# Configure paths for BLAS/LAPACK discovery
if(SLICOT_USE_VCPKG)
  message(STATUS "Using vcpkg for BLAS/LAPACK dependencies")
  # Check if vcpkg toolchain is already loaded (via CMakePresets.json or command line)
  if(DEFINED CMAKE_TOOLCHAIN_FILE AND CMAKE_TOOLCHAIN_FILE MATCHES "vcpkg.cmake")
    message(STATUS "vcpkg toolchain file already configured: ${CMAKE_TOOLCHAIN_FILE}")
  elseif(DEFINED ENV{VCPKG_ROOT})
    # Set toolchain file if not already set
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "vcpkg toolchain file")
    message(STATUS "Setting vcpkg toolchain file: ${CMAKE_TOOLCHAIN_FILE}")
  else()
    message(WARNING "VCPKG_ROOT environment variable not defined. Set CMAKE_TOOLCHAIN_FILE manually or use CMakePresets.json")
  endif()
else()
  message(STATUS "Using system-installed or user-specified BLAS/LAPACK")
  # Allow user to specify custom BLAS/LAPACK paths
  if(DEFINED BLAS_LAPACK_DIR)
    set(CMAKE_PREFIX_PATH "${BLAS_LAPACK_DIR}" ${CMAKE_PREFIX_PATH})
    message(STATUS "Adding user-specified BLAS/LAPACK directory to CMAKE_PREFIX_PATH: ${BLAS_LAPACK_DIR}")
  endif()
    
  # Special handling for Intel MKL - Only for Windows and Linux, not for macOS
  if(BLA_VENDOR MATCHES "Intel" AND DEFINED ENV{MKLROOT} AND NOT APPLE)
    message(STATUS "Intel MKL requested and found at: $ENV{MKLROOT}")
    # For Intel compilers with MKL, directly set the libraries to avoid issues with FindBLAS
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" OR CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
      if(WIN32)
        # Check multiple possible paths for the MKL library
        set(MKL_PATHS 
            "$ENV{MKLROOT}\\lib\\intel64\\mkl_rt.lib"
            "$ENV{MKLROOT}\\lib\\mkl_rt.lib"
            "$ENV{MKLROOT}\\lib\\intel64_win\\mkl_rt.lib")
        
        set(MKL_LIB_PATH "")
        foreach(path ${MKL_PATHS})
          if(EXISTS "${path}")
            set(MKL_LIB_PATH "${path}")
            break()
          endif()
        endforeach()
        
        # If MKL found, use it, otherwise fall back to vcpkg
        if(MKL_LIB_PATH)
          set(BLAS_LIBRARIES "${MKL_LIB_PATH}")
          set(LAPACK_LIBRARIES "${MKL_LIB_PATH}")
          message(STATUS "Using Intel MKL on Windows: ${BLAS_LIBRARIES}")
        else()
          # Fallback to using vcpkg if MKL library is not found
          message(WARNING "Intel MKL library not found in any standard locations. Falling back to vcpkg libraries.")
          set(SLICOT_USE_VCPKG ON)
        endif()
      else()
        # Linux paths
        set(BLAS_LIBRARIES "$ENV{MKLROOT}/lib/libmkl_rt.so")
        set(LAPACK_LIBRARIES "$ENV{MKLROOT}/lib/libmkl_rt.so")
        message(STATUS "Using Intel MKL on Linux: ${BLAS_LIBRARIES}")
      endif()
      # Skip the standard find_package for BLAS/LAPACK since we've set the libraries manually
      if(NOT SLICOT_USE_VCPKG)
        set(BLAS_FOUND TRUE)
        set(LAPACK_FOUND TRUE)
      endif()
    else()
      # For other compilers, just add the MKL path to the search path
      set(CMAKE_PREFIX_PATH "$ENV{MKLROOT}" ${CMAKE_PREFIX_PATH})
    endif()
  # Common locations for Intel oneMKL (when not using Intel compilers or not specifically requesting MKL)
  # Skip MKL on macOS as we only support Apple Silicon
  elseif(DEFINED ENV{MKLROOT} AND NOT APPLE)
    set(CMAKE_PREFIX_PATH "$ENV{MKLROOT}" ${CMAKE_PREFIX_PATH})
    message(STATUS "Found Intel MKL installation at: $ENV{MKLROOT}")
  endif()

  # Force Apple Accelerate framework on macOS
  if(APPLE)
    message(STATUS "On macOS, using Apple's Accelerate framework")
    set(BLA_VENDOR "Apple")
  endif()

  # Common locations for AMD AOCL
  if(DEFINED ENV{AOCL_ROOT})
    set(CMAKE_PREFIX_PATH "$ENV{AOCL_ROOT}" ${CMAKE_PREFIX_PATH})
    message(STATUS "Found AMD AOCL installation at: $ENV{AOCL_ROOT}")
  endif()
endif()

# Configure MKL to use LP64 interface with sequential threading
# This prevents thread contention and ensures 32-bit integer compatibility
if(BLA_VENDOR MATCHES "Intel")
    # Set MKL environment variables for consistent behavior
    set(ENV{MKL_INTERFACE_LAYER} "LP64")
    set(ENV{MKL_THREADING_LAYER} "SEQUENTIAL")
    message(STATUS "MKL configured: LP64 interface, SEQUENTIAL threading")
endif()

# Control BLAS/LAPACK library type to match SLICOT build type
if(SLICOT_BUILD_SHARED_LIBS)
    set(BLA_STATIC OFF)
    message(STATUS "Searching for shared BLAS/LAPACK libraries")
else()
    set(BLA_STATIC ON)
    message(STATUS "Searching for static BLAS/LAPACK libraries")
endif()

# Try to find BLAS and LAPACK if not already manually set
if(NOT BLAS_FOUND)
  find_package(BLAS REQUIRED)
endif()

if(NOT LAPACK_FOUND)
  find_package(LAPACK REQUIRED)
endif()

# Validate BLAS/LAPACK use 32-bit integers (LP64 interface)
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.22")
    if(DEFINED BLA_SIZEOF_INTEGER)
        if(NOT BLA_SIZEOF_INTEGER EQUAL 4)
            message(FATAL_ERROR
                "BLAS/LAPACK use ${BLA_SIZEOF_INTEGER}-byte integers, but SLICOT requires "
                "4-byte integers (LP64 interface). Please use LP64 BLAS/LAPACK libraries, "
                "not ILP64 variants. For Intel MKL, use mkl_intel_lp64 not mkl_intel_ilp64.")
        else()
            message(STATUS "BLAS/LAPACK integer size validated: ${BLA_SIZEOF_INTEGER} bytes (LP64)")
        endif()
    else()
        message(STATUS "BLAS/LAPACK integer size check not available (CMake < 3.22 or vendor doesn't report)")
    endif()
endif()

# Store the found libraries/targets for linking.
if(TARGET BLAS::blas AND TARGET LAPACK::lapack)
    set(SLICOT_LINK_TARGETS BLAS::blas LAPACK::lapack)
    message(STATUS "Using BLAS/LAPACK targets: ${SLICOT_LINK_TARGETS}")
elseif(BLAS_LIBRARIES AND LAPACK_LIBRARIES)
    set(SLICOT_LINK_TARGETS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    message(STATUS "Using BLAS/LAPACK libraries: ${SLICOT_LINK_TARGETS}")
else()
    message(FATAL_ERROR "Could not determine BLAS/LAPACK link information after find_package.")
endif()

# Compiler-specific settings are now handled through CMake presets
# Global flags can be overridden by CMAKE_Fortran_FLAGS set in presets
message(STATUS "Using compiler flags from preset configuration")

# Add subdirectories
if(SLICOT_LINK_AUX_WRAPPERS) # Conditionally include src_aux based on the option
  message(STATUS "Auxiliary wrappers requested, adding src_aux.")
  add_subdirectory(src_aux)
else()
  message(STATUS "Auxiliary wrappers not requested.")
endif()
add_subdirectory(src)

# Add C wrapper subdirectory if requested
if(SLICOT_BUILD_C_WRAPPERS)
    message(STATUS "C wrappers requested, adding src_c_wrapper.")
    add_subdirectory(src_c_wrapper)
else()
    message(STATUS "C wrappers not requested.")
endif()

# Link object libraries after they have been defined
# We can rely on the add_subdirectory() order to ensure targets exist when we get here

# Link C wrapper objects if enabled and target exists
if(SLICOT_BUILD_C_WRAPPERS AND TARGET slicot_c_wrapper_obj)
    message(STATUS "Linking slicot library with C wrapper objects from slicot_c_wrapper_obj.")
    target_link_libraries(slicot PRIVATE $<TARGET_OBJECTS:slicot_c_wrapper_obj>)
    # Properly propagate include directories from the C wrapper object library
    target_include_directories(slicot PUBLIC 
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    )
endif()

# Link auxiliary objects if enabled and target exists
if(SLICOT_LINK_AUX_WRAPPERS AND TARGET lpkaux)
    message(STATUS "Linking slicot library with auxiliary wrappers from lpkaux.")
    target_link_libraries(slicot PRIVATE $<TARGET_OBJECTS:lpkaux>)
endif()

# Only build examples if requested
if(SLICOT_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Add tests subdirectory if requested
if(SLICOT_BUILD_TESTS)
    if(SLICOT_BUILD_C_WRAPPERS) # Only add tests if C wrappers are built
        message(STATUS "Building C wrapper tests.")
        add_subdirectory(tests)
    else()
        message(WARNING "SLICOT_BUILD_TESTS is ON, but SLICOT_BUILD_C_WRAPPERS is OFF. Skipping tests.")
    endif()
endif()

# Install the export targets
install(EXPORT SLICOTTargets
    FILE SLICOTTargets.cmake
    NAMESPACE SLICOT::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SLICOT)

# Create and install the config file
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/SLICOTConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SLICOTConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SLICOT)

# Create and install the version file
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SLICOTConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

# Install the config files
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/SLICOTConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SLICOTConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SLICOT)

# Install header files
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h"
)

# Conditionally install example documentation
if(SLICOT_BUILD_EXAMPLES)
    # The actual installation command is inside examples/CMakeLists.txt
    # but we might want to mention the component here if needed.
    # install(DIRECTORY examples/ DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples ...)
endif()

# Print build configuration summary
slicot_print_build_summary()
