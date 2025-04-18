cmake_minimum_required(VERSION 3.15)

# ---
# NOTE: The SLICOT C wrapper generator always places Fortran CHARACTER string length arguments at the end of the argument list.
# This is portable for all compilers except Intel Fortran (ifort/ifx), which by default passes string lengths immediately after each string.
# The CMake configuration below sets the Intel Fortran flag -assume noestring_param_xlate to ensure string lengths are passed at the end,
# making the wrappers portable across compilers. Do not change this unless you know what you are doing.
# ---

# Project definition
project(SLICOT VERSION 5.9.0 LANGUAGES Fortran C)

# Options
option(SLICOT_BUILD_SHARED_LIBS "Build shared libraries" ON) # Default to shared for DLLs
option(SLICOT_BUILD_EXAMPLES "Build the SLICOT example programs and tests" OFF)
option(SLICOT_USE_ILP64 "Build with 64-bit integers (ILP64)" OFF) # Add ILP64 option
option(SLICOT_LINK_AUX_WRAPPERS "Link auxiliary LAPACK wrappers (for incomplete LAPACK)" ON) # Add option for aux wrappers
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

# Add ILP64 flags if requested
if(SLICOT_USE_ILP64)
  message(STATUS "ILP64 build requested.")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(ILP64_FLAGS "-fdefault-integer-8")
    message(STATUS "Using gfortran ILP64 flag: ${ILP64_FLAGS}")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" OR CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    set(ILP64_FLAGS "-i8")
    message(STATUS "Using Intel Fortran ILP64 flag: ${ILP64_FLAGS}")
  else()
    message(WARNING "ILP64 requested, but compiler flags are unknown for ${CMAKE_Fortran_COMPILER_ID}. Manual configuration may be required.")
    set(ILP64_FLAGS "")
  endif()
  # Add flags globally - apply to slicot and examples
  add_compile_options(${ILP64_FLAGS})
endif()

# Include standard modules
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

# Find BLAS and LAPACK - BLA_VENDOR might be set by presets
# Default to searching all vendors if not set
if(NOT DEFINED BLA_VENDOR)
  set(BLA_VENDOR "All")
endif()
message(STATUS "Searching for BLAS/LAPACK vendor: ${BLA_VENDOR}")

# Configure paths for BLAS/LAPACK discovery
if(SLICOT_USE_VCPKG)
  message(STATUS "Using vcpkg for BLAS/LAPACK dependencies")
  # Provide hints to standard FindBLAS/FindLAPACK using vcpkg install path
  if(DEFINED ENV{VCPKG_ROOT} AND IS_DIRECTORY "$ENV{VCPKG_ROOT}/installed/${VCPKG_TARGET_TRIPLET}")
      set(CMAKE_PREFIX_PATH "$ENV{VCPKG_ROOT}/installed/${VCPKG_TARGET_TRIPLET}" ${CMAKE_PREFIX_PATH})
      list(APPEND CMAKE_SYSTEM_PREFIX_PATH "$ENV{VCPKG_ROOT}/installed/${VCPKG_TARGET_TRIPLET}")
      message(STATUS "Adding vcpkg install prefix to CMAKE_PREFIX_PATH: $ENV{VCPKG_ROOT}/installed/${VCPKG_TARGET_TRIPLET}")
  else()
      message(WARNING "VCPKG_ROOT environment variable not defined or install directory not found. FindBLAS/LAPACK might fail.")
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

# Try to find BLAS and LAPACK if not already manually set
if(NOT BLAS_FOUND)
  find_package(BLAS REQUIRED)
endif()

if(NOT LAPACK_FOUND)
  find_package(LAPACK REQUIRED)
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

# Compilation settings
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # Add flags for GNU Fortran compiler
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    # Ensure lowercase symbols with trailing underscore
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-underscoring -fno-second-underscore")
    # Override with explicit lowercase and single underscore
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcase-lower -fsecond-underscore")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel" OR CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    # Intel Fortran compiler flags
    # String handling for Intel - ensures string lengths at end of argument list
    if(WIN32)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /iface:nomixed_str_len_arg")
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -nomixed-str-len-arg")
    endif()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
    # Flang compiler flags
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcase-lower -fno-underscoring")
    # Re-add single trailing underscore
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fsecond-underscore")
endif()

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

# Only build examples if requested
if(SLICOT_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Conditionally enable testing
if(SLICOT_BUILD_EXAMPLES)
    enable_testing()
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

# Conditionally install example documentation
if(SLICOT_BUILD_EXAMPLES)
    # The actual installation command is inside examples/CMakeLists.txt
    # but we might want to mention the component here if needed.
    # install(DIRECTORY examples/ DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples ...)
endif()

message(STATUS "  Build Examples: ${SLICOT_BUILD_EXAMPLES}")
