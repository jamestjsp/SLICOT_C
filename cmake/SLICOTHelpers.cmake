# SLICOT CMake Helper Functions
# This file contains utility functions for SLICOT build configuration

# Function to detect and configure optimal BLAS/LAPACK based on system
function(slicot_configure_blas_lapack)
    # Set BLA_VENDOR based on system if not already set
    if(NOT DEFINED BLA_VENDOR OR BLA_VENDOR STREQUAL "")
        if(APPLE)
            set(BLA_VENDOR "Apple" PARENT_SCOPE)
            message(STATUS "Auto-detected Apple Accelerate framework for BLAS/LAPACK")
        elseif(DEFINED ENV{MKLROOT})
            if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
                set(BLA_VENDOR "Intel10_64lp" PARENT_SCOPE)
            else()
                set(BLA_VENDOR "Intel10_64lp_seq" PARENT_SCOPE)
            endif()
            message(STATUS "Auto-detected Intel MKL for BLAS/LAPACK")
        else()
            set(BLA_VENDOR "Generic" PARENT_SCOPE)
            message(STATUS "Using generic BLAS/LAPACK detection")
        endif()
    endif()
endfunction()

# Function to print build configuration summary
function(slicot_print_build_summary)
    message(STATUS "")
    message(STATUS "SLICOT Build Configuration Summary:")
    message(STATUS "==================================")
    message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
    message(STATUS "  Install Prefix: ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "  Build Examples: ${SLICOT_BUILD_EXAMPLES}")
    message(STATUS "  Build Tests: ${SLICOT_BUILD_TESTS}")
    message(STATUS "  Build C Wrappers: ${SLICOT_BUILD_C_WRAPPERS}")
    message(STATUS "  Build Shared Libraries: ${SLICOT_BUILD_SHARED_LIBS}")
    message(STATUS "  Link Aux Wrappers: ${SLICOT_LINK_AUX_WRAPPERS}")
    message(STATUS "  Use ILP64: ${SLICOT_USE_ILP64}")
    message(STATUS "  Use vcpkg: ${SLICOT_USE_VCPKG}")
    
    if(DEFINED BLA_VENDOR)
        message(STATUS "  BLAS/LAPACK Vendor: ${BLA_VENDOR}")
    endif()
    
    message(STATUS "  Compilers:")
    message(STATUS "    C: ${CMAKE_C_COMPILER}")
    message(STATUS "    C++: ${CMAKE_CXX_COMPILER}")
    message(STATUS "    Fortran: ${CMAKE_Fortran_COMPILER}")
    message(STATUS "")
endfunction()

# Function to apply additional target-specific compiler flags for SLICOT
# Note: Platform and compiler-specific flags are now primarily handled via CMake presets
function(slicot_configure_compiler_flags target)
    # Apply any additional target-specific flags if needed
    # Most compiler flags are now set globally via CMAKE_Fortran_FLAGS in presets
    # This function is kept for future extensibility
    
    # Example: Add debug-specific flags that might not be in presets
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
            target_compile_options(${target} PRIVATE
                $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
                $<$<COMPILE_LANGUAGE:Fortran>:-Wextra>
            )
        endif()
    endif()
    
    message(STATUS "Applied additional compiler flags to target: ${target}")
endfunction()