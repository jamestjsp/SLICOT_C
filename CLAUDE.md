# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SLICOT (Subroutine Library In COntrol Theory) is a numerical library for control theoretical computations implemented in Fortran 77. This repository provides:
- Original SLICOT Fortran routines in `src/`
- C wrapper layer in `src_c_wrapper/` for C/C++ integration
- GoogleTest-based test suite in `tests/`
- Example programs in `examples/`
- Both CMake (modern, recommended) and legacy Makefile build systems

## Build System Commands

### CMake Build (Recommended)
```bash
# Configure using preset (choose appropriate for your platform)
cmake --preset macos-x64-debug     # For macOS Debug
cmake --preset macos-x64-release   # For macOS Release
cmake --preset linux-x64-debug     # For Linux Debug  
cmake --preset windows-x64-debug   # For Windows Debug

# Build
cmake --build --preset macos-x64-debug-build

# Run tests
ctest --preset macos-x64-debug-test

# Install
cmake --build --preset macos-x64-debug-build --target install
```

### Legacy Makefile Build
```bash
# Unix-like systems
make -f makefile_Unix

# Windows with nmake
nmake -f makefile
```

### Key CMake Options
- `SLICOT_BUILD_EXAMPLES=ON/OFF` - Build example programs (default: OFF)
- `SLICOT_BUILD_TESTS=ON/OFF` - Build GoogleTest test suite (default: ON)
- `SLICOT_BUILD_C_WRAPPERS=ON/OFF` - Build C wrapper layer (default: ON)
- `SLICOT_BUILD_SHARED_LIBS=ON/OFF` - Build shared libraries (default: ON)
- `SLICOT_USE_ILP64=ON/OFF` - Use 64-bit integers (default: OFF)

## Architecture

### Core Components
- **`src/`** - Original SLICOT Fortran 77 routines (500+ files)
- **`src_c_wrapper/`** - C wrapper implementations calling Fortran via F77_FUNC macro
- **`src_aux/`** - Auxiliary LAPACK compatibility routines
- **`include/`** - C header files for wrapper functions
- **`tests/`** - GoogleTest-based unit tests for C wrappers

### Library Organization (SLICOT chapters)
- **A** - Analysis Routines
- **B** - Benchmark and Test Problems  
- **D** - Data Analysis
- **F** - Filtering
- **I** - Identification
- **M** - Mathematical Routines
- **N** - Nonlinear Systems
- **S** - Synthesis Routines
- **T** - Transformation Routines
- **U** - Utility Routines

### C Wrapper Pattern
C wrappers follow consistent patterns:
- Use `F77_FUNC(fortran_name, FORTRAN_NAME)` for Fortran interop
- Handle memory allocation/deallocation
- Transpose matrices between C (row-major) and Fortran (column-major) ordering
- Perform input validation and error checking
- Use `SLICOT_EXPORT` macros for Windows DLL compatibility

### Test Structure
- Each C wrapper has corresponding `*_test.cpp` file
- Tests use GoogleTest framework
- Test data often loaded from CSV files in `tests/data/`
- Reference results compared against known outputs from HTML documentation

## Development Workflow

### Adding New Functionality
1. Fortran routines go in `src/` following SLICOT naming conventions
2. C wrapper goes in `src_c_wrapper/` with corresponding header in `include/`
3. Add sources to `cmake/slicot_sources.cmake`
4. Create unit test in `tests/` with test data
5. Update build files if needed

### Running Single Test
```bash
# Build first
cmake --build --preset macos-x64-debug-build

# Run specific test
./build/macos-x64-debug/tests/slicot_c_wrapper_tests --gtest_filter="AB01MD*"
```

### Common Build Issues
- Ensure BLAS/LAPACK are available (Intel MKL, OpenBLAS, Apple Accelerate, etc.)
- For Intel Fortran: CMake automatically sets `-nomixed-str-len-arg` flag
- For gfortran: CMake sets `-fallow-argument-mismatch -fno-second-underscore`
- Use `SLICOT_USE_VCPKG=ON` if using vcpkg for dependencies

### Test Data Format
- Fortran uses column-major storage, HTML examples show row-major presentation
- When parsing test data from HTML docs, read sequentially but account for potential transpose
- CSV test data files maintain consistency with Fortran READ statements
- See `tests/Readme.md` for detailed data parsing guidelines

## Dependencies

- **Required**: Fortran compiler (gfortran, Intel ifort/ifx), BLAS, LAPACK
- **Build**: CMake 3.15+, C/C++ compiler 
- **Testing**: GoogleTest (automatically fetched by CMake)
- **Optional**: Python 3.x for some test verification scripts

## Key Files to Understand
- `CMakeLists.txt` - Main build configuration
- `CMakePresets.json` - Platform-specific build presets  
- `cmake/slicot_sources.cmake` - Centralized source file lists
- `include/slicot_f77.h` - Fortran interop definitions
- `include/slicot_utils.h` - Common C wrapper utilities
- `tests/test_utils.h` - Test helper functions

## Known Limitations
- Some test cases have known failures documented in `tests/Known_Issues.txt`
- ILP64 support is experimental
- Windows shared library builds require careful symbol export handling