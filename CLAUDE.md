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
- **Internal Workspace Allocation**: Wrappers **must allocate and free workspace arrays** (`iwork`, `dwork`, etc.) internally using `malloc`/`free`. Do not expect callers to provide workspace
- **Row-Major Handling**: Accept `row_major` flag. When true, allocate temporary column-major buffers (`_cm`), transpose C input matrices to column-major before Fortran calls, then transpose results back to row-major after Fortran returns (if `info == 0`)
- **Leading Dimensions**: For Fortran calls, `LDA` always refers to number of **rows**. For row-major C input, `LDA` refers to number of **columns** (wrapper must calculate Fortran LDA internally)
- **Zero-Dimension Validation**: Wrapper validation must align with underlying Fortran routine behavior. Many SLICOT routines handle zero dimensions as valid quick-exit cases (`INFO=0`). Don't flag as C-level errors if Fortran accepts them
- **NULL Pointer Handling**: For zero-sized arrays, `NULL` pointers are often acceptable if Fortran routine doesn't reference the array. Check SLICOT documentation for specific behavior
- Use `SLICOT_EXPORT` macros for Windows DLL compatibility
- Use `CHECK_ALLOC` macro after `malloc` and implement `goto cleanup` for error handling

### Test Structure
- Each C wrapper has corresponding `*_test.cpp` file
- Tests use GoogleTest framework with `ColMajor` and `RowMajor` fixtures
- **Test Data Sources**: For small datasets (~10 samples), embed data directly in test fixture. For larger datasets, create CSV file in `tests/data/` and use `load_test_data_from_csv` utility
- **CSV Data Loading**: Use correct function signature:
  ```cpp
  bool load_test_data_from_csv(
      const std::string& filepath,
      const std::vector<std::string>& input_cols,
      const std::vector<std::string>& output_cols,
      std::vector<double>& u,
      std::vector<double>& y,
      int& num_samples);
  ```
- **Numerical Tolerance**: Use realistic tolerances for `EXPECT_NEAR`. Start with `1e-3` but adjust to `5e-3` or higher if algorithm shows natural numerical variations. Trust computed results over documentation if consistently different
- **Zero-Dimension Testing**: Include specific test cases for zero-dimension inputs (e.g., `N=0`, `M=0`, `NSMP=0`) to verify wrapper validation aligns with Fortran behavior
- Reference results compared against known outputs from HTML documentation and verified against independent simulations when possible

## Development Workflow

### Adding New Functionality
1. **Analyze Documentation**: Review HTML docs, Fortran examples (`examples/*.f`, `*.dat`, `*.res`), and SLICOT standards (`rep96-1.pdf`) for parameter constraints, workspace formulas, and zero-dimension behavior
2. **Parameter Validation Design**: Ensure C wrapper validation aligns with (and is not stricter than) underlying Fortran routine behavior, especially for edge cases like `N=0`
3. Fortran routines go in `src/` following SLICOT naming conventions
4. C wrapper goes in `src_c_wrapper/` with corresponding header in `include/` (use lowercase filenames)
5. Add sources to `cmake/slicot_sources.cmake`
6. Create unit test in `tests/` with test data, including zero-dimension test cases
7. Update build files if needed

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

## Helper Functions Reference

### Matrix Transpose Functions (slicot_utils.h)
```c
// Convert from C row-major to Fortran column-major (for Fortran calls)
void slicot_transpose_to_fortran_with_ld(const void *src, void *dest, int rows, int cols,
                                        int ld_src, int ld_dest, size_t elem_size);

// Convert from Fortran column-major to C row-major (after Fortran calls)  
void slicot_transpose_to_c_with_ld(const void *src, void *dest, int rows, int cols,
                                  int ld_src, int ld_dest, size_t elem_size);

// Basic transpose functions (compatible default leading dimensions)
void slicot_transpose_to_fortran(const void *src, void *dest, int rows, int cols, size_t elem_size);
void slicot_transpose_to_c(const void *src, void *dest, int rows, int cols, size_t elem_size);
```

### Test Data Loading (test_utils.h)
```cpp
// CSV data loading function - use exact signature
bool load_test_data_from_csv(
    const std::string& filepath,
    const std::vector<std::string>& input_cols,
    const std::vector<std::string>& output_cols,
    std::vector<double>& u,
    std::vector<double>& y,
    int& num_samples);
```

### Essential Macros (slicot_utils.h)
```c
#define CHECK_ALLOC(ptr) /* Memory allocation checking - use after malloc */
#define SLICOT_MEMORY_ERROR -1010 /* Error code for allocation failures */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define SLICOT_EXPORT /* Platform-specific DLL export/import handling */
```

### **CRITICAL**: Function Name Corrections
**DO NOT USE these function names (they do not exist):**
- `slicot_transpose_from_fortran_with_ld()` ❌
- `load_csv_data()` or `load_test_csv()` ❌

**ALWAYS USE these correct function names:**
- `slicot_transpose_to_fortran_with_ld()` ✅ (C row-major → Fortran column-major)
- `slicot_transpose_to_c_with_ld()` ✅ (Fortran column-major → C row-major)
- `load_test_data_from_csv()` ✅

## Known Limitations
- Some test cases have known failures documented in `tests/Known_Issues.txt`
- ILP64 support is experimental
- Windows shared library builds require careful symbol export handling