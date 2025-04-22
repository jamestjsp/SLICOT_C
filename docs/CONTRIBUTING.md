# SLICOT_C Contribution Guidelines

This document outlines the coding patterns, conventions, and best practices to follow when contributing to the SLICOT_C wrapper library.

## 0. Project Directory Structure

The SLICOT-Reference project has the following directory structure to help you navigate the codebase:

```
SLICOT-Reference/
├── build/                   # Build output directories for different configurations
├── cmake/                   # CMake configuration files
├── doc/                     # SLICOT function documentation (HTML)
├── docs/                    # Project documentation files
├── include/                 # Public header files
│   └── *.h                  # C wrapper headers
├── src/                     # Original SLICOT Fortran source code
├── src_aux/                 # Auxiliary source files
├── src_c_wrapper/           # C wrapper implementation files
│   └── *.c                  # C wrapper source files
└── tests/                   # Test files
    ├── *_test.cpp           # Test files for individual SLICOT routines
    ├── CMakeLists.txt       # Test build configuration
    └── Readme.md            # Test documentation

Key Files:
- CMakeLists.txt             # Main build configuration
- README.md                  # Project overview
- CONTRIBUTING.md            # This file
```

### File Naming Conventions

- **Fortran source files**: Named according to original SLICOT routine names (e.g., `ab05od.f`)
- **C wrapper headers**: Named with lowercase routine names (e.g., `ab05od.h`)
- **C wrapper source**: Match header names (e.g., `ab05od.c`)
- **Test files**: Named with function name followed by `_test.cpp` (e.g., `ab05od_test.cpp`)

## 1. Overall Code Structure

### Wrapper File Organization
Each wrapper function should be organized in the following sequence:
1. File documentation header comment
2. Required includes
3. Declaration of the external Fortran function using `F77_FUNC` macro
4. C wrapper function definition with proper parameter validation
5. Row-major to column-major conversion (if needed)
6. Workspace allocation
7. Fortran function call
8. Result conversion back to row-major format (if needed)
9. Memory cleanup

### Naming Conventions
- C wrapper functions should be prefixed with `slicot_`
- Variable names should follow a consistent pattern:
  - Use `*_cm` suffix for column-major copies of matrices
  - Use `*_f` suffix for Fortran version of parameters (e.g., `lda_f` for Fortran leading dimension)
  - Use `*_ptr` suffix for pointers that will be passed to Fortran functions

## 2. Memory Management Patterns

### Allocation Patterns
```c
size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;
if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
```

### Cleanup Pattern
Always use the cleanup pattern with a label:
```c
cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(iwork);  // Safe even if NULL
    free(a_cm); free(b_cm); free(c_cm); free(d_cm);
    return info;
```

### Array Size Calculation
When calculating array sizes, always handle edge cases where dimensions might be zero:
```c
size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;
size_t b_size = (size_t)n * m; if (n == 0 || m == 0) b_size = 0;
```

## 3. Matrix Format Conversion

### Row-Major to Column-Major Conversion
For input matrices, convert from row-major to column-major format:
```c
if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
```

### Column-Major to Row-Major Conversion
For output matrices, convert back from column-major to row-major format:
```c
if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
```

### Leading Dimension Handling
When copying matrices with different leading dimensions:
```c
slicot_transpose_to_c_with_ld(a_cm, a, nr_val, nr_val, n, lda, elem_size);
```

## 4. Error Handling

### Input Parameter Validation
```c
if (n < 0) { info = -2; goto cleanup; }
if (char_upper != 'N' && char_upper != 'I') { info = -1; goto cleanup; }
```

### Memory Allocation Checks
Always use the `CHECK_ALLOC` macro after memory allocations:
```c
dwork = (double*)malloc((size_t)ldwork * sizeof(double));
CHECK_ALLOC(dwork);
```

## 5. Fortran Interface Patterns

### Character Parameters
Convert character parameters to uppercase and set the hidden length argument:
```c
char type_upper = toupper(type);
const int type_len = 1;  // Fortran expects 1-based length for strings
```

### Leading Dimensions
Calculate Fortran leading dimensions with consideration for zero-sized arrays:
```c
int lda_f = (a_rows > 0) ? a_rows : 1;
```

### Workspace Query
For workspace optimization, perform a workspace query:
```c
ldwork = -1;  // Query mode
F77_FUNC(routine, ROUTINE)(..., &dwork_query, &ldwork, &info, ...);
ldwork = (int)dwork_query;
ldwork = MAX(ldwork, min_ldwork);  // Ensure minimum size
```

### Workspace Allocation Strategy

#### Complete Workspace Query Pattern (Recommended)
For routines supporting workspace query, follow this complete pattern:

```c
// 1. Allocate IWORK with the required size from documentation
iwork_size = MAX(1, 2 * m);  // Ensure minimum size 1
iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
CHECK_ALLOC(iwork);

// 2. Set up for workspace query
double dwork_query[1];
ldwork = -1;  // Query mode

// 3. Use adjusted leading dimensions for query if in row-major mode
int lda_q = row_major ? MAX(1, n) : lda;
int ldb_q = row_major ? MAX(1, n) : ldb;
int ldc_q = row_major ? MAX(1, m) : ldc;
int ldd_q = row_major ? MAX(1, m) : ldd;

// 4. Perform workspace query
F77_FUNC(routine, ROUTINE)(..., &lda_q, ..., &ldb_q, ..., 
                           &ldc_q, ..., &ldd_q, ...,
                           iwork, dwork_query, &ldwork, &info);

// 5. Handle special return codes that can occur during query
if (info != 0 && info != expected_warning_code) {
    // Query failed for reasons other than expected warnings
    goto cleanup;
}
info = 0;  // Reset info after successful query

// 6. Get optimal workspace size and check against minimum
ldwork = (int)dwork_query[0];
ldwork = MAX(ldwork, MAX(1, 4 * m));  // Minimum from documentation

// 7. Allocate the workspace with optimal size
dwork = (double*)malloc((size_t)ldwork * sizeof(double));
CHECK_ALLOC(dwork);
```

#### Notes on Workspace Query
- Use temporary leading dimension variables for the query when in row-major mode
- For some routines, the query might return special warning codes (e.g., singularity warnings)
- Always reset `info` to 0 after a successful query, even if it contained warnings
- Even with optimal sizes from query, always enforce the minimum size from documentation
- The workspace query should work with the original arrays even in row-major mode since the query doesn't modify or use array values

#### Alternative Approach (No Query Available)
If the routine doesn't support workspace query:

```c
// Calculate minimum workspace size from SLICOT documentation
int min_ldwork = /* formula from documentation */;

// If documentation mentions "for better performance", allocate extra space
// A common approach is to double the minimum size
ldwork = 2 * min_ldwork;  // Double for better performance
ldwork = MAX(min_ldwork, ldwork);  // Ensure we have at least the minimum

// Allocate workspace
dwork = (double*)malloc((size_t)ldwork * sizeof(double));
CHECK_ALLOC(dwork);
```

#### Integer Workspace (IWORK)
For integer workspaces:

```c
iwork_size = /* formula from documentation */;
iwork_size = MAX(1, iwork_size);  // Ensure at least size 1 for valid allocation
iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
CHECK_ALLOC(iwork);
```

## 6. Documentation Standards

### File Headers
```c
/**
 * @file filename.c
 * @brief C wrapper implementation for SLICOT routine XXXXX
 *
 * This file provides a C wrapper implementation for the SLICOT routine XXXXX,
 * which [brief description of what the routine does].
 */
```

### Function Documentation
All public functions should include documentation in the header file:
```c
/**
 * @brief [Brief description]
 *
 * [Detailed description of what the function does]
 *
 * @param param1 [Description of parameter]
 * @param param2 [Description of parameter]
 * @return [Description of return value]
 */
```

## 7. Matrix Dimension Handling

### Leading Dimension Validation
Validate leading dimensions according to storage order:
```c
if (row_major) {
    // For row-major C, LDA is the number of columns
    if (lda < min_lda_rm_cols) { info = -7; goto cleanup; }
} else {
    // For column-major C, LDA is the number of rows
    if (lda < min_lda_f) { info = -7; goto cleanup; }
}
```

### Identity Matrix Initialization
When an identity matrix is needed:
```c
set_identity(n, matrix, ld, row_major);  // row_major=0 for column-major
```

## 8. Symmetric Matrix Handling

### Symmetric Matrix Transposition
For symmetric matrices:
```c
// Row-major symmetric to column-major full
slicot_transpose_symmetric_to_fortran(src, dest, n, 'U', elem_size);

// Column-major full to row-major symmetric
slicot_transpose_symmetric_to_c(src, dest, n, 'U', elem_size);
```

## 9. Optional Parameters

### NULL Pointers for Optional Parameters
Handle optional parameters with:
```c
(param_needed ? param_ptr : NULL)
```

## 10. Test Development Guidelines

### Test Organization Pattern

Each wrapper function should have at least two test suites:

```cpp
// Column-major test fixture (original Fortran ordering)
class FunctionNameTestColMajor : public ::testing::Test {
protected:
    // Define common variables and expected values from documentation
    int N = 7, M = 2, P = 3;
    char DICO = 'C';
    double TOL = 1.0e-10;
    int ROW_MAJOR = 0;  // Column-major
    double check_tol = 1e-5;  // Tolerance for result verification
    
    // Expected results from documentation
    double expected_value = 7.93948;
    int expected_info = 0;
    
    // Input matrices in column-major order (as they appear in SLICOT docs)
    std::vector<double> A_in = {
        // Column 1
        -0.04165, -5.2100, 0.0000, 0.5450, 0.0000, 0.0000, 0.0000,
        // Column 2
        0.0000, -12.500, 3.3300, 0.0000, 0.0000, 0.0000, 0.0000,
        // ...more columns...
    };
};

// Row-major test fixture (inherit from column-major)
class FunctionNameTestRowMajor : public FunctionNameTestColMajor {
public:
    FunctionNameTestRowMajor() {
        ROW_MAJOR = 1;  // Override for row-major tests
    }
};
```

### Handling Input Data from Documentation

When reading example data from SLICOT HTML documentation:

1. **Identify Matrix Dimensions**: Extract dimensions N, M, P from the documentation.

2. **Parse Matrices from Program Data Section**:
   - SLICOT examples typically present matrices row by row for readability
   - For column-major tests, store them in column-major order
   - For row-major tests, use `slicot_transpose_to_c` to convert

3. **Extract Expected Results from Program Results Section**:
   - Get expected return values, modified matrices, and any warnings/flags

Example for parsing input matrices from documentation:

```cpp
// Column-major test case
TEST_F(FunctionNameTestColMajor, DocExample) {
    // Define leading dimensions for column-major format
    int LDA = N;  // >= N
    int LDB = N;  // >= N
    int LDC = P;  // >= P
    int LDD = P;  // >= P
    
    // Copy input matrices (they may be modified by the function)
    std::vector<double> A = A_in;
    std::vector<double> B = B_in;
    
    // Define output variables
    double result;
    int info;
    
    // Call wrapper function
    result = slicot_functionname(...parameters..., ROW_MAJOR);
    
    // Verify results
    ASSERT_EQ(info, expected_info);
    EXPECT_NEAR(result, expected_value, check_tol);
    // Check modified matrices if applicable
}

// Row-major test case
TEST_F(FunctionNameTestRowMajor, DocExample) {
    // Define leading dimensions for row-major format
    int LDA = N;  // Cols for A (NxN)
    int LDB = M;  // Cols for B (NxM)
    int LDC = N;  // Cols for C (PxN)
    
    // Convert column-major inputs to row-major
    std::vector<double> A(N * LDA);
    std::vector<double> B(N * LDB);
    
    slicot_transpose_to_c(A_in.data(), A.data(), N, N, sizeof(double));
    slicot_transpose_to_c(B_in.data(), B.data(), N, M, sizeof(double));
    
    // Call and verify as in column-major case
}
```

### Matrix Testing Patterns

1. **Basic Input Matrix Validation**:
   ```cpp
   ASSERT_EQ(A.size(), N * LDA);  // Ensure matrix is correctly sized
   ```

2. **Matrix Element Comparison**:
   ```cpp
   // Column-major comparison
   for (int j = 0; j < expected_NR; ++j) { // Col
       for (int i = 0; i < expected_NR; ++i) { // Row
           EXPECT_NEAR(A[i + j*LDA], expected_A[i + j*expected_NR], check_tol)
               << "A[" << i << "," << j << "] mismatch";
       }
   }
   
   // Row-major comparison
   for (int i = 0; i < expected_NR; ++i) { // Row
       for (int j = 0; j < expected_NR; ++j) { // Col
           EXPECT_NEAR(A[i*LDA + j], expected_A_rm[i*expected_NR + j], check_tol)
               << "A[" << i << "," << j << "] mismatch";
       }
   }
   ```

3. **Converting Expected Results for Row-Major Tests**:
   ```cpp
   // Convert expected column-major results to row-major for comparison
   std::vector<double> expected_A_rm(expected_NR * expected_NR);
   slicot_transpose_to_c(expected_A.data(), expected_A_rm.data(), expected_NR, expected_NR, sizeof(double));
   ```

4. **Helper Function for Matrix Debugging**:
   ```cpp
   void printMatrixD(const std::string& name, const double* data, int rows, int cols, 
                     int ld, bool rowMajor) {
       std::cout << name << " (" << rows << "x" << cols << "):\n";
       for (int i = 0; i < rows; ++i) {
           for (int j = 0; j < cols; ++j) {
               double val = rowMajor ? data[i*ld + j] : data[i + j*ld];
               std::cout << std::setw(9) << std::fixed << std::setprecision(4) 
                         << val << " ";
           }
           std::cout << "\n";
       }
   }
   ```

### Reading SLICOT Documentation Examples

Refer to the [documentation parsing guide](/Users/josephj/Workspace/SLICOT_C/tests/Readme.md) for details on how to interpret matrix data from SLICOT HTML documentation.

Key points:
1. Matrix data is typically presented row by row in the documentation
2. The Fortran code uses column-major storage, so be careful when parsing
3. Use the matrix dimensions from the documentation to guide parsing
4. For input matrices, ensure they match the dimensions specified in the documentation
5. For output matrices, compare with the expected results in the documentation

## 11. Additional Utility Functions and Macros

### Symmetric Matrix Utilities

Beyond the basic symmetric matrix transposition functions, these utilities can simplify symmetric matrix handling:

```c
// Copy only specified triangle (upper/lower) from a symmetric matrix to a full matrix
slicot_copy_symmetric_part(src, dest, n, 'U', ld, elem_size);
```

### Complex Number Handling

For routines that work with complex numbers:

```c
// Define complex variables
slicot_complex_double z;

// Access real part of complex number
double real_part = SLICOT_COMPLEX_REAL(z);
```

### Printing Matrices for Debugging

During development, this helper function can be included in test files.
**Note:** For printing matrices during development or debugging, use the `printMatrixD` function provided in `slicot_utils.h`. Do not redefine this function in your test files; simply include the header and call `printMatrixD` as needed.

### Fortran Name Mangling

Use the `F77_FUNC` macro for platform-independent Fortran function name mangling:

```c
// For function name "ab01nd" that is "AB01ND" in uppercase
extern void F77_FUNC(ab01nd, AB01ND)(...);
```

This macro handles different compiler conventions:
- On Windows with Intel Fortran: expands to uppercase name (AB01ND)
- On Windows with MinGW/gfortran: expands to lowercase with underscore (ab01nd_)
- On Linux/macOS systems: expands to lowercase with underscore (ab01nd_)

### DLL Export/Import Macros

When developing functions that need to be exported from a DLL:

```c
// In the implementation file
SLICOT_C_WRAPPER_API
void my_exported_function(...)
```

The `SLICOT_C_WRAPPER_API` macro handles platform-specific export declarations:
- On Windows: `__declspec(dllexport)` when building, `__declspec(dllimport)` when using
- On Linux/macOS with GCC/Clang: `__attribute__((visibility("default")))`
- When building static library: empty definition

### Error Code Constants

Use predefined error codes in your functions:

```c
if (!ptr) {
    info = SLICOT_MEMORY_ERROR;  // Standard memory allocation error (-1010)
    goto cleanup;
}
```

### Expected Return Codes and Warnings

Handle expected return codes and warnings from SLICOT routines:

```c
// Routine-specific warnings that should not terminate execution
if (info != 0 && info != expected_warning_code) {
    goto cleanup;
}

// Numerical singularity might be reported as a specific positive code
if (info == singularity_warning) {
    // Handle or report warning but continue
    info = 0;  // Clear the warning if it's acceptable
}
```

### Working with In-place Transposition

For square matrices that need to be transposed in-place:

```c
// Returns 0 on success, -1 if matrix is not square
int result = slicot_transpose_inplace(matrix, n, n, sizeof(double));
```

## 12. Handling LAPACK Dependencies

Many SLICOT routines depend on LAPACK functions that must be properly declared and used in the C wrappers.

### Declaring LAPACK Function Dependencies

When a SLICOT routine requires calling a LAPACK function (often mentioned in the documentation), declare the external function using the F77_FUNC macro:

```c
/* 
 * Declare LAPACK function dependencies 
 * (mentioned in SLICOT documentation)
 */
extern void F77_FUNC(dorgqr, DORGQR)(
    const int* m,      // INTEGER M
    const int* n,      // INTEGER N
    const int* k,      // INTEGER K
    double* a,         // DOUBLE PRECISION A(LDA,*)
    const int* lda,    // INTEGER LDA
    const double* tau, // DOUBLE PRECISION TAU(*)
    double* work,      // DOUBLE PRECISION WORK(*)
    const int* lwork,  // INTEGER LWORK
    int* info          // INTEGER INFO
);
```

### Identifying LAPACK Dependencies

Look for LAPACK dependencies in:
1. The "Further Comments" section of SLICOT documentation
2. Places where the documentation mentions forming matrices from factorized results
3. Statements like "can be obtained by calling the LAPACK library routine XXX"

Common LAPACK dependencies include:
- `DORGQR`/`ZUNGQR`: Generate orthogonal/unitary matrix from elementary reflectors
- `DORMQR`/`ZUNMQR`: Apply orthogonal/unitary matrix from QR factorization
- `DGESVD`/`ZGESVD`: Singular value decomposition
- `DGEQRF`/`ZGEQRF`: QR factorization

### Using LAPACK Functions in Wrappers

When calling LAPACK functions from wrappers:

```c
/* If JOBZ='F', we need to form the complete orthogonal matrix using DORGQR */
if (info == 0 && jobz_upper == 'F' && n > 0) {
    /* Construct the orthogonal matrix from elementary reflectors */
    int dorgqr_info = 0;
    
    if (row_major) {
        /* Call DORGQR to form the orthogonal matrix in column-major format */
        F77_FUNC(dorgqr, DORGQR)(&n, &n, &n, z_cm, &ldz_f, tau, dwork, &ldwork, &dorgqr_info);
        if (dorgqr_info != 0) {
            info = dorgqr_info;
            goto cleanup;
        }
    } else {
        /* Call DORGQR to form the orthogonal matrix directly */
        F77_FUNC(dorgqr, DORGQR)(&n, &n, &n, z, &ldz, tau, dwork, &ldwork, &dorgqr_info);
        if (dorgqr_info != 0) {
            info = dorgqr_info;
            goto cleanup;
        }
    }
}
```

### Workspace Considerations for LAPACK Functions

When a LAPACK function will be used, ensure the workspace is sized properly:

```c
/* According to AB01ND documentation, LDWORK >= MAX(1, N, 3*M) is required. */
/* For JOBZ='F', DORGQR will need additional workspace, we'll reuse this array */
ldwork = MAX(1, MAX(n, 3*m));
if (jobz_upper == 'F' && n > 0) {
    /* Need extra space for DORGQR: optimal LDWORK >= N */
    ldwork = MAX(ldwork, n);
}
ldwork = ldwork*2; // Double the size for better performance
```

### Testing LAPACK-Dependent Features

When testing functions with LAPACK dependencies:
1. Explicitly test the cases that trigger LAPACK function calls (e.g., JOBZ='F')
2. Verify the results match expected outputs from the documentation
3. Consider adding specific tests for the LAPACK-dependent functionality

## 13. Guidelines for AI Coding Agents

When using AI coding agents to implement SLICOT C wrappers and tests, following these guidelines can improve efficiency and accuracy:

### Efficient Information Gathering

1. **HTML Documentation Analysis**:
   - Analyze the corresponding HTML documentation (in the `doc/` directory) first
   - Extract key information in this order:
     1. Purpose and mathematical function of the routine
     2. Input/output parameters with dimensions
     3. Workspace requirements
     4. Example data/results for testing
     5. Special cases and error conditions

2. **Pattern Matching**:
   - Identify similar routines that have already been implemented
   - Look for wrappers with similar parameter structures or functionality
   - For example, routines in the same family (AB05xx, AB09xx, etc.) often follow similar patterns

3. **Fortran Source Analysis**:
   - If the wrapper implementation is complex, examine the original Fortran source (in `src/`)
   - Focus on array dimensions, workspace usage, and return code handling

### Implementation Strategy

1. **Wrapper Implementation Sequence**:
   - Follow the established structure sequence exactly (see Section 1)
   - Keep the parameter list as close as possible to the original Fortran function
   - Validate all input parameters before allocating resources

2. **Matrix Format Conversion Strategy**:
   - For large matrices, verify dimensions before allocating memory
   - When implementing new routines, follow the example of existing wrappers for similar routines
   - For any new helper functions, add detailed comments explaining the matrix layout conversions

3. **Test Design**:
   - Create tests for both column-major and row-major formats
   - Include at least one test with the example data from the HTML documentation
   - Add boundary condition tests (zero dimensions, minimum workspace, etc.)

### Troubleshooting Common Issues

1. **Memory Management Issues**:
   - Double-check allocation/deallocation patterns for consistency
   - Ensure every allocated resource has a corresponding free in the cleanup section
   - Watch for special cases where zero dimensions might cause issues

2. **Matrix Transposition Errors**:
   - Verify the dimensions used in transposition functions
   - Ensure leading dimensions are handled correctly for both input and output
   - For symmetric matrices, use the specific symmetric handling utilities

3. **Testing Discrepancies**:
   - If test results don't match documentation examples:
     1. Verify the matrix data was parsed correctly from the HTML docs
     2. Check if the original source has any undocumented behaviors
     3. Consider numerical precision issues (use appropriate check_tol values)

### Documentation and Commenting

1. **AI-Generated Comments Quality**:
   - Document parameter validation logic clearly
   - Explain any deviations from the standard wrapper pattern
   - Include references to the SLICOT documentation when implementing special cases

2. **Test Case Documentation**:
   - Note the source of test values (e.g., "from AB05OD.html example")
   - Document expected behavior for negative tests
   - Explain any transformations applied to the test data

By following these guidelines, AI coding agents can more efficiently and accurately implement SLICOT C wrappers and tests, maintaining consistency with the project's coding standards and patterns.

## 14. Code Template Examples

### Standard Wrapper Template

```c
/**
 * @file xxxxx.c
 * @brief C wrapper implementation for SLICOT routine XXXXX
 *
 * This file provides a C wrapper implementation for the SLICOT routine XXXXX,
 * which [brief description of what the routine does].
 */

#include <stdlib.h>
#include <ctype.h>  // For toupper
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "xxxxx.h"

// Include necessary SLICOT utility headers
#include "slicot_utils.h" // For MAX, CHECK_ALLOC, etc.
#include "slicot_f77.h"   // For F77_FUNC macro

// Declare the external Fortran routine
extern void F77_FUNC(xxxxx, XXXXX)(
    // Fortran function parameters with types
    // ...
);

/* C wrapper function definition */
SLICOT_C_WRAPPER_API
int slicot_xxxxx(/* C function parameters */)
{
    // Variable declarations
    int info = 0;
    
    // Input parameter validation
    
    // Memory allocation for column-major copies (if needed)
    
    // Perform row-major to column-major conversion (if needed)
    
    // Workspace allocation
    
    // Call the Fortran routine
    
    // Convert results back to row-major format (if needed)
    
cleanup:
    /* --- Cleanup --- */
    // Free all allocated memory
    
    return info;
}
```

### Standard Test Template

```cpp
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

// Include the header for the function being tested
#include "xxxxx.h"

// Define a fixture for column-major tests
class XXXXXTestColMajor : public ::testing::Test {
protected:
    // Define common test variables
    // Input parameters, expected output values, etc.
    
    // Column-major input matrices and expected outputs
};

// Row-major test fixture (inherit from column-major)
class XXXXXTestRowMajor : public XXXXXTestColMajor {
public:
    XXXXXTestRowMajor() {
        ROW_MAJOR = 1;  // Override for row-major tests
        // Convert matrices from column-major to row-major
    }
};

// Test case for column-major format
TEST_F(XXXXXTestColMajor, DocExample) {
    // Set up test parameters
    // Call the wrapper function
    // Verify results
}

// Test case for row-major format
TEST_F(XXXXXTestRowMajor, DocExample) {
    // Set up test parameters
    // Call the wrapper function
    // Verify results  
}

// Additional test cases for error conditions, special cases, etc.
```
