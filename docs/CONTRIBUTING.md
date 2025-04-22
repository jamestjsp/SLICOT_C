<!-- filepath: .\SLICOT_C\docs\CONTRIBUTING.md -->
# SLICOT_C Contribution Guidelines

This document outlines the coding patterns, conventions, and best practices to follow when contributing to the SLICOT_C wrapper library.

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Implementation Workflow](#2-implementation-workflow)
3. [Code Structure and Patterns](#3-code-structure-and-patterns)
4. [Memory Management](#4-memory-management)
5. [Matrix Format Handling](#5-matrix-format-handling)
6. [Error Handling](#6-error-handling)
7. [Fortran Interface](#7-fortran-interface)
8. [Documentation Standards](#8-documentation-standards)
9. [Test Development](#9-test-development)
10. [Test Data Extraction](#10-test-data-extraction)
11. [Common Error Cases](#11-common-error-cases)
12. [Utility Functions](#12-utility-functions)
13. [LAPACK Dependencies](#13-lapack-dependencies)
14. [Troubleshooting Guide](#14-troubleshooting-guide)
15. [Code Templates](#15-code-templates)

## 1. Project Overview

### 1.1 Directory Structure

```
SLICOT-Reference/
├── build/                   # Build output directories for different configurations
├── cmake/                   # CMake configuration files
├── doc/                     # SLICOT function documentation (HTML)
├── docs/                    # Project documentation files
├── examples/                # Example data files and test programs
├── include/                 # Public header files (*.h)
├── src/                     # Original SLICOT Fortran source code
├── src_aux/                 # Auxiliary source files
├── src_c_wrapper/           # C wrapper implementation files (*.c)
└── tests/                   # Test files (*_test.cpp)
```

### 1.2 File Naming Conventions

- **Fortran source files**: Named according to original SLICOT routine names (e.g., `ab05od.f`)
- **C wrapper headers**: Named with lowercase routine names (e.g., `ab05od.h`)
- **C wrapper source**: Match header names (e.g., `ab05od.c`)
- **Test files**: Named with function name followed by `_test.cpp` (e.g., `ab05od_test.cpp`)

## 2. Implementation Workflow

Follow this sequence when implementing a new wrapper:

1. Analyze the HTML documentation (e.g., `doc/AB05OD.html`)
2. Extract key information: parameters, dimensions, workspace requirements
3. Identify similar existing wrappers (same function family)
4. Create C wrapper following the template pattern
5. Create test cases using example data from documentation or example files

## 3. Code Structure and Patterns

### 3.1 Wrapper File Organization

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

### 3.2 Naming Conventions

- C wrapper functions should be prefixed with `slicot_`
- Variable names should follow a consistent pattern:
  - Use `*_cm` suffix for column-major copies of matrices
  - Use `*_f` suffix for Fortran version of parameters (e.g., `lda_f` for Fortran leading dimension)
  - Use `*_ptr` suffix for pointers that will be passed to Fortran functions

## 4. Memory Management

### 4.1 Allocation Patterns

```c
// Calculate size with zero dimension handling
size_t a_size = (size_t)m * n; if (m == 0 || n == 0) a_size = 0;

// Allocate memory if needed
if (a_size > 0) {
    a_cm = (double*)malloc(a_size * sizeof(double));
    CHECK_ALLOC(a_cm);
}
```

### 4.2 Cleanup Pattern

Always use the cleanup pattern with a label:

```c
cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(iwork);  // Safe even if NULL
    free(a_cm); free(b_cm); free(c_cm); free(d_cm);
    return info;
```

### 4.3 Array Size Calculation

When calculating array sizes, always handle edge cases where dimensions might be zero:

```c
size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;
size_t b_size = (size_t)n * m; if (n == 0 || m == 0) b_size = 0;
```

## 5. Matrix Format Handling

### 5.1 Row-Major to Column-Major Conversion

For input matrices, convert from row-major to column-major format:

```c
if (row_major && a_size > 0) {
    slicot_transpose_to_fortran(a, a_cm, m, n, sizeof(double));
}
```

### 5.2 Column-Major to Row-Major Conversion

For output matrices, convert back from column-major to row-major format:

```c
if (row_major && a_size > 0 && info == 0) {
    slicot_transpose_to_c(a_cm, a, m, n, sizeof(double));
}
```

### 5.3 Leading Dimension Handling

When copying matrices with different leading dimensions:

```c
slicot_transpose_to_c_with_ld(a_cm, a, nr_val, nr_val, n, lda, sizeof(double));
```

### 5.4 Symmetric Matrix Handling

For symmetric matrices:

```c
// Row-major symmetric to column-major full
slicot_transpose_symmetric_to_fortran(src, dest, n, 'U', elem_size);

// Column-major full to row-major symmetric
slicot_transpose_symmetric_to_c(src, dest, n, 'U', elem_size);
```

## 6. Error Handling

### 6.1 Input Parameter Validation

```c
// Check dimensions
if (n < 0) { info = -2; goto cleanup; }
if (m < 0) { info = -3; goto cleanup; }

// Check option parameters
if (toupper(dico) != 'C' && toupper(dico) != 'D') { info = -1; goto cleanup; }

// Check leading dimensions
if (row_major) {
    if (lda < n) { info = -6; goto cleanup; }  // Columns in row-major
} else {
    if (lda < m) { info = -6; goto cleanup; }  // Rows in column-major
}
```

### 6.2 Memory Allocation Checks

Always use the `CHECK_ALLOC` macro after memory allocations:

```c
dwork = (double*)malloc((size_t)ldwork * sizeof(double));
CHECK_ALLOC(dwork);
```

### 6.3 Error Code Constants

Use predefined error codes in your functions:

```c
if (!ptr) {
    info = SLICOT_MEMORY_ERROR;  // Standard memory allocation error (-1010)
    goto cleanup;
}
```

### 6.4 Expected Return Codes and Warnings

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

## 7. Fortran Interface

### 7.1 Character Parameters

Convert character parameters to uppercase and set the hidden length argument:

```c
char dico_upper = toupper(dico);
const int dico_len = 1;  // Hidden length for Fortran strings

// Pass in the Fortran call
F77_FUNC(routine, ROUTINE)(&dico_upper, ..., dico_len);
```

### 7.2 Leading Dimensions

Calculate Fortran leading dimensions with consideration for zero-sized arrays:

```c
int lda_f = (a_rows > 0) ? a_rows : 1;
```

### 7.3 Workspace Query

For workspace optimization, perform a workspace query:

```c
// Set up for workspace query
double dwork_query[1];
ldwork = -1;  // Query mode

// Perform query
F77_FUNC(routine, ROUTINE)(..., dwork_query, &ldwork, &info);

// Get optimal workspace size
ldwork = (int)dwork_query[0];
ldwork = MAX(ldwork, min_ldwork);  // Ensure minimum size

// Allocate workspace
dwork = (double*)malloc((size_t)ldwork * sizeof(double));
CHECK_ALLOC(dwork);
```

### 7.4 Complete Workspace Query Pattern

For routines supporting workspace query:

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

### 7.5 Alternative Approach (No Query Available)

If the routine doesn't support workspace query:

```c
// Calculate minimum workspace size from SLICOT documentation
int min_ldwork = /* formula from documentation */;

// If documentation mentions "for better performance", allocate extra space
ldwork = 2 * min_ldwork;  // Double for better performance
ldwork = MAX(min_ldwork, ldwork);  // Ensure we have at least the minimum

// Allocate workspace
dwork = (double*)malloc((size_t)ldwork * sizeof(double));
CHECK_ALLOC(dwork);
```

### 7.6 Integer Workspace (IWORK)

For integer workspaces:

```c
iwork_size = /* formula from documentation */;
iwork_size = MAX(1, iwork_size);  // Ensure at least size 1 for valid allocation
iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
CHECK_ALLOC(iwork);
```

### 7.7 Fortran Name Mangling

Use the `F77_FUNC` macro for platform-independent Fortran function name mangling:

```c
// For function name "ab01nd" that is "AB01ND" in uppercase
extern void F77_FUNC(ab01nd, AB01ND)(...);
```

This macro handles different compiler conventions:
- On Windows with Intel Fortran: expands to uppercase name (AB01ND)
- On Windows with MinGW/gfortran: expands to lowercase with underscore (ab01nd_)
- On Linux/macOS systems: expands to lowercase with underscore (ab01nd_)

## 8. Documentation Standards

### 8.1 File Headers

```c
/**
 * @file filename.c
 * @brief C wrapper implementation for SLICOT routine XXXXX
 *
 * This file provides a C wrapper implementation for the SLICOT routine XXXXX,
 * which [brief description of what the routine does].
 */
```

### 8.2 Function Documentation

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

## 9. Test Development

### 9.1 Test Organization Pattern

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
    
    // Input matrices in column-major order
    std::vector<double> A_in = {
        // Column 1
        -0.04165, -5.2100, 0.0000, 0.5450, 0.0000, 0.0000, 0.0000,
        // Column 2
        0.0000, -12.500, 3.3300, 0.0000, 0.0000, 0.0000, 0.0000,
        // ...more columns...
    };
    
    // Expected results
    double expected_value = 7.93948;
    int expected_info = 0;
};

// Row-major test fixture (inherit from column-major)
class FunctionNameTestRowMajor : public FunctionNameTestColMajor {
public:
    FunctionNameTestRowMajor() {
        ROW_MAJOR = 1;  // Override for row-major tests
    }
};
```

### 9.2 Column-Major Test Case

```cpp
TEST_F(FunctionNameTestColMajor, DocExample) {
    // Define leading dimensions for column-major format
    int LDA = N;  // >= N
    int LDB = N;  // >= N
    int LDC = P;  // >= P
    int LDD = P;  // >= P
    
    // Copy input matrices (they may be modified by the function)
    std::vector<double> A = A_in;
    std::vector<double> B = B_in;
    
    // Call wrapper function
    int info = slicot_functionname(...parameters..., ROW_MAJOR);
    
    // Verify results
    ASSERT_EQ(info, expected_info);
    
    // Check matrix elements
    for (int j = 0; j < N; ++j) { // Col
        for (int i = 0; i < N; ++i) { // Row
            EXPECT_NEAR(A[i + j*LDA], A_expected[i + j*N], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
}
```

### 9.3 Row-Major Test Case

```cpp
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
    
    // Call wrapper function
    int info = slicot_functionname(...parameters..., ROW_MAJOR);
    
    // Verify results
    ASSERT_EQ(info, expected_info);
    
    // Check matrix elements in row-major format
    for (int i = 0; i < N; ++i) { // Row
        for (int j = 0; j < N; ++j) { // Col
            EXPECT_NEAR(A[i*LDA + j], A_expected_rm[i*N + j], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
}
```

### 9.4 Required Test Types

1. **Documentation example test**: Using exact data from HTML docs or example files
2. **Input validation tests**: Error cases with invalid parameters
3. **Edge cases**: Zero dimensions, minimum/maximum values
4. **Optional**: Large matrix tests for performance

### 9.5 Inherited Test Fixture Initialization

When a row-major test fixture inherits from column-major:

- Explicitly initialize all inherited dimensions and variables
- Set combined dimensions (e.g., N = N1 + N2) in the constructor before allocating any vectors
- Always use the expected values from the documentation examples
- Avoid transposing vectors in constructors if they will be populated later in test methods

## 10. Test Data Extraction

### 10.1 Using Examples Directory for Test Data

The **examples** directory contains three important file types for creating test cases:

1. **Data Files (*.dat)**: 
   - Contain the exact input data in a well-structured format
   - Data is organized row-wise for readability

2. **Result Files (*.res)**:
   - Contain the expected output values
   - Results are presented row by row

3. **Test Programs (T*.f)**:
   - Fortran test programs showing how the input data is read
   - Demonstrates how to interpret and format the results

### 10.2 Reading Data from Example Files

For a routine like AB05ND:

1. **Identifying Matrix Dimensions**:
   - Parse the dimension line in the .dat file (e.g., `3 2 2 3` for N1, M1, P1, N2)

2. **Reading Input Matrices**:
   - Read matrices row by row from the .dat file
   - Convert to column-major format for test fixtures

3. **Reading Expected Results**:
   - Parse the result matrix from the .res file
   - Convert to column-major format for comparison

### 10.3 Converting Between Storage Formats

Convert row-major (as presented in documentation) to column-major (for Fortran):

```cpp
// For a matrix presented row-by-row in documentation:
//   1.0   0.0  -1.0
//   0.0  -1.0   1.0
//   1.0   1.0   2.0

// Store in column-major format for tests:
std::vector<double> A = {
    // Column 1
    1.0, 0.0, 1.0,
    // Column 2
    0.0, -1.0, 1.0,
    // Column 3
    -1.0, 1.0, 2.0
};
```

### 10.4 Important Notes on Matrix Representation

1. **Verify Matrix Dimensions**: Always check that your parsed data has the correct dimensions.

2. **Row-wise vs. Column-wise Presentation**: 
   - In documentation, data is typically presented row-wise for readability
   - For test fixtures, convert to column-major format for Fortran compatibility
   - For row-major tests, convert from column-major to row-major

3. **Inconsistencies in Data Presentation**:
   - Sometimes matrices may be stored transposed in the documentation
   - If dimensions don't match expectations, check the Fortran READ statements

4. **Handling Edge Cases**:
   - For zero-dimension cases, ensure vectors are properly sized
   - For non-square matrices, be especially careful with row/column ordering
   - Verify leading dimension parameters are handled correctly

## 11. Common Error Cases

### 11.1 Input Validation Errors

| Error Type | Check | Code |
|------------|-------|------|
| Negative dimensions | `if (n < 0)` | `info = -N` |
| Invalid option | `if (toupper(job) != 'B' && toupper(job) != 'F')` | `info = -N` |
| Small leading dimension | `if (lda < required_min)` | `info = -N` |
| NULL pointer for required array | `if (a == NULL)` | `info = -N` |

### 11.2 Runtime Errors

| Error Type | Fortran Info | Handling |
|------------|--------------|----------|
| Memory allocation | `SLICOT_MEMORY_ERROR` | `-1010` |
| Workspace too small | Typically positive | Allocate larger workspace |
| Algorithmic singularity | Typically positive | May be acceptable warning |

## 12. Utility Functions

### 12.1 Matrix Dimension Handling

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

### 12.2 Identity Matrix Initialization

```c
set_identity(n, matrix, ld, row_major);  // row_major=0 for column-major
```

### 12.3 Symmetric Matrix Utilities

```c
// Copy only specified triangle (upper/lower) from a symmetric matrix to a full matrix
slicot_copy_symmetric_part(src, dest, n, 'U', ld, elem_size);
```

### 12.4 Complex Number Handling

```c
// Define complex variables
slicot_complex_double z;

// Access real part of complex number
double real_part = SLICOT_COMPLEX_REAL(z);
```

### 12.5 DLL Export/Import Macros

```c
// In the implementation file
SLICOT_C_WRAPPER_API
void my_exported_function(...)
```

The `SLICOT_C_WRAPPER_API` macro handles platform-specific export declarations:
- On Windows: `__declspec(dllexport)` when building, `__declspec(dllimport)` when using
- On Linux/macOS with GCC/Clang: `__attribute__((visibility("default")))`
- When building static library: empty definition

### 12.6 In-place Transposition

```c
// Returns 0 on success, -1 if matrix is not square
int result = slicot_transpose_inplace(matrix, n, n, sizeof(double));
```

## 13. LAPACK Dependencies

### 13.1 Declaring LAPACK Function Dependencies

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

### 13.2 Identifying LAPACK Dependencies

Look for LAPACK dependencies in:
1. The "Further Comments" section of SLICOT documentation
2. Places where the documentation mentions forming matrices from factorized results
3. Statements like "can be obtained by calling the LAPACK library routine XXX"

Common LAPACK dependencies include:
- `DORGQR`/`ZUNGQR`: Generate orthogonal/unitary matrix from elementary reflectors
- `DORMQR`/`ZUNMQR`: Apply orthogonal/unitary matrix from QR factorization
- `DGESVD`/`ZGESVD`: Singular value decomposition
- `DGEQRF`/`ZGEQRF`: QR factorization

### 13.3 Using LAPACK Functions in Wrappers

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

## 14. Troubleshooting Guide

### 14.1 Matrix Dimensions Mismatch

- Double-check row vs column interpretation
- Verify leading dimension handling
- Check if transposition is needed

### 14.2 Memory Errors

- Verify zero dimension handling
- Ensure all memory is freed in cleanup block
- Check allocation size calculations
- Watch for vector::_M_default_append errors, indicating access beyond allocated memory

### 14.3 Test Failures

- Compare with documentation examples
- Check row-major/column-major conversion
- Verify matrix element indexing
- Consider numerical precision issues
- Try increasing tolerance for numerical functions (e.g., check_tol = 3.5 vs 1e-4)

### 14.4 Parameter Validation

- Review error code assignments
- Check if parameter validation matches Fortran routine
- Ensure character parameters are converted to uppercase

### 14.5 Fixture Initialization Issues

- When fixtures inherit from other fixtures, ensure parent variables are properly initialized
- Be careful with variables set in SetUp() that are needed in constructors
- For edge cases with zero dimensions, ensure proper allocation of dummy arrays

## 15. Code Templates

### 15.1 Standard Wrapper Template

```c
/**
 * @file function_name.c
 * @brief C wrapper for SLICOT routine FUNCTION_NAME
 *
 * [Description of what the routine does]
 */

#include <stdlib.h>
#include <ctype.h>
#include <stddef.h>

#include "function_name.h"
#include "slicot_utils.h"
#include "slicot_f77.h"

// Declare external Fortran function
extern void F77_FUNC(function_name, FUNCTION_NAME)(
    // Fortran function parameters with types
    // ...
);

/* C wrapper function definition */
SLICOT_C_WRAPPER_API
int slicot_function_name(/* C function parameters */)
{
    // 1. Variable declarations
    int info = 0;
    
    // 2. Input parameter validation
    if (n < 0) { info = -2; goto cleanup; }
    if (toupper(job) != 'B' && toupper(job) != 'F') { info = -1; goto cleanup; }
    
    // 3. Memory allocation for column-major copies
    size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;
    double* a_cm = NULL;
    if (a_size > 0) { 
        a_cm = (double*)malloc(a_size * sizeof(double)); 
        CHECK_ALLOC(a_cm); 
    }
    
    // 4. Row-major to column-major conversion
    if (row_major && a_size > 0) {
        slicot_transpose_to_fortran(a, a_cm, n, n, sizeof(double));
    }
    
    // 5. Workspace allocation
    int ldwork = MAX(1, 3*n);
    double* dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);
    
    // 6. Prepare Fortran parameters
    char job_upper = toupper(job);
    const int job_len = 1;
    int lda_f = (n > 0) ? n : 1;
    
    // 7. Call Fortran function
    F77_FUNC(function_name, FUNCTION_NAME)(
        &job_upper, &n,
        (row_major && a_size > 0) ? a_cm : a,
        &lda_f,
        dwork, &ldwork, &info,
        job_len
    );
    
    // 8. Convert results back to row-major
    if (row_major && a_size > 0 && info == 0) {
        slicot_transpose_to_c(a_cm, a, n, n, sizeof(double));
    }
    
cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(a_cm);
    return info;
}
```

### 15.2 Standard Test Template

```cpp
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "function_name.h"

// Column-major test fixture
class FunctionNameTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    int N = 3;
    char PARAM = 'X';
    int ROW_MAJOR = 0;  // Column-major
    double check_tol = 1e-4;
    
    // Input matrices from documentation (column-major order)
    std::vector<double> A_in = {
        /* Column 1 */ 1.0, 0.0, 1.0,
        /* Column 2 */ 0.0, -1.0, 1.0,
        /* Column 3 */ -1.0, 1.0, 2.0
    };
    
    // Expected outputs
    std::vector<double> A_expected = {
        /* Expected output matrix in column-major order */
    };
    int expected_info = 0;
};

// Row-major test fixture
class FunctionNameTestRowMajor : public FunctionNameTestColMajor {
public:
    FunctionNameTestRowMajor() {
        ROW_MAJOR = 1;  // Row-major
        
        // Convert test data to row-major
        A_rm.resize(N * N);
        A_expected_rm.resize(N * N);
        
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                A_rm[i*N + j] = A_in[i + j*N];
                A_expected_rm[i*N + j] = A_expected[i + j*N];
            }
        }
    }
    
    std::vector<double> A_rm, A_expected_rm;
};

// Column-major test
TEST_F(FunctionNameTestColMajor, DocExample) {
    // Set up
    int LDA = N;
    std::vector<double> A = A_in;
    
    // Call function
    int info = slicot_function_name(PARAM, N, A.data(), LDA, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(A[i + j*LDA], A_expected[i + j*N], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
}

// Row-major test
TEST_F(FunctionNameTestRowMajor, DocExample) {
    // Set up
    int LDA = N;
    std::vector<double> A = A_rm;
    
    // Call function
    int info = slicot_function_name(PARAM, N, A.data(), LDA, ROW_MAJOR);
    
    // Verify
    ASSERT_EQ(info, expected_info);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            EXPECT_NEAR(A[i*LDA + j], A_expected_rm[i*N + j], check_tol)
                << "A[" << i << "," << j << "] mismatch";
        }
    }
}
```
