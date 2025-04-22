# SLICOT_C Contribution Guidelines for AI Coding Agents

This document is structured specifically to help AI coding agents efficiently implement C wrappers and tests for the SLICOT library.

## 1. Project Structure

```
SLICOT-Reference/
├── doc/                     # SLICOT function documentation (HTML)
├── include/                 # Public header files (*.h)
├── src/                     # Original SLICOT Fortran source code
├── src_c_wrapper/           # C wrapper implementation files (*.c)
└── tests/                   # Test files (*_test.cpp)
```

## 2. Implementation Workflow

Follow this sequence when implementing a new wrapper:

1. Analyze the HTML documentation (e.g., `doc/AB05OD.html`)
2. Extract key information: parameters, dimensions, workspace requirements
3. Identify similar existing wrappers (same function family)
4. Create C wrapper following the template pattern
5. Create test cases using example data from documentation

## 3. C Wrapper Template

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

## 4. Test Case Template

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

## 5. Key Patterns for AI Implementation

### 5.1 Memory Allocation

```c
// 1. Calculate size with zero dimension handling
size_t a_size = (size_t)m * n; if (m == 0 || n == 0) a_size = 0;

// 2. Allocate memory if needed
if (a_size > 0) {
    a_cm = (double*)malloc(a_size * sizeof(double));
    CHECK_ALLOC(a_cm);
}

// 3. Use conditional allocation checks in processing
if (a_size > 0 && row_major) {
    slicot_transpose_to_fortran(a, a_cm, m, n, sizeof(double));
}
```

### 5.2 Parameter Validation

```c
// Check dimensions
if (n < 0) { info = -2; goto cleanup; }
if (m < 0) { info = -3; goto cleanup; }

// Check option parameters
if (char_dico != 'C' && char_dico != 'D') { info = -1; goto cleanup; }

// Check leading dimensions
if (row_major) {
    if (lda < n) { info = -6; goto cleanup; }  // Columns in row-major
} else {
    if (lda < m) { info = -6; goto cleanup; }  // Rows in column-major
}
```

### 5.3 Matrix Format Conversion

```c
// Row-major to column-major (before Fortran call)
if (row_major && a_size > 0) {
    slicot_transpose_to_fortran(a, a_cm, m, n, sizeof(double));
}

// Column-major to row-major (after Fortran call)
if (row_major && a_size > 0 && info == 0) {
    slicot_transpose_to_c(a_cm, a, m, n, sizeof(double));
}

// With different leading dimensions
if (row_major && a_size > 0 && info == 0) {
    slicot_transpose_to_c_with_ld(a_cm, a, m, n, m, lda, sizeof(double));
}
```

### 5.4 Workspace Handling

```c
// Basic workspace allocation
int ldwork = MAX(1, 3*n);
double* dwork = (double*)malloc((size_t)ldwork * sizeof(double));
CHECK_ALLOC(dwork);

// Optimal workspace query
double dwork_query;
ldwork = -1;  // Query mode
F77_FUNC(routine, ROUTINE)(..., &dwork_query, &ldwork, &info, ...);
ldwork = (int)dwork_query;
ldwork = MAX(ldwork, min_ldwork);  // Ensure minimum size
dwork = (double*)malloc((size_t)ldwork * sizeof(double));
CHECK_ALLOC(dwork);
```

### 5.5 Character Parameter Handling

```c
// Convert to uppercase and handle hidden length
char dico_upper = toupper(dico);
const int dico_len = 1;  // Hidden length for Fortran strings

// Pass in the Fortran call
F77_FUNC(routine, ROUTINE)(&dico_upper, ..., dico_len);
```

### 5.6 Cleanup Pattern

```c
cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(iwork);
    free(a_cm); free(b_cm); free(c_cm); free(d_cm);
    return info;
```

## 6. Data Extraction from HTML Documentation

### 6.1 Locating Example Data

1. Find the "Program Data" section in the HTML file
2. Extract matrices in row-major textual format
3. Convert to column-major for test fixtures
4. Find the "Program Results" section for expected outputs

### 6.2 Matrix Data Representation

For column-major test data:
```cpp
// Matrix data from documentation (converted to column-major)
std::vector<double> A_in = {
    // Column 1
    1.0, 0.0, 1.0,
    // Column 2
    0.0, -1.0, 1.0,
    // Column 3
    -1.0, 1.0, 2.0
};
```

For row-major conversions in tests:
```cpp
// Convert matrix from column-major to row-major
for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
        A_rm[i*N + j] = A_in[i + j*N];
    }
}
```

## 7. Common Error Cases to Handle

### 7.1 Input Validation Errors

| Error Type | Check | Code |
|------------|-------|------|
| Negative dimensions | `if (n < 0)` | `info = -N` |
| Invalid option | `if (toupper(job) != 'B' && toupper(job) != 'F')` | `info = -N` |
| Small leading dimension | `if (lda < required_min)` | `info = -N` |
| NULL pointer for required array | `if (a == NULL)` | `info = -N` |

### 7.2 Runtime Errors

| Error Type | Fortran Info | Handling |
|------------|--------------|----------|
| Memory allocation | `SLICOT_MEMORY_ERROR` | `-1010` |
| Workspace too small | Typically positive | Allocate larger workspace |
| Algorithmic singularity | Typically positive | May be acceptable warning |

## 8. Test Case Structure Guidelines

1. **Base Test Fixtures**:
   - One fixture for column-major tests (original Fortran format)
   - One derived fixture for row-major tests with conversion methods

2. **Required Test Types**:
   - Documentation example test (exact data from HTML docs)
   - Input validation tests (error cases)
   - Edge cases (zero dimensions)
   - Optional: Large matrix tests

3. **Verification Pattern**:
   - First check return code is as expected
   - Then verify matrix elements against expected values
   - Use appropriate tolerance for floating point comparisons

## 9. AI-Specific Implementation Tips

1. **Analyze similar wrappers first**: 
   - Look for wrappers in the same family (e.g., AB05xx functions)
   - Copy patterns for parameter validation, workspace handling

2. **Reuse existing structures**: 
   - Fortran call declarations
   - Parameter validation sequences
   - Matrix conversion patterns

3. **Generate complete tests**:
   - Implement both column-major and row-major tests
   - Test all significant code paths
   - Ensure boundary conditions are tested

4. **Documentation priority**:
   - Focus on parameter descriptions in comments
   - Clearly document matrix dimension requirements
   - Note any special cases or error conditions

5. **Log your reasoning**:
   - Document why specific implementation choices were made
   - Note any deviations from the usual patterns
   - Reference similar wrappers that influenced implementation

## 10. Troubleshooting Guide

1. **Matrix dimensions mismatch**:
   - Double-check row vs column interpretation
   - Verify leading dimension handling
   - Check if transposition is needed

2. **Memory errors**:
   - Verify zero dimension handling
   - Ensure all memory is freed in cleanup block
   - Check allocation size calculations

3. **Test failures**:
   - Compare with documentation examples
   - Check row-major/column-major conversion
   - Verify matrix element indexing
   - Consider numerical precision issues

4. **Parameter validation**:
   - Review error code assignments
   - Check if parameter validation matches Fortran routine
   - Ensure character parameters are converted to uppercase

By following these structured guidelines, AI coding agents can efficiently implement consistent, high-quality C wrappers and tests for the SLICOT library.
