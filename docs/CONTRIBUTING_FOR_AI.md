# SLICOT_C Contribution Guidelines for AI Coding Agents

This document is structured specifically to help AI coding agents efficiently implement C wrappers and tests for the SLICOT library.

## 1. Project Structure

```
SLICOT-Reference/
├── doc/                     # SLICOT function documentation (HTML)
├── examples/                # Example data files (*.dat, *.res) and Fortran example programs (*.f)
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
SLICOT_C_WRAPPER_API // Use this macro for DLL export/import handling
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

## 6. Finding and Extracting Test Data

### 6.1 Hierarchical Approach to Finding Test Data

Follow this order when searching for test data sources:

1. **First Check the Examples Directory**
   - Look for matching *.dat and *.res files in the examples directory:
     ```
     examples/[FUNCTION_NAME].dat  # Input test data
     examples/[FUNCTION_NAME].res  # Expected output
     examples/T[FUNCTION_NAME].f   # Fortran test program
     ```
   - Example: For function AB07ND, look for AB07ND.dat, AB07ND.res, and TAB07ND.f

2. **If No Related Examples, Check HTML Documentation**
   - Look in the doc directory for the function's HTML documentation:
     ```
     doc/[FUNCTION_NAME].html  # Function documentation
     ```
   - Parse examples in HTML documentation, paying attention to:
     - Example matrices provided in formatted tables
     - Input parameters and their descriptions
     - Expected outputs described in the text
   - Extract numerical values directly from code samples or examples


3. **If No Examples in Documentation, Look at Fortran Source**
   - Examine the original Fortran source code in the src directory:
     ```
     src/[FUNCTION_NAME].f  # Original Fortran implementation
     ```
   - Extract test cases from:
     - Documentation comments at the beginning of the file
     - Parameter validation and range checks
     - Special case handling for specific values

4. **For Last Resort: Generate Synthetic Test Data**
   - Create composite test data by using inputs from related functions:
     - Identify functions with similar parameter structures (e.g., AB05MD for AB05ND)
     - Extract relevant matrices of appropriate dimensions from their example data
     - Adapt matrices to match the specific requirements of your target function
   
   - Create matrices with specific mathematical properties relevant to the function:
     - For state-space functions: controllable/observable canonical forms
     - For stability functions: matrices with known eigenvalue distributions
     - For filtering functions: stable system matrices with known frequency responses
     - For decomposition functions: matrices with specific rank or condition number

   - When adapting matrices from related examples:
     - Preserve structural properties relevant to your function
     - Ensure dimensions are compatible with your function's requirements
     - Adjust values to ensure mathematical validity (e.g., matrix invertibility if required)
     - Scale values appropriately to avoid numerical issues

   - Run progressive validation:
     - Start with simple cases (identity, diagonal, or block matrices)
     - Gradually introduce more complex structures
     - Use the function with synthetic inputs in a "learning mode"
     - Capture the actual output for future test validation
     - Cross-verify with related functions where mathematical relationships exist

   - Document your synthetic test approach:
     - Clearly note the origin of any adapted matrices
     - Document any mathematical properties designed into test matrices
     - Note which aspects of function behavior are being tested
     - Clearly indicate these tests validate consistency rather than correctness against known standard results
     - Include formulas or theoretical expectations where applicable

### 6.2 Using Examples Directory for Test Data

The **examples** directory contains important files that should be used as the primary source for creating test cases:

1. **Data Files (*.dat)**: 
   - Contain the exact input data for each function in a well-structured format
   - Each value is precisely as required by the function
   - Data is organized row-wise and can be directly parsed

2. **Result Files (*.res)**:
   - Contain the expected output values from running the function
   - These are the exact values you should verify against in tests
   - Results are presented row by row in a clear format

3. **Test Programs (T*.f)**:
   - Fortran test programs showing how the input data is read
   - Shows how the SLICOT function is called with the inputs
   - Demonstrates how to interpret and format the results

Using these example files gives you a more direct and accurate source for test data than parsing the HTML documentation.

### 6.2 Parsing SLICOT Example Data

SLICOT routines are written in Fortran, which stores multi-dimensional arrays in **column-major order**. This means elements of a column are stored contiguously in memory. However, standard mathematical notation and many other programming languages (like C/C++) typically use **row-major order**, where elements of a row are contiguous.

The input and output data in the example *.dat and *.res files is typically presented **row-wise for readability**.

#### Key Steps for Using Example Files

When creating test cases:

1. **Locate the Example Files**:
   - Look for the appropriate *.dat file (e.g., `AB05ND.dat`) for input data
   - Look for the corresponding *.res file (e.g., `AB05ND.res`) for expected outputs
   - Review the Fortran test program (e.g., `TAB05ND.f`) to understand data processing

2. **Parse the Input Data**:
   - Read the values in the *.dat file sequentially
   - Values are organized by matrix, row by row
   - The first lines usually contain dimension parameters (N, M, P, etc.)

3. **Parse the Expected Output**:
   - Read the values in the *.res file sequentially for each output matrix
   - These are the expected values your test should verify against

#### Example: Reading from *.dat and *.res Files

Consider AB05ND as an example:

**1. Reading the input data (AB05ND.dat):**

```
AB05ND EXAMPLE PROGRAM DATA
   3     2     2     3
   1.0   0.0  -1.0
   0.0  -1.0   1.0
   1.0   1.0   2.0
   1.0   1.0   0.0
   2.0   0.0   1.0
   3.0  -2.0   1.0
   0.0   1.0   0.0
   1.0   0.0
   0.0   1.0
  -3.0   0.0   0.0
   1.0   0.0   1.0
   0.0  -1.0   2.0
   0.0  -1.0   0.0
   1.0   0.0   2.0
   1.0   1.0   0.0
   1.0   1.0  -1.0
   1.0   1.0
   0.0   1.0
```

The first line is a header.
The second line contains key dimensions: N1=3, M1=2, P1=2, N2=3.
Then we have matrices in this order (row by row):
- A1 (3×3): 3 rows of 3 values each
- B1 (3×2): 3 rows of 2 values each (note the data layout!)
- C1 (2×3): 2 rows of 3 values each
- D1 (2×2): 2 rows of 2 values each
...and so on for A2, B2, C2, D2

**2. Reading the expected output (AB05ND.res):**

```
AB05ND EXAMPLE PROGRAM RESULTS

 The state transition matrix of the connected system is
  -0.5000  -0.2500  -1.5000  -1.2500  -1.2500   0.7500
  -1.5000  -0.2500   0.5000  -0.2500  -0.2500  -0.2500
   1.0000   0.5000   2.0000  -0.5000  -0.5000   0.5000
   0.0000   0.5000   0.0000  -3.5000  -0.5000   0.5000
  -1.5000   1.2500  -0.5000   1.2500   0.2500   1.2500
   0.0000   1.0000   0.0000  -1.0000  -2.0000   3.0000
```

This shows the expected matrix A, presented row by row, with dimensions 6×6.

**3. Converting for test cases:**

For column-major format in C/C++ tests:

```cpp
// Matrix A1 data from AB05ND.dat (converted to column-major)
std::vector<double> A1 = {
    // Column 1
    1.0, 0.0, 1.0,
    // Column 2
    0.0, -1.0, 1.0,
    // Column 3
    -1.0, 1.0, 2.0
};

// Expected output matrix A from AB05ND.res (converted to column-major)
std::vector<double> A_expected = {
    // Column 1
    -0.5000, -1.5000, 1.0000, 0.0000, -1.5000, 0.0000,
    // Column 2
    -0.2500, -0.2500, 0.5000, 0.5000, 1.2500, 1.0000,
    // And so on...
};
```

For row-major tests, transpose from column-major to row-major:

```cpp
// Convert matrix from column-major to row-major
for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
        A_rm[i*N + j] = A1[i + j*N];
    }
}
```

### 6.2 Important Notes on Matrix Representation

1. **Verify Matrix Dimensions**: Always check that your parsed data has the correct dimensions. For a matrix of size MxN, you should have M*N elements.

2. **Row-wise vs. Column-wise Presentation**: 
   - In the HTML documentation, data is typically presented row-wise
   - For test fixtures, you need to convert this to column-major format
   - For row-major tests, you'll need to convert from column-major to row-major

3. **Inconsistencies in Data Presentation**: 
   - Sometimes, the data presentation might be inconsistent
   - In some examples, matrices (especially B matrices) may be stored transposed in the documentation
   - If the dimensions don't match what you expect, double-check the Fortran READ statements in the example code
   - Cross-reference with the Fortran code patterns like `((B(I,J), J=1,M), I=1,N)` which reads row-wise

4. **Handling Edge Cases**:
   - For zero-dimension cases, ensure your vectors are properly sized (potentially to 0)
   - For non-square matrices, be especially careful with row/column ordering
   - Verify that leading dimension parameters are handled correctly

### 6.3 Testing Zero-Dimension Cases

When testing SLICOT functions with zero-dimension inputs (e.g., N=0, M=0, P=0), extra attention is needed:

1. **Leading Dimensions for Zero-Dimension Matrices**:
   - Always follow the function's minimum requirements for leading dimensions
   - For column-major format, typically:
     - State matrices (A): `LDA >= max(1, N)`
     - Input matrices (B): `LDB >= max(1, N)` 
     - Output matrices (C): `LDC >= max(1, P)`
     - Feedthrough matrices (D): `LDD >= max(1, P)`
   - Even for zero-sized matrices, these requirements must be met

2. **Array Allocation**:
   - Allocate at least one element (never use nullptr) even for zero-sized matrices
   - For combined systems (e.g., N = N1 + N2), ensure output arrays are large enough
   ```cpp
   std::vector<double> A_zero(MAX(1, N1+N2), 0.0);
   ```

3. **Matrix Format Consistency**:
   - For zero-dimension tests, explicitly set matrix format (e.g., 0 for column-major)
   - Don't rely on class variables for this critical parameter
   - Use explicit constants to make intent clear:
   ```cpp
   slicot_function('N', n1_zero, m, p, n2,
                   /* matrices... */,
                   0); // Explicitly use column-major (0)
   ```

4. **Special Handling for Combined Systems**:
   - For functions dealing with combined systems (like AB05PD where N = N1 + N2):
     - When N1=0, set output leading dimensions based on N2
     - When N2=0, set output leading dimensions based on N1
     - Be explicit about each case rather than using general variables

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

4. **Inherited Test Fixture Initialization**:
   - When a row-major test fixture inherits from column-major, explicitly initialize all inherited dimensions and variables
   - Set combined dimensions (e.g., N = N1 + N2) in the constructor before allocating any vectors
   - Always use the expected values from the documentation examples, not a "capture and verify" approach
   - Avoid transposing vectors in constructors if they will be populated later in test methods

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

6. **Use DLL Export Macro**:
   - Always prefix the C wrapper function definition with `SLICOT_C_WRAPPER_API`.
   - This macro handles platform-specific details for exporting functions from shared libraries (DLLs on Windows). See `CONTRIBUTING.md` Section 12.5 for details.

## 10. Troubleshooting Guide

1. **Matrix dimensions mismatch**:
   - Double-check row vs column interpretation
   - Verify leading dimension handling
   - Check if transposition is needed

2. **Memory errors**:
   - Verify zero dimension handling
   - Ensure all memory is freed in cleanup block
   - Check allocation size calculations
   - Watch for vector::_M_default_append errors, indicating access beyond allocated memory

3. **Test failures**:
   - Compare with documentation examples
   - Check row-major/column-major conversion
   - Verify matrix element indexing
   - Consider numerical precision issues
   - Try increasing tolerance for numerical functions (e.g., check_tol = 3.5 vs 1e-4)
   - Verify the right test pattern: sometimes exact expected values won't match implementation

4. **Parameter validation**:
   - Review error code assignments
   - Check if parameter validation matches Fortran routine
   - Ensure character parameters are converted to uppercase

5. **Fixture initialization issues**:
   - When fixtures inherit from other fixtures, ensure parent variables are properly initialized
   - Always initialize combined dimensions (e.g., `N = N1 + N2`) at the beginning of the derived class constructor before using them for vector resizing to avoid "vector::_M_default_append" errors
   - Be careful with variables set in SetUp() that are needed in constructors - these won't be available during constructor execution
   - For edge cases with zero dimensions, ensure proper allocation of dummy arrays with required dimensions

6. **Test design patterns**:
   - Always use expected values from documentation rather than a "capture and verify" approach
   - Parse example data from HTML documentation carefully according to the Readme.md guidance
   - For tests with singular matrices or potential instability, adapt expected error codes as documented
   - For edge cases where N1=0 or N2=0, always initialize all necessary matrices with minimum valid dimensions

By following these structured guidelines, AI coding agents can efficiently implement consistent, high-quality C wrappers and tests for the SLICOT library.
