# SLICOT_C Contribution Guidelines for AI Coding Agents

**Version 9**

This document provides structured guidelines for AI coding agents to efficiently implement C wrappers and corresponding C++ tests for the SLICOT library, focusing on consistency, clarity, and reducing misinterpretations.

## 0. Key Principles (TL;DR for AI Agents)

### Goal
Create a C wrapper function (e.g., `slicot_function_name`) that calls the corresponding SLICOT Fortran routine (`F77_FUNC(function_name, FUNCTION_NAME)`).

### Workspace
The C wrapper **must allocate and free workspace arrays** (`iwork`, `dwork`, etc.) internally using `malloc`/`free`. Calculate sizes based on formulas in SLICOT documentation (see `rep96-1.pdf`, Sec 2.2.10). For routines that support it (check SLICOT routine documentation), query the Fortran routine for optimal workspace size first (e.g., returned in `DWORK(1)`), then fall back to documented formulas if the query fails or is not supported. **Do not expect the caller to provide workspace.**

### Row-Major Handling
The wrapper must accept a `row_major` flag. If `row_major` is true (1), the wrapper must:

1. Allocate temporary column-major buffers (`_cm`) for input/output matrices passed by the C caller
2. Transpose row-major C input matrices into the temporary column-major buffers before calling Fortran
3. Pass the column-major buffers and corresponding Fortran-style leading dimensions (number of rows) to the Fortran routine
4. If the Fortran routine modifies inputs or computes outputs, transpose results from column-major buffers back to the original C row-major arrays after the Fortran call (if `info == 0`)

### Leading Dimensions (LD)
- **Fortran (and Column-Major C):** `LDA`, `LDB`, etc., always refer to the number of **rows** allocated for the matrix in memory
- **Row-Major C:** `LDA`, `LDB`, etc., passed to the C wrapper refer to the number of **columns** allocated for the matrix in memory. The wrapper must calculate the required Fortran LDA (rows) internally

### Parameter Validation and Zero Dimensions
- Validate inputs rigorously in the C wrapper before calling Fortran
- **Crucially**, for zero-dimension cases (e.g., `N=0`, `M=0`), the C wrapper's validation must align with the underlying SLICOT Fortran routine's behavior. Many SLICOT routines handle zero dimensions as valid inputs, often performing a quick exit with `INFO=0` (see `rep96-1.pdf`, Sec 2.3.10). The C wrapper should not erroneously flag such cases as errors if the Fortran routine would accept them
- Pay close attention to the SLICOT documentation for each routine to understand if `NULL` pointers are acceptable for arrays that become zero-sized (e.g., an L-by-N matrix where N=0). If `NULL` is passed, ensure corresponding Fortran leading dimensions are still set to valid values (e.g., 1) if the Fortran interface expects them, even if the pointer itself is `NULL`

### Testing
- Use GTest fixtures (`ColMajor`, `RowMajor`)
- For small datasets (e.g., ~10 samples), embed data directly in the test fixture
- For larger datasets, create a CSV file in the `tests/data/` directory and use the `load_test_data_from_csv` utility
- **CSV Data Loading**: Use the correct function signature:
  ```cpp
  bool load_test_data_from_csv(
      const std::string& filepath,
      const std::vector<std::string>& input_cols,
      const std::vector<std::string>& output_cols,
      std::vector<double>& u,
      std::vector<double>& y,
      int& num_samples);
  ```
- **CSV Data Format**: The `load_test_data_from_csv` utility reads the specified columns (e.g., "U1", "U2") and loads them into separate vectors
- **Data Rearrangement**: For multi-column time-series data loaded via CSV, the test fixture must handle the loaded vectors appropriately for the expected matrix storage format
- The CSV header names **MUST** exactly match the `input_columns`/`output_columns` specified in the test fixture
- `SetUp` loads/defines data, performs necessary rearrangements, updates `NSMP` (number of samples), calculates LDs, and sizes output vectors
- Tests call the C wrapper, `ASSERT_EQ` the info code, and `EXPECT_NEAR` the numerical results
- **Numerical Tolerance**: Use realistic tolerances for numerical comparisons. Start with `1e-3` but adjust to `5e-3` or higher if the algorithm shows natural numerical variations. The computed results should be very close to expected values but may have small differences due to compiler optimizations, different numerical libraries, or natural precision variations in iterative algorithms
- **Verify Expected Results**: Double-check expected results (`.res` files, documentation examples, original Fortran example programs `examples/*.f` as per `rep96-1.pdf` Ch. 4 & App. C) against independent simulations (e.g., using Python control libraries) if possible, as documentation examples may contain errors. Use the verified results in the test fixture. **Important**: If computed results are consistently close to each other but slightly different from documentation values, trust the computed results and adjust expected values accordingly
- Include specific test cases for zero-dimension inputs to ensure correct handling by the wrapper and alignment with Fortran routine behavior

### Additional Requirements
- **Filenames**: All C source and header files must use lowercase naming (e.g., `ab05od.c`, not `AB05OD.c`)
- **Error Handling**: Use `CHECK_ALLOC` after `malloc`. Implement `goto cleanup` for error exits. Free all allocated memory in the cleanup block. Return appropriate info codes (negative for wrapper errors, Fortran info otherwise)

## 1. Project Structure

```
SLICOT_C/
├── benchmark_data/   # Performance benchmarking data
├── build/            # Build artifacts (generated)
├── cmake/            # CMake build configuration files
├── doc/              # SLICOT function documentation (HTML)
├── docs/             # Project documentation
├── examples/         # Original example data (*.dat, *.res) and Fortran programs (*.f)
├── include/          # Public C header files for wrappers (*.h, lowercase naming)
├── src/              # Original SLICOT Fortran source code
├── src_aux/          # Auxiliary source files
├── src_c_wrapper/    # C wrapper implementation files (*.c, lowercase naming)
└── tests/            # C++ Test files (*_test.cpp, test_utils.h/cpp)
    └── data/          # Test data files (*.csv) derived from examples
```

## 2. Implementation Workflow

Follow this sequence precisely:

### 2.1 Analyze Documentation (Crucial for Edge Cases)
- Review the function's HTML documentation (e.g., `doc/AB05OD.html`)
- Study the corresponding Fortran example files (`examples/*.dat`, `examples/*.res`, `examples/T*.f`)
- Consult the general SLICOT standards (`rep96-1.pdf`), especially sections on:
  - User Interface Standards (Sec 2.2), particularly problem dimensions (2.2.6) and leading dimensions (2.2.7)
  - Programming Standards (Sec 2.3), particularly error checking for input parameters and zero dimensions (2.3.10)
  - Workspace formulas and query mechanisms (Sec 2.2.10)
- Pay extremely close attention to how the specific Fortran routine is documented to handle zero-dimension inputs (e.g., `N=0`, `M=0`). Understand if it's a valid quick-exit case (often `INFO=0`) or an error. This dictates the C wrapper's validation logic
- Determine if `NULL` pointers are acceptable for arrays that become effectively zero-sized due to other dimension parameters

### 2.2 Prepare/Verify Test Data
**Option A (For larger datasets):**
- Locate the `.dat` file in `examples/` (and the corresponding `.f` example program and `.res` results file). These are primary sources for understanding data structure and expected outputs
- Create a corresponding CSV file in `tests/data/`. Use lowercase for the filename
- The first row **MUST** be a header. Choose clear, unique names
- Copy numerical data from `.dat` into subsequent CSV rows

**Option B (Preferred for small datasets ~10 samples):** 
Define test data directly in the C++ test fixture.

**Verification:** Crucially, verify expected output data (from `.res` file or documentation) against an independent source (e.g., Python simulation). If discrepancies arise, trust verified simulation results and document this.

### 2.3 Implementation Steps
1. **Extract Key Information**: Note parameters, types, dimensions, LD rules, workspace formulas, and specific behaviors for zero-dimension inputs from the Fortran documentation
2. **Identify Similar Wrappers**: Find wrappers in `src_c_wrapper/` for similar routines as a reference, but always prioritize the specific documentation for the target routine
3. **Create C Wrapper**: Implement the `.c` file using the Internal Workspace Allocation template (Section 3)
   - Ensure correct `row_major` handling
   - Calculate workspace based on formulas/queries
   - Implement C-level parameter validation that is consistent with, and not overly restrictive compared to, the documented behavior of the Fortran routine, especially for zero-dimension cases
4. **Create Header File**: Implement the `.h` file in `include/` using the template (Section 3.1). Document parameters clearly, especially `row_major` handling and conditions for `NULL` pointers
5. **Create Test Cases**: Implement the `_test.cpp` file using the GTest template (Section 4)
   - If using CSV, define `input_columns` and `output_columns` to exactly match CSV headers. Use `load_test_data_from_csv`
   - **Data Rearrangement**: For multi-column time-series data from CSV, explicitly rearrange to Fortran column-major order in `SetUp`
   - If embedding data, define it preferably in Fortran column-major order
   - Implement tests for column-major, row-major, and parameter validation
   - **Crucially**, add specific test cases for various zero-dimension scenarios (e.g., `N=0`, `M=0`, `NSMP=0`, combinations thereof) to verify the wrapper's validation logic and alignment with expected Fortran behavior

## 3. C Wrapper Template (Internal Workspace Allocation)

```c
/**
 * @file function_name.c
 * @brief C wrapper for SLICOT routine FUNCTION_NAME.
 * @details [Description of what the routine does, its purpose, and key algorithms.]
 * Workspace (IWORK, DWORK, etc.) is allocated internally by this wrapper.
 * Input/output matrix format is handled via the row_major parameter.
 * C-level validation for zero dimensions aims to align with Fortran routine behavior.
 */

#include <stdlib.h> // For malloc, free
#include <ctype.h>  // For toupper
#include <stddef.h> // For size_t
#include <math.h>   // For MAX/MIN if needed (often provided by slicot_utils.h)
#include <stdio.h>  // For error logging (optional)

#include "function_name.h" // Public header for this wrapper
#include "slicot_utils.h"  // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX/MIN, transpose functions etc.
#include "slicot_f77.h"    // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */
extern void F77_FUNC(function_name, FUNCTION_NAME)(
    /* Fortran function parameters with const for inputs */
    // const char* job, const int* n, ...,
    /* Workspace arrays passed to Fortran */
    int* iwork, double* dwork, const int* ldwork, // Example for DWORK, add others if needed
    int* info
    /* Hidden string lengths if any */
    // int job_len
);

/* C wrapper function definition */
SLICOT_EXPORT // Macro for DLL export/import handling
int slicot_function_name(/* C function parameters, excluding workspace */
                         // const char* job, int n, ..., double* a, int lda, ..., int row_major
                        )
{
    // 1. Variable declarations
    int info = 0;
    int local_iwarn = 0; // If routine has IWARN
    char job_upper; // Example for char parameter

    int *iwork = NULL;
    double *dwork = NULL;
    int liwork = 0;
    int ldwork_alloc = 0; // Renamed from ldwork to avoid conflict with Fortran arg name

    double* a_cm = NULL;
    // ... declare other _cm pointers as needed ...

    const double* a_ptr = a; // Pointer to pass to Fortran
    // ... declare other _ptr pointers ...

    int lda_f, ldb_f; // Fortran leading dimensions
    // ... declare other _f leading dimensions ...

    size_t a_size = 0; // Size in elements for allocation
    // ... declare other _size variables ...

    const int job_len = 1; // Example for char job

    // 2. Input parameter validation
    // **Consult SLICOT docs for each parameter's constraints and behavior with zero dimensions.**
    job_upper = toupper(*job); // Example: dereference if char* job
    if (job_upper != 'B' && job_upper != 'F') { info = -1 /* Map to job's arg index */; goto cleanup; }
    if (n < 0) { info = -2 /* Map to n's arg index */; goto cleanup; }
    // ... other scalar validations ...

    // Example: Matrix A (N x N_COLS_A)
    // If N=0, A might not be referenced by Fortran, or a NULL is fine.
    // If N>0, A is usually required.
    if (n > 0) { // A is relevant
        if (a == NULL) { info = -arg_idx_A; goto cleanup; }
        if (row_major) { // C LDA is columns
            if (lda < N_COLS_A) { info = -arg_idx_LDA; goto cleanup; }
        } else { // Fortran LDA is rows
            if (lda < MAX(1,n)) { info = -arg_idx_LDA; goto cleanup; }
        }
    } else { // N == 0
         // If A is not NULL when N=0, LDA must still be >= 1.
         // If A is NULL when N=0, this is often okay if Fortran handles N=0.
         if (a != NULL && lda < 1) {info = -arg_idx_LDA; goto cleanup;}
    }
    // ... similar detailed validation for other matrices (B, C, D, U, Y) ...
    // **Key: If a dimension (like N or M) is zero, arrays dependent on that dimension
    // might be allowed to be NULL by the Fortran routine. The C validation
    // should permit this if the Fortran routine does (check docs!).**

    // If any C-level validation fails before even considering Fortran's own checks
    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation
    // (Same logic as in CONTRIBUTING_FOR_AI.md Version 7, using documented formulas,
    // ensuring MAX(1, ...) or MAX(2, ...) for workspace array sizes as per SLICOT docs)
    // Example for dwork:
    // long long n_ll = n; ... (calculate ldw_formula_ll based on params) ...
    // ldwork_alloc = (int)MAX(MIN_DWORK_SIZE_FROM_DOC, ldw_formula_ll);
    // if (ldwork_alloc > 0) { dwork = ...; CHECK_ALLOC(dwork); } else { dwork = NULL; }

    // 4. Memory allocation for column-major copies (if row_major)
    // (Same logic as Version 7, calculate a_size, b_size etc. considering zero dimensions)
    // if (n > 0 && N_COLS_A > 0) a_size = (size_t)n * N_COLS_A; else a_size = 0;
    // if (row_major) {
    //    if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
    // }

    // 5. Prepare Fortran parameters and perform conversions
    // (Same logic as Version 7, setting _ptr and _f variables, handling NULLs for zero-size)
    // Fortran LDs (lda_f, ldb_f, etc.) must be >= 1 if the corresponding _ptr is not NULL.
    // If _ptr is NULL (e.g. for a zero-sized array), lda_f can be 1.

    // Example for A:
    // lda_f = (n == 0) ? 1 : MAX(1, n); // Fortran LDA is rows
    // if (row_major) {
    //    if (a_size > 0) { slicot_transpose_to_fortran_with_ld(a, a_cm, n, N_COLS_A, lda, lda_f, sizeof(double)); a_ptr = a_cm; }
    //    else { a_ptr = NULL; } // if N=0 or N_COLS_A=0
    // } else { // Column-major C
    //    if (a_size == 0) a_ptr = NULL;
    //    lda_f = lda; // Use C LDA (which is rows)
    // }
    // if (a_ptr != NULL && lda_f < 1) lda_f = 1; // Final safety for non-NULL pointers

    // 7. Call Fortran function
    F77_FUNC(function_name, FUNCTION_NAME)(
        &job_upper, &n, // ... other scalar params ...
        a_ptr, &lda_f,
        // ... other matrix/vector params ...
        iwork, dwork, &ldwork_alloc, // Pass internally allocated workspace
        &info // Fortran routine's info
        // &local_iwarn, // if applicable
        // job_len // Pass hidden length if applicable
    );
    // The Fortran routine's 'info' now holds the result.

    // 8. Convert results back to row-major (if needed)
    // (Same logic as Version 7)

cleanup:
    free(iwork);
    free(dwork);
    // ... free other workspace arrays ...
    if (row_major) {
        free(a_cm);
        // ... free other _cm arrays ...
    }

    // The 'info' variable now contains the status from the Fortran call,
    // or from C-level validation / memory allocation errors.
    return info;
}
```

## 3.1 Header File Documentation (include/function_name.h)

```c
/**
 * @file function_name.h
 * @brief Header for C wrapper of SLICOT routine FUNCTION_NAME.
 */

#ifndef SLICOT_WRAPPER_FUNCTION_NAME_H // Use unique guard
#define SLICOT_WRAPPER_FUNCTION_NAME_H

#include "slicot_utils.h" // Provides SLICOT_EXPORT macro

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief [Brief description of what the function does].
 * @details [More detailed explanation, algorithms used, purpose].
 * This is a C wrapper for the SLICOT Fortran routine FUNCTION_NAME.
 * **Workspace is allocated internally.**
 * **Note on Zero Dimensions:** Behavior for zero dimensions (e.g., N=0)
 * aligns with the underlying Fortran routine. Consult SLICOT documentation
 * for FUNCTION_NAME regarding specific handling and whether NULL pointers
 * are acceptable for zero-sized arrays.
 *
 * @param job [in] Description of job parameter (e.g., char: 'B', 'F', etc.).
 * @param n [in] Description of dimension n. If N=0, certain arrays may not be referenced or can be NULL.
 * @param a [in/out/in,out] Description of matrix A. Dimensions (e.g., n x n_cols_a).
 * Stored column-wise if row_major=0, row-wise if row_major=1.
 * Can often be NULL if n=0 or n_cols_a=0, check specific routine docs.
 * @param lda [in] Leading dimension of the C array storing A.
 * If row_major=0 (column-major), lda >= max(1, number of rows of A).
 * If row_major=1 (row-major), lda >= max(1, number of columns of A).
 * Must be >= 1 if 'a' is not NULL.
 * @param[out] x Output parameter x description.
 * @param row_major [in] Specifies matrix storage for input/output matrices like A, B, etc.:
 * 0 for column-major (Fortran default),
 * 1 for row-major (C default).
 *
 * @return info Error indicator:
 * = 0: successful exit
 * < 0: if info = -i, the i-th argument had an illegal value (wrapper or Fortran validation)
 * > 0: Fortran routine specific error (see SLICOT documentation for FUNCTION_NAME)
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.
 */
SLICOT_EXPORT
int slicot_function_name(/* C function parameters matching .c file, excluding workspace */
                         // const char* job, int n, double* a, int lda, ..., int* x, int row_major
                        );

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_WRAPPER_FUNCTION_NAME_H */
```

## 4. Test Case Template (Using CSV Loader)

```cpp
#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept> // For std::runtime_error
#include <algorithm> // For std::max

#include "function_name.h" // Include the wrapper header
#include "slicot_utils.h"  // For transpose functions if needed for setup/verification
#include "test_utils.h"    // For load_test_data_from_csv
#include "test_config.h"   // Include CMake-generated configuration for data path

// Use the TEST_DATA_DIR macro defined in test_config.h
const std::string DATA_FILE_PATH_PREFIX = TEST_DATA_DIR; // Allows appending filename

// --- Column-Major Test Fixture ---
class FunctionNameTestColMajor : public ::testing::Test {
protected:
    // Test parameters (set based on .dat/.res file and Fortran example program)
    int N = 3;
    // ... (rest of fixture setup) ...

    std::string csv_filename = "function_name.csv"; // Default, can be overridden

    void SetUp() override {
        // ... (initial setup) ...

        // Load data using correct signature
        std::vector<std::string> input_columns = {"U1", "U2"};
        std::vector<std::string> output_columns = {"Y1"};
        
        try {
            std::vector<double> loaded_u, loaded_y;
            int loaded_samples = 0;
            
            bool success = load_test_data_from_csv(DATA_FILE_PATH_PREFIX + csv_filename, 
                                                  input_columns, output_columns,
                                                  loaded_u, loaded_y, loaded_samples);
            
            if (success && loaded_samples > 0) {
                // Update NSMP based on loaded data
                NSMP = loaded_samples;
                // ... resize and copy data ...
            } else {
                throw std::runtime_error("Failed to load CSV data");
            }
        } catch (const std::exception& e) {
            // Fallback to embedded data
            // ... fallback data setup ...
        }
    }
};

class FunctionNameTestRowMajor : public FunctionNameTestColMajor {
protected:
    void SetUp() override {
        FunctionNameTestColMajor::SetUp();
        
        // Convert using correct function signatures - note the function name
        slicot_transpose_to_c_with_ld(A.data(), A_rm.data(), N, N_COLS_A, 
                                     fortran_lda, row_major_lda, sizeof(double));
        // ... other conversions ...
    }
};
```

## 5. Interpreting SLICOT Fortran Documentation for Edge Cases

This section is critical for developing robust C wrappers.

### 5.1 Zero Dimensions (N=0, M=0, L=0, NSMP=0, etc.)

- SLICOT routines often have specific, documented behavior for zero dimensions. Many will perform a "quick exit" with `INFO = 0` (see `rep96-1.pdf`, Sec 2.3.10: "If zero dimensions are encountered, which lead to immediate termination, the subroutine must set INFO = 0 and return...")
- Your C wrapper's validation logic must not be stricter than the Fortran routine. If the Fortran routine gracefully handles `N=0` by returning `INFO=0`, your C wrapper should allow `N=0` to pass through to the Fortran layer, rather than flagging it as a C-level error (unless `N<0`, which is always an error)
- Check the "Arguments" section of the SLICOT HTML documentation carefully. For example, if `N=0`, an array like `A(LDA,N)` becomes `A(LDA,0)`. The documentation might state "If N=0, this array is not referenced."

### 5.2 NULL Pointers for Zero-Sized Arrays

- If a dimension parameter (e.g., `N`) causes an array to become zero-sized (e.g., A becomes N-by-0 or 0-by-0), the corresponding C pointer argument (`a`) can often be `NULL`
- Verify this by:
  - Checking if the SLICOT documentation says the array "is not referenced" under those dimensional conditions
  - Testing: Pass `NULL` in your GTest cases for these scenarios
- If the Fortran routine does expect a non-NULL pointer even for a zero-sized array (less common, but possible for some interfaces or older code), the C wrapper might need to pass a pointer to a dummy 1-element static array. However, aim to pass `NULL` if the Fortran routine doesn't reference the data
- **Transpose Function Names**: When converting between row-major and column-major formats, use the correct function names from `slicot_utils.h`:
  - `slicot_transpose_to_fortran_with_ld()` - convert from C row-major to Fortran column-major
  - `slicot_transpose_to_c_with_ld()` - convert from Fortran column-major to C row-major
  - **Do NOT use** `slicot_transpose_from_fortran_with_ld()` as this function name does not exist in the actual utility header

### 5.3 Leading Dimensions (LDs) for NULL or Zero-Sized Arrays

- Even if an array pointer is `NULL` (because the array is zero-sized), the corresponding Fortran leading dimension argument (e.g., `LDA_F`) passed to the Fortran call should still be a valid positive integer, typically `MAX(1, relevant_dimension_for_ld)`
- For example, if `A_PTR` is `NULL` because `N=0`, `LDA_F` should still be passed as 1 (or `MAX(1,N)` which evaluates to 1). This is because the Fortran subroutine signature expects an integer, and some compilers/linkers might have issues with `LD=0` even if the array isn't accessed
- The SLICOT standard (`rep96-1.pdf`, Sec 2.3.10) implies non-positive leading dimensions are errors if the array is to be used. For unreferenced zero-sized arrays, `LD=1` is a safe default

### 5.4 Workspace for Zero Dimensions

Workspace calculation formulas might simplify or evaluate to small values (e.g., `MAX(1, N)` becomes 1 if `N=0`). Ensure the allocated workspace is at least the minimum specified by the documentation (often `MAX(1, ...)` or `MAX(2, ...)`). If a formula legitimately yields 0 and the routine needs no workspace for that case, passing `NULL` for the workspace pointer and 0 for its size might be acceptable if the Fortran routine handles it. When in doubt, allocate a minimal valid workspace (e.g., 1 element).

## 6. Helper Function Reference

To prevent hallucination and ensure correct usage, here are the exact function signatures available in the project:

### 6.1 Matrix Transpose Functions (from slicot_utils.h)

```c
// Basic transpose functions (compatible default leading dimensions)
void slicot_transpose_to_fortran(const void *src, void *dest, int rows, int cols, size_t elem_size);
void slicot_transpose_to_c(const void *src, void *dest, int rows, int cols, size_t elem_size);

// Transpose functions with custom leading dimensions
void slicot_transpose_to_fortran_with_ld(const void *src, void *dest, int rows, int cols,
                                        int ld_src, int ld_dest, size_t elem_size);
void slicot_transpose_to_c_with_ld(const void *src, void *dest, int rows, int cols,
                                  int ld_src, int ld_dest, size_t elem_size);

// In-place transpose for square matrices
int slicot_transpose_inplace(void *matrix, int rows, int cols, size_t elem_size);

// Symmetric matrix handling
void slicot_copy_symmetric_part(const void *src, void *dest, int n, char uplo, int ld, size_t elem_size);
void slicot_transpose_symmetric_to_fortran(const void *src, void *dest, int n, char uplo, size_t elem_size);
void slicot_transpose_symmetric_to_c(const void *src, void *dest, int n, char uplo, size_t elem_size);
```

### 6.2 Utility Functions (from slicot_utils.h)

```c
// Matrix initialization
void set_identity(int n, double* mat, int ld, int row_major);

// Debug/printing utilities
void printMatrixD(const char* name, const double* data, int rows, int cols, int ld, int rowMajor);
```

### 6.3 Test Data Loading (from test_utils.h)

```cpp
// CSV data loading function
bool load_test_data_from_csv(
    const std::string& filepath,
    const std::vector<std::string>& input_cols,
    const std::vector<std::string>& output_cols,
    std::vector<double>& u,
    std::vector<double>& y,
    int& num_samples);
```

### 6.4 Essential Macros (from slicot_utils.h)

```c
// Memory allocation checking
#define CHECK_ALLOC(ptr) \
    do { \
        if ((ptr) == NULL) { \
            info = SLICOT_MEMORY_ERROR; \
            goto cleanup; \
        } \
    } while (0)

// Error codes
#define SLICOT_MEMORY_ERROR -1010

// Math utilities
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

// Export/import macros
#define SLICOT_EXPORT // Platform-specific DLL export/import handling

// Complex number handling
#define SLICOT_COMPLEX_REAL(z) // Gets real part of complex number
```

### 6.5 **IMPORTANT**: Function Name Corrections

**DO NOT USE these function names (they do not exist):**
- `slicot_transpose_from_fortran_with_ld()` ❌
- `slicot_transpose_from_c_with_ld()` ❌
- `load_csv_data()` ❌
- `load_test_csv()` ❌

**ALWAYS USE these correct function names:**
- `slicot_transpose_to_fortran_with_ld()` ✅
- `slicot_transpose_to_c_with_ld()` ✅
- `load_test_data_from_csv()` ✅

### 6.6 Common Usage Patterns

**Converting from row-major to column-major (for Fortran calls):**
```c
slicot_transpose_to_fortran_with_ld(src_row_major, dest_col_major, rows, cols, 
                                   row_major_ld, fortran_ld, sizeof(double));
```

**Converting from column-major to row-major (after Fortran calls):**
```c
slicot_transpose_to_c_with_ld(src_col_major, dest_row_major, rows, cols, 
                             fortran_ld, row_major_ld, sizeof(double));
```

**Loading test data from CSV:**
```cpp
std::vector<std::string> input_columns = {"U1", "U2"};
std::vector<std::string> output_columns = {"Y1"};
std::vector<double> loaded_u, loaded_y;
int loaded_samples = 0;

bool success = load_test_data_from_csv(filepath, input_columns, output_columns,
                                      loaded_u, loaded_y, loaded_samples);
```

By meticulously checking the specific SLICOT routine's documentation and applying these principles, the C wrapper can more accurately mirror the Fortran routine's intended behavior, especially for edge cases.

## 7. Building the Project and Running Tests

This section provides instructions for AI agents and developers on how to build the project and run tests.

### 7.1 CMake Presets Overview

The project uses CMake presets (defined in `CMakePresets.json`) to standardize build configurations across different platforms and compiler toolchains. These presets define common configurations for Windows (Intel and MinGW toolchains) and macOS.

### 7.2 Configuring the Project

To configure the project using a preset:

```bash
# List available configuration presets
cmake --list-presets

# Configure using a specific preset (example: macOS debug)
cmake --preset macos-x64-debug
```

Available configuration presets include:
- `windows-x64-debug-intel`: Windows with Intel compilers (Debug)
- `windows-x64-release-intel`: Windows with Intel compilers (Release)
- `windows-x64-debug-mingw`: Windows with MinGW toolchain (Debug)
- `windows-x64-release-mingw`: Windows with MinGW toolchain (Release)
- `macos-x64-debug`: macOS with clang/gfortran toolchain (Debug)
- `macos-x64-release`: macOS with clang/gfortran toolchain (Release)

### 7.3 Building the Project

After configuration, build the project using the corresponding build preset:

```bash
# List available build presets
cmake --list-presets=build

# Build using a specific preset (example: macOS debug)
cmake --build --preset macos-x64-debug-build
```

Build presets are named with the `-build` suffix (e.g., `macos-x64-debug-build`).

### 7.4 Running Tests with CTest

Once the project is built, you can run the tests using CTest:

```bash
# Run all tests for a specific configuration
cd build/macos-x64-debug  # Navigate to the build directory
ctest
```

#### 7.4.1 Running Specific Tests

To run specific tests that match a pattern:

```bash
# Run all tests containing "function_name" in their name
ctest -R function_name
```

#### 7.4.2 Verbose Test Output

For more detailed test output, use the `--verbose` flag:

```bash
# Run tests with detailed output
ctest --verbose

# Combine with -R to run specific tests with detailed output
ctest -R function_name --verbose
```

#### 7.4.3 Running Tests via CMake Presets

Alternatively, you can run tests using the test presets:

```bash
# List available test presets
cmake --list-presets=test

# Run tests using a specific preset (example: macOS debug)
cmake --build --preset macos-x64-debug-test
```

Test presets are named with the `-test` suffix (e.g., `macos-x64-debug-test`).

### 7.5 Build and Test Examples

Complete workflow example for macOS:

```bash
# Configure
cmake --preset macos-x64-debug

# Build
cmake --build --preset macos-x64-debug-build

# Run all tests
cd build/macos-x64-debug
ctest

# Run specific tests with detailed output
ctest -R ab01od --verbose
```

AI agents implementing C wrappers should ensure all tests pass for both column-major and row-major storage formats and test edge cases thoroughly.