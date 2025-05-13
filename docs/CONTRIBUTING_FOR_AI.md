# **SLICOT_C Contribution Guidelines for AI Coding Agents**

**Version 6**

This document provides structured guidelines for AI coding agents to efficiently implement C wrappers and corresponding C++ tests for the SLICOT library, focusing on consistency, clarity, and reducing misinterpretations.

## **0. Key Principles (TL;DR for AI Agents)**

* **Goal:** Create a C wrapper function (e.g., `slicot_function_name`) that calls the corresponding SLICOT Fortran routine (`F77_FUNC(function_name, FUNCTION_NAME)`).  
* **Workspace:** The C wrapper **must allocate and free workspace arrays (`iwork`, `dwork`) internally** using `malloc`/`free`. Calculate sizes based on formulas in SLICOT documentation. Do **not** expect the caller to provide workspace.  
* **Row-Major Handling:** The wrapper must accept a `row_major` flag. If `row_major` is true (1), the wrapper must:  
  * Allocate temporary column-major buffers (`_cm`) for input/output matrices passed by the C caller.  
  * Transpose row-major C input matrices into the temporary column-major buffers before calling Fortran.  
  * Pass the *column-major buffers* and corresponding *Fortran-style leading dimensions* (number of **rows**) to the Fortran routine.  
  * If the Fortran routine modifies inputs or computes outputs, transpose results from column-major buffers back to the original C row-major arrays after the Fortran call (if `info == 0`).  
* **Leading Dimensions (LD):**  
  * **Fortran (and Column-Major C):** `LDA`, `LDB`, etc., always refer to the number of **rows** allocated for the matrix in memory.  
  * **Row-Major C:** `LDA`, `LDB`, etc., passed to the C wrapper refer to the number of **columns** allocated for the matrix in memory. The wrapper must calculate the required Fortran LDA (rows) internally.  
* **Testing:**  
  * Use GTest fixtures (`ColMajor`, `RowMajor`).  
  * **For small datasets** (e.g., ~10 samples), **embed data directly** in the test fixture.
  * **For larger datasets**, create a CSV file in the `tests/data/` directory and use the `load_test_data_from_csv` utility.
  * **CSV Data Format:** The `load_test_data_from_csv` utility reads the specified columns (e.g., "U1", "U2") and concatenates them into a flat vector (e.g., `[U1(1..NY), U2(1..NY), ...]`).
  * **Data Rearrangement:** For multi-column time-series data loaded via CSV (like input `U` or expected output `Y`), the test fixture **must explicitly rearrange** the loaded flat vector into the **true Fortran column-major order** (`[val(1,1), ..., val(rows,1), val(1,2), ..., val(rows,2), ...]`) before passing it to the wrapper or using it for comparison.
  * The **CSV header names MUST exactly match** the `input_columns`/`output_columns` specified in the test fixture.
  * `SetUp` loads/defines data, **performs necessary rearrangements**, updates `NSMP` (number of samples), calculates LDs, and sizes output vectors.
  * Tests call the C wrapper, `ASSERT_EQ` the `info` code, and `EXPECT_NEAR` the numerical results.
  * **Verify Expected Results:** Double-check expected results (`.res` files, documentation examples) against independent simulations (e.g., using Python control libraries) if possible, as documentation examples may contain errors. Use the verified results in the test fixture.
* **Filenames:** All C source and header files **must use lowercase naming** (e.g., `ab05od.c`, not `AB05OD.c`).
* **Error Handling:** Use `CHECK_ALLOC` after `malloc`. Implement `goto cleanup` for error exits. Free all allocated memory in the `cleanup` block. Return appropriate `info` codes (negative for wrapper errors, Fortran `info` otherwise).

## **1. Project Structure**
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
## **2. Implementation Workflow**

Follow this sequence precisely:

1. **Analyze Documentation:** Review the function's HTML documentation (e.g., `doc/AB05OD.html`) and the corresponding Fortran example files (`examples/*.dat`, `examples/*.res`, `examples/T*.f`). Understand parameters, dimensions, constraints, and **workspace formulas**.  
2. **Prepare/Verify Test Data:**
   * **Option A (For larger datasets):**
     * Locate the `.dat` file in `examples/`.
     * Create a corresponding CSV file in `tests/data/` (e.g., `tests/data/function_name.csv`). Use lowercase for the filename.
     * The **first row MUST be a header**. Choose clear, unique names (e.g., "U1", "Y1", "A11", "A12"). **These names are critical** for test data loading.
     * Copy the numerical data from the `.dat` file into subsequent rows of the CSV.
   * **Option B (Preferred for small datasets ~10 samples):** Define the test data directly within the C++ test fixture as `std::vector<double>`.
   * **Verification:** **Crucially, verify the expected output data** (from `.res` file or documentation) against an independent source (e.g., Python simulation using control or harold libraries). **If there are discrepancies between expected outputs and simulation results, trust the simulation results** and document this in the test code. Document the source of verified data if it differs from the SLICOT example.
3. **Extract Key Information:** Note parameters, types, dimensions, LD rules, and **exact workspace formulas**.  
4. **Identify Similar Wrappers:** Find wrappers in `src_c_wrapper/` for similar routines (e.g., AB05xx) as a reference.  
5. **Create C Wrapper:** Implement the `.c` file using the **Internal Workspace Allocation** template (Section 3). Ensure it correctly handles `row_major` and calculates workspace based on formulas.  
6. **Create Header File:** Implement the `.h` file in `include/` using the template (Section 3.1). Document parameters clearly, especially `row_major` handling.  
7. **Create Test Cases:** Implement the `_test.cpp` file using the GTest template (Section 4).  
   * **If using CSV:** Define `input_columns` and `output_columns` in the fixture to **exactly match the header names** created in the CSV file. Use `load_test_data_from_csv` in `SetUp`.
   * **Data Rearrangement (Important!):** If loading multi-column time-series data (e.g., inputs `U(M, NY)`, expected outputs `Y(P, NY)`) via CSV, **explicitly rearrange** the data loaded by `load_test_data_from_csv` from its concatenated format (e.g., `[U1(1..NY), U2(1..NY)]`) into the **true Fortran column-major order** (e.g., `[U(1,1), U(2,1), ..., U(M,1), U(1,2), ...]`) within the `SetUp` method. This rearranged data should be used for comparisons and passed to the wrapper (for inputs).
   * **If embedding data (preferred for small datasets):** Define the input and **verified** expected output data directly in the test fixture, preferably already in Fortran column-major order.
   * Implement tests for column-major, row-major, and parameter validation.

## **3. C Wrapper Template (Internal Workspace Allocation)**

```c
/**  
 * @file function_name.c  
 * @brief C wrapper for SLICOT routine FUNCTION_NAME.  
 * @details [Description of what the routine does, its purpose, and key algorithms.]  
 * Workspace (IWORK, DWORK) is allocated internally by this wrapper.  
 * Input/output matrix format is handled via the row_major parameter.  
 */

#include <stdlib.h> // For malloc, free  
#include <ctype.h>  // For toupper  
#include <stddef.h> // For size_t  
#include <math.h>   // For MAX/MIN if needed (often provided by slicot_utils.h)  
#include <stdio.h>  // For error logging (optional)

#include "function_name.h" // Public header for this wrapper  
#include "slicot_utils.h"  // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX/MIN, transpose functions etc.  
#include "slicot_f77.h"    // Provides F77_FUNC macro for Fortran name mangling

* External Fortran routine declaration */  
extern void F77_FUNC(function_name, FUNCTION_NAME)(  
    * Fortran function parameters with const for inputs */  
    // const char* job, const int* n, ...,  
    * Workspace arrays passed by Fortran */  
    int* iwork, double* dwork, const int* ldwork,  
    int* info  
    * Hidden string lengths if any */  
    // int job_len  
);

* C wrapper function definition */  
SLICOT_EXPORT // Macro for DLL export/import handling  
int slicot_function_name(* C function parameters, excluding workspace */  
                         // const char* job, int n, ..., double* a, int lda, ..., int row_major  
                        )  
{  
    // 1. Variable declarations  
    int info = 0;          // Return status code  
    int *iwork = NULL;     // Internal integer workspace pointer  
    double *dwork = NULL;  // Internal double workspace pointer  
    int liwork = 0;        // Calculated size for iwork  
    int ldwork = 0;        // Calculated size for dwork

    // Pointers for column-major copies (if row_major is used)  
    double* a_cm = NULL;  
    // ... declare other _cm pointers as needed ...

    // 2. Input parameter validation (Check BEFORE allocating memory)  
    // Check scalar parameters first  
    if (n < 0) { info = -2; goto cleanup; } // Example: Check dimension n  
    char job_upper = toupper(job); // Convert character options once  
    if (job_upper != 'B' && job_upper != 'F') { info = -1; goto cleanup; } // Example: Check option job

    // Check pointers for required arrays (handle optional arrays carefully)  
    // Example: If 'a' is required when n > 0  
    if (a == NULL && n > 0) { info = -6 * Map to correct Fortran argument index */; goto cleanup; }

    // Check leading dimensions based on row_major flag  
    int min_lda_f = MAX(1, n); // Minimum Fortran LDA (number of ROWS)  
    if (row_major) {  
        // C LDA is number of COLUMNS  
        if (n > 0 && lda < n) { info = -6; goto cleanup; } // Check minimum columns for A  
    } else {  
        // C LDA is number of rows (Fortran style)  
        if (lda < min_lda_f) { info = -6; goto cleanup; } // Check minimum rows for A  
    }  
    // ... validate other pointers and leading dimensions ...

    // Exit if any validation failed before allocating memory  
    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation  
    // STEP 1: Check if the routine supports workspace query
    // Look in the documentation for: "On exit, if INFO = 0, DWORK(1) returns the optimal value of LDWORK"
    // If that text exists, the routine supports workspace query
    
    // METHOD A: FOR ROUTINES SUPPORTING WORKSPACE QUERY
    // Use workspace query to get optimal workspace size
    ldwork = -1; // Signal workspace query
    double dwork_query[1];
    
    // Call Fortran routine with workspace query flag (ldwork = -1)
    F77_FUNC(function_name, FUNCTION_NAME)(
        // ...pass parameters, with NULL for arrays that have size formulas...
        NULL, NULL, dwork_query, &ldwork, &info
    );
    
    // If query successful, use returned optimal workspace size
    if (info == 0 || info == -19) { // -19 is often returned for workspace query
        ldwork = (int)dwork_query[0];
        // Ensure minimum size based on documentation formula 
        int min_ldwork = MAX(1, /* formula for ldwork based on n, m, p etc. */);
        ldwork = MAX(ldwork, min_ldwork);
    } else {
        // Fallback to formula if query fails
        ldwork = MAX(1, /* formula for ldwork based on n, m, p etc. */);
    }
    
    // METHOD B: FOR ROUTINES NOT SUPPORTING WORKSPACE QUERY
    // Skip the query and directly use the formula from documentation
    ldwork = MAX(1, /* formula for ldwork based on n, m, p etc. */);
    
    // HANDLE ZERO DIMENSIONS EXPLICITLY
    // Formula may not work correctly for n=0, m=0, p=0
    if (n == 0) { 
        ldwork = 1; // Minimum workspace size for n=0
    }
    
    // Allocate workspace
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));  
    CHECK_ALLOC(dwork);

    // Handle integer workspace similarly with zero dimensions
    liwork = MAX(1, /* formula for liwork */);
    if (n == 0 && liwork_depends_on_n) {
        // Only allocate if positive size needed (n>0)
        iwork = NULL;
    } else {
        iwork = (int*)malloc((size_t)liwork * sizeof(int));  
        CHECK_ALLOC(iwork);
    }

    // 4. Memory allocation for column-major copies (if row_major)  
    size_t a_size = (size_t)n * n; if (n == 0) a_size = 0; // Example for matrix A (n x n)  
    // ... calculate other sizes ...  
    if (row_major) {  
        if (a_size > 0) {  
            a_cm = (double*)malloc(a_size * sizeof(double));  
            CHECK_ALLOC(a_cm);  
        }  
        // ... allocate other _cm arrays ...  
    }

    // 5. Prepare Fortran parameters and perform conversions  
    double* a_ptr = a; // Default: point to original C data  
    int lda_f = lda;   // Default: use original C LDA (which is rows if col-major, cols if row-major)

    if (row_major) {  
        // Set Fortran LDA (number of ROWS)  
        lda_f = MAX(1, n);  
        // Point to column-major copy if allocated, else NULL  
        if (a_size > 0) {  
            slicot_transpose_to_fortran(a, a_cm, n, n, sizeof(double)); // Convert input A  
            a_ptr = a_cm;  
        } else {  
            a_ptr = NULL; // Pass NULL if size is 0  
        }  
        // ... prepare other pointers and Fortran LDs (ROWS) for row_major ...  
    } else {  
        // Column-major C: Ensure NULL is passed if size is 0  
        if (a_size == 0) a_ptr = NULL;  
        // ... prepare other pointers for column_major (pass NULL if size 0)...  
        // lda_f is already correctly set to the C LDA (rows)  
    }

    // Prepare other Fortran parameters (e.g., hidden string lengths)  
    const int job_len = 1; // Example

    // 7. Call Fortran function  
    F77_FUNC(function_name, FUNCTION_NAME)(  
        \&job_upper, \&n,  
        a_ptr, \&lda_f, // Pass potentially converted pointer and Fortran LDA (ROWS)  
        // ... other parameters ...  
        iwork, dwork, \&ldwork, // Pass internally allocated workspace  
        \&info,  
        job_len // Pass hidden length if applicable  
    );

    // 8. Convert results back to row-major (if needed)  
    // This is crucial if Fortran modifies input arrays (e.g., EQUIL='S')  
    // or if output arrays need conversion.  
    if (row_major && info == 0) {  
        // Example: Input 'a' was modified and needs copying back  
        // if (job_upper == 'S' && a_size > 0) {  
        //     slicot_transpose_to_c(a_cm, a, n, n, sizeof(double));  
        // }  
        // Example: Output 'x' (p x q) was written to x_cm (Fortran ldx_f rows)  
        //          C array 'x' has ldx columns.  
        // if (x_size > 0) {  
        //     slicot_transpose_to_c_with_ld(x_cm, x, p, q, ldx_f, ldx, sizeof(double));  
        // }  
    }

cleanup:  
    /* --- Cleanup --- */  
    // Free internally allocated workspace FIRST  
    free(iwork); // Safe to call free(NULL)  
    free(dwork); // Safe to call free(NULL)  
    // Free temporary column-major arrays if allocated  
    free(a_cm);  
    // ... free other _cm arrays ...

    // Check if info was set by CHECK_ALLOC during workspace/copy allocation  
    if (info == SLICOT_MEMORY_ERROR) {  
       fprintf(stderr, "Error: Memory allocation failed in slicot_%s.\\n", "function_name");  
    }  
    // Return the info code (either from validation, CHECK_ALLOC, or Fortran)  
    return info;  
}
```

### **3.1 Header File Documentation (include/function_name.h)**

```cpp
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
 *  
 * @param job [Description of job parameter ('B', 'F', etc.)].  
 * @param n [Description of dimension n].  
 * @param a [in/out/in,out] Description of matrix A. Dimensions (n x n).  
 * Stored column-wise if row_major=0, row-wise if row_major=1.  
 * @param lda Leading dimension of the C array storing A.  
 * If row_major=0 (column-major), lda >= max(1, n) (number of rows).  
 * If row_major=1 (row-major), lda >= max(1, n) (number of columns).  
 * @param[out] x Output parameter x description.  
 * @param row_major Specifies matrix storage for input/output matrices like A, B, etc.:  
 * 0 for column-major (Fortran default),  
 * 1 for row-major (C default).  
 * **Note:** Some Fortran output arrays might only be returned in column-major format  
 * even if `row_major=1`. Check wrapper implementation/docs.  
 *  
 * @return info Error indicator:  
 * = 0: successful exit  
 * < 0: if info = -i, the i-th argument had an illegal value (wrapper or Fortran validation)  
 * > 0: Fortran routine specific error (see SLICOT documentation for FUNCTION_NAME)  
 * = SLICOT_MEMORY_ERROR (-1010): internal memory allocation failed.  
 */  
SLICOT_EXPORT  
int slicot_function_name(* C function parameters matching .c file, excluding workspace */  
                         // char job, int n, double* a, int lda, ..., int* x, int row_major  
                        );

#ifdef __cplusplus  
}  
#endif

#endif * SLICOT_WRAPPER_FUNCTION_NAME_H */
```

## **4. Test Case Template (Using CSV Loader)**

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
const std::string DATA_FILE_PATH = TEST_DATA_DIR "function_name.csv";

// --- Column-Major Test Fixture ---  
class FunctionNameTestColMajor : public ::testing::Test {  
protected:  
    // Test parameters (set based on .dat/.res file)  
    int N = 3;  
    int M = 1;  
    int P = 1;  
    int NSMP = 0; // Number of samples - **determined by loader in SetUp**  
    char PARAM = 'X';  
    // ... other scalar parameters ...

    // Column names to load from CSV (**MUST match CSV header exactly**)  
    std::vector<std::string> input_columns = {"U1"}; // Example  
    std::vector<std::string> output_columns = {"Y1"}; // Example

    // Verification tolerance  
    double check_tol = 1e-4;

    // Input/Output data vectors (column-major)  
    std::vector<double> U; // Loaded inputs from CSV  
    std::vector<double> Y; // Loaded outputs from CSV  
    std::vector<double> A; // Example input matrix (if needed, initialized in SetUp)  
    // ... other input/output vectors ...

    // Expected results (from .res file or manual calculation)  
    std::vector<double> A_expected;  
    // ... other expected output vectors ...  
    int expected_info = 0;  
    // ... other expected scalar outputs ...

    // Result variables (initialized in test body)  
    int info_result = -999; // Initialize to indicate not run  
    // ... other result variables ...

    // Leading dimensions (**calculated in SetUp AFTER NSMP is known**)  
    int LDA = 1;  
    int LDU = 1;  
    int LDY = 1;  
    // ... other leading dimensions ...

    // SetUp method: Load data, initialize inputs, size outputs  
    void SetUp() override {  
        // --- Load Time Series Data (if applicable) ---  
        int samples_loaded = 0;  
        // Ensure fixture M/P match the number of columns requested  
        ASSERT_EQ(input_columns.size(), M) << "Fixture M != input_columns size";  
        ASSERT_EQ(output_columns.size(), P) << "Fixture P != output_columns size";

        try {  
            // Use DATA_FILE_PATH which now uses TEST_DATA_DIR
            bool success = load_test_data_from_csv(  
                DATA_FILE_PATH, input_columns, output_columns,  
                U, Y, samples_loaded // U, Y are populated here  
            );  
            ASSERT_TRUE(success) << "CSV loading reported failure for " << DATA_FILE_PATH;  
            ASSERT_GT(samples_loaded, 0) << "No data samples loaded from CSV: " << DATA_FILE_PATH;  
            // **CRITICAL: Update NSMP based on actual data loaded**  
            NSMP = samples_loaded;

        } catch (const std::runtime_error& e) {  
            FAIL() << "CSV data loading failed: " << e.what();  
        } catch (...) {  
            FAIL() << "Caught unknown exception during CSV data loading.";  
        }

        // --- Initialize other inputs / Size outputs ---  
        // Example: Initialize input matrix A (column-major)  
        A = { * values for A */ };  
        A_expected = { * expected values for A (if modified) or other output */ };

        // Size output arrays based on parameters and NSMP  
        // Example: Output vector X (size N)  
        // X_out.resize(N);

        // **Calculate Leading Dimensions AFTER NSMP is set** (Column Major)  
        LDA = std::max(1, N); // Fortran/Col-major LDA is rows  
        LDU = (M > 0) ? std::max(1, NSMP) : 1; // Fortran/Col-major LDU is rows  
        LDY = (P > 0) ? std::max(1, NSMP) : 1; // Fortran/Col-major LDY is rows  
        // ... calculate other LDs ...  
    }  
};

// --- Row-Major Test Fixture ---  
class FunctionNameTestRowMajor : public FunctionNameTestColMajor {  
protected:  
    // Input data vectors in row-major format  
    std::vector<double> U_rm;  
    std::vector<double> Y_rm;  
    std::vector<double> A_rm; // Example input matrix  
    // ... other row-major input vectors ...

    // Expected output vectors in row-major format (if needed for comparison)  
    std::vector<double> A_expected_rm;  
     // ... other row-major expected vectors ...

    void SetUp() override {  
        // --- Load Column-Major Data First ---  
        int samples_loaded = 0;  
        ASSERT_EQ(input_columns.size(), M) << "Fixture M != input_columns size";  
        ASSERT_EQ(output_columns.size(), P) << "Fixture P != output_columns size";

        std::vector<double> U_col, Y_col; // Temporary column-major storage  
        try {  
             // Use DATA_FILE_PATH which now uses TEST_DATA_DIR
             bool success = load_test_data_from_csv(  
                DATA_FILE_PATH, input_columns, output_columns,  
                U_col, Y_col, samples_loaded  
            );  
            ASSERT_TRUE(success) << "CSV loading reported failure for " << DATA_FILE_PATH;  
            ASSERT_GT(samples_loaded, 0) << "No data samples loaded from CSV: " << DATA_FILE_PATH;  
            // **CRITICAL: Update NSMP based on actual data loaded**  
            NSMP = samples_loaded; // Update base class NSMP as well

        } catch (const std::runtime_error& e) {  
            FAIL() << "CSV data loading failed: " << e.what();  
        } catch (...) {  
            FAIL() << "Caught unknown exception during CSV data loading.";  
        }

        // --- Convert Loaded Data to Row-Major ---  
        U_rm.resize(U_col.size());  
        Y_rm.resize(Y_col.size());  
        // Use U_col, Y_col as source for transpose  
        if (M > 0 && !U_col.empty()) slicot_transpose_to_c(U_col.data(), U_rm.data(), NSMP, M, sizeof(double));  
        if (P > 0 && !Y_col.empty()) slicot_transpose_to_c(Y_col.data(), Y_rm.data(), NSMP, P, sizeof(double));

        // --- Initialize/Convert other inputs ---  
        // Example: Initialize input matrix A (row-major)  
        // A_rm = { * row-major values for A */ };  
        // Example: Convert expected results if needed for comparison  
        // A_expected_rm.resize(A_expected.size());  
        // if (!A_expected.empty()) slicot_transpose_to_c(A_expected.data(), A_expected_rm.data(), N, N, sizeof(double));

        // --- Size outputs (same as column-major, wrapper handles output format) ---  
        // Example: Output vector X (size N)  
        // X_out.resize(N);

        // **Calculate Leading Dimensions AFTER NSMP is set** (Row Major C requires cols)  
        LDA = std::max(1, N); // Cols of A  
        LDU = (M > 0) ? std::max(1, M) : 1; // Cols of U  
        LDY = (P > 0) ? std::max(1, P) : 1; // Cols of Y  
         // ... calculate other LDs for C row-major arrays ...  
    }  
};

// --- Test Cases ---

// Test: Documentation Example (Column-Major)  
TEST_F(FunctionNameTestColMajor, DocExample) {  
    // Call C wrapper function (workspace handled internally)  
    info_result = slicot_function_name(PARAM, N, M, P, NSMP, // Use NSMP from SetUp  
                                       A.data(), LDA, // Pass potentially initialized A  
                                       (M > 0 ? U.data() : nullptr), LDU, // Pass loaded U  
                                       (P > 0 ? Y.data() : nullptr), LDY, // Pass loaded Y  
                                       * other args... */,  
                                       0 * row_major = false */);

    // Verify return code  
    ASSERT_EQ(info_result, expected_info);

    // Verify output parameters and arrays against expected column-major values  
    // Example: Check output vector X_out against X_expected  
    // ASSERT_EQ(X_out.size(), X_expected.size());  
    // for (size_t i = 0; i < X_expected.size(); \++i) {  
    //    EXPECT_NEAR(X_out[i], X_expected[i], check_tol) << "Mismatch at index " << i;  
    // }  
}

// Test: Documentation Example (Row-Major)  
TEST_F(FunctionNameTestRowMajor, DocExample) {  
    // Call C wrapper function (workspace handled internally)  
    info_result = slicot_function_name(PARAM, N, M, P, NSMP, // Use NSMP from SetUp  
                                       A_rm.data(), LDA, // Pass row-major A and its LDA (cols)  
                                       (M > 0 ? U_rm.data() : nullptr), LDU, // Pass row-major U and its LDU (cols)  
                                       (P > 0 ? Y_rm.data() : nullptr), LDY, // Pass row-major Y and its LDY (cols)  
                                       * other args... */,  
                                       1 * row_major = true */);

    // Verify return code  
    ASSERT_EQ(info_result, expected_info);

    // Verify output parameters and arrays.  
    // **CRITICAL**: Wrapper output arrays (like X_out) are likely still column-major.  
    // Input arrays modified in place (like A_rm if EQUIL='S') will be row-major.  
    // Compare accordingly, potentially transposing actual or expected results for verification.  
    // Example: Compare modified A_rm against row-major A_expected_rm  
    // for (size_t i = 0; i < A_expected_rm.size(); \++i) {  
    //     EXPECT_NEAR(A_rm[i], A_expected_rm[i], check_tol);  
    // }  
}

// Test: Parameter Validation (specific checks for the wrapper)  
TEST_F(FunctionNameTestColMajor, ParameterValidation) {  
    // This test checks if the C wrapper correctly validates inputs.  
    // It does NOT need valid computational data from CSV, just dummy placeholders.  
    std::vector<double> dummy_A(1); // Minimal allocation  
    int n_out = 0; // Dummy output scalar

    // Test invalid N  
    info_result = slicot_function_name(PARAM, -1, M, P, 10, * Use dummy NSMP */  
                                       dummy_A.data(), 1, * Dummy LDA */  
                                       nullptr, 1, nullptr, 1, * Dummy U, Y */  
                                       * other args... */, 0);  
    EXPECT_EQ(info_result, -2 * Expected error code for N */);

    // Test invalid PARAM  
    info_result = slicot_function_name('Z', N, M, P, 10, dummy_A.data(), N, * ... */, 0);  
    EXPECT_EQ(info_result, -1 * Expected error code for PARAM */);

    // Test invalid LDA (Column-Major)  
    info_result = slicot_function_name(PARAM, N, M, P, 10, dummy_A.data(), 0 * LDA=0 */, * ... */, 0);  
    EXPECT_EQ(info_result, -6 * Expected error code for LDA */);

    // Test NULL pointer for required array A (assuming N > 0)  
    if (N > 0) {  
        info_result = slicot_function_name(PARAM, N, M, P, 10, nullptr * A=NULL */, N, * ... */, 0);  
        EXPECT_EQ(info_result, -6 * Expected error code for A */);  
    }

    // Add more specific validation checks for other parameters...  
}
```

## **5. Key Additions to Contributing Guidelines for Zero-Dimension Cases**

### **5.1 Memory Allocation (Column-Major Copies)**

```c
// --- Inside the C wrapper function ---

// 1. Declare pointers for column-major copies  
double *a_cm = NULL, *b_cm = NULL;

// 2. Calculate size with zero dimension handling (after input validation)  
size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;  
size_t b_size = (size_t)n * m; if (n == 0 || m == 0) b_size = 0;  
if (row_major) {  
    if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }  
    if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }  
}

// 3. Set pointers for Fortran call (a_ptr, b_ptr) and Fortran LDs (lda_f, ldb_f)  
double* a_ptr = a; // Default to original C pointer  
double* b_ptr = b; // Default to original C pointer  
int lda_f = lda;   // Default to original C LDA  
int ldb_f = ldb;   // Default to original C LDB

if (row_major) {  
    // Fortran expects column-major, so LD is number of rows  
    lda_f = MAX(1, n); // Fortran LDA = rows  
    ldb_f = MAX(1, n); // Fortran LDB = rows  
    // Point to the allocated column-major copy if size > 0, otherwise pass NULL  
    a_ptr = (a_size > 0) ? a_cm : NULL;  
    b_ptr = (b_size > 0) ? b_cm : NULL;  
} else {  
    // Column-major C: LDs are already rows. Pass NULL if size is 0.  
    if (a_size == 0) a_ptr = NULL;  
    if (b_size == 0) b_ptr = NULL;  
}

// 4. Perform conversion before Fortran call (copy data into _cm buffers)  
if (row_major) {  
    if (a_size > 0) {  
        slicot_transpose_to_fortran(a, a_cm, n, n, sizeof(double));  
    }  
    if (b_size > 0) {  
        slicot_transpose_to_fortran(b, b_cm, n, m, sizeof(double));  
    }  
}

// --- Fortran Call using a_ptr, lda_f, b_ptr, ldb_f ---

// --- Cleanup Section ---  
cleanup:  
    free(a_cm); // Safe to call free on NULL  
    free(b_cm);  
    // ... other frees ...  
    return info;
```

### **5.2 Parameter Validation**

```c
// --- Inside the C wrapper function ---  
// **Perform ALL validation BEFORE allocating any memory (workspace or copies)**

// Check scalar parameters first  
if (n < 0) { info = -2; goto cleanup; } // Map to Fortran index for N  
if (m < 0) { info = -3; goto cleanup; } // Map to Fortran index for M  
char job_upper = toupper(job); // Convert once  
if (job_upper != 'B' && job_upper != 'F') { info = -1; goto cleanup; } // Map to Fortran index for JOB

// Check required pointers (handle optional arrays based on dimensions/flags)  
// Example: A is required if N > 0  
if (a == NULL && n > 0) { info = -6; goto cleanup; } // Map to Fortran index for A  
// Example: B is required only if N > 0, M > 0 and JOB = 'B'  
if (b == NULL && n > 0 && m > 0 && job_upper == 'B') { info = -8; goto cleanup; } // Map to Fortran index for B

// Check leading dimensions based on row_major flag and dimensions  
int min_lda_f = MAX(1, n); // Minimum Fortran LDA (rows) for A  
int min_ldb_f = MAX(1, n); // Minimum Fortran LDB (rows) for B  
if (row_major) {  
    // C LDA is number of columns  
    if (n > 0 && lda < n) { info = -6; goto cleanup; } // Check minimum columns for A  
    if (m > 0 && ldb < m) { info = -8; goto cleanup; } // Check minimum columns for B  
} else {  
    // C LDA is number of rows  
    if (lda < min_lda_f) { info = -6; goto cleanup; } // Check minimum rows for A  
    if (ldb < min_ldb_f) { info = -8; goto cleanup; } // Check minimum rows for B  
}

// If any validation fails, jump to cleanup *before* allocations  
if (info != 0) { goto cleanup; }

// --- Proceed to Workspace Allocation ---
```

### **5.3 Matrix Format Conversion**

```c
// --- Inside the C wrapper function ---

// 1. Allocate column-major temporary storage (if row_major)
double *a_cm = NULL, *x_cm = NULL; // a is input (n x n), x is output (p x q)  
size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;  
size_t x_size = (size_t)p * q; if (p == 0 || q == 0) x_size = 0;  
if (row_major) {  
    if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }  
    // Allocate buffer for Fortran to write output into if row_major  
    if (x_size > 0) { x_cm = (double*)malloc(x_size * sizeof(double)); CHECK_ALLOC(x_cm); }  
}

// 2. Set Fortran pointers and LDs  
double* a_ptr = a;  
double* x_ptr = x; // C pointer to the output array  
int lda_f = lda;   // C LDA  
int ldx_f = ldx;   // C LDX

if (row_major) {  
    // Fortran needs pointers to column-major buffers and LDs as ROWS  
    lda_f = MAX(1, n); // Fortran LDA = rows  
    ldx_f = MAX(1, p); // Fortran LDX = rows  
    a_ptr = (a_size > 0) ? a_cm : NULL;  
    x_ptr = (x_size > 0) ? x_cm : NULL; // Fortran writes output to temp buffer  
} else {  
    // Column-major C: Pass original pointers (or NULL if size is 0) and C LDs (rows)  
    if (a_size == 0) a_ptr = NULL;  
    if (x_size == 0) x_ptr = NULL;  
}

// 3. Convert C row-major INPUT 'a' to column-major 'a_cm' (before Fortran call)  
if (row_major && a_size > 0) { 
    // IMPORTANT: Use the _with_ld variant when matrices have leading dimensions
    // different from their actual dimensions (common case)
    slicot_transpose_to_fortran_with_ld(a, a_cm, n, n, lda, lda_f, sizeof(double));
    
    // Use the basic transpose function only when the matrix is tightly packed
    // slicot_transpose_to_fortran(a, a_cm, n, n, sizeof(double));
}

// --- Call Fortran routine with a_ptr, lda_f, x_ptr, ldx_f ---  
// Fortran computes output and writes to x_ptr (which is x_cm if row_major, or x if col-major)

// 4. Convert Fortran column-major OUTPUT 'x_cm' back to C row-major 'x' (after Fortran call)  
if (row_major && info == 0) {
    // Convert output matrix 'x' (p x q) from column-major x_cm to row-major x
    // Handle leading dimensions properly with _with_ld variant 
    if (x_size > 0) {  
        slicot_transpose_to_c_with_ld(x_cm, x, p, q, ldx_f, ldx, sizeof(double));  
    }  
}

// --- Cleanup Section ---  
cleanup:  
    free(a_cm);  
    free(x_cm); // Free the temporary output buffer if allocated  
    return info;
```

### **5.4 Handling Zero Dimensions**

Proper handling of zero-dimension cases is **critical** for wrapper robustness:

```c
// --- Calculate sizes with zero dimension handling ---
size_t a_size = (size_t)n * n;
size_t b_size = (size_t)n * m;

// These checks ensure a_size = 0 if n = 0, preventing allocation
if (n == 0) a_size = 0;
if (n == 0 || m == 0) b_size = 0;

// --- Allocate only if positive size ---
if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }

// --- Workspace calculation for n = 0 case ---
int liwork = 0;
if (job_upper == 'N') {
    liwork = MAX(1, n);  // Formula for positive n
    if (n > 0) {
        // Only allocate if n > 0
        iwork = (int*)malloc((size_t)liwork * sizeof(int));
        CHECK_ALLOC(iwork);
    } else {
        iwork = NULL;  // For n=0, pass NULL to Fortran
    }
}

// --- Set pointers properly for zero dimensions ---
double* a_ptr = (a_size > 0) ? a : NULL;  // Pass NULL if size is 0
```

When validating parameters for zero dimensions:

```c
// --- Check leading dimensions only when dimensions > 0 ---
if (row_major) {
    if (n > 0 && lda < n) { info = -6; goto cleanup; }  // Only check if n > 0
    if (n > 0 && m > 0 && ldb < m) { info = -8; goto cleanup; }  // Check if n > 0 AND m > 0
}
```

Important to test explicitly with zero-dimension tests:
```cpp
// Test the special case where N=0
TEST_F(FunctionNameTestColMajor, ZeroDimension) {
    int zero_n = 0;
    // Call wrapper with N=0 and check that it succeeds
    info_result = slicot_function_name(PARAM, zero_n, M, P, 
                                      nullptr, 1, nullptr, 1, nullptr, 1,
                                      0 /* row_major = false */);
    EXPECT_EQ(info_result, 0);
}
```