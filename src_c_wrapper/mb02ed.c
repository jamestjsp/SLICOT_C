/**
 * @file mb02ed.c
 * @brief C wrapper implementation for SLICOT routine MB02ED
 *
 * This file provides a C wrapper implementation for the SLICOT routine MB02ED,
 * which solves T*X = B or X*T = B for a symmetric positive definite
 * block Toeplitz matrix T.
 */

#include <stdlib.h>
#include <ctype.h>
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "mb02ed.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * Note T and B are input/output.
 */
extern void F77_FUNC(mb02ed, MB02ED)(
    const char* typet,      // CHARACTER*1 TYPET
    const int* k,           // INTEGER K
    const int* n,           // INTEGER N
    const int* nrhs,        // INTEGER NRHS
    double* t,              // DOUBLE PRECISION T(LDT,*) (in/out)
    const int* ldt,         // INTEGER LDT
    double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
    const int* ldb,         // INTEGER LDB
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info,              // INTEGER INFO (output)
    int typet_len           // Hidden length
);


/* C wrapper function definition */
SLICOT_C_WRAPPER_API
int slicot_mb02ed(char typet, int k, int n, int nrhs,
                  double* t, int ldt, double* b, int ldb,
                  int row_major)
{
    /* Local variables */
    int info = 0;
    int ldwork = -1; /* Use -1 for workspace query */
    double dwork_query;
    double* dwork = NULL;
    const int typet_len = 1;

    char typet_upper = toupper(typet);

    /* Pointers for column-major copies if needed */
    double *t_cm = NULL, *b_cm = NULL;

    /* --- Input Parameter Validation --- */

    if (k < 0) { info = -2; goto cleanup; }
    if (n < 0) { info = -3; goto cleanup; }
    if (nrhs < 0) { info = -4; goto cleanup; }
    if (typet_upper != 'R' && typet_upper != 'C') { info = -1; goto cleanup; }

    // Determine dimensions based on TYPET
    int t_fort_rows, t_fort_cols, b_fort_rows, b_fort_cols;
    int t_c_rows, t_c_cols, b_c_rows, b_c_cols;
    int min_ldt_f, min_ldb_f;

    if (typet_upper == 'R') {
        t_fort_rows = k; t_fort_cols = n * k;
        b_fort_rows = nrhs; b_fort_cols = n * k;
        t_c_rows = k; t_c_cols = n * k; // C row-major interpretation
        b_c_rows = nrhs; b_c_cols = n * k; // C row-major interpretation
        min_ldt_f = MAX(1, k);
        min_ldb_f = MAX(1, nrhs);
    } else { // TYPET = 'C'
        t_fort_rows = n * k; t_fort_cols = k;
        b_fort_rows = n * k; b_fort_cols = nrhs;
        t_c_rows = n * k; t_c_cols = k; // C row-major interpretation
        b_c_rows = n * k; b_c_cols = nrhs; // C row-major interpretation
        min_ldt_f = MAX(1, n * k);
        min_ldb_f = MAX(1, n * k);
    }


    // Check leading dimensions based on storage order
    if (row_major) {
        // For row-major C, LD is the number of columns
        int min_ldt_rm_cols = t_c_cols;
        int min_ldb_rm_cols = b_c_cols;
        if (ldt < min_ldt_rm_cols) { info = -6; goto cleanup; }
        if (ldb < min_ldb_rm_cols) { info = -8; goto cleanup; }
    } else {
        // For column-major C, LD is the number of rows (Fortran style)
        if (ldt < min_ldt_f) { info = -6; goto cleanup; }
        if (ldb < min_ldb_f) { info = -8; goto cleanup; }
    }

    /* --- Workspace Allocation --- */

    // Allocate DWORK based on query
    ldwork = -1; // Query mode
    F77_FUNC(mb02ed, MB02ED)(&typet_upper, &k, &n, &nrhs, t, &ldt, b, &ldb,
                             &dwork_query, &ldwork, &info, typet_len);

    // Allow INFO = -10 from query (indicates LDWORK too small, returns minimum)
    if (info < 0 && info != -10) { goto cleanup; }
    info = 0; // Reset info after query

    // Get the required dwork size from query result
    ldwork = (int)dwork_query;
    // Check against minimum documented size: MAX(1, N*K*K + (N+2)*K)
    int min_ldwork = MAX(1, n * k * k + (n + 2) * k);
    ldwork = MAX(ldwork, min_ldwork);

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

    /* --- Prepare Arrays and Call Fortran Routine --- */
    size_t elem_size = sizeof(double);

    if (row_major) {
        /* --- Row-Major Case --- */

        /* Allocate memory for column-major copies */
        size_t t_size = (size_t)t_c_rows * t_c_cols; if (t_c_rows == 0 || t_c_cols == 0) t_size = 0;
        size_t b_size = (size_t)b_c_rows * b_c_cols; if (b_c_rows == 0 || b_c_cols == 0) b_size = 0;

        if (t_size > 0) { t_cm = (double*)malloc(t_size * elem_size); CHECK_ALLOC(t_cm); }
        if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }

        /* Transpose C (row-major) inputs to Fortran (column-major) copies */
        if (t_size > 0) slicot_transpose_to_fortran(t, t_cm, t_c_rows, t_c_cols, elem_size);
        if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_c_rows, b_c_cols, elem_size);

        /* Fortran leading dimensions */
        int ldt_f = (t_c_rows > 0) ? t_c_rows : 1;
        int ldb_f = (b_c_rows > 0) ? b_c_rows : 1;

        /* Call the Fortran routine */
        F77_FUNC(mb02ed, MB02ED)(&typet_upper, &k, &n, &nrhs,
                                 t_cm, &ldt_f,           // Pass CM T
                                 b_cm, &ldb_f,           // Pass CM B
                                 dwork, &ldwork, &info,
                                 typet_len);

        /* Copy back results from column-major temps to original row-major arrays */
        if (info == 0) {
             if (t_size > 0) slicot_transpose_to_c(t_cm, t, t_c_rows, t_c_cols, elem_size);
             if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_c_rows, b_c_cols, elem_size);
        }
        /* Column-major copies will be freed in cleanup */

    } else {
        /* --- Column-Major Case --- */

        /* Call the Fortran routine directly with user-provided arrays */
        F77_FUNC(mb02ed, MB02ED)(&typet_upper, &k, &n, &nrhs,
                                 t, &ldt,                // Pass original T
                                 b, &ldb,                // Pass original B
                                 dwork, &ldwork, &info,
                                 typet_len);
        // T and B are modified in place by the Fortran call.
    }

cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(t_cm);
    free(b_cm);

    return info;
}
