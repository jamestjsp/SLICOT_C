/**
 * @file mb05nd.c
 * @brief C wrapper implementation for SLICOT routine MB05ND
 *
 * This file provides a C wrapper implementation for the SLICOT routine MB05ND,
 * which computes the matrix exponential exp(A*delta) and its integral.
 */

#include <stdlib.h>
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "mb05nd.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * Note A is input only. EX, EXINT are output.
 */
extern void F77_FUNC(mb05nd, MB05ND)(
    const int* n,           // INTEGER N
    const double* delta,    // DOUBLE PRECISION DELTA
    const double* a,        // DOUBLE PRECISION A(LDA,*)
    const int* lda,         // INTEGER LDA
    double* ex,             // DOUBLE PRECISION EX(LDEX,*) (output)
    const int* ldex,        // INTEGER LDEX
    double* exint,          // DOUBLE PRECISION EXINT(LDEXIN,*) (output)
    const int* ldexin,      // INTEGER LDEXIN
    const double* tol,      // DOUBLE PRECISION TOL
    int* iwork,             // INTEGER IWORK(*)
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info               // INTEGER INFO (output)
);


/* C wrapper function definition */
SLICOT_C_WRAPPER_API
int slicot_mb05nd(int n, double delta, const double* a, int lda,
                  double* ex, int ldex, double* exint, int ldexin,
                  double tol, int row_major)
{
    /* Local variables */
    int info = 0;
    int ldwork = -1; /* Use -1 for workspace query */
    double dwork_query;
    double* dwork = NULL;
    int* iwork = NULL;
    int iwork_size = 0;

    /* Pointers for column-major copies if needed */
    double *a_cm = NULL, *ex_cm = NULL, *exint_cm = NULL;

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -1; goto cleanup; }
    // No check for DELTA needed based on docs (except potential overflow handled by INFO=N+1)
    // No check for TOL needed based on docs

    // Check leading dimensions based on storage order
    int min_lda_f = MAX(1, n);
    int min_ldex_f = MAX(1, n);
    int min_ldexin_f = MAX(1, n);

    if (row_major) {
        // For row-major C, LDA is the number of columns
        int min_lda_rm_cols = n;
        int min_ldex_rm_cols = n;
        int min_ldexin_rm_cols = n;
        if (lda < min_lda_rm_cols) { info = -4; goto cleanup; }
        if (ldex < min_ldex_rm_cols) { info = -6; goto cleanup; }
        if (ldexin < min_ldexin_rm_cols) { info = -8; goto cleanup; }
    } else {
        // For column-major C, LDA is the number of rows (Fortran style)
        if (lda < min_lda_f) { info = -4; goto cleanup; }
        if (ldex < min_ldex_f) { info = -6; goto cleanup; }
        if (ldexin < min_ldexin_f) { info = -8; goto cleanup; }
    }

    /* --- Workspace Allocation --- */

    // Allocate IWORK (size N)
    iwork_size = MAX(1, n);
    iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
    CHECK_ALLOC(iwork);

    // Allocate DWORK based on query
    ldwork = -1; // Query mode
    F77_FUNC(mb05nd, MB05ND)(&n, &delta, a, &lda, ex, &ldex, exint, &ldexin,
                             &tol, iwork, &dwork_query, &ldwork, &info);

    if (info < 0) { goto cleanup; } // Query failed due to invalid argument
    info = 0; // Reset info after query

    // Get the required dwork size from query result
    ldwork = (int)dwork_query;
    // Check against minimum documented size: MAX(1, N*(N+1))
    int min_ldwork = MAX(1, n * (n + 1));
    ldwork = MAX(ldwork, min_ldwork);

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

    /* --- Prepare Arrays and Call Fortran Routine --- */
    size_t elem_size = sizeof(double);

    if (row_major) {
        /* --- Row-Major Case --- */

        /* Allocate memory for column-major copies */
        size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
        size_t ex_rows = n; size_t ex_cols = n; size_t ex_size = ex_rows * ex_cols;
        size_t exint_rows = n; size_t exint_cols = n; size_t exint_size = exint_rows * exint_cols;

        if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
        if (ex_size > 0) { ex_cm = (double*)malloc(ex_size * elem_size); CHECK_ALLOC(ex_cm); }
        if (exint_size > 0) { exint_cm = (double*)malloc(exint_size * elem_size); CHECK_ALLOC(exint_cm); }

        /* Transpose C (row-major) input A to Fortran (column-major) copy */
        if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);

        /* Fortran leading dimensions */
        int lda_f = (a_rows > 0) ? a_rows : 1;
        int ldex_f = (ex_rows > 0) ? ex_rows : 1;
        int ldexin_f = (exint_rows > 0) ? exint_rows : 1;

        /* Call the Fortran routine */
        F77_FUNC(mb05nd, MB05ND)(&n, &delta,
                                 a_cm, &lda_f,           // Pass CM A
                                 ex_cm, &ldex_f,         // Pass CM EX (out)
                                 exint_cm, &ldexin_f,    // Pass CM EXINT (out)
                                 &tol, iwork, dwork, &ldwork, &info);

        /* Copy back results from column-major temps to original row-major arrays */
        if (info == 0 || info == (n + 1)) { // Copy back even if INFO = N+1 (overflow warning)
            if (ex_size > 0) slicot_transpose_to_c(ex_cm, ex, ex_rows, ex_cols, elem_size);
            if (exint_size > 0) slicot_transpose_to_c(exint_cm, exint, exint_rows, exint_cols, elem_size);
        }
        /* Column-major copies will be freed in cleanup */

    } else {
        /* --- Column-Major Case --- */

        /* Call the Fortran routine directly with user-provided arrays */
        F77_FUNC(mb05nd, MB05ND)(&n, &delta,
                                 a, &lda,                // Pass original A
                                 ex, &ldex,              // Pass original EX
                                 exint, &ldexin,         // Pass original EXINT
                                 &tol, iwork, dwork, &ldwork, &info);
        // EX and EXINT are modified in place.
    }

cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(iwork);
    free(a_cm);
    free(ex_cm);
    free(exint_cm);

    return info;
}
