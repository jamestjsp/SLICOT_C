/**
 * @file mb05nd.c
 * @brief C wrapper implementation for SLICOT routine MB05ND
 *
 * This file provides a C wrapper implementation for the SLICOT routine MB05ND,
 * which computes the matrix exponential exp(A*delta) and its integral.
 * This version is refactored to align with the ab01nd.c pattern.
 */

#include <stdlib.h>
#include <stddef.h> // For size_t
#include <ctype.h>  // For toupper (though not used here, kept for pattern consistency)
#include <string.h> // For memcpy/memset if needed (usually handled by utils)

// Include the header file for this wrapper
#include "mb05nd.h" // Assuming mb05nd.h exists
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
SLICOT_EXPORT
int slicot_mb05nd(int n, double delta, const double* a, int lda,
                  double* ex, int ldex, double* exint, int ldexin,
                  double tol, int row_major)
{
    /* --- Local Variables --- */
    int info = 0;           // SLICOT routine return status
    int ldwork_actual;      // Workspace size for dwork
    double* dwork_allocated_buffer = NULL;   // Double workspace array
    int* iwork_allocated_buffer = NULL;      // Integer workspace array
    int liwork_actual_size = 0;              // Size of iwork

    /* Pointers for column-major copies (if row_major is true) */
    double *a_cm = NULL, *ex_cm = NULL, *exint_cm = NULL;

    /* Fortran-compatible leading dimensions */
    int lda_f, ldex_f, ldexin_f;

    /* Array dimensions and sizes */
    size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
    size_t ex_rows = n; size_t ex_cols = n; size_t ex_size = ex_rows * ex_cols;
    size_t exint_rows = n; size_t exint_cols = n; size_t exint_size = exint_rows * exint_cols;
    const size_t elem_size = sizeof(double); // Size of a double element

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -1; goto cleanup; }
    // No check for DELTA needed based on docs (except potential overflow handled by INFO=N+1)
    // No check for TOL needed based on docs

    // Check leading dimensions based on storage order
    int min_lda_f = MAX(1, n);
    int min_ldex_f = MAX(1, n);
    int min_ldexin_f = MAX(1, n);

    if (row_major) {
        // For row-major C, LD is the number of columns
        int min_lda_rm_cols = n;
        int min_ldex_rm_cols = n;
        int min_ldexin_rm_cols = n;
        if (lda < min_lda_rm_cols) { info = -4; goto cleanup; }
        if (ldex < min_ldex_rm_cols) { info = -6; goto cleanup; }
        if (ldexin < min_ldexin_rm_cols) { info = -8; goto cleanup; }
    } else {
        // For column-major C, LD is the number of rows (Fortran style)
        if (lda < min_lda_f) { info = -4; goto cleanup; }
        if (ldex < min_ldex_f) { info = -6; goto cleanup; }
        if (ldexin < min_ldexin_f) { info = -8; goto cleanup; }
    }

    /* --- Workspace Allocation --- */

    // LIWORK: Fixed size based on documentation
    liwork_actual_size = MAX(1, n);
    iwork_allocated_buffer = (int*)malloc((size_t)liwork_actual_size * sizeof(int));
    CHECK_ALLOC(iwork_allocated_buffer);

    // LDWORK: Use the formula provided in the documentation
    if (n == 0) {
        ldwork_actual = 1;
    } else {
        ldwork_actual = MAX(1, n * (n + 1)); // Minimum required size
        ldwork_actual = MAX(ldwork_actual, 2 * n * n); // Optimum size
    }

    dwork_allocated_buffer = (double*)malloc((size_t)ldwork_actual * elem_size);
    CHECK_ALLOC(dwork_allocated_buffer);

    /* --- Prepare Arrays and Call Fortran Routine --- */

    if (row_major) {
        /* --- Row-Major Case --- */

        /* Allocate memory for column-major copies */
        if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
        if (ex_size > 0) { ex_cm = (double*)malloc(ex_size * elem_size); CHECK_ALLOC(ex_cm); }
        if (exint_size > 0) { exint_cm = (double*)malloc(exint_size * elem_size); CHECK_ALLOC(exint_cm); }

        /* Transpose C (row-major) input A to Fortran (column-major) copy */
        if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);

        /* Fortran leading dimensions */
        lda_f = (a_rows > 0) ? a_rows : 1;
        ldex_f = (ex_rows > 0) ? ex_rows : 1;
        ldexin_f = (exint_rows > 0) ? exint_rows : 1;

        /* Call the Fortran routine with column-major data */
        F77_FUNC(mb05nd, MB05ND)(&n, &delta,
                                 a_cm, &lda_f,           // Pass CM A (in)
                                 ex_cm, &ldex_f,         // Pass CM EX (out)
                                 exint_cm, &ldexin_f,    // Pass CM EXINT (out)
                                 &tol, iwork_allocated_buffer, dwork_allocated_buffer, &ldwork_actual, &info);

        /* Copy back results from column-major temps to original row-major arrays */
        if (info == 0 || info == (n + 1)) { // Copy back even if INFO = N+1 (overflow warning)
            if (ex_size > 0) slicot_transpose_to_c(ex_cm, ex, ex_rows, ex_cols, elem_size);
            if (exint_size > 0) slicot_transpose_to_c(exint_cm, exint, exint_rows, exint_cols, elem_size);
        }

    } else {
        /* --- Column-Major Case --- */

        /* Fortran leading dimensions are the same as C */
        lda_f = lda;
        ldex_f = ldex;
        ldexin_f = ldexin;

        /* Call the Fortran routine directly with user-provided arrays */
        F77_FUNC(mb05nd, MB05ND)(&n, &delta,
                                 a, &lda_f,              // Pass original A
                                 ex, &ldex_f,            // Pass original EX
                                 exint, &ldexin_f,       // Pass original EXINT
                                 &tol, iwork_allocated_buffer, dwork_allocated_buffer, &ldwork_actual, &info);
        // EX and EXINT are modified in place.
    }

cleanup:
    /* --- Cleanup --- */
    // Free allocated memory (safe to free NULL pointers)
    free(dwork_allocated_buffer);
    free(iwork_allocated_buffer);
    free(a_cm);
    free(ex_cm);
    free(exint_cm);

    /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
    return info;
}
