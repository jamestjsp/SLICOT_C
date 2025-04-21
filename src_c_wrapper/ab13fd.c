/**
 * @file ab13fd.c
 * @brief C wrapper implementation for SLICOT routine AB13FD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB13FD,
 * which computes the distance from a real matrix A to the nearest
 * complex matrix with an eigenvalue on the imaginary axis, using SVD.
 * Refactored to align with ab01nd.c structure.
 */

#include <stdlib.h>
#include <complex.h> // For creal with C99 complex types
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "ab13fd.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_complex_double, SLICOT_COMPLEX_REAL
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * Note A is input only (const). BETA, OMEGA are output.
 * Uses COMPLEX*16 workspace CWORK.
 */
extern void F77_FUNC(ab13fd, AB13FD)(
    const int* n,           // INTEGER N
    const double* a,        // DOUBLE PRECISION A(LDA,*)
    const int* lda,         // INTEGER LDA
    double* beta,           // DOUBLE PRECISION BETA (output)
    double* omega,          // DOUBLE PRECISION OMEGA (output)
    const double* tol,      // DOUBLE PRECISION TOL
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    slicot_complex_double* cwork, // COMPLEX*16 CWORK(*)
    const int* lcwork,      // INTEGER LCWORK
    int* info               // INTEGER INFO (output)
);


/* C wrapper function definition */
SLICOT_C_WRAPPER_API
int slicot_ab13fd(int n, const double* a, int lda,
                  double* beta, double* omega, double tol,
                  int row_major)
{
    /* Local variables */
    int info = 0;
    int ldwork = -1; /* Use -1 for workspace query */
    int lcwork = -1; /* Use -1 for workspace query */
    double dwork_query;
    slicot_complex_double cwork_query;
    double* dwork = NULL;
    slicot_complex_double* cwork = NULL;

    /* Pointers for column-major copies if needed */
    double *a_cm = NULL;

    /* Pointers to pass to Fortran */
    const double *a_ptr;
    int lda_f;

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -1; goto cleanup; }
    // TOL check: Fortran routine uses eps if tol < eps.
    // No explicit check needed here unless we want to enforce tol >= 0.
    if (tol < 0.0) { info = -6; goto cleanup; }

    // Check leading dimensions based on storage order
    int min_lda_f = MAX(1, n);

    if (row_major) {
        // For row-major C, LDA is the number of columns
        int min_lda_rm_cols = n;
        if (n > 0 && lda < min_lda_rm_cols) { info = -3; goto cleanup; }
    } else {
        // For column-major C, LDA is the number of rows (Fortran style)
        if (lda < min_lda_f) { info = -3; goto cleanup; }
    }

    /* --- Workspace Allocation --- */

    // Query DWORK and CWORK sizes
    ldwork = -1; // Query mode
    lcwork = -1; // Query mode
    // Use dummy LDs for query if dimensions are 0
    int lda_q = row_major ? MAX(1, n) : lda;

    F77_FUNC(ab13fd, AB13FD)(&n,
                             NULL, &lda_q,          // NULL array for query
                             beta, omega, &tol,
                             &dwork_query, &ldwork,
                             &cwork_query, &lcwork, &info);

    if (info < 0) { goto cleanup; } // Query failed due to invalid argument
    info = 0; // Reset info after query

    // Get the required workspace sizes from query results
    ldwork = (int)dwork_query;
    lcwork = (int)SLICOT_COMPLEX_REAL(cwork_query); // Use macro to get real part

    // Check against minimum documented sizes
    int min_ldwork = MAX(1, 3 * n * (n + 2));
    int min_lcwork = MAX(1, n * (n + 3));
    ldwork = MAX(ldwork, min_ldwork);
    lcwork = MAX(lcwork, min_lcwork);

    // Allocate workspaces
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);
    cwork = (slicot_complex_double*)malloc((size_t)lcwork * sizeof(slicot_complex_double));
    CHECK_ALLOC(cwork);

    /* --- Prepare Arrays and Call Fortran Routine --- */
    size_t elem_size = sizeof(double);

    if (row_major) {
        /* --- Row-Major Case --- */

        /* Allocate memory for column-major copy of A */
        size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
        if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }

        /* Transpose C (row-major) input A to Fortran (column-major) copy */
        if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);

        /* Fortran leading dimension */
        lda_f = MAX(1, a_rows);

        /* Set pointers for Fortran call */
        a_ptr = a_cm;

    } else {
        /* --- Column-Major Case --- */
        lda_f = lda;
        a_ptr = a;
    }

    /* Call the computational routine */
    F77_FUNC(ab13fd, AB13FD)(&n,
                             a_ptr, &lda_f,           // Pass A ptr
                             beta, omega, &tol,     // Pass output pointers, address of tol
                             dwork, &ldwork, cwork, &lcwork, &info); // Pass workspaces

    /* No copy-back needed for A as it's input only */
    /* BETA and OMEGA are filled directly */

cleanup:
    /* --- Cleanup --- */
    free(cwork);
    free(dwork);
    free(a_cm); // Safe even if NULL

    return info;
}
