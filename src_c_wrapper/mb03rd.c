/**
 * @file mb03rd.c
 * @brief C wrapper implementation for SLICOT routine MB03RD
 *
 * This file provides a C wrapper implementation for the SLICOT routine MB03RD,
 * which reduces a matrix A in real Schur form to block-diagonal form
 * using well-conditioned non-orthogonal similarity transformations.
 */

#include <stdlib.h>
#include <ctype.h>
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "mb03rd.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * Note A is input/output. X is input/output if JOBX='U'.
 * NBLCKS, BLSIZE, WR, WI are output.
 */
extern void F77_FUNC(mb03rd, MB03RD)(
    const char* jobx,       // CHARACTER*1 JOBX
    const char* sort,       // CHARACTER*1 SORT
    const int* n,           // INTEGER N
    const double* pmax,     // DOUBLE PRECISION PMAX
    double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
    const int* lda,         // INTEGER LDA
    double* x,              // DOUBLE PRECISION X(LDX,*) (in/out)
    const int* ldx,         // INTEGER LDX
    int* nblcks,            // INTEGER NBLCKS (output)
    int* blsize,            // INTEGER BLSIZE(*) (output)
    double* wr,             // DOUBLE PRECISION WR(*) (output)
    double* wi,             // DOUBLE PRECISION WI(*) (output)
    const double* tol,      // DOUBLE PRECISION TOL
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    int* info,              // INTEGER INFO (output)
    int jobx_len,           // Hidden length
    int sort_len            // Hidden length
);


/* C wrapper function definition */
int slicot_mb03rd(char jobx, char sort, int n, double pmax,
                  double* a, int lda, double* x, int ldx,
                  int* nblcks, int* blsize, double* wr, double* wi,
                  double tol, int row_major)
{
    /* Local variables */
    int info = 0;
    double* dwork = NULL; // Workspace
    int dwork_size = 0;
    const int jobx_len = 1, sort_len = 1;

    char jobx_upper = toupper(jobx);
    char sort_upper = toupper(sort);

    /* Pointers for column-major copies if needed */
    double *a_cm = NULL;
    double *x_cm = NULL; // Needed if row_major and jobx == 'U'

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -3; goto cleanup; }
    if (pmax < 1.0) { info = -4; goto cleanup; }
    if (jobx_upper != 'N' && jobx_upper != 'U') { info = -1; goto cleanup; }
    if (sort_upper != 'N' && sort_upper != 'S' && sort_upper != 'C' && sort_upper != 'B') {
        info = -2; goto cleanup;
    }

    // Check leading dimensions based on storage order and JOBX
    int min_lda_f = MAX(1, n);
    int min_ldx_f = (jobx_upper == 'U') ? MAX(1, n) : 1;

    if (row_major) {
        // For row-major C, LDA is the number of columns
        int min_lda_rm_cols = n;
        int min_ldx_rm_cols = (jobx_upper == 'U') ? n : 1;
        if (lda < min_lda_rm_cols) { info = -6; goto cleanup; }
        if (ldx < min_ldx_rm_cols) { info = -8; goto cleanup; } // Check even if JOBX='N'? Fortran requires >=1
    } else {
        // For column-major C, LDA is the number of rows (Fortran style)
        if (lda < min_lda_f) { info = -6; goto cleanup; }
        if (ldx < min_ldx_f) { info = -8; goto cleanup; }
    }

    /* --- Workspace Allocation --- */

    // Allocate DWORK (size N) - No query needed
    dwork_size = MAX(1, n); // Ensure minimum size 1
    dwork = (double*)malloc((size_t)dwork_size * sizeof(double));
    CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

    /* --- Prepare Arrays and Call Fortran Routine --- */
    size_t elem_size = sizeof(double);

    if (row_major) {
        /* --- Row-Major Case --- */

        /* Allocate memory for column-major copies */
        size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
        size_t x_rows = n; size_t x_cols = n; size_t x_size = x_rows * x_cols;

        if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
        if (jobx_upper == 'U' && x_size > 0) {
            x_cm = (double*)malloc(x_size * elem_size); CHECK_ALLOC(x_cm);
        }

        /* Transpose C (row-major) inputs to Fortran (column-major) copies */
        if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
        if (jobx_upper == 'U' && x_size > 0) {
            slicot_transpose_to_fortran(x, x_cm, x_rows, x_cols, elem_size);
        }

        /* Fortran leading dimensions */
        int lda_f = (a_rows > 0) ? a_rows : 1;
        int ldx_f = (x_rows > 0) ? x_rows : 1;

        /* Call the Fortran routine */
        F77_FUNC(mb03rd, MB03RD)(&jobx_upper, &sort_upper, &n, &pmax,
                                 a_cm, &lda_f,           // Pass CM A
                                 (jobx_upper == 'U' ? x_cm : NULL), &ldx_f, // Pass CM X or NULL
                                 nblcks, blsize, wr, wi,
                                 &tol, dwork, &info,
                                 jobx_len, sort_len);

        /* Copy back results from column-major temps to original row-major arrays */
        if (info == 0) {
            // Copy back modified A
            if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
            // Copy back modified X if accumulated
            if (jobx_upper == 'U' && x_size > 0) {
                slicot_transpose_to_c(x_cm, x, x_rows, x_cols, elem_size);
            }
            // NBLCKS, BLSIZE, WR, WI are filled directly.
        }
        /* Column-major copies will be freed in cleanup */

    } else {
        /* --- Column-Major Case --- */

        /* Call the Fortran routine directly with user-provided arrays */
        F77_FUNC(mb03rd, MB03RD)(&jobx_upper, &sort_upper, &n, &pmax,
                                 a, &lda,                // Pass original A
                                 (jobx_upper == 'U' ? x : NULL), &ldx, // Pass original X or NULL
                                 nblcks, blsize, wr, wi,
                                 &tol, dwork, &info,
                                 jobx_len, sort_len);
        // A, X (if JOBX='U'), NBLCKS, BLSIZE, WR, WI are modified in place.
    }

cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(a_cm);
    free(x_cm); // Safe even if NULL

    return info;
}
