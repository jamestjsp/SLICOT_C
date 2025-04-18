/**
 * @file mb05md.c
 * @brief C wrapper implementation for SLICOT routine MB05MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine MB05MD,
 * which computes the matrix exponential exp(A*delta) for a real
 * non-defective matrix A using eigenvalue decomposition.
 */

#include <stdlib.h>
#include <ctype.h>
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "mb05md.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * Note A is input/output. V, Y, VALR, VALI are output.
 */
extern void F77_FUNC(mb05md, MB05MD)(
    const char* balanc,     // CHARACTER*1 BALANC
    const int* n,           // INTEGER N
    const double* delta,    // DOUBLE PRECISION DELTA
    double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
    const int* lda,         // INTEGER LDA
    double* v,              // DOUBLE PRECISION V(LDV,*) (output)
    const int* ldv,         // INTEGER LDV
    double* y,              // DOUBLE PRECISION Y(LDY,*) (output)
    const int* ldy,         // INTEGER LDY
    double* valr,           // DOUBLE PRECISION VALR(*) (output)
    double* vali,           // DOUBLE PRECISION VALI(*) (output)
    int* iwork,             // INTEGER IWORK(*)
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info,              // INTEGER INFO (output)
    int balanc_len          // Hidden length
);


/* C wrapper function definition */
int slicot_mb05md(char balanc, int n, double delta,
                  double* a, int lda,
                  double* v, int ldv,
                  double* y, int ldy,
                  double* valr, double* vali,
                  int row_major)
{
    /* Local variables */
    int info = 0;
    int ldwork = -1; /* Use -1 for workspace query */
    double dwork_query;
    double* dwork = NULL;
    int* iwork = NULL;
    int iwork_size = 0;
    const int balanc_len = 1;

    char balanc_upper = toupper(balanc);

    /* Pointers for column-major copies if needed */
    double *a_cm = NULL, *v_cm = NULL, *y_cm = NULL;

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -2; goto cleanup; }
    if (balanc_upper != 'N' && balanc_upper != 'S') { info = -1; goto cleanup; }

    // Check leading dimensions based on storage order
    int min_lda_f = MAX(1, n);
    int min_ldv_f = MAX(1, n);
    int min_ldy_f = MAX(1, n);

    if (row_major) {
        // For row-major C, LDA is the number of columns
        int min_lda_rm_cols = n;
        int min_ldv_rm_cols = n;
        int min_ldy_rm_cols = n;
        if (lda < min_lda_rm_cols) { info = -5; goto cleanup; }
        if (ldv < min_ldv_rm_cols) { info = -7; goto cleanup; }
        if (ldy < min_ldy_rm_cols) { info = -9; goto cleanup; }
    } else {
        // For column-major C, LDA is the number of rows (Fortran style)
        if (lda < min_lda_f) { info = -5; goto cleanup; }
        if (ldv < min_ldv_f) { info = -7; goto cleanup; }
        if (ldy < min_ldy_f) { info = -9; goto cleanup; }
    }

    /* --- Workspace Allocation --- */

    // Allocate IWORK (size N)
    iwork_size = MAX(1, n);
    iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
    CHECK_ALLOC(iwork);

    // Allocate DWORK based on query
    ldwork = -1; // Query mode
    F77_FUNC(mb05md, MB05MD)(&balanc_upper, &n, &delta, a, &lda, v, &ldv, y, &ldy,
                             valr, vali, iwork, &dwork_query, &ldwork, &info,
                             balanc_len);

    if (info < 0) { goto cleanup; } // Query failed due to invalid argument
    info = 0; // Reset info after query

    // Get the required dwork size from query result
    ldwork = (int)dwork_query;
    // Check against minimum documented size: MAX(1, 4*N)
    int min_ldwork = MAX(1, 4 * n);
    ldwork = MAX(ldwork, min_ldwork);

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

    /* --- Prepare Arrays and Call Fortran Routine --- */
    size_t elem_size = sizeof(double);

    if (row_major) {
        /* --- Row-Major Case --- */

        /* Allocate memory for column-major copies */
        size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
        size_t v_rows = n; size_t v_cols = n; size_t v_size = v_rows * v_cols;
        size_t y_rows = n; size_t y_cols = n; size_t y_size = y_rows * y_cols;

        if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
        if (v_size > 0) { v_cm = (double*)malloc(v_size * elem_size); CHECK_ALLOC(v_cm); }
        if (y_size > 0) { y_cm = (double*)malloc(y_size * elem_size); CHECK_ALLOC(y_cm); }

        /* Transpose C (row-major) input A to Fortran (column-major) copy */
        if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);

        /* Fortran leading dimensions */
        int lda_f = (a_rows > 0) ? a_rows : 1;
        int ldv_f = (v_rows > 0) ? v_rows : 1;
        int ldy_f = (y_rows > 0) ? y_rows : 1;

        /* Call the Fortran routine */
        F77_FUNC(mb05md, MB05MD)(&balanc_upper, &n, &delta,
                                 a_cm, &lda_f,           // Pass CM A (in/out)
                                 v_cm, &ldv_f,           // Pass CM V (out)
                                 y_cm, &ldy_f,           // Pass CM Y (out)
                                 valr, vali, iwork, dwork, &ldwork, &info,
                                 balanc_len);

        /* Copy back results from column-major temps to original row-major arrays */
        if (info == 0 || info == (n + 1) || info == (n + 2)) { // Copy back even on warnings/errors N+1, N+2
            if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size); // Result exp(A*delta)
            if (v_size > 0) slicot_transpose_to_c(v_cm, v, v_rows, v_cols, elem_size); // Eigenvectors
            if (y_size > 0) slicot_transpose_to_c(y_cm, y, y_rows, y_cols, elem_size); // Intermediate Y
            // VALR, VALI are filled directly.
        }
        /* Column-major copies will be freed in cleanup */

    } else {
        /* --- Column-Major Case --- */

        /* Call the Fortran routine directly with user-provided arrays */
        F77_FUNC(mb05md, MB05MD)(&balanc_upper, &n, &delta,
                                 a, &lda,                // Pass original A
                                 v, &ldv,                // Pass original V
                                 y, &ldy,                // Pass original Y
                                 valr, vali, iwork, dwork, &ldwork, &info,
                                 balanc_len);
        // A, V, Y, VALR, VALI are modified in place.
    }

cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(iwork);
    free(a_cm);
    free(v_cm);
    free(y_cm);

    return info;
}
