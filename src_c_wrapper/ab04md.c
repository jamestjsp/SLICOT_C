/**
 * @file ab04md.c
 * @brief C wrapper implementation for SLICOT routine AB04MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB04MD,
 * which performs discrete-time <-> continuous-time system conversion using
 * a bilinear transformation.
 */

#include <stdlib.h>
#include <ctype.h>  // For toupper

// Include the header file for this wrapper
#include "ab04md.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * This handles potential name mangling issues between C and Fortran compilers.
 * All arguments are passed by reference (pointers).
 * The hidden string length argument for 'type' is explicitly included at the end.
 */
extern void F77_FUNC(ab04md, AB04MD)(
    const char* type,       // CHARACTER*1 TYPE
    const int* n,           // INTEGER N
    const int* m,           // INTEGER M
    const int* p,           // INTEGER P
    const double* alpha,    // DOUBLE PRECISION ALPHA
    const double* beta,     // DOUBLE PRECISION BETA
    double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
    const int* lda,         // INTEGER LDA
    double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
    const int* ldb,         // INTEGER LDB
    double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
    const int* ldc,         // INTEGER LDC
    double* d,              // DOUBLE PRECISION D(LDD,*) (in/out)
    const int* ldd,         // INTEGER LDD
    int* iwork,             // INTEGER IWORK(*)
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info,              // INTEGER INFO (output)
    int type_len            // Hidden length argument for type (integer)
);

/* C wrapper function definition */
SLICOT_C_WRAPPER_API
int slicot_ab04md(char type, int n, int m, int p,
                 double alpha, double beta,
                 double* a, int lda,
                 double* b, int ldb,
                 double* c, int ldc,
                 double* d, int ldd,
                 int row_major)
{
    /* Local variables */
    int info = 0;
    int ldwork = -1;  /* Start with workspace query */
    double dwork_query;
    double* dwork = NULL;
    int* iwork = NULL;
    int type_len = 1;  // Fortran expects 1-based length for strings

    char type_upper = toupper(type);

    /* Pointers for column-major copies if needed */
    double *a_cm = NULL;
    double *b_cm = NULL;
    double *c_cm = NULL;
    double *d_cm = NULL;

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -2; goto cleanup; }
    if (m < 0) { info = -3; goto cleanup; }
    if (p < 0) { info = -4; goto cleanup; }

    if (type_upper != 'D' && type_upper != 'C') {
        info = -1; goto cleanup;
    }
    
    // Check for non-zero alpha and beta
    if (alpha == 0.0) { info = -5; goto cleanup; }
    if (beta == 0.0) { info = -6; goto cleanup; }

    // Check leading dimensions based on storage order
    int min_lda_f = MAX(1, n);
    int min_ldb_f = MAX(1, n);
    int min_ldc_f = MAX(1, p);
    int min_ldd_f = MAX(1, p);

    if (row_major) {
        // For row-major C, LDA is the number of columns
        int min_lda_rm_cols = n;
        int min_ldb_rm_cols = m;
        int min_ldc_rm_cols = n;
        int min_ldd_rm_cols = m;
        if (lda < min_lda_rm_cols) { info = -8; goto cleanup; }
        if (ldb < min_ldb_rm_cols) { info = -10; goto cleanup; }
        if (ldc < min_ldc_rm_cols) { info = -12; goto cleanup; }
        if (ldd < min_ldd_rm_cols) { info = -14; goto cleanup; }
    } else {
        // For column-major C, LDA is the number of rows (Fortran style)
        if (lda < min_lda_f) { info = -8; goto cleanup; }
        if (ldb < min_ldb_f) { info = -10; goto cleanup; }
        if (ldc < min_ldc_f) { info = -12; goto cleanup; }
        if (ldd < min_ldd_f) { info = -14; goto cleanup; }
    }

    /* --- Workspace Allocation --- */

    // Allocate IWORK (size N)
    iwork = (int*)malloc((size_t)n * sizeof(int));
    CHECK_ALLOC(iwork);

    // Perform workspace query
    ldwork = -1; // Query mode
    F77_FUNC(ab04md, AB04MD)(&type_upper, &n, &m, &p, &alpha, &beta,
                            a, &lda, b, &ldb, c, &ldc, d, &ldd,
                            iwork, &dwork_query, &ldwork, &info,
                            type_len);

    if (info != 0) {
        // Query failed, use documented minimum size
        ldwork = MAX(1, n);
    } else {
        // Get optimal workspace size from query
        ldwork = (int)dwork_query;
        // Ensure minimum size from documentation
        ldwork = MAX(ldwork, MAX(1, n));
    }

    // Allocate DWORK with optimal/minimum size
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);

    /* --- Prepare Arrays and Call Fortran Routine --- */

    if (row_major) {
        /* --- Row-Major Case --- */

        /* Allocate memory for column-major copies */
        size_t a_rows = n; size_t a_cols = n; size_t a_size = a_rows * a_cols;
        size_t b_rows = n; size_t b_cols = m; size_t b_size = b_rows * b_cols;
        size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
        size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;

        // Allocate column-major copies - only if dimensions > 0
        if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
        if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
        if (c_size > 0) { c_cm = (double*)malloc(c_size * sizeof(double)); CHECK_ALLOC(c_cm); }
        if (d_size > 0) { d_cm = (double*)malloc(d_size * sizeof(double)); CHECK_ALLOC(d_cm); }

        /* Transpose C (row-major) inputs to Fortran (column-major) copies */
        if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, sizeof(double));
        if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, sizeof(double));
        if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, sizeof(double));
        if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, sizeof(double));

        /* Fortran leading dimensions */
        int lda_f = (a_rows > 0) ? a_rows : 1;
        int ldb_f = (b_rows > 0) ? b_rows : 1;
        int ldc_f = (c_rows > 0) ? c_rows : 1;
        int ldd_f = (d_rows > 0) ? d_rows : 1;

        /* Call the Fortran routine with column-major data */
        F77_FUNC(ab04md, AB04MD)(&type_upper, &n, &m, &p, &alpha, &beta,
                                a_cm, &lda_f, b_cm, &ldb_f, c_cm, &ldc_f, d_cm, &ldd_f,
                                iwork, dwork, &ldwork, &info, type_len);

        /* Copy back results from column-major temps to original row-major arrays */
        if (info == 0) { // Only copy back if successful
            if (a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, sizeof(double));
            if (b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, sizeof(double));
            if (c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows, c_cols, sizeof(double));
            if (d_size > 0) slicot_transpose_to_c(d_cm, d, d_rows, d_cols, sizeof(double));
        }
    } else {
        /* --- Column-Major Case --- */
        /* Call the Fortran routine directly with user-provided arrays */
        F77_FUNC(ab04md, AB04MD)(&type_upper, &n, &m, &p, &alpha, &beta,
                                a, &lda, b, &ldb, c, &ldc, d, &ldd,
                                iwork, dwork, &ldwork, &info, type_len);
    }

cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(iwork);
    free(a_cm);
    free(b_cm);
    free(c_cm);
    free(d_cm);

    return info;
}
