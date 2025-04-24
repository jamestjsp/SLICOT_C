/**
 * @file ab05md.c
 * @brief C wrapper implementation for SLICOT routine AB05MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB05MD,
 * which computes the state-space model (A,B,C,D) for the cascaded
 * inter-connection of two systems.
 */

#include <stdlib.h>
#include <ctype.h>  // For toupper

// Include the header file for this wrapper
#include "ab05md.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * This handles potential name mangling issues between C and Fortran compilers.
 * All arguments are passed by reference (pointers).
 * The hidden string length arguments are explicitly included at the end.
 */
extern void F77_FUNC(ab05md, AB05MD)(
    const char* uplo,       // CHARACTER*1 UPLO
    const char* over,       // CHARACTER*1 OVER
    const int* n1,          // INTEGER N1
    const int* m1,          // INTEGER M1
    const int* p1,          // INTEGER P1
    const int* n2,          // INTEGER N2
    const int* p2,          // INTEGER P2
    const double* a1,       // DOUBLE PRECISION A1(LDA1,*)
    const int* lda1,        // INTEGER LDA1
    const double* b1,       // DOUBLE PRECISION B1(LDB1,*)
    const int* ldb1,        // INTEGER LDB1
    const double* c1,       // DOUBLE PRECISION C1(LDC1,*)
    const int* ldc1,        // INTEGER LDC1
    const double* d1,       // DOUBLE PRECISION D1(LDD1,*)
    const int* ldd1,        // INTEGER LDD1
    const double* a2,       // DOUBLE PRECISION A2(LDA2,*)
    const int* lda2,        // INTEGER LDA2
    const double* b2,       // DOUBLE PRECISION B2(LDB2,*)
    const int* ldb2,        // INTEGER LDB2
    const double* c2,       // DOUBLE PRECISION C2(LDC2,*)
    const int* ldc2,        // INTEGER LDC2
    const double* d2,       // DOUBLE PRECISION D2(LDD2,*)
    const int* ldd2,        // INTEGER LDD2
    int* n,                 // INTEGER N (output)
    double* a,              // DOUBLE PRECISION A(LDA,*) (output)
    const int* lda,         // INTEGER LDA
    double* b,              // DOUBLE PRECISION B(LDB,*) (output)
    const int* ldb,         // INTEGER LDB
    double* c,              // DOUBLE PRECISION C(LDC,*) (output)
    const int* ldc,         // INTEGER LDC
    double* d,              // DOUBLE PRECISION D(LDD,*) (output)
    const int* ldd,         // INTEGER LDD
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info,              // INTEGER INFO (output)
    int uplo_len,           // Hidden length for uplo
    int over_len            // Hidden length for over
);


/* C wrapper function definition */
SLICOT_EXPORT
int slicot_ab05md(char uplo, char over,
                 int n1, int m1, int p1, int n2, int p2,
                 const double* a1, int lda1, const double* b1, int ldb1,
                 const double* c1, int ldc1, const double* d1, int ldd1,
                 const double* a2, int lda2, const double* b2, int ldb2,
                 const double* c2, int ldc2, const double* d2, int ldd2,
                 int* n,
                 double* a, int lda, double* b, int ldb,
                 double* c, int ldc, double* d, int ldd,
                 int row_major)
{
    /* Local variables */
    int info = 0;
    double* dwork = NULL;
    int ldwork = 1;  // Default for OVER='N'
    int n_calc = n1 + n2;  // Calculated total state dimension
    const int uplo_len = 1;
    const int over_len = 1;

    char uplo_upper = toupper(uplo);
    char over_upper = toupper(over);

    /* Pointers for column-major copies if needed */
    double *a1_cm = NULL, *b1_cm = NULL, *c1_cm = NULL, *d1_cm = NULL;
    double *a2_cm = NULL, *b2_cm = NULL, *c2_cm = NULL, *d2_cm = NULL;
    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;

    /* --- Input Parameter Validation --- */

    // Check input dimensions
    if (n1 < 0) { info = -3; goto cleanup; }
    if (m1 < 0) { info = -4; goto cleanup; }
    if (p1 < 0) { info = -5; goto cleanup; }
    if (n2 < 0) { info = -6; goto cleanup; }
    if (p2 < 0) { info = -7; goto cleanup; }

    // Check UPLO and OVER parameters
    if (uplo_upper != 'U' && uplo_upper != 'L') { 
        info = -1; goto cleanup; 
    }
    if (over_upper != 'N' && over_upper != 'O') { 
        info = -2; goto cleanup; 
    }

    // Set output N
    if (n) { *n = n_calc; }

    // Check leading dimensions based on storage order
    // First, minimum values for Fortran-style (column-major) arrays
    int min_lda1_f = MAX(1, n1);
    int min_ldb1_f = MAX(1, n1);
    int min_ldc1_f = (n1 > 0) ? MAX(1, p1) : 1;
    int min_ldd1_f = MAX(1, p1);
    int min_lda2_f = MAX(1, n2);
    int min_ldb2_f = MAX(1, n2);
    int min_ldc2_f = (n2 > 0) ? MAX(1, p2) : 1;
    int min_ldd2_f = MAX(1, p2);
    int min_lda_f = MAX(1, n_calc);
    int min_ldb_f = MAX(1, n_calc);
    int min_ldc_f = (n_calc > 0) ? MAX(1, p2) : 1;
    int min_ldd_f = MAX(1, p2);

    if (row_major) {
        // For row-major C arrays, check number of columns
        if (lda1 < n1) { info = -9; goto cleanup; }
        if (ldb1 < m1) { info = -11; goto cleanup; }
        if (ldc1 < n1) { info = -13; goto cleanup; }
        if (ldd1 < m1) { info = -15; goto cleanup; }
        if (lda2 < n2) { info = -17; goto cleanup; }
        if (ldb2 < p1) { info = -19; goto cleanup; }
        if (ldc2 < n2) { info = -21; goto cleanup; }
        if (ldd2 < p1) { info = -23; goto cleanup; }
        // Check output array dimensions
        if (lda < n_calc) { info = -26; goto cleanup; }
        if (ldb < m1) { info = -28; goto cleanup; }
        if (ldc < n_calc) { info = -30; goto cleanup; }
        if (ldd < m1) { info = -32; goto cleanup; }
    } else {
        // For column-major C arrays, check standard Fortran leading dimensions
        if (lda1 < min_lda1_f) { info = -9; goto cleanup; }
        if (ldb1 < min_ldb1_f) { info = -11; goto cleanup; }
        if (ldc1 < min_ldc1_f) { info = -13; goto cleanup; }
        if (ldd1 < min_ldd1_f) { info = -15; goto cleanup; }
        if (lda2 < min_lda2_f) { info = -17; goto cleanup; }
        if (ldb2 < min_ldb2_f) { info = -19; goto cleanup; }
        if (ldc2 < min_ldc2_f) { info = -21; goto cleanup; }
        if (ldd2 < min_ldd2_f) { info = -23; goto cleanup; }
        // Check output array dimensions
        if (lda < min_lda_f) { info = -26; goto cleanup; }
        if (ldb < min_ldb_f) { info = -28; goto cleanup; }
        if (ldc < min_ldc_f) { info = -30; goto cleanup; }
        if (ldd < min_ldd_f) { info = -32; goto cleanup; }
    }

    /* --- Workspace Allocation --- */
    // Workspace is needed if OVER='O'
    if (over_upper == 'O') {
        // According to docs: LDWORK >= MAX(1, P1*MAX(N1, M1, N2, P2))
        int max_dim = MAX(MAX(n1, m1), MAX(n2, p2));
        ldwork = MAX(1, p1 * max_dim);
        
        dwork = (double*)malloc((size_t)ldwork * sizeof(double));
        CHECK_ALLOC(dwork);
    }

    /* --- Prepare Arrays and Call Fortran Routine --- */
    if (row_major) {
        /* --- Row-Major Case --- */
        
        /* Allocate memory for column-major copies */
        size_t a1_size = (size_t)n1 * n1;
        size_t b1_size = (size_t)n1 * m1;
        size_t c1_size = (size_t)p1 * n1;
        size_t d1_size = (size_t)p1 * m1;
        size_t a2_size = (size_t)n2 * n2;
        size_t b2_size = (size_t)n2 * p1;
        size_t c2_size = (size_t)p2 * n2;
        size_t d2_size = (size_t)p2 * p1;
        
        // Allocate and transpose input arrays
        if (a1_size > 0) { a1_cm = (double*)malloc(a1_size * sizeof(double)); CHECK_ALLOC(a1_cm); }
        if (b1_size > 0) { b1_cm = (double*)malloc(b1_size * sizeof(double)); CHECK_ALLOC(b1_cm); }
        if (c1_size > 0) { c1_cm = (double*)malloc(c1_size * sizeof(double)); CHECK_ALLOC(c1_cm); }
        if (d1_size > 0) { d1_cm = (double*)malloc(d1_size * sizeof(double)); CHECK_ALLOC(d1_cm); }
        if (a2_size > 0) { a2_cm = (double*)malloc(a2_size * sizeof(double)); CHECK_ALLOC(a2_cm); }
        if (b2_size > 0) { b2_cm = (double*)malloc(b2_size * sizeof(double)); CHECK_ALLOC(b2_cm); }
        if (c2_size > 0) { c2_cm = (double*)malloc(c2_size * sizeof(double)); CHECK_ALLOC(c2_cm); }
        if (d2_size > 0) { d2_cm = (double*)malloc(d2_size * sizeof(double)); CHECK_ALLOC(d2_cm); }
        
        // Allocate memory for column-major outputs
        size_t a_size = (size_t)n_calc * n_calc;
        size_t b_size = (size_t)n_calc * m1;
        size_t c_size = (size_t)p2 * n_calc;
        size_t d_size = (size_t)p2 * m1;
        
        if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
        if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
        if (c_size > 0) { c_cm = (double*)malloc(c_size * sizeof(double)); CHECK_ALLOC(c_cm); }
        if (d_size > 0) { d_cm = (double*)malloc(d_size * sizeof(double)); CHECK_ALLOC(d_cm); }

        /* Transpose row-major inputs to column-major */
        if (a1_size > 0) slicot_transpose_to_fortran(a1, a1_cm, n1, n1, sizeof(double));
        if (b1_size > 0) slicot_transpose_to_fortran(b1, b1_cm, n1, m1, sizeof(double));
        if (c1_size > 0) slicot_transpose_to_fortran(c1, c1_cm, p1, n1, sizeof(double));
        if (d1_size > 0) slicot_transpose_to_fortran(d1, d1_cm, p1, m1, sizeof(double));
        if (a2_size > 0) slicot_transpose_to_fortran(a2, a2_cm, n2, n2, sizeof(double));
        if (b2_size > 0) slicot_transpose_to_fortran(b2, b2_cm, n2, p1, sizeof(double));
        if (c2_size > 0) slicot_transpose_to_fortran(c2, c2_cm, p2, n2, sizeof(double));
        if (d2_size > 0) slicot_transpose_to_fortran(d2, d2_cm, p2, p1, sizeof(double));

        /* Fortran leading dimensions */
        int lda1_f = (n1 > 0) ? n1 : 1;
        int ldb1_f = (n1 > 0) ? n1 : 1;
        int ldc1_f = (p1 > 0) ? p1 : 1;
        int ldd1_f = (p1 > 0) ? p1 : 1;
        int lda2_f = (n2 > 0) ? n2 : 1;
        int ldb2_f = (n2 > 0) ? n2 : 1;
        int ldc2_f = (p2 > 0) ? p2 : 1;
        int ldd2_f = (p2 > 0) ? p2 : 1;
        int lda_f = (n_calc > 0) ? n_calc : 1;
        int ldb_f = (n_calc > 0) ? n_calc : 1;
        int ldc_f = (p2 > 0) ? p2 : 1;
        int ldd_f = (p2 > 0) ? p2 : 1;

        /* Call the Fortran routine */
        F77_FUNC(ab05md, AB05MD)(&uplo_upper, &over_upper,
                                &n1, &m1, &p1, &n2, &p2,
                                a1_cm, &lda1_f, b1_cm, &ldb1_f,
                                c1_cm, &ldc1_f, d1_cm, &ldd1_f,
                                a2_cm, &lda2_f, b2_cm, &ldb2_f,
                                c2_cm, &ldc2_f, d2_cm, &ldd2_f,
                                &n_calc,
                                a_cm, &lda_f, b_cm, &ldb_f,
                                c_cm, &ldc_f, d_cm, &ldd_f,
                                dwork, &ldwork, &info,
                                uplo_len, over_len);

        /* Copy results back to row-major arrays */
        if (info == 0) {
            if (a_size > 0) slicot_transpose_to_c(a_cm, a, n_calc, n_calc, sizeof(double));
            if (b_size > 0) slicot_transpose_to_c(b_cm, b, n_calc, m1, sizeof(double));
            if (c_size > 0) slicot_transpose_to_c(c_cm, c, p2, n_calc, sizeof(double));
            if (d_size > 0) slicot_transpose_to_c(d_cm, d, p2, m1, sizeof(double));
        }
    } else {
        /* --- Column-Major Case --- */
        
        /* Call Fortran routine directly */
        F77_FUNC(ab05md, AB05MD)(&uplo_upper, &over_upper,
                                &n1, &m1, &p1, &n2, &p2,
                                a1, &lda1, b1, &ldb1,
                                c1, &ldc1, d1, &ldd1,
                                a2, &lda2, b2, &ldb2,
                                c2, &ldc2, d2, &ldd2,
                                &n_calc,
                                a, &lda, b, &ldb,
                                c, &ldc, d, &ldd,
                                dwork, &ldwork, &info,
                                uplo_len, over_len);
    }

cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(a1_cm); free(b1_cm); free(c1_cm); free(d1_cm);
    free(a2_cm); free(b2_cm); free(c2_cm); free(d2_cm);
    free(a_cm); free(b_cm); free(c_cm); free(d_cm);

    return info;
}
