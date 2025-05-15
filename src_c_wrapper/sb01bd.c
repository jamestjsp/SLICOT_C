/**
 * @file sb01bd.c
 * @brief C wrapper implementation for SLICOT routine SB01BD
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB01BD,
 * which determines the state feedback matrix F for a given system (A,B)
 * such that the closed-loop state matrix A+B*F has specified eigenvalues.
 */

#include <stdlib.h>
#include <ctype.h>
#include <stddef.h> // For size_t
#include <stdio.h>  // For fprintf, stderr

// Include the header file for this wrapper
#include "sb01bd.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * A, WR, WI are input/output. B is input.
 * NFP, NAP, NUP, F, Z, IWARN are output.
 */
extern void F77_FUNC(sb01bd, SB01BD)(
    const char* dico,       // CHARACTER*1 DICO
    const int* n,           // INTEGER N
    const int* m,           // INTEGER M
    const int* np,          // INTEGER NP
    const double* alpha,    // DOUBLE PRECISION ALPHA
    double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
    const int* lda,         // INTEGER LDA
    const double* b,        // DOUBLE PRECISION B(LDB,*) (in)
    const int* ldb,         // INTEGER LDB
    double* wr,             // DOUBLE PRECISION WR(*) (in/out)
    double* wi,             // DOUBLE PRECISION WI(*) (in/out)
    int* nfp,               // INTEGER NFP (output)
    int* nap,               // INTEGER NAP (output)
    int* nup,               // INTEGER NUP (output)
    double* f,              // DOUBLE PRECISION F(LDF,*) (output)
    const int* ldf,         // INTEGER LDF
    double* z,              // DOUBLE PRECISION Z(LDZ,*) (output)
    const int* ldz,         // INTEGER LDZ
    const double* tol,      // DOUBLE PRECISION TOL
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* iwarn,             // INTEGER IWARN (output)
    int* info,              // INTEGER INFO (output)
    int dico_len            // Hidden length
);


/* C wrapper function definition */
SLICOT_EXPORT
int slicot_sb01bd(char dico, int n, int m, int np, double alpha,
                  double* a, int lda, const double* b, int ldb,
                  double* wr, double* wi,
                  int* nfp, int* nap, int* nup,
                  double* f, int ldf, double* z, int ldz,
                  double tol, int* iwarn, int row_major)
{
    /* Local variables */
    int info = 0;
    double* dwork = NULL;
    int ldwork = 0;
    // No iwork needed for this routine
    const int dico_len = 1;

    char dico_upper = toupper(dico);

    /* Pointers for column-major copies if needed */
    double *a_cm = NULL;
    double *b_cm = NULL; // b is const input, so b_cm is a copy
    double *f_cm = NULL;
    double *z_cm = NULL;

    double* a_ptr;
    const double* b_ptr; 
    double* f_ptr;
    double* z_ptr;

    int lda_f, ldb_f, ldf_f, ldz_f;

    /* --- Input Parameter Validation --- */
    /* Perform ALL validation BEFORE allocating any memory */

    if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
    if (n < 0) { info = -2; goto cleanup; }
    if (m < 0) { info = -3; goto cleanup; }
    if (np < 0 || (n > 0 && np > n) || (n == 0 && np > 0) ) { info = -4; goto cleanup; }
    if (dico_upper == 'D' && alpha < 0.0) { info = -5; goto cleanup; } 

    // Check required pointers only if dimensions are non-zero
    if (n > 0 && a == NULL) { info = -6; goto cleanup; }
    if (n > 0 && m > 0 && b == NULL) { info = -8; goto cleanup; }
    if (np > 0 && (wr == NULL || wi == NULL)) { info = -10; goto cleanup; } // Assuming WR is 10th, WI is 11th
    if (n > 0 && m > 0 && f == NULL) { info = -15; goto cleanup; } // F is 15th
    if (n > 0 && z == NULL) { info = -17; goto cleanup; } // Z is 17th
    if (iwarn == NULL) { info = -21; goto cleanup; } // IWARN is 21st
    if (nfp == NULL || nap == NULL || nup == NULL) { info = -12; goto cleanup; } // NFP is 12th

    // Check leading dimensions based on storage order and dimensions
    int min_lda_f_rows = MAX(1, n);
    int min_ldb_f_rows = MAX(1, n);
    int min_ldf_f_rows = MAX(1, m);
    int min_ldz_f_rows = MAX(1, n);

    if (row_major) {
        // For row-major C, LDA is the number of columns
        if (n > 0 && lda < n) { info = -7; goto cleanup; }
        if (m > 0 && ldb < m) { info = -9; goto cleanup; } // B is N x M, C LDB is cols of B
        if (n > 0 && ldf < n) { info = -16; goto cleanup; } // F is M x N, C LDF is cols of F
        if (n > 0 && ldz < n) { info = -18; goto cleanup; }
    } else {
        // For column-major C, LDA is the number of rows (Fortran style)
        if (lda < min_lda_f_rows) { info = -7; goto cleanup; }
        if (ldb < min_ldb_f_rows) { info = -9; goto cleanup; }
        if (ldf < min_ldf_f_rows) { info = -16; goto cleanup; }
        if (ldz < min_ldz_f_rows) { info = -18; goto cleanup; }
    }
    
    if (info != 0) { goto cleanup; }

    /* --- Workspace allocation (Method B: Formula from documentation) --- */
    // LDWORK >= MAX( 1,5*M,5*N,2*N+4*M ).
    ldwork = 1; // Start with MAX(1, ...)
    if (m > 0) {
        ldwork = MAX(ldwork, 5 * m);
        if (n > 0) { // Only consider terms with N if N > 0
            ldwork = MAX(ldwork, 2 * n + 4 * m);
        }
    }
    if (n > 0) {
        ldwork = MAX(ldwork, 5 * n);
    }
    // Ensure ldwork is at least 1 even for zero dimensions if formula resulted in 0
    ldwork = MAX(1, ldwork);

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);

    /* --- Prepare arrays for column-major format if using row-major --- */
    size_t elem_size = sizeof(double);
    size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;
    size_t b_size = (size_t)n * m; if (n == 0 || m == 0) b_size = 0;
    size_t f_size = (size_t)m * n; if (m == 0 || n == 0) f_size = 0;
    size_t z_size = (size_t)n * n; if (n == 0) z_size = 0;

    if (row_major) {
        /* Allocate memory for column-major copies */
        if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
        if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
        if (f_size > 0) { f_cm = (double*)malloc(f_size * elem_size); CHECK_ALLOC(f_cm); }
        if (z_size > 0) { z_cm = (double*)malloc(z_size * elem_size); CHECK_ALLOC(z_cm); }

        /* Fortran leading dimensions (number of rows) */
        lda_f = MAX(1, n);
        ldb_f = MAX(1, n);
        ldf_f = MAX(1, m);
        ldz_f = MAX(1, n);

        /* Transpose C (row-major) inputs to Fortran (column-major) copies */
        if (a_size > 0) slicot_transpose_to_fortran_with_ld(a, a_cm, n, n, lda, lda_f, elem_size);
        if (b_size > 0) slicot_transpose_to_fortran_with_ld(b, b_cm, n, m, ldb, ldb_f, elem_size);

        /* Set pointers to _cm buffers */
        a_ptr = (a_size > 0) ? a_cm : NULL;
        b_ptr = (b_size > 0) ? b_cm : NULL;
        f_ptr = (f_size > 0) ? f_cm : NULL;
        z_ptr = (z_size > 0) ? z_cm : NULL;

    } else {
        /* Column-major case - use original arrays and C LDs */
        lda_f = lda;
        ldb_f = ldb;
        ldf_f = ldf;
        ldz_f = ldz;

        a_ptr = (a_size > 0) ? a : NULL;
        b_ptr = (b_size > 0) ? b : NULL;
        f_ptr = (f_size > 0) ? f : NULL;
        z_ptr = (z_size > 0) ? z : NULL;
    }
     // WR, WI are 1D, no transpose needed, pass directly.
     // For zero-sized arrays, pointers should be NULL.
    double* wr_ptr = (np > 0) ? wr : NULL;
    double* wi_ptr = (np > 0) ? wi : NULL;

    /* --- Call the computational routine --- */
    F77_FUNC(sb01bd, SB01BD)(&dico_upper, &n, &m, &np, &alpha,
                             a_ptr, &lda_f,
                             b_ptr, &ldb_f,
                             wr_ptr, wi_ptr,
                             nfp, nap, nup,
                             f_ptr, &ldf_f,
                             z_ptr, &ldz_f,
                             &tol, dwork, &ldwork, iwarn, &info,
                             dico_len);

    /* --- Copy results back to row-major format if needed --- */
    if (row_major && (info == 0 || info == 3 || info == 4)) { 
        if (a_size > 0 && a_cm != NULL && a != NULL) slicot_transpose_to_c_with_ld(a_cm, a, n, n, lda_f, lda, elem_size);
        if (f_size > 0 && f_cm != NULL && f != NULL) slicot_transpose_to_c_with_ld(f_cm, f, m, n, ldf_f, ldf, elem_size);
        if (z_size > 0 && z_cm != NULL && z != NULL) slicot_transpose_to_c_with_ld(z_cm, z, n, n, ldz_f, ldz, elem_size);
    }

cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(a_cm);
    free(b_cm); 
    free(f_cm);
    free(z_cm);

    if (info == SLICOT_MEMORY_ERROR) {
       fprintf(stderr, "Error: Memory allocation failed in slicot_sb01bd.\n");
    }
    return info;
}