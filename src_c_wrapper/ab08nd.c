/**
 * @file ab08nd.c
 * @brief C wrapper implementation for SLICOT routine AB08ND
 *
 * This file provides a C implementation of the SLICOT routine AB08ND wrapper,
 * which constructs a regular pencil for a given system and computes its
 * invariant zeros and Kronecker indices.
 * Refactored to align with ab01nd.c structure.
 */

#include <stdlib.h>
#include <math.h>   // For isnan, isinf (optional checks)
#include <ctype.h>  // For toupper
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "ab08nd.h"
// Include necessary SLICOT utility headers
// Ensure slicot_utils.h provides: MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR,
// slicot_transpose_to_fortran, slicot_transpose_to_c, slicot_transpose_to_c_with_ld
#include "slicot_utils.h"
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external FORTRAN routine using the F77_FUNC macro.
 * Explicitly includes the hidden 'int' argument for the length of 'equil'.
 */
extern void F77_FUNC(ab08nd, AB08ND)(
    const char* equil, const int* n, const int* m, const int* p,
    double* a, const int* lda, double* b, const int* ldb,
    double* c, const int* ldc, double* d, const int* ldd,
    int* nu, int* rank, int* dinfz, int* nkror, int* nkrol,
    int* infz, int* kronr, int* kronl,
    double* af, const int* ldaf, double* bf, const int* ldbf,
    const double* tol, int* iwork, double* dwork, const int* ldwork,
    int* info,
    int equil_len /* Hidden length argument for equil */);

/* C wrapper for AB08ND */
SLICOT_EXPORT
int slicot_ab08nd(char equil, int n, int m, int p,
                  double* a, int lda,
                  double* b, int ldb,
                  double* c, int ldc,
                  double* d, int ldd,
                  int* nu, int* rank,
                  int* dinfz, int* nkror, int* nkrol,
                  int* infz, int* kronr, int* kronl,
                  double* af, int ldaf,
                  double* bf, int ldbf,
                  double tol, int row_major)
{
    /* Local variables */
    int info = 0;
    int ldwork = -1; /* Use -1 for workspace query */
    double *dwork = NULL;
    int *iwork = NULL;
    int iwork_size = 0;
    const int equil_len = 1; // Fortran expects 1-based length for strings

    char equil_upper = toupper(equil);

    /* Pointers for column-major copies if needed */
    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
    double *af_cm = NULL, *bf_cm = NULL;

    /* Pointers to pass to Fortran */
    double *a_ptr, *b_ptr, *c_ptr, *d_ptr;
    double *af_ptr, *bf_ptr;
    int lda_f, ldb_f, ldc_f, ldd_f;
    int ldaf_f, ldbf_f;

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -2; goto cleanup; }
    if (m < 0) { info = -3; goto cleanup; }
    if (p < 0) { info = -4; goto cleanup; }
    if (equil_upper != 'S' && equil_upper != 'N') { info = -1; goto cleanup; }
    // Optional: Check TOL range if necessary
    // if (tol < 0.0) { info = -23; goto cleanup; }

    // Determine minimum required Fortran leading dimensions
    int min_lda_f = MAX(1, n);
    int min_ldb_f = MAX(1, n);
    int min_ldc_f = MAX(1, p);
    int min_ldd_f = MAX(1, p);
    int min_ldaf_f = MAX(1, n + m); // Fortran routine requires LDAF >= N+M (rows)
    int min_ldbf_f = MAX(1, n + p); // Fortran routine requires LDBF >= N+P (rows)

    if (row_major) {
        // For row-major C, LDA/LDB/LDC/LDD are the number of columns
        // LDAF/LDBF are also number of columns for the C arrays af/bf
        int min_lda_rm_cols = n;
        int min_ldb_rm_cols = m;
        int min_ldc_rm_cols = n;
        int min_ldd_rm_cols = m;
        // Check C array dimensions (number of columns)
        if (n > 0 && lda < min_lda_rm_cols) { info = -6; goto cleanup; }
        if (n > 0 && ldb < min_ldb_rm_cols) { info = -8; goto cleanup; }
        if (p > 0 && ldc < min_ldc_rm_cols) { info = -10; goto cleanup; }
        if (p > 0 && ldd < min_ldd_rm_cols) { info = -12; goto cleanup; }
        // Check if C ldaf/ldbf are sufficient to act as Fortran LDAF/LDBF (rows)
        if (ldaf < min_ldaf_f) { info = -20; goto cleanup; } // C ldaf (cols) must be >= Fortran LDAF (rows) requirement
        if (ldbf < min_ldbf_f) { info = -22; goto cleanup; } // C ldbf (cols) must be >= Fortran LDBF (rows) requirement
    } else {
        // For column-major C, LDA/LDB/LDC/LDD/LDAF/LDBF are the number of rows (Fortran style)
        if (lda < min_lda_f) { info = -6; goto cleanup; }
        if (ldb < min_ldb_f) { info = -8; goto cleanup; }
        if (ldc < min_ldc_f) { info = -10; goto cleanup; }
        if (ldd < min_ldd_f) { info = -12; goto cleanup; }
        if (ldaf < min_ldaf_f) { info = -20; goto cleanup; }
        if (ldbf < min_ldbf_f) { info = -22; goto cleanup; }
    }

    /* --- Workspace Allocation --- */

    // Determine iwork size from AB08ND documentation: MAX(M,P)
    iwork_size = MAX(1, MAX(m, p)); // Ensure minimum size 1
    iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
    CHECK_ALLOC(iwork);

    // Call Fortran routine for DWORK workspace query
    ldwork = -1; // Query mode
    double dwork_query[1]; /* Use array instead of scalar for workspace query */
    
    // Use dummy LDs for query if dimensions are 0
    int lda_q = row_major ? MAX(1, n) : lda;
    int ldb_q = row_major ? MAX(1, n) : ldb;
    int ldc_q = row_major ? MAX(1, p) : ldc;
    int ldd_q = row_major ? MAX(1, p) : ldd;
    int ldaf_q = row_major ? ldaf : ldaf; // Use C ldaf/ldbf for Fortran LDAF/LDBF
    int ldbf_q = row_major ? ldbf : ldbf;

    F77_FUNC(ab08nd, AB08ND)(&equil_upper, &n, &m, &p,
                             NULL, &lda_q, NULL, &ldb_q, // NULL arrays for query
                             NULL, &ldc_q, NULL, &ldd_q,
                             nu, rank, dinfz, nkror, nkrol, // Pass original pointers
                             NULL, NULL, NULL, // NULL integer arrays
                             NULL, &ldaf_q, NULL, &ldbf_q, // NULL output arrays
                             &tol, iwork, dwork_query, &ldwork, &info, // ldwork = -1
                             equil_len /* Explicit length for equil */);

    if (info != 0) {
        // Query failed, likely due to invalid N, M, P, or LDs passed to query
        goto cleanup;
    }

    // Get the required dwork size from query result
    ldwork = (int)dwork_query[0];
    // Ensure minimum size based on documentation formula
    int min_ldwork_doc = 1;
    if (n > 0 || m > 0 || p > 0) {
        min_ldwork_doc = MAX(1, MIN(p, m) + MAX(3 * m - 1, n));
        min_ldwork_doc = MAX(min_ldwork_doc, MIN(p, n) + MAX(3 * p - 1, MAX(n + p, n + m)));
        min_ldwork_doc = MAX(min_ldwork_doc, MIN(m, n) + MAX(3 * m - 1, n + m));
    }
    ldwork = MAX(ldwork, min_ldwork_doc);

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork); // Use CHECK_ALLOC macro

    /* --- Prepare Arrays and Call Fortran Routine --- */
    size_t elem_size = sizeof(double);

    if (row_major) {
        /* --- Row-Major Case --- */

        /* Allocate memory for column-major copies */
        // Dimensions for temporary Fortran-style matrices
        size_t a_rows_f = n; size_t a_cols_f = n; size_t a_size = a_rows_f * a_cols_f;
        size_t b_rows_f = n; size_t b_cols_f = m; size_t b_size = b_rows_f * b_cols_f;
        size_t c_rows_f = p; size_t c_cols_f = n; size_t c_size = c_rows_f * c_cols_f;
        size_t d_rows_f = p; size_t d_cols_f = m; size_t d_size = d_rows_f * d_cols_f;
        // Fortran AF is LDAF_f x (N+MIN(P,M)), BF is LDBF_f x (N+M)
        // Use C ldaf/ldbf as Fortran LDAF_f/LDBF_f (rows)
        ldaf_f = ldaf; // Fortran LDAF (rows) = C ldaf (cols)
        ldbf_f = ldbf; // Fortran LDBF (rows) = C ldbf (cols)
        size_t af_rows_f = ldaf_f; size_t af_cols_f = n + MIN(p, m); size_t af_size = af_rows_f * af_cols_f;
        size_t bf_rows_f = ldbf_f; size_t bf_cols_f = n + m; size_t bf_size = bf_rows_f * bf_cols_f;

        // Allocate only if size > 0
        if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
        if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
        if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
        if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
        if (af_size > 0) { af_cm = (double*)malloc(af_size * elem_size); CHECK_ALLOC(af_cm); }
        if (bf_size > 0) { bf_cm = (double*)malloc(bf_size * elem_size); CHECK_ALLOC(bf_cm); }

        /* Convert row-major inputs to column-major */
        // Pass C dimensions (rows, cols) to transpose function
        if (a_size > 0) slicot_transpose_to_fortran(a, a_cm, n, n, elem_size); // C A is n x n (lda cols)
        if (b_size > 0) slicot_transpose_to_fortran(b, b_cm, n, m, elem_size); // C B is n x m (ldb cols)
        if (c_size > 0) slicot_transpose_to_fortran(c, c_cm, p, n, elem_size); // C C is p x n (ldc cols)
        if (d_size > 0) slicot_transpose_to_fortran(d, d_cm, p, m, elem_size); // C D is p x m (ldd cols)

        /* Fortran leading dimensions for the call (rows) */
        lda_f = MAX(1, a_rows_f);
        ldb_f = MAX(1, b_rows_f);
        ldc_f = MAX(1, c_rows_f);
        ldd_f = MAX(1, d_rows_f);
        // ldaf_f and ldbf_f already set above

        /* Set pointers for Fortran call */
        a_ptr = a_cm; b_ptr = b_cm; c_ptr = c_cm; d_ptr = d_cm;
        af_ptr = af_cm; bf_ptr = bf_cm;

    } else {
        /* --- Column-Major Case --- */
        // Fortran LDs are the same as C LDs
        lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
        ldaf_f = ldaf; ldbf_f = ldbf;
        // Pointers point directly to user data
        a_ptr = a; b_ptr = b; c_ptr = c; d_ptr = d;
        af_ptr = af; bf_ptr = bf;
    }

    /* Call the computational routine */
    F77_FUNC(ab08nd, AB08ND)(&equil_upper, &n, &m, &p,
                             a_ptr, &lda_f, b_ptr, &ldb_f,
                             c_ptr, &ldc_f, d_ptr, &ldd_f,
                             nu, rank, dinfz, nkror, nkrol,
                             infz, kronr, kronl,
                             af_ptr, &ldaf_f, bf_ptr, &ldbf_f,
                             &tol, iwork, dwork, &ldwork, &info,
                             equil_len);

    /* Copy back results if row_major */
    if (row_major && info == 0) {
         int nu_val = *nu;
         // Check if C caller allocated enough columns for the nu x nu result
         // C ldaf/ldbf are number of columns
         if (ldaf < nu_val) { info = -20; goto cleanup; } // Invalid C LDAF (cols) for computed NU
         if (ldbf < nu_val) { info = -22; goto cleanup; } // Invalid C LDBF (cols) for computed NU

         // Copy back AF and BF (only the NU x NU part)
         if (nu_val > 0) {
             // We need to correctly handle the leading dimensions during transpose
             // For AF: source is column-major with ldaf_f rows, destination is row-major with ldaf columns
             if (af_cm) {
                 // Use our new helper function that handles different leading dimensions
                 slicot_transpose_to_c_with_ld(af_cm, af, nu_val, nu_val, ldaf_f, ldaf, elem_size);
             }
             
             // For BF: source is column-major with ldbf_f rows, destination is row-major with ldbf columns
             if (bf_cm) {
                 slicot_transpose_to_c_with_ld(bf_cm, bf, nu_val, nu_val, ldbf_f, ldbf, elem_size);
             }
         }
         
         /* Update original input matrices if changed by the routine (EQUIL='S') */
         if (equil_upper == 'S') {
             // Dimensions of the original C matrices
             size_t a_rows_c = n; size_t a_cols_c = n; size_t a_size = a_rows_c * a_cols_c;
             size_t b_rows_c = n; size_t b_cols_c = m; size_t b_size = b_rows_c * b_cols_c;
             size_t c_rows_c = p; size_t c_cols_c = n; size_t c_size = c_rows_c * c_cols_c;
             size_t d_rows_c = p; size_t d_cols_c = m; size_t d_size = d_rows_c * d_cols_c;

             // Transpose back from temporary column-major storage
             if (a_cm && a_size > 0) slicot_transpose_to_c(a_cm, a, a_rows_c, a_cols_c, elem_size);
             if (b_cm && b_size > 0) slicot_transpose_to_c(b_cm, b, b_rows_c, b_cols_c, elem_size);
             if (c_cm && c_size > 0) slicot_transpose_to_c(c_cm, c, c_rows_c, c_cols_c, elem_size);
             if (d_cm && d_size > 0) slicot_transpose_to_c(d_cm, d, d_rows_c, d_cols_c, elem_size);
         }
    }

cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(iwork);
    // Free temporary column-major arrays if allocated
    free(a_cm);
    free(b_cm);
    free(c_cm);
    free(d_cm);
    free(af_cm);
    free(bf_cm);

    return info;
}
