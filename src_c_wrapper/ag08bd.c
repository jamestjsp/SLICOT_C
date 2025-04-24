/**
 * @file ag08bd.c
 * @brief C wrapper implementation for SLICOT routine AG08BD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AG08BD,
 * which computes the zeros and Kronecker structure of a descriptor
 * system pencil S(lambda) = [A-lambda*E, B; C, D].
 * Refactored to align with ab01nd.c structure.
 */

#include <stdlib.h>
#include <ctype.h>  // For toupper
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "ag08bd.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * Note A, E, B, C are input/output. D is input.
 * Many integer outputs.
 */
extern void F77_FUNC(ag08bd, AG08BD)(
    const char* equil,      // CHARACTER*1 EQUIL
    const int* l,           // INTEGER L
    const int* n,           // INTEGER N
    const int* m,           // INTEGER M
    const int* p,           // INTEGER P
    double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
    const int* lda,         // INTEGER LDA
    double* e,              // DOUBLE PRECISION E(LDE,*) (in/out)
    const int* lde,         // INTEGER LDE
    double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
    const int* ldb,         // INTEGER LDB
    double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
    const int* ldc,         // INTEGER LDC
    const double* d,        // DOUBLE PRECISION D(LDD,*)
    const int* ldd,         // INTEGER LDD
    int* nfz,               // INTEGER NFZ (output)
    int* nrank,             // INTEGER NRANK (output)
    int* niz,               // INTEGER NIZ (output)
    int* dinfz,             // INTEGER DINFZ (output)
    int* nkror,             // INTEGER NKROR (output)
    int* ninfe,             // INTEGER NINFE (output)
    int* nkrol,             // INTEGER NKROL (output)
    int* infz,              // INTEGER INFZ(*) (output)
    int* kronr,             // INTEGER KRONR(*) (output)
    int* infe,              // INTEGER INFE(*) (output)
    int* kronl,             // INTEGER KRONL(*) (output)
    const double* tol,      // DOUBLE PRECISION TOL
    int* iwork,             // INTEGER IWORK(*)
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info,              // INTEGER INFO (output)
    int equil_len           // Hidden length
);


/* C wrapper function definition */
SLICOT_EXPORT
int slicot_ag08bd(char equil, int l, int n, int m, int p,
                  double* a, int lda, double* e, int lde,
                  double* b, int ldb, double* c, int ldc,
                  const double* d, int ldd,
                  int* nfz, int* nrank, int* niz, int* dinfz,
                  int* nkror, int* ninfe, int* nkrol,
                  int* infz, int* kronr, int* infe, int* kronl,
                  double tol, int row_major)
{
    /* Local variables */
    int info = 0;
    int ldwork = -1; /* Use -1 for workspace query */
    double dwork_query;
    double* dwork = NULL;
    int* iwork = NULL;
    int iwork_size = 0;
    const int equil_len = 1;

    char equil_upper = toupper(equil);

    /* Pointers for column-major copies if needed */
    double *a_cm = NULL, *e_cm = NULL, *b_cm = NULL, *c_cm = NULL;
    double *d_cm = NULL; // D is const input, but need temp if row_major

    /* Pointers to pass to Fortran */
    double *a_ptr, *e_ptr, *b_ptr, *c_ptr;
    const double *d_ptr;
    int lda_f, lde_f, ldb_f, ldc_f, ldd_f;


    /* --- Input Parameter Validation --- */

    if (l < 0) { info = -2; goto cleanup; }
    if (n < 0) { info = -3; goto cleanup; }
    if (m < 0) { info = -4; goto cleanup; }
    if (p < 0) { info = -5; goto cleanup; }
    if (equil_upper != 'S' && equil_upper != 'N') { info = -1; goto cleanup; }
    if (tol >= 1.0 || tol < 0.0) { info = -26; goto cleanup; } // TOL must be [0, 1)

    // Check leading dimensions based on storage order
    int min_lda_f = MAX(1, l); int min_lde_f = MAX(1, l);
    int min_ldb_f = (m > 0) ? MAX(1, l) : 1;
    int min_ldc_f = MAX(1, p); int min_ldd_f = MAX(1, p);

    if (row_major) {
        // For row-major C, LDA is the number of columns
        int min_lda_rm_cols = n; int min_lde_rm_cols = n;
        int min_ldb_rm_cols = m;
        int min_ldc_rm_cols = n; int min_ldd_rm_cols = m;
        if (l > 0 && lda < min_lda_rm_cols) { info = -7; goto cleanup; }
        if (l > 0 && lde < min_lde_rm_cols) { info = -9; goto cleanup; }
        if (l > 0 && ldb < min_ldb_rm_cols) { info = -11; goto cleanup; }
        if (p > 0 && ldc < min_ldc_rm_cols) { info = -13; goto cleanup; }
        if (p > 0 && ldd < min_ldd_rm_cols) { info = -15; goto cleanup; }
    } else {
        // For column-major C, LDA is the number of rows (Fortran style)
        if (lda < min_lda_f) { info = -7; goto cleanup; }
        if (lde < min_lde_f) { info = -9; goto cleanup; }
        if (ldb < min_ldb_f) { info = -11; goto cleanup; }
        if (ldc < min_ldc_f) { info = -13; goto cleanup; }
        if (ldd < min_ldd_f) { info = -15; goto cleanup; }
    }

    /* --- Workspace Allocation --- */

    // Allocate IWORK (size N + max(1,M))
    iwork_size = MAX(1, n + MAX(1, m));
    iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
    CHECK_ALLOC(iwork);

    // Allocate DWORK based on query
    ldwork = -1; // Query mode
    // Use dummy LDs for query if dimensions are 0
    int lda_q = row_major ? MAX(1, l) : lda;
    int lde_q = row_major ? MAX(1, l) : lde;
    int ldb_q = row_major ? MAX(1, l) : ldb;
    int ldc_q = row_major ? MAX(1, p) : ldc;
    int ldd_q = row_major ? MAX(1, p) : ldd;

    F77_FUNC(ag08bd, AG08BD)(&equil_upper, &l, &n, &m, &p,
                             NULL, &lda_q, NULL, &lde_q, NULL, &ldb_q, NULL, &ldc_q, NULL, &ldd_q, // NULL arrays
                             nfz, nrank, niz, dinfz, nkror, ninfe, nkrol,
                             NULL, NULL, NULL, NULL, // NULL integer arrays
                             &tol, iwork, &dwork_query, &ldwork, &info,
                             equil_len);

    if (info < 0) { goto cleanup; } // Query failed due to invalid argument
    info = 0; // Reset info after query

    // Get the required dwork size from query result
    ldwork = (int)dwork_query;
    // Check against minimum documented size
    int ldw = MAX(1, 5 * MAX(l + p, m + n));
    // Simplified check, full formula complex: MAX(ldw, (L+P+M+N)*(M+N))
    ldw = MAX(ldw, (l + p + m + n) * (m + n));
    int min_ldwork = ldw;
    if (equil_upper == 'S') {
        min_ldwork = MAX(4 * (l + n), ldw);
    }
    ldwork = MAX(ldwork, min_ldwork);

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

    /* --- Prepare Arrays and Call Fortran Routine --- */
    size_t elem_size = sizeof(double);

    if (row_major) {
        /* --- Row-Major Case --- */

        /* Allocate memory for column-major copies */
        size_t a_rows = l; size_t a_cols = n; size_t a_size = a_rows * a_cols;
        size_t e_rows = l; size_t e_cols = n; size_t e_size = e_rows * e_cols;
        size_t b_rows = l; size_t b_cols = m; size_t b_size = b_rows * b_cols;
        size_t c_rows = p; size_t c_cols = n; size_t c_size = c_rows * c_cols;
        size_t d_rows = p; size_t d_cols = m; size_t d_size = d_rows * d_cols;

        if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
        if (e_size > 0) { e_cm = (double*)malloc(e_size * elem_size); CHECK_ALLOC(e_cm); }
        if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
        if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
        if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); } // D is input only

        /* Transpose C (row-major) inputs to Fortran (column-major) copies */
        if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
        if (e_cm) slicot_transpose_to_fortran(e, e_cm, e_rows, e_cols, elem_size);
        if (b_cm) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
        if (c_cm) slicot_transpose_to_fortran(c, c_cm, c_rows, c_cols, elem_size);
        if (d_cm) slicot_transpose_to_fortran(d, d_cm, d_rows, d_cols, elem_size);

        /* Fortran leading dimensions */
        lda_f = MAX(1, a_rows);
        lde_f = MAX(1, e_rows);
        ldb_f = MAX(1, b_rows);
        ldc_f = MAX(1, c_rows);
        ldd_f = MAX(1, d_rows);

        /* Set pointers for Fortran call */
        a_ptr = a_cm; e_ptr = e_cm; b_ptr = b_cm; c_ptr = c_cm;
        d_ptr = d_cm; // Pass pointer to temp CM copy of D

    } else {
        /* --- Column-Major Case --- */
        lda_f = lda; lde_f = lde; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
        a_ptr = a; e_ptr = e; b_ptr = b; c_ptr = c;
        d_ptr = d; // Pass original pointer
    }

    /* Call the computational routine */
    F77_FUNC(ag08bd, AG08BD)(&equil_upper, &l, &n, &m, &p,
                             a_ptr, &lda_f, e_ptr, &lde_f, // Pass A, E ptrs
                             b_ptr, &ldb_f, c_ptr, &ldc_f, // Pass B, C ptrs
                             d_ptr, &ldd_f,              // Pass D ptr
                             nfz, nrank, niz, dinfz, nkror, ninfe, nkrol,
                             infz, kronr, infe, kronl,
                             &tol, iwork, dwork, &ldwork, &info,
                             equil_len);

    /* Copy back results from column-major temps to original row-major arrays */
    if (row_major && info == 0) {
        int nfz_val = *nfz;
        // Copy back modified A and E (reduced pencil Af, Ef)
        if (nfz_val >= 0) { // Allow NFZ=0
            size_t a_rows = l; size_t a_cols = n; size_t a_size = a_rows * a_cols;
            size_t e_rows = l; size_t e_cols = n; size_t e_size = e_rows * e_cols;
            if (a_cm && a_size > 0) slicot_transpose_to_c(a_cm, a, nfz_val, nfz_val, elem_size);
            if (e_cm && e_size > 0) slicot_transpose_to_c(e_cm, e, nfz_val, nfz_val, elem_size);
        }
        // B and C are overwritten with useless info, no need to copy back.
        // D is input only, no copy back needed.
        // Integer output arrays are filled directly.
    }

cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(iwork);
    free(a_cm);
    free(e_cm);
    free(b_cm);
    free(c_cm);
    free(d_cm); // Safe even if NULL

    return info;
}
