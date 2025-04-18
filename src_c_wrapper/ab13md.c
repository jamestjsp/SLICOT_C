/**
 * @file ab13md.c
 * @brief C wrapper implementation for SLICOT routine AB13MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB13MD,
 * which computes an upper bound on the structured singular value for a
 * square complex matrix Z with a given block uncertainty structure.
 */

#include <stdlib.h>
#include <ctype.h>
#include <stddef.h> // For size_t
#include <complex.h> // For creal

// Include the header file for this wrapper
#include "ab13md.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_complex_double
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * Note use of slicot_complex_double for Z and ZWORK.
 * X is input/output. BOUND, D, G are output.
 */
extern void F77_FUNC(ab13md, AB13MD)(
    const char* fact,       // CHARACTER*1 FACT
    const int* n,           // INTEGER N
    const slicot_complex_double* z, // COMPLEX*16 Z(LDZ,*)
    const int* ldz,         // INTEGER LDZ
    const int* m,           // INTEGER M
    const int* nblock,      // INTEGER NBLOCK(*)
    const int* itype,       // INTEGER ITYPE(*)
    double* x,              // DOUBLE PRECISION X(*) (in/out)
    double* bound,          // DOUBLE PRECISION BOUND (output)
    double* d,              // DOUBLE PRECISION D(*) (output)
    double* g,              // DOUBLE PRECISION G(*) (output)
    int* iwork,             // INTEGER IWORK(*)
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    slicot_complex_double* zwork, // COMPLEX*16 ZWORK(*)
    const int* lzwork,      // INTEGER LZWORK
    int* info,              // INTEGER INFO (output)
    int fact_len            // Hidden length
);


/* C wrapper function definition */
SLICOT_C_WRAPPER_API
int slicot_ab13md(char fact, int n, const slicot_complex_double* z, int ldz,
                  int m, const int* nblock, const int* itype,
                  double* x, double* bound, double* d, double* g,
                  int row_major)
{
    /* Local variables */
    int info = 0;
    int ldwork = -1; /* Use -1 for workspace query */
    int lzwork = -1; /* Use -1 for workspace query */
    double dwork_query;
    slicot_complex_double zwork_query;
    double* dwork = NULL;
    slicot_complex_double* zwork = NULL;
    int* iwork = NULL;
    int iwork_size = 0;
    const int fact_len = 1;

    char fact_upper = toupper(fact);

    /* Pointers for column-major copies if needed */
    slicot_complex_double *z_cm = NULL;

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -2; goto cleanup; }
    if (m < 1) { info = -5; goto cleanup; } // M must be >= 1
    if (fact_upper != 'F' && fact_upper != 'N') { info = -1; goto cleanup; }

    // Check leading dimensions based on storage order
    int min_ldz_f = MAX(1, n);

    if (row_major) {
        // For row-major C, LDZ is the number of columns
        int min_ldz_rm_cols = n;
        if (ldz < min_ldz_rm_cols) { info = -4; goto cleanup; }
    } else {
        // For column-major C, LDZ is the number of rows (Fortran style)
        if (ldz < min_ldz_f) { info = -4; goto cleanup; }
    }
    // Further validation (block sizes positive, sum=N, itype=1/2, real block size=1)
    // is handled by the Fortran routine (INFO = 1, 2, 3, 4).

    /* --- Workspace Allocation --- */

    // Allocate IWORK (size MAX(4*M-2, N))
    iwork_size = MAX(1, MAX(4 * m - 2, n));
    iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
    CHECK_ALLOC(iwork);

    // Query DWORK and ZWORK sizes
    ldwork = -1; // Query mode
    lzwork = -1; // Query mode
    F77_FUNC(ab13md, AB13MD)(&fact_upper, &n, z, &ldz, &m, nblock, itype,
                             x, bound, d, g, iwork,
                             &dwork_query, &ldwork,
                             &zwork_query, &lzwork, &info,
                             fact_len);

    if (info < 0) { goto cleanup; } // Query failed due to invalid argument
    info = 0; // Reset info after query

    // Get the required workspace sizes from query results
    ldwork = (int)dwork_query;
    lzwork = (int)SLICOT_COMPLEX_REAL(zwork_query); // Use macro to get real part

    // Check against minimum documented sizes
    int min_ldwork = 1;
    if (n > 0 && m > 0) { // Avoid potential negative sizes if n=0 or m=0
        min_ldwork = 2*n*n*m - n*n + 9*m*m + n*m + 11*n + 33*m - 11;
    }
    min_ldwork = MAX(1, min_ldwork);

    int min_lzwork = 1;
     if (n > 0 || m > 0) { // Avoid potential negative sizes if n=0 and m=0
        min_lzwork = 6*n*n*m + 12*n*n + 6*m + 6*n - 3;
    }
    min_lzwork = MAX(1, min_lzwork);

    ldwork = MAX(ldwork, min_ldwork);
    lzwork = MAX(lzwork, min_lzwork);

    // Allocate workspaces
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);
    zwork = (slicot_complex_double*)malloc((size_t)lzwork * sizeof(slicot_complex_double));
    CHECK_ALLOC(zwork);

    /* --- Prepare Arrays and Call Fortran Routine --- */
    size_t elem_size = sizeof(slicot_complex_double);

    if (row_major) {
        /* --- Row-Major Case --- */

        /* Allocate memory for column-major copy of Z */
        size_t z_rows = n; size_t z_cols = n; size_t z_size = z_rows * z_cols;
        if (z_size > 0) { z_cm = (slicot_complex_double*)malloc(z_size * elem_size); CHECK_ALLOC(z_cm); }

        /* Transpose C (row-major) input Z to Fortran (column-major) copy */
        if (z_size > 0) slicot_transpose_to_fortran(z, z_cm, z_rows, z_cols, elem_size);

        /* Fortran leading dimension */
        int ldz_f = (z_rows > 0) ? z_rows : 1;

        /* Call the Fortran routine */
        F77_FUNC(ab13md, AB13MD)(&fact_upper, &n,
                                 z_cm, &ldz_f,          // Pass CM Z
                                 &m, nblock, itype,
                                 x, bound, d, g, iwork,
                                 dwork, &ldwork, zwork, &lzwork, &info,
                                 fact_len);

        /* No copy-back needed for Z as it's input only */
        /* X, BOUND, D, G are modified directly */
        /* Temp array z_cm will be freed in cleanup */

    } else {
        /* --- Column-Major Case --- */

        /* Call the Fortran routine directly with user-provided arrays */
        F77_FUNC(ab13md, AB13MD)(&fact_upper, &n,
                                 z, &ldz,                // Pass original Z
                                 &m, nblock, itype,
                                 x, bound, d, g, iwork,
                                 dwork, &ldwork, zwork, &lzwork, &info,
                                 fact_len);
        // X, BOUND, D, G are modified directly.
    }

cleanup:
    /* --- Cleanup --- */
    free(zwork);
    free(dwork);
    free(iwork);
    free(z_cm); // Safe even if NULL

    return info;
}
