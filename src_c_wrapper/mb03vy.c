/**
 * @file mb03vy.c
 * @brief C wrapper implementation for SLICOT routine MB03VY
 *
 * This file provides a C wrapper implementation for the SLICOT routine MB03VY,
 * which generates the orthogonal matrices Q_j from the elementary
 * reflectors computed by MB03VD.
 * NOTE: This wrapper assumes the 3D array 'a' is passed in column-major order
 * and the 'row_major' flag is IGNORED for 'a'.
 */

#include <stdlib.h>
#include <stddef.h> // For size_t

// Include the header file for this wrapper
#include "mb03vy.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, MIN, CHECK_ALLOC, SLICOT_MEMORY_ERROR
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * A is input/output. TAU is input.
 */
extern void F77_FUNC(mb03vy, MB03VY)(
    const int* n,           // INTEGER N
    const int* p,           // INTEGER P
    const int* ilo,         // INTEGER ILO
    const int* ihi,         // INTEGER IHI
    double* a,              // DOUBLE PRECISION A(LDA1,LDA2,P) (in/out)
    const int* lda1,        // INTEGER LDA1
    const int* lda2,        // INTEGER LDA2
    const double* tau,      // DOUBLE PRECISION TAU(LDTAU,P)
    const int* ldtau,       // INTEGER LDTAU
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info               // INTEGER INFO (output)
);


/* C wrapper function definition */
SLICOT_C_WRAPPER_API
int slicot_mb03vy(int n, int p, int ilo, int ihi,
                  double* a, int lda1, int lda2,
                  const double* tau, int ldtau,
                  int row_major) // row_major ignored for 'a'
{
    /* Local variables */
    int info = 0;
    int ldwork = -1; /* Use -1 for workspace query */
    double dwork_query;
    double* dwork = NULL;
    // No _cm pointers needed as 'a' is assumed column-major and row_major is ignored.
    // No iwork needed for this routine.

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -1; goto cleanup; }
    if (p < 1) { info = -2; goto cleanup; }
    if (ilo < 1 || ilo > MAX(1, n)) { info = -3; goto cleanup; }
    if (ihi < MIN(ilo, n) || ihi > n) { info = -4; goto cleanup; }
    // NOTE: LDA1, LDA2, LDTAU checks assume 'a' and 'tau' are column-major.
    if (lda1 < MAX(1, n)) { info = -6; goto cleanup; }
    if (lda2 < MAX(1, n)) { info = -7; goto cleanup; }
    if (ldtau < MAX(1, n - 1)) { info = -9; goto cleanup; }
    // The 'row_major' flag is intentionally ignored for validation here due to the 3D array 'a'.

    /* --- Workspace Allocation --- */

    // Perform workspace query for DWORK
    ldwork = -1; // Query mode
    F77_FUNC(mb03vy, MB03VY)(&n, &p, &ilo, &ihi, a, &lda1, &lda2, tau, &ldtau,
                             &dwork_query, &ldwork, &info);

    if (info < 0 && info != -11) { goto cleanup; } // Query failed due to invalid argument (allow INFO=-11 from query)
    info = 0; // Reset info after query

    // Get the required dwork size from query result
    ldwork = (int)dwork_query;
    // Check against minimum documented size: MAX(1, N)
    int min_ldwork = MAX(1, n);
    ldwork = MAX(ldwork, min_ldwork);

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

    /* --- Call the computational routine --- */

    // NOTE: Assuming 'a' is already in column-major (Fortran) layout.
    // No transposition is performed for the 3D array 'a'.
    // 'row_major' flag is ignored. 'tau' is input-only, no copy-back needed.
    F77_FUNC(mb03vy, MB03VY)(&n, &p, &ilo, &ihi,
                             a, &lda1, &lda2,
                             tau, &ldtau,
                             dwork, &ldwork, &info);
    // A is modified in place.

    /* --- Copy results back to row-major format if needed --- */
    // No copy-back needed as 'a' is assumed column-major and modified in-place.

cleanup:
    /* --- Cleanup --- */
    free(dwork);

    /* Return the info code from the Fortran routine or SLICOT_MEMORY_ERROR */
    return info;
}