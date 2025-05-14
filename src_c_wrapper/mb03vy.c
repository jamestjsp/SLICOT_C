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
SLICOT_EXPORT
int slicot_mb03vy(int n, int p, int ilo, int ihi,
                  double* a, int lda1, int lda2,
                  const double* tau, int ldtau,
                  int row_major) // row_major ignored for 'a'
{
    /* Local variables */
    int info = 0;
    int ldwork_val; // Renamed to avoid conflict with ldwork parameter if it existed
    double* dwork = NULL;
    // No _cm pointers needed as 'a' is assumed column-major and row_major is ignored.
    // No iwork needed for this routine.

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -1; goto cleanup; } // N
    if (p < 1) { info = -2; goto cleanup; } // P
    // ILO (arg 3)
    if (n == 0) { // Special case for N=0
        if (ilo != 1) { info = -3; goto cleanup; }
    } else { // N > 0
        if (ilo < 1 || ilo > n) { info = -3; goto cleanup; }
    }
    // IHI (arg 4)
    if (n == 0) { // Special case for N=0
         if (ihi != 0) { info = -4; goto cleanup; }
    } else { // N > 0
        if (ihi < MIN(ilo, n) || ihi > n) { info = -4; goto cleanup; }
    }
    
    // A (arg 5) - NULL check
    if (n > 0 && p > 0 && a == NULL) { info = -5; goto cleanup; }

    // LDA1 (arg 6)
    if (lda1 < MAX(1, n)) { info = -6; goto cleanup; }
    // LDA2 (arg 7)
    if (lda2 < MAX(1, n)) { info = -7; goto cleanup; }

    // TAU (arg 8) - NULL check
    // TAU is dimensioned (LDTAU, P). It's required if P > 0.
    // N-1 part of LDTAU matters if N > 1 for actual reflector data.
    if (p > 0 && tau == NULL) { info = -8; goto cleanup; }
    
    // LDTAU (arg 9)
    if (ldtau < MAX(1, (n > 0 ? n - 1 : 0) )) { info = -9; goto cleanup; }
    // The 'row_major' flag is intentionally ignored for validation here due to the 3D array 'a'.

    /* --- Workspace Allocation --- */

    // Allocate a small DWORK array for the workspace query
    double dwork_temp_query[1];
    int ldwork_query_flag = -1; // Query mode
    
    // Perform workspace query
    // For query, pass potentially NULL 'a' and 'tau' if they were NULL and N=0 (already handled by validation if N>0)
    // Or if the routine is robust to NULL for query when dimensions are zero.
    // However, if N>0, 'a' and 'tau' have been validated to be non-NULL.
    double* a_for_query = (n > 0 && p > 0) ? a : NULL; // Use actual 'a' if valid, else NULL if N=0 or P=0 (though P>=1)
    const double* tau_for_query = (p > 0) ? tau : NULL; // Use actual 'tau' if P>=1
    int lda1_for_query = lda1;
    int lda2_for_query = lda2;
    int ldtau_for_query = ldtau;
    if (n == 0) { // If N=0, LDs must be >=1 for Fortran
        lda1_for_query = MAX(1, lda1_for_query);
        lda2_for_query = MAX(1, lda2_for_query);
        ldtau_for_query = MAX(1, ldtau_for_query);
    }


    F77_FUNC(mb03vy, MB03VY)(&n, &p, &ilo, &ihi, 
                             a_for_query, &lda1_for_query, &lda2_for_query, 
                             tau_for_query, &ldtau_for_query,
                             dwork_temp_query, &ldwork_query_flag, &info);

    if (info < 0 && info != SLICOT_MEMORY_ERROR) { // If query itself returns error (e.g. -1 to -9 due to other params)
        goto cleanup; 
    }
    info = 0; // Reset info if query was "successful" in terms of not erroring on other params
    
    // Get the optimal workspace size from the query result
    ldwork_val = (int)dwork_temp_query[0];
    
    // Check against minimum documented size: MAX(1, N)
    int min_ldwork_formula = MAX(1, n);
    ldwork_val = MAX(ldwork_val, min_ldwork_formula);

    // Allocate the workspace
    if (ldwork_val > 0) {
        dwork = (double*)malloc((size_t)ldwork_val * sizeof(double));
        CHECK_ALLOC(dwork); // Sets info to SLICOT_MEMORY_ERROR and jumps to cleanup on failure
    } else { // Should not happen as ldwork_val is MAX(..., 1)
        dwork = NULL; 
    }


    /* --- Call the computational routine --- */
    // 'a' and 'tau' are used directly as they are assumed column-major.
    // If N=0, pointers a and tau might be NULL but should not be dereferenced by Fortran.
    // If N>0, they've been checked to be non-NULL.
    F77_FUNC(mb03vy, MB03VY)(&n, &p, &ilo, &ihi,
                             a, &lda1, &lda2, // Pass original 'a'
                             tau, &ldtau,     // Pass original 'tau'
                             dwork, &ldwork_val, &info);
    // A is modified in place.

cleanup:
    /* --- Cleanup --- */
    free(dwork);

    return info;
}
