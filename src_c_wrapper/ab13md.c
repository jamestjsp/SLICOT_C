/**
 * @file ab13md.c
 * @brief C wrapper implementation for SLICOT routine AB13MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB13MD,
 * which computes an upper bound on the structured singular value for a
 * square complex matrix Z with a given block uncertainty structure.
 * This version skips workspace queries and uses slightly more generous,
 * formula-based workspace sizes derived from documentation.
 */

#include <stdlib.h>
#include <ctype.h>  // For toupper
#include <stddef.h> // For size_t
#include <string.h> // For memset

// Include the header file for this wrapper
#include "ab13md.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" 
#include "slicot_f77.h"   

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
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
SLICOT_EXPORT
int slicot_ab13md(char fact, int n, const slicot_complex_double* z, int ldz,
                  int m, const int* nblock, const int* itype,
                  double* x, double* bound, double* d, double* g,
                  int row_major)
{
    /* Local variables */
    int info = 0;
    int ldwork_val = 1; 
    int lzwork_val = 1; 
    double* dwork = NULL;
    slicot_complex_double* zwork = NULL;
    int* iwork = NULL;
    int iwork_size = 0;
    const int fact_len = 1;

    char fact_upper = toupper(fact);

    /* Pointers for column-major copies if needed */
    slicot_complex_double *z_cm = NULL;

    /* Pointers to pass to Fortran */
    const slicot_complex_double *z_ptr;
    int ldz_f;

    /* --- Input Parameter Validation (Essential C-level checks) --- */

    if (fact_upper != 'F' && fact_upper != 'N') { info = -1; goto cleanup; } 
    if (n < 0) { info = -2; goto cleanup; } 
    if (n > 0 && z == NULL) { info = -3; goto cleanup; }
    
    int min_ldz_f_val = MAX(1, n);
    if (row_major) {
        if (n > 0 && ldz < n) { info = -4; goto cleanup; } 
    } else {
        if (n > 0 && ldz < min_ldz_f_val) { info = -4; goto cleanup; } 
    }

    if (m < 1) { info = -5; goto cleanup; } 
    if (m > 0 && nblock == NULL) { info = -6; goto cleanup; }
    if (m > 0 && itype == NULL) { info = -7; goto cleanup; }
    
    if (m > 0 && nblock != NULL && itype != NULL) { 
        int first_block_size = nblock[0]; 
        if (fact_upper == 'F' || (n > 0 && first_block_size != n)) { 
            int num_real_blocks = 0;
            for (int i_mr = 0; i_mr < m; ++i_mr) {
                if (itype[i_mr] == 1) num_real_blocks++;
            }
            int x_dim_calculated = m + num_real_blocks - 1;
            if (x_dim_calculated > 0 && x == NULL) {
                info = -8; goto cleanup;
            }
        }
    } else if (m > 0 && (nblock == NULL || itype == NULL)) {
        if (fact_upper == 'F' && x == NULL) { info = -8; goto cleanup; }
    }

    if (bound == NULL) { info = -9; goto cleanup; } 
    if (n > 0 && d == NULL) { info = -10; goto cleanup; }
    if (n > 0 && g == NULL) { info = -11; goto cleanup; }

    /* --- Workspace Allocation --- */

    // IWORK: size MAX(1, MAX(4*M-2, N))
    iwork_size = MAX(1, MAX(4 * m - 2, n)); 
    iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
    if (iwork == NULL && iwork_size > 0) { info = SLICOT_MEMORY_ERROR; goto cleanup;} // CHECK_ALLOC behavior
    // Initialize iwork
    if (iwork) memset(iwork, 0, (size_t)iwork_size * sizeof(int));


    // DWORK and ZWORK sizes determined by formula (no query, slightly more generous)
    if (n == 0) { // M must be >= 1
        ldwork_val = 9*m*m + 33*m - 11;
        lzwork_val = 6*m - 3;
    } else { // N > 0 and M >= 1
        // Minimum documented formula
        long long ldwork_min_formula = (long long)2*n*n*m - (long long)n*n + (long long)9*m*m + 
                                       (long long)n*m + (long long)11*n + (long long)33*m - 11;
        // Add simpler part of "best performance" recommendation
        ldwork_val = (int)(ldwork_min_formula + (long long)5*n);

        long long lzwork_min_formula = (long long)6*n*n*m + (long long)12*n*n + (long long)6*m + 
                                       (long long)6*n - 3;
        // Add simpler part of "best performance" recommendation
        lzwork_val = (int)(lzwork_min_formula + (long long)3*n);
    }
    // Ensure at least 1 for all workspace arrays
    ldwork_val = MAX(1, ldwork_val); 
    lzwork_val = MAX(1, lzwork_val);

    dwork = (double*)malloc((size_t)ldwork_val * sizeof(double));
    if (dwork == NULL && ldwork_val > 0) { info = SLICOT_MEMORY_ERROR; goto cleanup;} // CHECK_ALLOC behavior
    
    zwork = (slicot_complex_double*)malloc((size_t)lzwork_val * sizeof(slicot_complex_double));
    if (zwork == NULL && lzwork_val > 0) { info = SLICOT_MEMORY_ERROR; goto cleanup;} // CHECK_ALLOC behavior
    
    // Initialize dwork and zwork
    if (dwork) memset(dwork, 0, (size_t)ldwork_val * sizeof(double));
    if (zwork) memset(zwork, 0, (size_t)lzwork_val * sizeof(slicot_complex_double));


    /* --- Prepare Z Array for Fortran --- */
    size_t z_elem_size = sizeof(slicot_complex_double);
    if (row_major && n > 0) {
        ldz_f = MAX(1, n); 
        size_t z_matrix_elements = (size_t)n * n;
        z_cm = (slicot_complex_double*)malloc(z_matrix_elements * z_elem_size);
        if (z_cm == NULL && z_matrix_elements > 0) { info = SLICOT_MEMORY_ERROR; goto cleanup;} // CHECK_ALLOC behavior
        if (z_cm) { // Ensure z_cm is not NULL before transposing
           slicot_transpose_to_fortran_with_ld(z, z_cm, n, n, ldz, ldz_f, z_elem_size);
        }
        z_ptr = z_cm;
    } else { 
        ldz_f = (n > 0) ? ldz : 1; 
        z_ptr = (n > 0 && z != NULL) ? z : NULL; 
    }
    if (n == 0) ldz_f = 1; 

    /* Call the computational Fortran routine */
    F77_FUNC(ab13md, AB13MD)(&fact_upper, &n,
                             z_ptr, &ldz_f,
                             &m, nblock, itype,
                             x, bound, d, g, iwork,
                             dwork, &ldwork_val, zwork, &lzwork_val, &info,
                             fact_len);

cleanup:
    /* --- Cleanup --- */
    free(zwork);
    free(dwork);
    free(iwork);
    if (z_cm) free(z_cm); 

    return info;
}
