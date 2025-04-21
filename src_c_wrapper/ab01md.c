/**
 * @file ab01md.c
 * @brief C wrapper implementation for SLICOT routine AB01MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB01MD,
 * which finds a controllable realization for a linear time-invariant
 * single-input system using orthogonal transformations.
 */

#include <stdlib.h>
#include <string.h> // For memset
#include <ctype.h>  // For toupper

// Include the header file for this wrapper
#include "ab01md.h"
// Include necessary SLICOT utility headers
#include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines
#include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

/*
 * Declare the external Fortran routine using the F77_FUNC macro.
 * This handles potential name mangling issues between C and Fortran compilers.
 * All arguments are passed by reference (pointers).
 * The hidden string length argument for 'jobz' is explicitly included at the end.
 */
extern void F77_FUNC(ab01md, AB01MD)(
    const char* jobz,       // CHARACTER*1 JOBZ
    const int* n,           // INTEGER N
    double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
    const int* lda,         // INTEGER LDA
    double* b,              // DOUBLE PRECISION B(*) (in/out)
    int* ncont,             // INTEGER NCONT (output)
    double* z,              // DOUBLE PRECISION Z(LDZ,*) (output)
    const int* ldz,         // INTEGER LDZ
    double* tau,            // DOUBLE PRECISION TAU(*) (output)
    const double* tol,      // DOUBLE PRECISION TOL
    double* dwork,          // DOUBLE PRECISION DWORK(*)
    const int* ldwork,      // INTEGER LDWORK
    int* info,              // INTEGER INFO (output)
    int jobz_len            // Hidden length argument for jobz (integer)
);

/* Declaration of LAPACK DORGQR routine for forming orthogonal matrix from elementary reflectors */
extern void F77_FUNC(dorgqr, DORGQR)(
    const int* m,          // INTEGER M (number of rows)
    const int* n,          // INTEGER N (number of columns)
    const int* k,          // INTEGER K (number of elementary reflectors)
    double* a,             // DOUBLE PRECISION A(LDA,*) (in/out)
    const int* lda,        // INTEGER LDA
    const double* tau,     // DOUBLE PRECISION TAU(*) (input)
    double* work,          // DOUBLE PRECISION WORK(*) (workspace)
    const int* lwork,      // INTEGER LWORK
    int* info              // INTEGER INFO (output)
);

/* C wrapper function definition */
SLICOT_C_WRAPPER_API
int slicot_ab01md(char jobz, int n,
                  double* a, int lda,
                  double* b,
                  int* ncont,
                  double* z, int ldz,
                  double* tau, double tol, int row_major)
{
    /* Local variables */
    int info = 0;
    int ldwork = 0;
    double* dwork = NULL;
    int jobz_len = 1; // Fortran expects 1-based length for strings

    char jobz_upper = toupper(jobz);

    /* Pointers for column-major copies if needed */
    double *a_cm = NULL;
    double *z_cm = NULL; // Needed if row_major and jobz != 'N'

    /* --- Input Parameter Validation --- */

    if (n < 0) { info = -2; goto cleanup; }

    if (jobz_upper != 'N' && jobz_upper != 'F' && jobz_upper != 'I') {
        info = -1; goto cleanup;
    }

    // Check leading dimensions based on storage order and JOBZ
    int min_lda_f = MAX(1, n);
    int min_ldz_f = (jobz_upper == 'N') ? 1 : MAX(1, n);

    if (row_major) {
        // For row-major C, LDA/LDZ is the number of columns
        int min_lda_rm_cols = n;
        int min_ldz_rm_cols = (jobz_upper == 'N') ? 1 : n; // Need n columns for Z
        if (lda < min_lda_rm_cols) { info = -4; goto cleanup; }
        if (jobz_upper != 'N' && ldz < min_ldz_rm_cols) { info = -8; goto cleanup; }
    } else {
        // For column-major C, LDA/LDZ is the number of rows (Fortran style)
        if (lda < min_lda_f) { info = -4; goto cleanup; }
        if (ldz < min_ldz_f) { info = -8; goto cleanup; } // Check even if JOBZ='N' as per docs
    }

    /* --- Prepare arrays for column-major format if using row-major --- */
    int lda_f, ldz_f;
    double *a_ptr, *z_ptr;
    
    if (row_major) {
        /* Allocate and convert A from row-major to column-major */
        if (n > 0) {
            a_cm = (double*)malloc((size_t)n * n * sizeof(double));
            CHECK_ALLOC(a_cm);
            slicot_transpose_to_fortran(a, a_cm, n, n, sizeof(double));
        }
        
        /* Allocate Z if needed and initialize to identity if JOBZ='I' */
        if (jobz_upper != 'N' && n > 0) {
            z_cm = (double*)malloc((size_t)n * n * sizeof(double));
            CHECK_ALLOC(z_cm);
            
            if (jobz_upper == 'I') {
                set_identity(n, z_cm, n, 0); // 0 for column-major
            }
        }
        
        lda_f = n; // Column-major LDA for Fortran
        ldz_f = (jobz_upper == 'N') ? 1 : n; // Column-major LDZ for Fortran
        a_ptr = a_cm;
        z_ptr = (jobz_upper == 'N') ? NULL : z_cm;
    } else {
        /* Column-major case - use original arrays */
        lda_f = lda;
        ldz_f = ldz;
        a_ptr = a;
        z_ptr = (jobz_upper == 'N') ? NULL : z;
        
        /* Initialize Z to identity if JOBZ='I' in column-major format */
        if (jobz_upper == 'I' && n > 0) {
            set_identity(n, z, ldz, 0); // 0 for column-major
        }
    }

    /* --- Workspace allocation --- */
    /* According to AB01MD documentation, LDWORK >= MAX(1,N) is required.
     * For JOBZ='F', DORGQR will need additional workspace, we'll reuse this array 
     */
    ldwork = MAX(1, n*2);
    
    /* For optimal performance with DORGQR when JOBZ='F', allocate more space */
    if (jobz_upper == 'F' && n > 0) {
        /* Need extra space for DORGQR: optimal LDWORK >= N */
        ldwork = ldwork*2; // Double the size for better performance
    }
    
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);
    
    /* --- Call the computational routine --- */
    F77_FUNC(ab01md, AB01MD)(&jobz_upper, &n, 
                           a_ptr, &lda_f, 
                           b, ncont, 
                           z_ptr, &ldz_f, 
                           tau, &tol, 
                           dwork, &ldwork, &info, 
                           jobz_len);
    
    /* If JOBZ='F', we need to form the complete orthogonal matrix using DORGQR */
    if (info == 0 && jobz_upper == 'F' && n > 0) {
        /* Construct the orthogonal matrix from elementary reflectors */
        int dorgqr_info = 0;
        
        if (row_major) {
            /* Call DORGQR to form the orthogonal matrix in column-major format */
            F77_FUNC(dorgqr, DORGQR)(&n, &n, &n, z_cm, &ldz_f, tau, dwork, &ldwork, &dorgqr_info);
            if (dorgqr_info != 0) {
                info = dorgqr_info;
                goto cleanup;
            }
        } else {
            /* Call DORGQR to form the orthogonal matrix directly */
            F77_FUNC(dorgqr, DORGQR)(&n, &n, &n, z, &ldz, tau, dwork, &ldwork, &dorgqr_info);
            if (dorgqr_info != 0) {
                info = dorgqr_info;
                goto cleanup;
            }
        }
    }
    
    /* --- Copy results back to row-major format if needed --- */
    if (row_major && info == 0 && n > 0) {
        slicot_transpose_to_c(a_cm, a, n, n, sizeof(double));
        
        if (jobz_upper != 'N') {
            slicot_transpose_to_c(z_cm, z, n, n, sizeof(double));
        }
    }
    
cleanup:
    /* --- Cleanup --- */
    free(dwork);
    free(a_cm);
    free(z_cm);
    
    return info;
}
