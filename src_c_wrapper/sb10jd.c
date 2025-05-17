/**
 * @file sb10jd.c
 * @brief C wrapper for SLICOT routine SB10JD.
 * @details Converts a descriptor state-space system into regular state-space form.
 * Matrices A, B, C, D, E are modified by this wrapper to store results or
 * intermediate values as per SLICOT routine behavior.
 * Workspace (DWORK) is allocated internally.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h> // For malloc, free
#include <string.h> // For memcpy
#include <math.h>   // For fmax (from C standard library)
#include "sb10jd.h" // Public header for this wrapper
#include "slicot_utils.h"  // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions etc.
#include "slicot_f77.h"    // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */
extern void F77_FUNC(sb10jd, SB10JD)(
    const int* n_fortran, const int* m_fortran, const int* np_fortran,
    double* a_fortran, const int* lda_fortran,
    double* b_fortran, const int* ldb_fortran,
    double* c_fortran, const int* ldc_fortran,
    double* d_fortran, const int* ldd_fortran,
    double* e_fortran, const int* lde_fortran,
    int* nsys_fortran,
    double* dwork_fortran, const int* ldwork_fortran,
    int* info_fortran
);

/* C wrapper function definition */
SLICOT_EXPORT
int slicot_sb10jd(int n_param, int m_param, int np_param,
                  double* a, int lda,
                  double* b, int ldb,
                  double* c, int ldc,
                  double* d, int ldd,
                  double* e, int lde,
                  int* nsys,
                  int row_major)
{
    // 1. Variable declarations
    int info = 0;
    double *dwork = NULL;
    int ldwork = 0;

    // Pointers for column-major copies for input/output matrices
    // A, B, C, D, E are input/output.
    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL, *e_cm = NULL;

    // Fortran-style leading dimensions
    int lda_f, ldb_f, ldc_f, ldd_f, lde_f;

    // 2. Input parameter validation
    if (n_param < 0) { info = -1; goto cleanup; }
    if (m_param < 0) { info = -2; goto cleanup; }
    if (np_param < 0) { info = -3; goto cleanup; }

    // Check pointers for required arrays (A, E are N x N; B is N x M; C is NP x N; D is NP x M)
    if (a == NULL && n_param > 0) { info = -4; goto cleanup; }
    if (b == NULL && n_param > 0 && m_param > 0) { info = -6; goto cleanup; }
    if (c == NULL && np_param > 0 && n_param > 0) { info = -8; goto cleanup; }
    if (d == NULL && np_param > 0 && m_param > 0) { info = -10; goto cleanup; }
    if (e == NULL && n_param > 0) { info = -11; goto cleanup; }
    if (nsys == NULL) { info = -13; goto cleanup; }


    // Check leading dimensions
    int min_lda_f = MAX(1, n_param);
    int min_ldb_f = MAX(1, n_param);
    int min_ldc_f = MAX(1, np_param);
    int min_ldd_f = MAX(1, np_param);
    int min_lde_f = MAX(1, n_param);

    if (row_major) { // C LDA is number of columns
        if (n_param > 0 && lda < n_param) { info = -5; goto cleanup; } // A is NxN
        if (n_param > 0 && m_param > 0 && ldb < m_param) { info = -7; goto cleanup; } // B is NxM
        if (np_param > 0 && n_param > 0 && ldc < n_param) { info = -9; goto cleanup; } // C is NPxN
        if (np_param > 0 && m_param > 0 && ldd < m_param) { info = -10; goto cleanup; } // D is NPxM
        if (n_param > 0 && lde < n_param) { info = -12; goto cleanup; } // E is NxN
    } else { // Column-major C (Fortran-style LDs)
        if (n_param > 0 && lda < min_lda_f) { info = -5; goto cleanup; }
        if (n_param > 0 && m_param > 0 && ldb < min_ldb_f) { info = -7; goto cleanup; }
        if (np_param > 0 && n_param > 0 && ldc < min_ldc_f) { info = -9; goto cleanup; }
        if (np_param > 0 && m_param > 0 && ldd < min_ldd_f) { info = -10; goto cleanup; }
        if (n_param > 0 && lde < min_lde_f) { info = -12; goto cleanup; }
    }
    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation
    // LDWORK >= max( 1, 2*N*N + 2*N + N*MAX( 5, N + M + NP ) )
    if (n_param == 0) {
        ldwork = 1;
    } else {
        ldwork = 2 * n_param * n_param + 2 * n_param + n_param * MAX(5, n_param + m_param + np_param);
        ldwork = MAX(1, ldwork);
    }
    
    // Workspace query (optional, SB10JD doc says DWORK(1) returns optimal on exit if INFO=0)
    // For SB10JD, let's try the query first, then fallback to formula if needed.
    int ldwork_query_val = -1;
    double dwork_query_result[1];
    int n_f_dummy = n_param, m_f_dummy = m_param, np_f_dummy = np_param;
    int lda_f_dummy = min_lda_f, ldb_f_dummy = min_ldb_f, ldc_f_dummy = min_ldc_f, ldd_f_dummy = min_ldd_f, lde_f_dummy = min_lde_f;
    int nsys_dummy; // nsys is output, not needed for query if other arrays are NULL

    // Call with NULL for matrices if allowed by Fortran for query, or minimal valid LDs
    // Since A,B,C,D,E are modified, we cannot pass NULL for query if they are accessed.
    // The doc says DWORK(1) on exit, implying a successful run.
    // It is safer to calculate LDWORK using the formula directly for SB10JD,
    // as the routine modifies its input arrays.
    // So, we will use the formula calculated above.

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);

    // 4. Memory allocation for column-major copies
    size_t a_size_elem = (size_t)n_param * n_param; if (n_param == 0) a_size_elem = 0;
    size_t b_size_elem = (size_t)n_param * m_param; if (n_param == 0 || m_param == 0) b_size_elem = 0;
    size_t c_size_elem = (size_t)np_param * n_param; if (np_param == 0 || n_param == 0) c_size_elem = 0;
    size_t d_size_elem = (size_t)np_param * m_param; if (np_param == 0 || m_param == 0) d_size_elem = 0;
    size_t e_size_elem = (size_t)n_param * n_param; if (n_param == 0) e_size_elem = 0;

    if (row_major) {
        if (a_size_elem > 0) { a_cm = (double*)malloc(a_size_elem * sizeof(double)); CHECK_ALLOC(a_cm); }
        if (b_size_elem > 0) { b_cm = (double*)malloc(b_size_elem * sizeof(double)); CHECK_ALLOC(b_cm); }
        if (c_size_elem > 0) { c_cm = (double*)malloc(c_size_elem * sizeof(double)); CHECK_ALLOC(c_cm); }
        if (d_size_elem > 0) { d_cm = (double*)malloc(d_size_elem * sizeof(double)); CHECK_ALLOC(d_cm); }
        if (e_size_elem > 0) { e_cm = (double*)malloc(e_size_elem * sizeof(double)); CHECK_ALLOC(e_cm); }
    }

    // 5. Prepare Fortran parameters and perform conversions
    double* a_ptr = a, *b_ptr = b, *c_ptr = c, *d_ptr = d, *e_ptr = e;

    lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd; lde_f = lde;

    if (row_major) {
        lda_f = MAX(1, n_param); ldb_f = MAX(1, n_param); ldc_f = MAX(1, np_param);
        ldd_f = MAX(1, np_param); lde_f = MAX(1, n_param);

        if (a_size_elem > 0) { slicot_transpose_to_fortran_with_ld(a, a_cm, n_param, n_param, lda, lda_f, sizeof(double)); a_ptr = a_cm; } else { a_ptr = NULL; }
        if (b_size_elem > 0) { slicot_transpose_to_fortran_with_ld(b, b_cm, n_param, m_param, ldb, ldb_f, sizeof(double)); b_ptr = b_cm; } else { b_ptr = NULL; }
        if (c_size_elem > 0) { slicot_transpose_to_fortran_with_ld(c, c_cm, np_param, n_param, ldc, ldc_f, sizeof(double)); c_ptr = c_cm; } else { c_ptr = NULL; }
        if (d_size_elem > 0) { slicot_transpose_to_fortran_with_ld(d, d_cm, np_param, m_param, ldd, ldd_f, sizeof(double)); d_ptr = d_cm; } else { d_ptr = NULL; }
        if (e_size_elem > 0) { slicot_transpose_to_fortran_with_ld(e, e_cm, n_param, n_param, lde, lde_f, sizeof(double)); e_ptr = e_cm; } else { e_ptr = NULL; }
    } else { // Column-major C: ensure NULL pointers for zero-sized arrays
        if (a_size_elem == 0) a_ptr = NULL;
        if (b_size_elem == 0) b_ptr = NULL;
        if (c_size_elem == 0) c_ptr = NULL;
        if (d_size_elem == 0) d_ptr = NULL;
        if (e_size_elem == 0) e_ptr = NULL;
    }

    int n_f_call = n_param, m_f_call = m_param, np_f_call = np_param;
    int ldwork_f_call = ldwork;
    int nsys_f_call; // Fortran will write to this

    // 7. Call Fortran function
    F77_FUNC(sb10jd, SB10JD)(&n_f_call, &m_f_call, &np_f_call,
                             a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f, e_ptr, &lde_f,
                             &nsys_f_call, dwork, &ldwork_f_call, &info);
    
    // Update C nsys from Fortran output
    if (nsys != NULL) {
        *nsys = nsys_f_call;
    }

    // 8. Convert results back to row-major (A, B, C, D are modified in place by Fortran)
    // E is modified but its output is not useful.
    if (row_major && info == 0) {
        // The output NSYS determines the actual dimensions of the modified A, B, C.
        // D is NPxM, its dimensions don't change based on NSYS.
        int current_nsys = (nsys != NULL) ? *nsys : 0;

        // A becomes NSYS x NSYS
        if (a_size_elem > 0 && current_nsys > 0 && a_cm != NULL) { // Check a_cm as it's the source
             slicot_transpose_to_c_with_ld(a_cm, a, current_nsys, current_nsys, lda_f, lda, sizeof(double));
        } else if (current_nsys == 0 && a != NULL && n_param > 0) { // If system becomes 0-order, clear original A
            // This might not be desired, depends on convention for 0-order output.
            // For now, we assume the caller expects the first current_nsys x current_nsys part.
        }

        // B becomes NSYS x M
        if (b_size_elem > 0 && current_nsys > 0 && m_param > 0 && b_cm != NULL) {
            slicot_transpose_to_c_with_ld(b_cm, b, current_nsys, m_param, ldb_f, ldb, sizeof(double));
        }

        // C becomes NP x NSYS
        if (c_size_elem > 0 && np_param > 0 && current_nsys > 0 && c_cm != NULL) {
            slicot_transpose_to_c_with_ld(c_cm, c, np_param, current_nsys, ldc_f, ldc, sizeof(double));
        }
        
        // D is NP x M, dimensions unchanged by NSYS, but content might be Dd.
        if (d_size_elem > 0 && d_cm != NULL) {
            slicot_transpose_to_c_with_ld(d_cm, d, np_param, m_param, ldd_f, ldd, sizeof(double));
        }
    }


cleanup:
    free(dwork);
    if (row_major) {
        free(a_cm); free(b_cm); free(c_cm); free(d_cm); free(e_cm);
    }
    
    if (info == SLICOT_MEMORY_ERROR) {
       // Memory allocation error was caught by CHECK_ALLOC
    }
    return info;
}
