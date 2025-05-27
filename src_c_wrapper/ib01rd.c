/**
 * @file ib01rd.c
 * @brief C wrapper for SLICOT routine IB01RD.
 * @details Estimates the initial state of a linear time-invariant (LTI)
 * discrete-time system, given the system matrices (A,B,C,D) and
 * the input and output trajectories of the system. Matrix A is
 * assumed to be in a real Schur form. Workspace (IWORK, DWORK)
 * is allocated internally by this wrapper. Input/output matrix
 * format is handled via the row_major parameter.
 */

#include <stdlib.h> // For malloc, free
#include <ctype.h>  // For toupper
#include <stddef.h> // For size_t
#include <math.h>   // For fabs (used by MAX/MIN if they involve floating point)
#include <stdio.h>  // For error logging (optional)

#include "ib01rd.h"       // Public header for this wrapper
#include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX/MIN, transpose functions etc.
#include "slicot_f77.h"   // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */
extern void F77_FUNC(ib01rd, IB01RD)(
    const char* job, const int* n, const int* m, const int* l, const int* nsmp,
    const double* a, const int* lda,
    const double* b, const int* ldb,
    const double* c, const int* ldc,
    const double* d, const int* ldd,
    const double* u, const int* ldu,
    const double* y, const int* ldy,
    double* x0,
    const double* tol,
    int* iwork,
    double* dwork, const int* ldwork,
    int* iwarn, int* info,
    int job_len);

/* C wrapper function definition */
SLICOT_EXPORT
int slicot_ib01rd(
    char job, int n, int m, int l, int nsmp,
    const double* a, int lda,
    const double* b, int ldb,
    const double* c, int ldc,
    const double* d, int ldd,
    const double* u, int ldu,
    const double* y, int ldy,
    double* x0,
    double tol,
    int* iwarn, /* Output warning indicator */
    double* dwork_info, /* Optional output: dwork[0]=opt_ldwork, dwork[1]=rcond */
    int row_major)
{
    // 1. Variable declarations
    int info = 0;
    int local_iwarn = 0;
    char job_upper;

    int *iwork = NULL;
    double *dwork = NULL;
    int liwork = 0;
    int ldwork_alloc = 0;

    double* a_cm = NULL;
    double* b_cm = NULL;
    double* c_cm = NULL;
    double* d_cm = NULL;
    double* u_cm = NULL;
    double* y_cm = NULL;

    const double* a_ptr = a;
    const double* b_ptr = b;
    const double* c_ptr = c;
    const double* d_ptr = d;
    const double* u_ptr = u;
    const double* y_ptr = y;

    int lda_f, ldb_f, ldc_f, ldd_f, ldu_f, ldy_f;
    size_t a_size = 0, b_size = 0, c_size = 0, d_size = 0, u_size = 0, y_size = 0;

    const int job_len = 1;

    // 2. Input parameter validation
    job_upper = toupper(job);
    if (job_upper != 'Z' && job_upper != 'N') { info = -1; goto cleanup; }
    if (n < 0) { info = -2; goto cleanup; }
    if (m < 0) { info = -3; goto cleanup; }
    if (l <= 0) { info = -4; goto cleanup; } // L must be > 0
    if (nsmp < n) { info = -5; goto cleanup; }

    // Validate matrix A and LDA
    if (n > 0) {
        if (a == NULL) { info = -6; goto cleanup; }
        if (row_major) { // C LDA is columns
            if (lda < n) { info = -7; goto cleanup; }
        } else { // Fortran LDA is rows
            if (lda < MAX(1,n)) { info = -7; goto cleanup; }
        }
    } else { // N == 0
         if (lda < 1 && a != NULL) {info = -7; goto cleanup;}
    }


    // Validate matrix B and LDB
    if (n > 0 && m > 0) {
        if (b == NULL) { info = -8; goto cleanup; }
        if (row_major) { // C LDB is columns
            if (ldb < m) { info = -9; goto cleanup; }
        } else { // Fortran LDB is rows
            if (ldb < n) { info = -9; goto cleanup; }
        }
    } else { // N == 0 or M == 0, B not strictly required or referenced
        if (ldb < 1 && b != NULL ) { info = -9; goto cleanup; }
    }

    // Validate matrix C and LDC
    // C is L-by-N. If N=0, C is L-by-0 (size 0), and c can be NULL.
    if (l > 0) { // L is always > 0 (checked earlier)
        if (n > 0 && c == NULL) { // C is required only if N > 0 (since L > 0)
            info = -10; goto cleanup;
        }
        // If c is not NULL (even if N=0), or if N > 0 (meaning c must be non-NULL), then LDC must be valid.
        if (c != NULL || n > 0) {
            if (row_major) { // C LDC is number of columns (N)
                if (n > 0 && ldc < n) { info = -11; goto cleanup; }
                // If N=0 and c is not NULL, LDC (cols) must be >=1 for a non-NULL ptr.
                if (n == 0 && c != NULL && ldc < 1 ) { info = -11; goto cleanup; }
            } else { // Fortran LDC is number of rows (L)
                if (ldc < l) { info = -11; goto cleanup; }
            }
        }
    }


    // Validate matrix D and LDD
    if (m > 0 && job_upper == 'N') {
        if (d == NULL) { info = -12; goto cleanup; }
        if (row_major) { // C LDD is columns (M)
            if (ldd < m) { info = -13; goto cleanup; }
        } else { // Fortran LDD is rows (L)
            if (ldd < l) { info = -13; goto cleanup; }
        }
    } else { // M == 0 or JOB == 'Z', D not strictly required or referenced
         if (ldd < 1 && d != NULL) { info = -13; goto cleanup; }
    }

    // Validate matrix U and LDU
    if (m > 0) {
        if (u == NULL && nsmp > 0) { info = -14; goto cleanup; }
        if (row_major) { // C LDU is columns (M)
            if (ldu < m) { info = -15; goto cleanup; }
        } else { // Fortran LDU is rows (NSMP)
            if (ldu < MAX(1,nsmp)) { info = -15; goto cleanup; }
        }
    } else { // M == 0
        if (ldu < 1 && u != NULL) { info = -15; goto cleanup; }
    }

    // Validate matrix Y and LDY
    // Y is NSMP-by-L. If NSMP=0, Y is 0-by-L (size 0), and y can be NULL.
    if (l > 0) { // L is always > 0
        if (nsmp > 0 && y == NULL ) { info = -16; goto cleanup; }
         if (y != NULL || nsmp > 0) { // Check LDY if y is provided or required
            if (row_major) { // C LDY is columns (L)
                if (ldy < l) { info = -17; goto cleanup; }
            } else { // Fortran LDY is rows (NSMP)
                if (ldy < MAX(1,nsmp)) { info = -17; goto cleanup; }
            }
        }
    }


    // Validate X0 (output, pointer must be valid if N > 0)
    if (n > 0 && x0 == NULL) { info = -18; goto cleanup; }

    // Validate TOL
    if (tol > 1.0) { info = -19; goto cleanup; }


    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation
    if (n > 0) {
        liwork = n;
        iwork = (int*)malloc((size_t)liwork * sizeof(int));
        CHECK_ALLOC(iwork);
    } else {
        liwork = 0; 
        iwork = NULL; 
    }

    long long t_ll = nsmp; 
    long long n_ll = n;
    long long l_ll = l;
    long long q_ll = n_ll * l_ll;

    long long ldw1_ll;
    if (n_ll == 0) { // Simplified formula for N=0 if applicable, or ensure it evaluates reasonably
        // Based on LDW1 structure, if N=0, terms with N become 0.
        // LDW1 = t*L*(0 + 1) + 0 + max(0, 0) = t*L
        // However, DWORK is always min 2.
        ldw1_ll = t_ll * l_ll;
    } else {
        ldw1_ll = t_ll * l_ll * (n_ll + 1) + 2 * n_ll + MAX(2 * n_ll * n_ll, 4 * n_ll);
    }
    ldwork_alloc = (int)MAX(2LL, ldw1_ll);


    if (ldwork_alloc > 0) {
        dwork = (double*)malloc((size_t)ldwork_alloc * sizeof(double));
        CHECK_ALLOC(dwork);
    } else { 
        dwork = NULL;
    }


    // 4. Memory allocation for column-major copies (if row_major)
    if (n > 0) a_size = (size_t)n * n; else a_size = 0;
    if (n > 0 && m > 0) b_size = (size_t)n * m; else b_size = 0;
    if (l > 0 && n > 0) c_size = (size_t)l * n; else c_size = 0;
    if (l > 0 && m > 0) d_size = (size_t)l * m; else d_size = 0;
    if (nsmp > 0 && m > 0) u_size = (size_t)nsmp * m; else u_size = 0;
    if (nsmp > 0 && l > 0) y_size = (size_t)nsmp * l; else y_size = 0;


    if (row_major) {
        if (a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
        if (b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
        if (c_size > 0) { c_cm = (double*)malloc(c_size * sizeof(double)); CHECK_ALLOC(c_cm); }
        if (d_size > 0 && job_upper == 'N') { d_cm = (double*)malloc(d_size * sizeof(double)); CHECK_ALLOC(d_cm); }
        if (u_size > 0) { u_cm = (double*)malloc(u_size * sizeof(double)); CHECK_ALLOC(u_cm); }
        if (y_size > 0) { y_cm = (double*)malloc(y_size * sizeof(double)); CHECK_ALLOC(y_cm); }
    }

    // 5. Prepare Fortran parameters and perform conversions
    lda_f = (n == 0) ? 1 : MAX(1, n); 
    ldb_f = (n == 0) ? 1 : MAX(1, n); 
    ldc_f = MAX(1, l);                
    ldd_f = MAX(1, l);                
    ldu_f = MAX(1, nsmp);             
    ldy_f = MAX(1, nsmp);             

    if (row_major) {
        if (a_size > 0) { slicot_transpose_to_fortran_with_ld(a, a_cm, n, n, lda, lda_f, sizeof(double)); a_ptr = a_cm; }
        else { a_ptr = NULL; }

        if (b_size > 0) { slicot_transpose_to_fortran_with_ld(b, b_cm, n, m, ldb, ldb_f, sizeof(double)); b_ptr = b_cm; }
        else { b_ptr = NULL; }

        if (c_size > 0) { slicot_transpose_to_fortran_with_ld(c, c_cm, l, n, ldc, ldc_f, sizeof(double)); c_ptr = c_cm; }
        else { c_ptr = NULL; } 

        if (d_size > 0 && job_upper == 'N') { slicot_transpose_to_fortran_with_ld(d, d_cm, l, m, ldd, ldd_f, sizeof(double)); d_ptr = d_cm; }
        else { d_ptr = NULL; } 

        if (u_size > 0) { slicot_transpose_to_fortran_with_ld(u, u_cm, nsmp, m, ldu, ldu_f, sizeof(double)); u_ptr = u_cm; }
        else { u_ptr = NULL; }

        if (y_size > 0) { slicot_transpose_to_fortran_with_ld(y, y_cm, nsmp, l, ldy, ldy_f, sizeof(double)); y_ptr = y_cm; }
        else { y_ptr = NULL; }

    } else { // Column-major C
        if (a_size == 0) a_ptr = NULL; lda_f = lda;
        if (b_size == 0) b_ptr = NULL; ldb_f = ldb;
        if (c_size == 0) c_ptr = NULL; ldc_f = ldc; 
        if (d_size == 0 || job_upper == 'Z') d_ptr = NULL; ldd_f = ldd; 
        if (u_size == 0) u_ptr = NULL; ldu_f = ldu;
        if (y_size == 0) y_ptr = NULL; ldy_f = ldy;
    }
    
    // Ensure Fortran LDs are at least 1 if the corresponding pointer is not NULL.
    // This is important for dummy arguments in Fortran.
    if (a_ptr != NULL && lda_f < 1) lda_f = 1;
    if (b_ptr != NULL && ldb_f < 1) ldb_f = 1;
    if (c_ptr != NULL && ldc_f < 1) ldc_f = 1;
    if (d_ptr != NULL && ldd_f < 1) ldd_f = 1;
    if (u_ptr != NULL && ldu_f < 1) ldu_f = 1;
    if (y_ptr != NULL && ldy_f < 1) ldy_f = 1;


    // 7. Call Fortran function
    F77_FUNC(ib01rd, IB01RD)(&job_upper, &n, &m, &l, &nsmp,
                             a_ptr, &lda_f,
                             b_ptr, &ldb_f,
                             c_ptr, &ldc_f,
                             d_ptr, &ldd_f,
                             u_ptr, &ldu_f,
                             y_ptr, &ldy_f,
                             x0, 
                             &tol,
                             iwork,
                             dwork, &ldwork_alloc,
                             &local_iwarn, &info,
                             job_len);

    if (iwarn != NULL) {
        *iwarn = local_iwarn;
    }

    if (dwork_info != NULL && dwork != NULL && ldwork_alloc >=2 && info == 0) {
        dwork_info[0] = dwork[0]; 
        dwork_info[1] = dwork[1]; 
    } else if (dwork_info != NULL && dwork != NULL && ldwork_alloc >=1 && info == -22) {
        dwork_info[0] = dwork[0]; 
        if (ldwork_alloc >=2) dwork_info[1] = 0.0; 
    } else if (dwork_info != NULL) {
        if (ldwork_alloc >=1) dwork_info[0] = 0.0;
        if (ldwork_alloc >=2) dwork_info[1] = 0.0;
    }


    // 8. Convert results back to row-major (if needed) - None for IB01RD

cleanup:
    free(iwork);
    free(dwork);
    if (row_major) {
        free(a_cm);
        free(b_cm);
        free(c_cm);
        free(d_cm);
        free(u_cm);
        free(y_cm);
    }

    return info;
}
