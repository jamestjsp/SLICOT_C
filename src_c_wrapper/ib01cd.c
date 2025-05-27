/**
 * @file ib01cd.c
 * @brief C wrapper for SLICOT routine IB01CD.
 * @details Estimating the initial state and system matrices B and D using A, B, and input-output data.
 * This routine estimates the initial state and, optionally, the system matrices B and D of a linear 
 * time-invariant discrete-time system. Workspace is allocated internally by this wrapper.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h>
#include <ctype.h>
#include <stddef.h>
#include <math.h>

#include "ib01cd.h"
#include "slicot_utils.h"
#include "slicot_f77.h"

/* External Fortran routine declaration */
extern void F77_FUNC(ib01cd, IB01CD)(
    const char* jobx0, const char* comuse, const char* job,
    const int* n, const int* m, const int* l, const int* nsmp,
    const double* a, const int* lda,
    double* b, const int* ldb,
    const double* c, const int* ldc,
    double* d, const int* ldd,
    const double* u, const int* ldu,
    const double* y, const int* ldy,
    double* x0,
    double* v, const int* ldv,
    const double* tol,
    int* iwork, double* dwork, const int* ldwork,
    int* iwarn, int* info,
    int jobx0_len, int comuse_len, int job_len
);

SLICOT_EXPORT
int slicot_ib01cd(const char* jobx0, const char* comuse, const char* job,
                  int n, int m, int l, int nsmp,
                  const double* a, int lda,
                  double* b, int ldb,
                  const double* c, int ldc,
                  double* d, int ldd,
                  const double* u, int ldu,
                  const double* y, int ldy,
                  double* x0,
                  double* v, int ldv,
                  double tol,
                  int* iwarn,
                  int row_major)
{
    int info = 0;
    int local_iwarn = 0;
    char jobx0_upper, comuse_upper, job_upper;

    int *iwork = NULL;
    double *dwork = NULL;
    int liwork = 0;
    int ldwork_alloc = 0;

    const double* a_cm = NULL;
    double* b_cm = NULL;
    const double* c_cm = NULL;
    double* d_cm = NULL;
    const double* u_cm = NULL;
    const double* y_cm = NULL;
    double* v_cm = NULL;

    const double* a_ptr = a;
    double* b_ptr = b;
    const double* c_ptr = c;
    double* d_ptr = d;
    const double* u_ptr = u;
    const double* y_ptr = y;
    double* v_ptr = v;

    int lda_f, ldb_f, ldc_f, ldd_f, ldu_f, ldy_f, ldv_f;

    size_t a_size = 0, b_size = 0, c_size = 0, d_size = 0;
    size_t u_size = 0, y_size = 0, v_size = 0;

    const int jobx0_len = 1, comuse_len = 1, job_len = 1;

    // Input parameter validation
    if (jobx0 == NULL) { info = -1; goto cleanup; }
    jobx0_upper = toupper(*jobx0);
    if (jobx0_upper != 'X' && jobx0_upper != 'N') { info = -1; goto cleanup; }

    if (comuse == NULL) { info = -2; goto cleanup; }
    comuse_upper = toupper(*comuse);
    if (comuse_upper != 'C' && comuse_upper != 'U' && comuse_upper != 'N') { info = -2; goto cleanup; }

    if (job == NULL) { info = -3; goto cleanup; }
    job_upper = toupper(*job);
    if (job_upper != 'B' && job_upper != 'D') { info = -3; goto cleanup; }

    if (n < 0) { info = -4; goto cleanup; }
    if (m < 0) { info = -5; goto cleanup; }
    if (l <= 0) { info = -6; goto cleanup; }
    if (nsmp < 0) { info = -7; goto cleanup; }

    // Matrix A validation (N x N)
    if (n > 0 && (jobx0_upper == 'X' || comuse_upper == 'C')) {
        if (a == NULL) { info = -8; goto cleanup; }
        if (row_major) {
            if (lda < n) { info = -9; goto cleanup; }
        } else {
            if (lda < MAX(1, n)) { info = -9; goto cleanup; }
        }
    } else if (a != NULL && lda < 1) { info = -9; goto cleanup; }

    // Matrix B validation (N x M)
    if (n > 0 && m > 0 && (comuse_upper == 'C' || (jobx0_upper == 'X' && comuse_upper == 'U'))) {
        if (b == NULL) { info = -10; goto cleanup; }
        if (row_major) {
            if (ldb < m) { info = -11; goto cleanup; }
        } else {
            if (ldb < MAX(1, n)) { info = -11; goto cleanup; }
        }
    } else if (b != NULL && ldb < 1) { info = -11; goto cleanup; }

    // Matrix C validation (L x N)
    if (n > 0 && (jobx0_upper == 'X' || comuse_upper == 'C')) {
        if (c == NULL) { info = -12; goto cleanup; }
        if (row_major) {
            if (ldc < n) { info = -13; goto cleanup; }
        } else {
            if (ldc < l) { info = -13; goto cleanup; }
        }
    } else if (c != NULL && ldc < 1) { info = -13; goto cleanup; }

    // Matrix D validation (L x M)
    if (m > 0 && comuse_upper == 'C' && job_upper == 'D') {
        if (d == NULL) { info = -14; goto cleanup; }
        if (row_major) {
            if (ldd < m) { info = -15; goto cleanup; }
        } else {
            if (ldd < l) { info = -15; goto cleanup; }
        }
    } else if (d != NULL && ldd < 1) { info = -15; goto cleanup; }

    // Matrix U validation (NSMP x M)
    if (m > 0 && nsmp > 0 && (comuse_upper == 'C' || (jobx0_upper == 'X' && comuse_upper == 'U'))) {
        if (u == NULL) { info = -16; goto cleanup; }
        if (row_major) {
            if (ldu < m) { info = -17; goto cleanup; }
        } else {
            if (ldu < MAX(1, nsmp)) { info = -17; goto cleanup; }
        }
    } else if (u != NULL && ldu < 1) { info = -17; goto cleanup; }

    // Matrix Y validation (NSMP x L)
    if (nsmp > 0 && (jobx0_upper == 'X' || comuse_upper == 'C')) {
        if (y == NULL) { info = -18; goto cleanup; }
        if (row_major) {
            if (ldy < l) { info = -19; goto cleanup; }
        } else {
            if (ldy < MAX(1, nsmp)) { info = -19; goto cleanup; }
        }
    } else if (y != NULL && ldy < 1) { info = -19; goto cleanup; }

    // X0 validation
    if (jobx0_upper == 'X' && n > 0 && x0 == NULL) { info = -20; goto cleanup; }

    // Matrix V validation (N x N)
    if (n > 0 && (jobx0_upper == 'X' || comuse_upper == 'C')) {
        if (v == NULL) { info = -21; goto cleanup; }
        if (row_major) {
            if (ldv < n) { info = -22; goto cleanup; }
        } else {
            if (ldv < MAX(1, n)) { info = -22; goto cleanup; }
        }
    } else if (v != NULL && ldv < 1) { info = -22; goto cleanup; }

    if (info != 0) { goto cleanup; }

    // Internal Workspace Allocation
    // Calculate LIWORK
    if (jobx0_upper == 'N' && comuse_upper != 'C') {
        liwork = 0;
    } else if (jobx0_upper == 'X' && comuse_upper != 'C') {
        liwork = n;
    } else if (comuse_upper == 'C') {
        int a_val = (jobx0_upper == 'N') ? 0 : n;
        if (job_upper == 'B') {
            liwork = n * m + a_val;
        } else { // job_upper == 'D'
            liwork = MAX(n * m + a_val, m);
        }
    }

    if (liwork > 0) {
        iwork = (int*)malloc(liwork * sizeof(int));
        CHECK_ALLOC(iwork);
    }

    // Calculate LDWORK
    if (jobx0_upper == 'N' && comuse_upper != 'C') {
        ldwork_alloc = 2;
    } else if (MAX(n, m) == 0) {
        ldwork_alloc = 2;
    } else {
        // Complex workspace calculation based on documentation
        int ldw1, ldw2, ldw3;
        int r, a_val, b_val, c_val, d_val, f_val, q_val;
        
        a_val = (jobx0_upper == 'N') ? 0 : n;
        r = n * m + a_val;
        b_val = (job_upper == 'B') ? 0 : l * m;
        c_val = (jobx0_upper == 'N') ? 0 : l * n;
        d_val = (jobx0_upper == 'N') ? 0 : (2 * n * n + n);
        
        if (job_upper == 'B' || m == 0) {
            f_val = 2 * r;
        } else {
            f_val = m + MAX(2 * r, m);
        }
        
        q_val = b_val + r * l;

        if (comuse_upper == 'C') {
            ldw1 = (m == 0 || job_upper == 'B') ? 2 : 3;
            
            int ldwa = nsmp * l * (r + 1) + MAX(n + MAX(d_val, f_val), 6 * r);
            ldw2 = ldwa;
            if (m > 0 && job_upper == 'D') {
                ldw2 = MAX(ldwa, nsmp * l * (r + 1) + 2 * m * m + 6 * m);
            }
            
            int ldwb = (b_val + r) * (r + 1) + MAX(q_val * (r + 1) + n * n * m + c_val + MAX(d_val, f_val), 6 * r);
            ldw3 = ldwb;
            if (m > 0 && job_upper == 'D') {
                ldw3 = MAX(ldwb, (b_val + r) * (r + 1) + 2 * m * m + 6 * m);
            }
        } else {
            ldw1 = 2;
            ldw2 = nsmp * l * (n + 1) + 2 * n + MAX(2 * n * n, 4 * n);
            ldw3 = n * (n + 1) + 2 * n + MAX(n * l * (n + 1) + 2 * n * n + l * n, 4 * n);
        }
        
        ldwork_alloc = ldw1 + n * (n + m + l) + MAX(5 * n, MAX(ldw1, MIN(ldw2, ldw3)));
    }

    ldwork_alloc = MAX(2, ldwork_alloc);
    if (ldwork_alloc > 0) {
        dwork = (double*)malloc(ldwork_alloc * sizeof(double));
        CHECK_ALLOC(dwork);
    }

    // Memory allocation for column-major copies (if row_major)
    if (n > 0) a_size = (size_t)n * n;
    if (n > 0 && m > 0) b_size = (size_t)n * m;
    if (l > 0 && n > 0) c_size = (size_t)l * n;
    if (l > 0 && m > 0) d_size = (size_t)l * m;
    if (nsmp > 0 && m > 0) u_size = (size_t)nsmp * m;
    if (nsmp > 0 && l > 0) y_size = (size_t)nsmp * l;
    if (n > 0) v_size = (size_t)n * n;

    if (row_major) {
        if (a_size > 0 && a != NULL) {
            a_cm = (double*)malloc(a_size * sizeof(double));
            CHECK_ALLOC(a_cm);
        }
        if (b_size > 0 && b != NULL) {
            b_cm = (double*)malloc(b_size * sizeof(double));
            CHECK_ALLOC(b_cm);
        }
        if (c_size > 0 && c != NULL) {
            c_cm = (double*)malloc(c_size * sizeof(double));
            CHECK_ALLOC(c_cm);
        }
        if (d_size > 0 && d != NULL) {
            d_cm = (double*)malloc(d_size * sizeof(double));
            CHECK_ALLOC(d_cm);
        }
        if (u_size > 0 && u != NULL) {
            u_cm = (double*)malloc(u_size * sizeof(double));
            CHECK_ALLOC(u_cm);
        }
        if (y_size > 0 && y != NULL) {
            y_cm = (double*)malloc(y_size * sizeof(double));
            CHECK_ALLOC(y_cm);
        }
        if (v_size > 0 && v != NULL) {
            v_cm = (double*)malloc(v_size * sizeof(double));
            CHECK_ALLOC(v_cm);
        }
    }

    // Prepare Fortran parameters and perform conversions
    lda_f = (n == 0) ? 1 : MAX(1, n);
    ldb_f = (n == 0) ? 1 : MAX(1, n);
    ldc_f = MAX(1, l);
    ldd_f = MAX(1, l);
    ldu_f = (nsmp == 0) ? 1 : MAX(1, nsmp);
    ldy_f = (nsmp == 0) ? 1 : MAX(1, nsmp);
    ldv_f = (n == 0) ? 1 : MAX(1, n);

    if (row_major) {
        if (a_size > 0 && a != NULL) {
            slicot_transpose_to_fortran_with_ld(a, (double*)a_cm, n, n, lda, lda_f, sizeof(double));
            a_ptr = a_cm;
        } else {
            a_ptr = a;
        }
        if (b_size > 0 && b != NULL) {
            slicot_transpose_to_fortran_with_ld(b, b_cm, n, m, ldb, ldb_f, sizeof(double));
            b_ptr = b_cm;
        } else {
            b_ptr = b;
        }
        if (c_size > 0 && c != NULL) {
            slicot_transpose_to_fortran_with_ld(c, (double*)c_cm, l, n, ldc, ldc_f, sizeof(double));
            c_ptr = c_cm;
        } else {
            c_ptr = c;
        }
        if (d_size > 0 && d != NULL) {
            slicot_transpose_to_fortran_with_ld(d, d_cm, l, m, ldd, ldd_f, sizeof(double));
            d_ptr = d_cm;
        } else {
            d_ptr = d;
        }
        if (u_size > 0 && u != NULL) {
            slicot_transpose_to_fortran_with_ld(u, (double*)u_cm, nsmp, m, ldu, ldu_f, sizeof(double));
            u_ptr = u_cm;
        } else {
            u_ptr = u;
        }
        if (y_size > 0 && y != NULL) {
            slicot_transpose_to_fortran_with_ld(y, (double*)y_cm, nsmp, l, ldy, ldy_f, sizeof(double));
            y_ptr = y_cm;
        } else {
            y_ptr = y;
        }
        if (v_size > 0 && v != NULL) {
            v_ptr = v_cm;
        } else {
            v_ptr = v;
        }
    } else {
        lda_f = (a != NULL) ? lda : 1;
        ldb_f = (b != NULL) ? ldb : 1;
        ldc_f = (c != NULL) ? ldc : 1;
        ldd_f = (d != NULL) ? ldd : 1;
        ldu_f = (u != NULL) ? ldu : 1;
        ldy_f = (y != NULL) ? ldy : 1;
        ldv_f = (v != NULL) ? ldv : 1;
    }

    // Call Fortran function
    F77_FUNC(ib01cd, IB01CD)(
        &jobx0_upper, &comuse_upper, &job_upper,
        &n, &m, &l, &nsmp,
        a_ptr, &lda_f,
        b_ptr, &ldb_f,
        c_ptr, &ldc_f,
        d_ptr, &ldd_f,
        u_ptr, &ldu_f,
        y_ptr, &ldy_f,
        x0,
        v_ptr, &ldv_f,
        &tol,
        iwork, dwork, &ldwork_alloc,
        &local_iwarn, &info,
        jobx0_len, comuse_len, job_len
    );

    if (iwarn != NULL) {
        *iwarn = local_iwarn;
    }

    // Convert results back to row-major (if needed)
    if (row_major && info == 0) {
        if (b_size > 0 && b != NULL && comuse_upper == 'C') {
            slicot_transpose_to_c_with_ld(b_ptr, b, n, m, ldb_f, ldb, sizeof(double));
        }
        if (d_size > 0 && d != NULL && comuse_upper == 'C' && job_upper == 'D') {
            slicot_transpose_to_c_with_ld(d_ptr, d, l, m, ldd_f, ldd, sizeof(double));
        }
        if (v_size > 0 && v != NULL) {
            slicot_transpose_to_c_with_ld(v_ptr, v, n, n, ldv_f, ldv, sizeof(double));
        }
    }

cleanup:
    free(iwork);
    free(dwork);
    if (row_major) {
        free((void*)a_cm);
        free(b_cm);
        free((void*)c_cm);
        free(d_cm);
        free((void*)u_cm);
        free((void*)y_cm);
        free(v_cm);
    }

    return info;
}
