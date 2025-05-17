/**
 * @file sb10hd.c
 * @brief C wrapper for SLICOT routine SB10HD.
 * @details Computes the matrices of the H2 optimal n-state controller
 * for a continuous-time system.
 * Workspace (IWORK, DWORK, BWORK) is allocated internally by this wrapper.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h> // For malloc, free
#include <string.h> // For memcpy, memset
#include <math.h>   // For fabs, fmax (from C standard library)
#include "sb10hd.h" // Public header for this wrapper
#include "slicot_utils.h"  // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions etc.
#include "slicot_f77.h"    // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */
extern void F77_FUNC(sb10hd, SB10HD)(
    const int* n_fortran, const int* m_fortran, const int* np_fortran, const int* ncon_fortran, const int* nmeas_fortran,
    const double* a_fortran, const int* lda_fortran,
    const double* b_fortran, const int* ldb_fortran,
    const double* c_fortran, const int* ldc_fortran,
    const double* d_fortran, const int* ldd_fortran,
    double* ak_fortran, const int* ldak_fortran,
    double* bk_fortran, const int* ldbk_fortran,
    double* ck_fortran, const int* ldck_fortran,
    double* dk_fortran, const int* lddk_fortran,
    double* rcond_fortran, const double* tol_fortran,
    int* iwork_fortran, double* dwork_fortran, const int* ldwork_fortran,
    int* bwork_fortran, /* Fortran LOGICAL, C int */
    int* info_fortran
);

/* C wrapper function definition */
SLICOT_EXPORT
int slicot_sb10hd(int n_param, int m_param, int np_param, int ncon_param, int nmeas_param,
                  double* a, int lda,
                  double* b, int ldb,
                  double* c, int ldc,
                  double* d, int ldd,
                  double* ak, int ldak,
                  double* bk, int ldbk,
                  double* ck, int ldck,
                  double* dk, int lddk,
                  double* rcond, double tol,
                  int row_major)
{
    // 1. Variable declarations
    int info = 0;
    int *iwork = NULL;
    double *dwork = NULL;
    int *bwork = NULL; 
    int liwork = 0;
    int ldwork = 0;
    int lbwork = 0;

    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
    double *ak_cm = NULL, *bk_cm = NULL, *ck_cm = NULL, *dk_cm = NULL;

    int lda_f, ldb_f, ldc_f, ldd_f;
    int ldak_f, ldbk_f, ldck_f, lddk_f;

    double tol_f = tol;

    // 2. Input parameter validation
    if (n_param < 0) { info = -1; goto cleanup; }
    if (m_param < 0) { info = -2; goto cleanup; }
    if (np_param < 0) { info = -3; goto cleanup; }
    if (ncon_param < 0 || ncon_param > m_param) { info = -4; goto cleanup; }
    if (nmeas_param < 0 || nmeas_param > np_param) { info = -5; goto cleanup; }
    
    // Constraints from documentation: NP-NMEAS >= NCON and M-NCON >= NMEAS
    if (np_param - nmeas_param < ncon_param) { info = -4; goto cleanup; } 
    if (m_param - ncon_param < nmeas_param) { info = -5; goto cleanup; }

    if (a == NULL && n_param > 0) { info = -6; goto cleanup; }
    if (b == NULL && n_param > 0 && m_param > 0) { info = -8; goto cleanup; }
    if (c == NULL && np_param > 0 && n_param > 0) { info = -10; goto cleanup; }
    if (d == NULL && np_param > 0 && m_param > 0) { info = -12; goto cleanup; }
    if (ak == NULL && n_param > 0) { info = -13; goto cleanup; }
    if (bk == NULL && n_param > 0 && nmeas_param > 0) { info = -15; goto cleanup; }
    if (ck == NULL && ncon_param > 0 && n_param > 0) { info = -17; goto cleanup; }
    if (dk == NULL && ncon_param > 0 && nmeas_param > 0) { info = -19; goto cleanup; }
    if (rcond == NULL) { info = -20; goto cleanup; }

    int min_lda_f = MAX(1, n_param);
    int min_ldb_f = MAX(1, n_param);
    int min_ldc_f = MAX(1, np_param);
    int min_ldd_f = MAX(1, np_param);
    int min_ldak_f = MAX(1, n_param);
    int min_ldbk_f = MAX(1, n_param);
    int min_ldck_f = MAX(1, ncon_param);
    int min_lddk_f = MAX(1, ncon_param);

    if (row_major) {
        if (n_param > 0 && lda < n_param) { info = -7; goto cleanup; }
        if (n_param > 0 && m_param > 0 && ldb < m_param) { info = -9; goto cleanup; }
        if (np_param > 0 && n_param > 0 && ldc < n_param) { info = -11; goto cleanup; }
        if (np_param > 0 && m_param > 0 && ldd < m_param) { info = -12; goto cleanup; } // D is NPxM, C_LDD is M cols
        if (n_param > 0 && ldak < n_param) { info = -14; goto cleanup; }
        if (n_param > 0 && nmeas_param > 0 && ldbk < nmeas_param) { info = -16; goto cleanup; }
        if (ncon_param > 0 && n_param > 0 && ldck < n_param) { info = -18; goto cleanup; }
        if (ncon_param > 0 && nmeas_param > 0 && lddk < nmeas_param) { info = -19; goto cleanup; } // DK is NCONxNMEAS, C_LDDK is NMEAS cols
    } else { // Column-major C
        if (n_param > 0 && lda < min_lda_f) { info = -7; goto cleanup; }
        if (n_param > 0 && m_param > 0 && ldb < min_ldb_f) { info = -9; goto cleanup; }
        if (np_param > 0 && n_param > 0 && ldc < min_ldc_f) { info = -11; goto cleanup; }
        if (np_param > 0 && m_param > 0 && ldd < min_ldd_f) { info = -12; goto cleanup; }
        if (n_param > 0 && ldak < min_ldak_f) { info = -14; goto cleanup; }
        if (n_param > 0 && nmeas_param > 0 && ldbk < min_ldbk_f) { info = -16; goto cleanup; }
        if (ncon_param > 0 && n_param > 0 && ldck < min_ldck_f) { info = -18; goto cleanup; }
        if (ncon_param > 0 && nmeas_param > 0 && lddk < min_lddk_f) { info = -19; goto cleanup; }
    }
    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation
    // LIWORK = max(2*N,N*N)
    liwork = MAX(1, MAX(2 * n_param, n_param * n_param));
    iwork = (int*)malloc((size_t)liwork * sizeof(int));
    CHECK_ALLOC(iwork);

    // BWORK size is 2*N
    lbwork = 2 * n_param;
    if (n_param > 0) {
        bwork = (int*)malloc((size_t)lbwork * sizeof(int));
        CHECK_ALLOC(bwork);
    } else {
        bwork = NULL;
    }

    // LDWORK calculation based on SLICOT documentation for SB10HD
    int N_ws = n_param;
    int M_ws = m_param;
    int NP_ws = np_param;
    int M2_ws = ncon_param;  // NCON
    int NP2_ws = nmeas_param; // NMEAS
    int M1_ws = M_ws - M2_ws;   // M1 in formula (M - NCON)
    int NP1_ws = NP_ws - NP2_ws; // NP1 in formula (NP - NMEAS)

    int term1 = M2_ws + NP1_ws*NP1_ws + MAX(NP1_ws*N_ws, MAX(3*M2_ws+NP1_ws,5*M2_ws));
    int term2 = NP2_ws + M1_ws*M1_ws + MAX(M1_ws*N_ws, MAX(3*NP2_ws+M1_ws,5*M1_ws));
    int term3 = N_ws*M2_ws;
    int term4 = NP2_ws*N_ws;
    int term5 = NP2_ws*M2_ws;
    int term6 = N_ws*(14*N_ws+12+M2_ws+NP2_ws)+5;
    
    ldwork = N_ws*M_ws + NP_ws*(N_ws+M_ws) + M2_ws*M2_ws + NP2_ws*NP2_ws +
             MAX(1, MAX(MAX(term1,term2), MAX(MAX(term3,term4), MAX(term5,term6))));
    ldwork = MAX(1, ldwork);

    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);

    // 4. Memory allocation for column-major copies
    size_t a_size_elem = (size_t)n_param * n_param; if (n_param == 0) a_size_elem = 0;
    size_t b_size_elem = (size_t)n_param * m_param; if (n_param == 0 || m_param == 0) b_size_elem = 0;
    size_t c_size_elem = (size_t)np_param * n_param; if (np_param == 0 || n_param == 0) c_size_elem = 0;
    size_t d_size_elem = (size_t)np_param * m_param; if (np_param == 0 || m_param == 0) d_size_elem = 0;
    size_t ak_size_elem = (size_t)n_param * n_param; if (n_param == 0) ak_size_elem = 0;
    size_t bk_size_elem = (size_t)n_param * nmeas_param; if (n_param == 0 || nmeas_param == 0) bk_size_elem = 0;
    size_t ck_size_elem = (size_t)ncon_param * n_param; if (ncon_param == 0 || n_param == 0) ck_size_elem = 0;
    size_t dk_size_elem = (size_t)ncon_param * nmeas_param; if (ncon_param == 0 || nmeas_param == 0) dk_size_elem = 0;

    if (row_major) {
        if (a_size_elem > 0) { a_cm = (double*)malloc(a_size_elem * sizeof(double)); CHECK_ALLOC(a_cm); }
        if (b_size_elem > 0) { b_cm = (double*)malloc(b_size_elem * sizeof(double)); CHECK_ALLOC(b_cm); }
        if (c_size_elem > 0) { c_cm = (double*)malloc(c_size_elem * sizeof(double)); CHECK_ALLOC(c_cm); }
        if (d_size_elem > 0) { d_cm = (double*)malloc(d_size_elem * sizeof(double)); CHECK_ALLOC(d_cm); }
        if (ak_size_elem > 0) { ak_cm = (double*)malloc(ak_size_elem * sizeof(double)); CHECK_ALLOC(ak_cm); }
        if (bk_size_elem > 0) { bk_cm = (double*)malloc(bk_size_elem * sizeof(double)); CHECK_ALLOC(bk_cm); }
        if (ck_size_elem > 0) { ck_cm = (double*)malloc(ck_size_elem * sizeof(double)); CHECK_ALLOC(ck_cm); }
        if (dk_size_elem > 0) { dk_cm = (double*)malloc(dk_size_elem * sizeof(double)); CHECK_ALLOC(dk_cm); }
    }

    // 5. Prepare Fortran parameters and perform conversions
    double* a_ptr = a, *b_ptr = b, *c_ptr = c, *d_ptr = d;
    double* ak_ptr = ak, *bk_ptr = bk, *ck_ptr = ck, *dk_ptr = dk;

    lda_f = lda; ldb_f = ldb; ldc_f = ldc; ldd_f = ldd;
    ldak_f = ldak; ldbk_f = ldbk; ldck_f = ldck; lddk_f = lddk;

    if (row_major) {
        lda_f = MAX(1, n_param); ldb_f = MAX(1, n_param); ldc_f = MAX(1, np_param); ldd_f = MAX(1, np_param);
        ldak_f = MAX(1, n_param); ldbk_f = MAX(1, n_param); ldck_f = MAX(1, ncon_param); lddk_f = MAX(1, ncon_param);

        if (a_size_elem > 0) { slicot_transpose_to_fortran_with_ld(a, a_cm, n_param, n_param, lda, lda_f, sizeof(double)); a_ptr = a_cm; } else { a_ptr = NULL; }
        if (b_size_elem > 0) { slicot_transpose_to_fortran_with_ld(b, b_cm, n_param, m_param, ldb, ldb_f, sizeof(double)); b_ptr = b_cm; } else { b_ptr = NULL; }
        if (c_size_elem > 0) { slicot_transpose_to_fortran_with_ld(c, c_cm, np_param, n_param, ldc, ldc_f, sizeof(double)); c_ptr = c_cm; } else { c_ptr = NULL; }
        if (d_size_elem > 0) { slicot_transpose_to_fortran_with_ld(d, d_cm, np_param, m_param, ldd, ldd_f, sizeof(double)); d_ptr = d_cm; } else { d_ptr = NULL; }
        
        ak_ptr = (ak_size_elem > 0) ? ak_cm : NULL;
        bk_ptr = (bk_size_elem > 0) ? bk_cm : NULL;
        ck_ptr = (ck_size_elem > 0) ? ck_cm : NULL;
        dk_ptr = (dk_size_elem > 0) ? dk_cm : NULL;
    } else {
        if (a_size_elem == 0) a_ptr = NULL; if (b_size_elem == 0) b_ptr = NULL;
        if (c_size_elem == 0) c_ptr = NULL; if (d_size_elem == 0) d_ptr = NULL;
        if (ak_size_elem == 0) ak_ptr = NULL; if (bk_size_elem == 0) bk_ptr = NULL;
        if (ck_size_elem == 0) ck_ptr = NULL; if (dk_size_elem == 0) dk_ptr = NULL;
    }

    int n_f_call = n_param, m_f_call = m_param, np_f_call = np_param;
    int ncon_f_call = ncon_param, nmeas_f_call = nmeas_param;
    int ldwork_f_call = ldwork;

    // 7. Call Fortran function
    F77_FUNC(sb10hd, SB10HD)(&n_f_call, &m_f_call, &np_f_call, &ncon_f_call, &nmeas_f_call,
                             a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                             ak_ptr, &ldak_f, bk_ptr, &ldbk_f, ck_ptr, &ldck_f, dk_ptr, &lddk_f,
                             rcond, &tol_f, iwork, dwork, &ldwork_f_call, bwork, &info);

    // 8. Convert results back to row-major
    if (row_major && info == 0) {
        if (ak_size_elem > 0) { slicot_transpose_to_c_with_ld(ak_cm, ak, n_param, n_param, ldak_f, ldak, sizeof(double)); }
        if (bk_size_elem > 0) { slicot_transpose_to_c_with_ld(bk_cm, bk, n_param, nmeas_param, ldbk_f, ldbk, sizeof(double)); }
        if (ck_size_elem > 0) { slicot_transpose_to_c_with_ld(ck_cm, ck, ncon_param, n_param, ldck_f, ldck, sizeof(double)); }
        if (dk_size_elem > 0) { slicot_transpose_to_c_with_ld(dk_cm, dk, ncon_param, nmeas_param, lddk_f, lddk, sizeof(double)); }
    }

cleanup:
    free(iwork);
    free(dwork);
    free(bwork);

    if (row_major) {
        free(a_cm); free(b_cm); free(c_cm); free(d_cm);
        free(ak_cm); free(bk_cm); free(ck_cm); free(dk_cm);
    }
    
    if (info == SLICOT_MEMORY_ERROR) {
       // Memory allocation error was caught by CHECK_ALLOC
    }
    return info;
}
