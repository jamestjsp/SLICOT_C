/**
 * @file sb10fd.c
 * @brief C wrapper for SLICOT routine SB10FD.
 * @details Computes the matrices of an H-infinity (sub)optimal n-state
 * controller for a continuous-time system.
 * Workspace (IWORK, DWORK, BWORK) is allocated internally by this wrapper.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h> // For malloc, free
#include <string.h> // For memcpy, memset
#include <math.h>   // For fabs, fmax
#include "sb10fd.h" // Public header for this wrapper
#include "slicot_utils.h"  // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions etc.
#include "slicot_f77.h"    // Provides F77_FUNC macro for Fortran name mangling

/* External Fortran routine declaration */
extern void F77_FUNC(sb10fd, SB10FD)(
    const int* n_fortran, const int* m_fortran, const int* np_fortran, const int* ncon_fortran, const int* nmeas_fortran,
    const double* gamma_fortran,
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
int slicot_sb10fd(int n_param, int m_param, int np_param, int ncon_param, int nmeas_param,
                  double gamma_val,
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
    int *bwork = NULL; // For Fortran LOGICAL array
    int liwork = 0;
    int ldwork = 0;
    int lbwork = 0;

    // Pointers for column-major copies
    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
    double *ak_cm = NULL, *bk_cm = NULL, *ck_cm = NULL, *dk_cm = NULL;

    // Fortran-style leading dimensions
    int lda_f, ldb_f, ldc_f, ldd_f;
    int ldak_f, ldbk_f, ldck_f, lddk_f;

    // Make a copy of gamma to pass by reference
    double gamma_f = gamma_val;
    double tol_f = tol;


    // 2. Input parameter validation (Check BEFORE allocating memory)
    if (n_param < 0) { info = -1; goto cleanup; }
    if (m_param < 0) { info = -2; goto cleanup; }
    if (np_param < 0) { info = -3; goto cleanup; }
    if (ncon_param < 0 || ncon_param > m_param) { info = -4; goto cleanup; }
    if (nmeas_param < 0 || nmeas_param > np_param) { info = -5; goto cleanup; }
    if (gamma_val < 0.0) { info = -6; goto cleanup; }

    // Additional constraints from documentation
    // NP-NMEAS >= NCON  => np_param - nmeas_param >= ncon_param
    if (np_param - nmeas_param < ncon_param) { info = -4; goto cleanup; }
    // M-NCON >= NMEAS   => m_param - ncon_param >= nmeas_param
    if (m_param - ncon_param < nmeas_param) { info = -5; goto cleanup; }

    // Check pointers for required arrays
    // For N=0, A,B,C,AK,BK,CK might be NULL if their dimensions involving N become zero.
    // D and DK depend on NP, M, NCON, NMEAS.
    if (a == NULL && n_param > 0) { info = -7; goto cleanup; }
    if (b == NULL && n_param > 0 && m_param > 0) { info = -9; goto cleanup; }
    if (c == NULL && np_param > 0 && n_param > 0) { info = -11; goto cleanup; }
    if (d == NULL && np_param > 0 && m_param > 0) { info = -13; goto cleanup; } // D is NPxM
    if (ak == NULL && n_param > 0) { info = -14; goto cleanup; }
    if (bk == NULL && n_param > 0 && nmeas_param > 0) { info = -16; goto cleanup; }
    if (ck == NULL && ncon_param > 0 && n_param > 0) { info = -18; goto cleanup; }
    if (dk == NULL && ncon_param > 0 && nmeas_param > 0) { info = -20; goto cleanup; } // DK is NCONxNMEAS
    if (rcond == NULL) { info = -21; goto cleanup; }


    // Check leading dimensions
    int min_lda_f = MAX(1, n_param);
    int min_ldb_f = MAX(1, n_param); // B is N x M, Fortran LDB is rows (N)
    int min_ldc_f = MAX(1, np_param); // C is NP x N, Fortran LDC is rows (NP)
    int min_ldd_f = MAX(1, np_param); // D is NP x M, Fortran LDD is rows (NP)
    int min_ldak_f = MAX(1, n_param);
    int min_ldbk_f = MAX(1, n_param); // BK is N x NMEAS, Fortran LDBK is rows (N)
    int min_ldck_f = MAX(1, ncon_param); // CK is NCON x N, Fortran LDCK is rows (NCON)
    int min_lddk_f = MAX(1, ncon_param); // DK is NCON x NMEAS, Fortran LDDK is rows (NCON)

    if (row_major) {
        // C LDA is number of columns
        if (n_param > 0 && lda < n_param) { info = -8; goto cleanup; } // A is NxN
        if (n_param > 0 && m_param > 0 && ldb < m_param) { info = -10; goto cleanup; } // B is NxM
        if (np_param > 0 && n_param > 0 && ldc < n_param) { info = -12; goto cleanup; } // C is NPxN
        if (np_param > 0 && m_param > 0 && ldd < m_param) { info = -13; goto cleanup; } // D is NPxM

        if (n_param > 0 && ldak < n_param) { info = -15; goto cleanup; } // AK is NxN
        if (n_param > 0 && nmeas_param > 0 && ldbk < nmeas_param) { info = -17; goto cleanup; } // BK is NxNMEAS
        if (ncon_param > 0 && n_param > 0 && ldck < n_param) { info = -19; goto cleanup; } // CK is NCONxN
        if (ncon_param > 0 && nmeas_param > 0 && lddk < nmeas_param) { info = -20; goto cleanup; } // DK is NCONxNMEAS
    } else { // Column-major C (Fortran-style LDs)
        if (n_param > 0 && lda < min_lda_f) { info = -8; goto cleanup; }
        if (n_param > 0 && m_param > 0 && ldb < min_ldb_f) { info = -10; goto cleanup; }
        if (np_param > 0 && n_param > 0 && ldc < min_ldc_f) { info = -12; goto cleanup; }
        if (np_param > 0 && m_param > 0 && ldd < min_ldd_f) { info = -13; goto cleanup; }

        if (n_param > 0 && ldak < min_ldak_f) { info = -15; goto cleanup; }
        if (n_param > 0 && nmeas_param > 0 && ldbk < min_ldbk_f) { info = -17; goto cleanup; }
        if (ncon_param > 0 && n_param > 0 && ldck < min_ldck_f) { info = -19; goto cleanup; }
        if (ncon_param > 0 && nmeas_param > 0 && lddk < min_lddk_f) { info = -20; goto cleanup; }
    }
    if (info != 0) { goto cleanup; }

    // 3. Internal Workspace Allocation
    // LIWORK = max(2*max(N,M-NCON,NP-NMEAS,NCON),N*N)
    int m1_ws = m_param - ncon_param;  // M1 in Fortran doc for workspace
    int np1_ws = np_param - nmeas_param; // NP1 in Fortran doc for workspace
    liwork = MAX(1, MAX(2 * MAX(n_param, MAX(m1_ws, MAX(np1_ws, ncon_param))), n_param * n_param));
    iwork = (int*)malloc((size_t)liwork * sizeof(int));
    CHECK_ALLOC(iwork);

    // BWORK size is 2*N
    lbwork = 2 * n_param;
    if (n_param > 0) { // Only allocate if N > 0
        bwork = (int*)malloc((size_t)lbwork * sizeof(int));
        CHECK_ALLOC(bwork);
    } else {
        bwork = NULL; 
    }
    
    // LDWORK calculation based on SLICOT documentation / Slycot wrapper (minimum size formula)
    // Renaming to avoid conflict with Fortran doc names if they were macros
    int N_ws = n_param;
    int M_ws = m_param;
    int NP_ws = np_param;
    int M2_ws = ncon_param;  // NCON
    int NP2_ws = nmeas_param; // NMEAS

    int m1_fmla = M_ws - M2_ws;   // M1 in formula (M - NCON)
    int np1_fmla = NP_ws - NP2_ws; // NP1 in formula (NP - NMEAS)
    
    int d1_fmla = np1_fmla - M2_ws;
    int d2_fmla = m1_fmla - NP2_ws;

    int lw1, lw2, lw3, lw4, lw5, lw6;

    lw1 = (N_ws + np1_fmla + 1)*(N_ws + M2_ws) + MAX(3*(N_ws + M2_ws) + N_ws + np1_fmla, 5*(N_ws + M2_ws));
    lw2 = (N_ws + NP2_ws)*(N_ws + m1_fmla + 1) + MAX(3*(N_ws + NP2_ws) + N_ws + m1_fmla, 5*(N_ws + NP2_ws));
    lw3 = M2_ws + np1_fmla*np1_fmla + MAX(np1_fmla*MAX(N_ws,m1_fmla), MAX(3*M2_ws+np1_fmla,5*M2_ws));
    lw4 = NP2_ws + m1_fmla*m1_fmla + MAX(MAX(N_ws,np1_fmla)*m1_fmla, MAX(3*NP2_ws+m1_fmla,5*NP2_ws));
    
    int term_lw5_1 = M_ws*M_ws + MAX(2*m1_fmla, 3*N_ws*N_ws + MAX(N_ws*M_ws, 10*N_ws*N_ws + 12*N_ws + 5));
    int term_lw5_2 = NP_ws*NP_ws + MAX(2*np1_fmla, 3*N_ws*N_ws + MAX(N_ws*NP_ws, 10*N_ws*N_ws + 12*N_ws + 5));
    lw5 = 2*N_ws*N_ws + N_ws*(M_ws+NP_ws) + MAX(1, MAX(term_lw5_1, term_lw5_2));

    int term_lw6_1 = d1_fmla*d1_fmla + MAX(2*d1_fmla, (d1_fmla+d2_fmla)*NP2_ws);
    int term_lw6_2 = d2_fmla*d2_fmla + MAX(2*d2_fmla, d2_fmla*M2_ws);
    int term_lw6_3 = N_ws*(2*NP2_ws + M2_ws) + MAX(2*N_ws*M2_ws, M2_ws*NP2_ws + MAX(M2_ws*M2_ws+3*M2_ws, NP2_ws*(2*NP2_ws + M2_ws + MAX(NP2_ws,N_ws))));
    lw6 = 2*N_ws*N_ws + N_ws*(M_ws+NP_ws) + MAX(1, M2_ws*NP2_ws + NP2_ws*NP2_ws + M2_ws*M2_ws + MAX(term_lw6_1, MAX(term_lw6_2, MAX(3*N_ws, term_lw6_3))));
    
    ldwork = N_ws*M_ws + NP_ws*(N_ws+M_ws) + M2_ws*M2_ws + NP2_ws*NP2_ws + MAX(1, MAX(lw1, MAX(lw2, MAX(lw3, MAX(lw4, MAX(lw5, lw6))))));
    ldwork = MAX(1, ldwork); // Ensure ldwork is at least 1, especially if all dimensions are 0.


    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);


    // 4. Memory allocation for column-major copies (if row_major)
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
        lda_f = MAX(1, n_param);    // Fortran LDA is rows for A (NxN)
        ldb_f = MAX(1, n_param);    // Fortran LDB is rows for B (NxM)
        ldc_f = MAX(1, np_param);   // Fortran LDC is rows for C (NPxN)
        ldd_f = MAX(1, np_param);   // Fortran LDD is rows for D (NPxM)
        ldak_f = MAX(1, n_param);   // Fortran LDAK is rows for AK (NxN)
        ldbk_f = MAX(1, n_param);   // Fortran LDBK is rows for BK (NxNMEAS)
        ldck_f = MAX(1, ncon_param);// Fortran LDCK is rows for CK (NCONxN)
        lddk_f = MAX(1, ncon_param);// Fortran LDDK is rows for DK (NCONxNMEAS)

        if (a_size_elem > 0) { slicot_transpose_to_fortran_with_ld(a, a_cm, n_param, n_param, lda, lda_f, sizeof(double)); a_ptr = a_cm; } else { a_ptr = NULL; }
        if (b_size_elem > 0) { slicot_transpose_to_fortran_with_ld(b, b_cm, n_param, m_param, ldb, ldb_f, sizeof(double)); b_ptr = b_cm; } else { b_ptr = NULL; }
        if (c_size_elem > 0) { slicot_transpose_to_fortran_with_ld(c, c_cm, np_param, n_param, ldc, ldc_f, sizeof(double)); c_ptr = c_cm; } else { c_ptr = NULL; }
        if (d_size_elem > 0) { slicot_transpose_to_fortran_with_ld(d, d_cm, np_param, m_param, ldd, ldd_f, sizeof(double)); d_ptr = d_cm; } else { d_ptr = NULL; }

        ak_ptr = (ak_size_elem > 0) ? ak_cm : NULL;
        bk_ptr = (bk_size_elem > 0) ? bk_cm : NULL;
        ck_ptr = (ck_size_elem > 0) ? ck_cm : NULL;
        dk_ptr = (dk_size_elem > 0) ? dk_cm : NULL;
    } else { 
        // For column major C, ensure NULL pointers are passed if dimensions are zero
        if (a_size_elem == 0) a_ptr = NULL;
        if (b_size_elem == 0) b_ptr = NULL;
        if (c_size_elem == 0) c_ptr = NULL;
        if (d_size_elem == 0) d_ptr = NULL;
        if (ak_size_elem == 0) ak_ptr = NULL;
        if (bk_size_elem == 0) bk_ptr = NULL;
        if (ck_size_elem == 0) ck_ptr = NULL;
        if (dk_size_elem == 0) dk_ptr = NULL;
    }

    // Prepare int params for Fortran call
    int n_f_call = n_param, m_f_call = m_param, np_f_call = np_param;
    int ncon_f_call = ncon_param, nmeas_f_call = nmeas_param;
    int ldwork_f_call = ldwork;


    // 7. Call Fortran function
    F77_FUNC(sb10fd, SB10FD)(&n_f_call, &m_f_call, &np_f_call, &ncon_f_call, &nmeas_f_call, &gamma_f,
                             a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                             ak_ptr, &ldak_f, bk_ptr, &ldbk_f, ck_ptr, &ldck_f, dk_ptr, &lddk_f,
                             rcond, &tol_f, iwork, dwork, &ldwork_f_call, bwork, &info);

    // 8. Convert results back to row-major (if needed)
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
       // fprintf(stderr, "Error: Memory allocation failed in slicot_sb10fd.\n");
    }
    return info;
}
