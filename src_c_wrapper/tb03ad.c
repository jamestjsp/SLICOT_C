/**
 * @file tb03ad.c
 * @brief C wrapper for SLICOT routine TB03AD.
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>

#include "tb03ad.h"
#include "slicot_utils.h"
#include "slicot_f77.h"

/* External Fortran routine declaration */
extern void F77_FUNC(tb03ad, TB03AD)(
    const char* leri, const char* equil,
    const int* n, const int* m, const int* p,
    double* a, const int* lda,
    double* b, const int* ldb,
    double* c, const int* ldc,
    const double* d, const int* ldd,
    int* nr,
    int* index,
    double* pcoeff, const int* ldpco1, const int* ldpco2,
    double* qcoeff, const int* ldqco1, const int* ldqco2,
    double* vcoeff, const int* ldvco1, const int* ldvco2,
    const double* tol,
    int* iwork, double* dwork, const int* ldwork,
    int* info,
    size_t leri_len, size_t equil_len
);

// Helper function to transpose a 3D matrix block from Fortran column-major to C row-major
static void transpose_3d_f2c(const double* f_arr_cm, double* c_arr_rm,
                             int d1, int d2, int d3,
                             int ld1_f, int ld2_f, /* Fortran leading dims of f_arr_cm */
                             int ld1_c, int ld2_c  /* C logical dims for indexing c_arr_rm (d1,d2,d3 used for iteration) */
                            ) {
    if (!f_arr_cm || !c_arr_rm || d1 == 0 || d2 == 0 || d3 == 0) return;
    for (int k = 0; k < d3; ++k) {     // Outermost loop in Fortran (slowest in memory for slice)
        for (int j = 0; j < d2; ++j) { // Middle loop in Fortran
            for (int i = 0; i < d1; ++i) { // Innermost loop in Fortran (fastest in memory for slice)
                // Fortran index: f_arr_cm[i + j*ld1_f + k*ld1_f*ld2_f]
                // C row-major index for c_arr_rm[i][j][k] where k is fastest:
                // c_arr_rm[i*ld2_c*d3 + j*d3 + k]
                // Assuming ld1_c = d1, ld2_c = d2 for the C output block
                c_arr_rm[i * ld2_c * d3 + j * d3 + k] = f_arr_cm[i + j * ld1_f + k * ld1_f * ld2_f];
            }
        }
    }
}

// Helper function to transpose a 3D matrix block from C row-major to Fortran column-major
static void transpose_3d_c2f(const double* c_arr_rm, double* f_arr_cm,
                             int d1, int d2, int d3,
                             int ld1_c, int ld2_c, /* C logical dims for indexing c_arr_rm */
                             int ld1_f, int ld2_f  /* Fortran leading dims of f_arr_cm */
                            ) {
    if (!c_arr_rm || !f_arr_cm || d1 == 0 || d2 == 0 || d3 == 0) return;
    for (int k = 0; k < d3; ++k) {
        for (int j = 0; j < d2; ++j) {
            for (int i = 0; i < d1; ++i) {
                f_arr_cm[i + j * ld1_f + k * ld1_f * ld2_f] = c_arr_rm[i * ld2_c * d3 + j * d3 + k];
            }
        }
    }
}


SLICOT_EXPORT
int slicot_tb03ad(
    char leri_param, char equil_param,
    int n_param, int m_param, int p_param,
    double* a_io, int lda_c,
    double* b_io, int ldb_c,
    double* c_io, int ldc_c,
    const double* d_in, int ldd_c,
    int* nr_out,
    int* index_out,
    double* pcoeff_out, int ldpco1_c, int ldpco2_c,
    double* qcoeff_out, int ldqco1_c, int ldqco2_c,
    double* vcoeff_out, int ldvco1_c, int ldvco2_c,
    double tol_param,
    int row_major)
{
    int info = 0;
    char leri_f = toupper(leri_param);
    char equil_f = toupper(equil_param);

    // Workspace
    int *iwork = NULL;
    double *dwork = NULL;
    int liwork, ldwork_val;

    // Column-major copies
    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
    double *pcoeff_cm = NULL, *qcoeff_cm = NULL, *vcoeff_cm = NULL;

    // Fortran leading dimensions
    int lda_f, ldb_f, ldc_f, ldd_f;
    int ldpco1_f, ldpco2_f, ldqco1_f, ldqco2_f, ldvco1_f, ldvco2_f;

    // Determine porm, porp based on LERI
    int porm, porp;
    if (leri_f == 'L') {
        porm = p_param;
        porp = m_param;
    } else if (leri_f == 'R') {
        porm = m_param;
        porp = p_param;
    } else {
        info = -1; goto cleanup; // Invalid LERI
    }

    // --- Parameter Validation ---
    if (equil_f != 'S' && equil_f != 'N') { info = -2; goto cleanup; }
    if (n_param < 0) { info = -3; goto cleanup; }
    if (m_param < 0) { info = -4; goto cleanup; }
    if (p_param < 0) { info = -5; goto cleanup; }

    // Validate pointers for essential outputs if dimensions > 0
    if (nr_out == NULL) { info = -11; goto cleanup; } // NR is 11th arg in Fortran list
    if (index_out == NULL && porm > 0) { info = -12; goto cleanup; }
    
    // Max dimensions for workspace arrays B, C, D
    int max_mp = MAX(1, MAX(m_param, p_param));

    // Fortran LDA, LDB, LDC, LDD
    lda_f = MAX(1, n_param);
    ldb_f = MAX(1, n_param); // B is N x max_mp in Fortran
    ldc_f = MAX(1, max_mp);  // C is max_mp x N in Fortran (LDC >= MAX(1,M,P))
    ldd_f = MAX(1, max_mp);  // D is max_mp x max_mp in Fortran (LDD >= MAX(1,M,P))

    // Validate C leading dimensions against Fortran requirements
    if (row_major) {
        if (n_param > 0 && lda_c < n_param) { info = -7; goto cleanup; }
        if (n_param > 0 && m_param > 0 && ldb_c < m_param) { info = -9; goto cleanup; } // ldb_c is cols of B data
        if (p_param > 0 && n_param > 0 && ldc_c < n_param) { info = -11; goto cleanup; } // ldc_c is cols of C data
        if (p_param > 0 && m_param > 0 && ldd_c < m_param) { info = -13; goto cleanup; } // ldd_c is cols of D data
    } else { // Column-major C
        if (n_param > 0 && lda_c < lda_f) { info = -7; goto cleanup; }
        if (n_param > 0 && ldb_c < ldb_f ) { info = -9; goto cleanup; } // ldb_c is rows of B
        if (p_param > 0 && ldc_c < p_param) { info = -11; goto cleanup; } // ldc_c is rows of C, but Fortran LDC needs to be max_mp
        if (ldc_c < ldc_f && n_param > 0) {info = -11; goto cleanup; } // Check against Fortran's LDC for C
        if (p_param > 0 && ldd_c < p_param) { info = -13; goto cleanup; } // ldd_c is rows of D
        if (ldd_c < ldd_f && (p_param > 0 || m_param > 0) ) {info = -13; goto cleanup; } // Check against Fortran's LDD for D
    }
    
    // Fortran leading dimensions for PCOEFF, QCOEFF, VCOEFF
    if (leri_f == 'L') { // porm = P, porp = M
        ldpco1_f = MAX(1, p_param); ldpco2_f = MAX(1, p_param);
        ldqco1_f = MAX(1, p_param); ldqco2_f = MAX(1, m_param);
        ldvco1_f = MAX(1, p_param); ldvco2_f = MAX(1, n_param); // Fortran LDVCO2 is N
    } else { // LERI = 'R', porm = M, porp = P
        ldpco1_f = MAX(1, m_param); ldpco2_f = MAX(1, m_param);
        // As per Fortran doc for TB03AD: LDQCO1/2 for LERI='R' are MAX(1,M,P)
        ldqco1_f = MAX(1, max_mp); ldqco2_f = MAX(1, max_mp);
        ldvco1_f = MAX(1, m_param); ldvco2_f = MAX(1, n_param); // Fortran LDVCO2 is N
    }

    // Validate C leading dimensions for coefficient arrays (user provides flat arrays)
    // User provides ldpco1_c, ldpco2_c etc. which are the actual dimensions of their data block
    if (pcoeff_out == NULL && porm > 0) { info = -15; goto cleanup; }
    if (qcoeff_out == NULL && porm > 0 && porp > 0) { info = -18; goto cleanup; }
    // VCOEFF can be NRx0, so NR can be 0.
    // If NR becomes 0, vcoeff_out might not be strictly needed if ldvco2_c is 0.
    // However, the third dimension is N+1, so if porm > 0, vcoeff_out should be non-NULL.
    if (vcoeff_out == NULL && porm > 0 && (n_param +1) > 0 ) { info = -21; goto cleanup; }


    // --- Workspace Allocation ---
    liwork = n_param + max_mp;
    liwork = MAX(1, liwork);
    iwork = (int*)malloc((size_t)liwork * sizeof(int));
    CHECK_ALLOC(iwork);

    int pm_ws = (leri_f == 'L') ? p_param : m_param;
    ldwork_val = MAX(1, MAX(n_param + MAX(n_param, MAX(3 * m_param, 3 * p_param)), pm_ws * (pm_ws + 2)));
    if (n_param == 0 && m_param == 0 && p_param == 0) ldwork_val = 1; // Ensure min 1 for zero N,M,P

    dwork = (double*)malloc((size_t)ldwork_val * sizeof(double));
    CHECK_ALLOC(dwork);

    // --- Memory for Column-Major Copies ---
    size_t a_nelems = (size_t)lda_f * n_param; if(n_param == 0) a_nelems = 0;
    size_t b_nelems_f = (size_t)ldb_f * max_mp; if(n_param == 0) b_nelems_f = 0; // B is N x max_mp in Fortran
    size_t c_nelems_f = (size_t)ldc_f * n_param; if(p_param == 0 && n_param == 0 && m_param == 0) c_nelems_f = 0; else if (p_param == 0 && n_param > 0) c_nelems_f = 0; // C is max_mp x N
    size_t d_nelems_f = (size_t)ldd_f * max_mp; if(p_param == 0 && m_param == 0) d_nelems_f = 0; // D is max_mp x max_mp

    size_t pcoeff_nelems_f = (size_t)ldpco1_f * ldpco2_f * (n_param + 1); if(porm == 0) pcoeff_nelems_f = 0;
    size_t qcoeff_nelems_f = (size_t)ldqco1_f * ldqco2_f * (n_param + 1); if(porm == 0 || porp == 0) qcoeff_nelems_f = 0;
    size_t vcoeff_nelems_f = (size_t)ldvco1_f * ldvco2_f * (n_param + 1); if(porm == 0) vcoeff_nelems_f = 0; // NR can be 0, ldvco2_f is N_param


    if (row_major) {
        if (a_nelems > 0 && a_io) { a_cm = (double*)malloc(a_nelems * sizeof(double)); CHECK_ALLOC(a_cm); }
        if (b_nelems_f > 0 && b_io) { b_cm = (double*)malloc(b_nelems_f * sizeof(double)); CHECK_ALLOC(b_cm); memset(b_cm, 0, b_nelems_f * sizeof(double));}
        if (c_nelems_f > 0 && c_io) { c_cm = (double*)malloc(c_nelems_f * sizeof(double)); CHECK_ALLOC(c_cm); memset(c_cm, 0, c_nelems_f * sizeof(double));}
        if (d_nelems_f > 0 && d_in) { d_cm = (double*)malloc(d_nelems_f * sizeof(double)); CHECK_ALLOC(d_cm); memset(d_cm, 0, d_nelems_f * sizeof(double));}
        
        if (pcoeff_nelems_f > 0 && pcoeff_out) { pcoeff_cm = (double*)malloc(pcoeff_nelems_f * sizeof(double)); CHECK_ALLOC(pcoeff_cm); }
        if (qcoeff_nelems_f > 0 && qcoeff_out) { qcoeff_cm = (double*)malloc(qcoeff_nelems_f * sizeof(double)); CHECK_ALLOC(qcoeff_cm); }
        if (vcoeff_nelems_f > 0 && vcoeff_out) { vcoeff_cm = (double*)malloc(vcoeff_nelems_f * sizeof(double)); CHECK_ALLOC(vcoeff_cm); }
    }

    // --- Prepare Fortran parameters and perform conversions ---
    double *a_ptr = a_io, *b_ptr = b_io, *c_ptr = c_io;
    const double *d_ptr = d_in;
    double *pcoeff_ptr = pcoeff_out, *qcoeff_ptr = qcoeff_out, *vcoeff_ptr = vcoeff_out;

    if (row_major) {
        // Inputs
        if (a_cm && a_io && n_param > 0) { slicot_transpose_to_fortran_with_ld(a_io, a_cm, n_param, n_param, lda_c, lda_f, sizeof(double)); }
        a_ptr = a_cm ? a_cm : NULL; // Use cm if allocated, else NULL (e.g. N=0)

        if (b_cm && b_io && n_param > 0 && m_param > 0) { // B is N x M from user
            slicot_transpose_to_fortran_with_ld(b_io, b_cm, n_param, m_param, ldb_c, ldb_f, sizeof(double));
        }
        b_ptr = b_cm ? b_cm : NULL;

        if (c_cm && c_io && p_param > 0 && n_param > 0) { // C is P x N from user
             slicot_transpose_to_fortran_with_ld(c_io, c_cm, p_param, n_param, ldc_c, ldc_f, sizeof(double));
        }
        c_ptr = c_cm ? c_cm : NULL;
        
        if (d_cm && d_in && p_param > 0 && m_param > 0) { // D is P x M from user
            slicot_transpose_to_fortran_with_ld(d_in, d_cm, p_param, m_param, ldd_c, ldd_f, sizeof(double));
        }
        d_ptr = d_cm ? d_cm : (double*)d_in; // d_in is const, cast if d_cm not used. If d_cm is NULL, d_ptr should be NULL if d_in is also null.

        // Outputs (Fortran will write to these _cm buffers)
        pcoeff_ptr = pcoeff_cm ? pcoeff_cm : NULL;
        qcoeff_ptr = qcoeff_cm ? qcoeff_cm : NULL;
        vcoeff_ptr = vcoeff_cm ? vcoeff_cm : NULL;

    } else { // Column-major C
        if (n_param == 0) { a_ptr = NULL; } // Fortran expects NULL for zero-dim arrays
        if (n_param == 0 || m_param == 0) { b_ptr = NULL; } // If B is effectively zero size for Fortran
        if (p_param == 0 || n_param == 0) { c_ptr = NULL; } // If C is effectively zero size
        if (p_param == 0 || m_param == 0) { d_ptr = NULL; }
        if (porm == 0) {pcoeff_ptr = NULL; qcoeff_ptr = NULL; vcoeff_ptr = NULL;}
        else if (porp == 0 && qcoeff_ptr) {qcoeff_ptr = NULL;}


        // For column major, ensure Fortran LDCs are met by user's C LDCs
        // lda_c, ldb_c, ldc_c, ldd_c are already Fortran-style (rows)
        // This was partially checked in validation.
        // For C, LDC is P, Fortran LDC is max_mp. If ldc_c < ldc_f, it's an issue.
    }
    
    int nr_val;

    // --- Call Fortran function ---
    F77_FUNC(tb03ad, TB03AD)(&leri_f, &equil_f, &n_param, &m_param, &p_param,
                             a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                             &nr_val, index_out,
                             pcoeff_ptr, &ldpco1_f, &ldpco2_f,
                             qcoeff_ptr, &ldqco1_f, &ldqco2_f,
                             vcoeff_ptr, &ldvco1_f, &ldvco2_f,
                             &tol_param, iwork, dwork, &ldwork_val, &info,
                             1, 1);

    *nr_out = nr_val;
    int current_nr = nr_val;
    int kpcoef_actual = 0;
    if (info == 0 && porm > 0) {
        for(int i=0; i<porm; ++i) {
            if (index_out[i] > kpcoef_actual) kpcoef_actual = index_out[i];
        }
        kpcoef_actual += 1;
    }
    if (porm == 0) kpcoef_actual = 0; // If porm is 0, index_out is not accessed.
    kpcoef_actual = MAX(0, kpcoef_actual); // Ensure non-negative


    // --- Convert results back to row-major if needed ---
    if (row_major && info == 0) {
        // A_io (NR x NR)
        if (a_cm && a_io && current_nr > 0) {
            slicot_transpose_to_c_with_ld(a_cm, a_io, current_nr, current_nr, lda_f, lda_c, sizeof(double));
        }
        // B_io (NR x M)
        if (b_cm && b_io && current_nr > 0 && m_param > 0) {
            slicot_transpose_to_c_with_ld(b_cm, b_io, current_nr, m_param, ldb_f, ldb_c, sizeof(double));
        }
        // C_io (P x NR)
        if (c_cm && c_io && p_param > 0 && current_nr > 0) {
            slicot_transpose_to_c_with_ld(c_cm, c_io, p_param, current_nr, ldc_f, ldc_c, sizeof(double));
        }

        // PCOEFF_out (porm x porm x kpcoef_actual)
        if (pcoeff_cm && pcoeff_out && porm > 0 && kpcoef_actual > 0) {
            transpose_3d_f2c(pcoeff_cm, pcoeff_out, porm, porm, kpcoef_actual, ldpco1_f, ldpco2_f, ldpco1_c, ldpco2_c);
        }

        // QCOEFF_out
        int q_d1, q_d2;
        if (leri_f == 'L') { q_d1 = p_param; q_d2 = m_param; } // P x M
        else { q_d1 = p_param; q_d2 = m_param; } // P x M (actual data block)
                                                 // Fortran QCOEFF(LDQCO1,LDQCO2,*)
                                                 // LERI='R': Q is P x M, LDQCO1_f=max_mp, LDQCO2_f=max_mp
                                                 // LERI='L': Q is P x M, LDQCO1_f=P, LDQCO2_f=M

        if (qcoeff_cm && qcoeff_out && q_d1 > 0 && q_d2 > 0 && kpcoef_actual > 0) {
             transpose_3d_f2c(qcoeff_cm, qcoeff_out, q_d1, q_d2, kpcoef_actual, ldqco1_f, ldqco2_f, ldqco1_c, ldqco2_c);
        }
        
        // VCOEFF_out (porm x NR x kpcoef_actual)
        // ldvco2_c is user's N, but actual data is NR columns
        if (vcoeff_cm && vcoeff_out && porm > 0 && current_nr > 0 && kpcoef_actual > 0) {
            transpose_3d_f2c(vcoeff_cm, vcoeff_out, porm, current_nr, kpcoef_actual, ldvco1_f, ldvco2_f, ldvco1_c, ldvco2_c);
        }
    }


cleanup:
    free(iwork);
    free(dwork);
    if (row_major) {
        free(a_cm); free(b_cm); free(c_cm); free(d_cm);
        free(pcoeff_cm); free(qcoeff_cm); free(vcoeff_cm);
    }
    
    if (info == SLICOT_MEMORY_ERROR && n_param > -1) { // Check n_param to avoid message during intentional error tests
       fprintf(stderr, "Error: Memory allocation failed in slicot_tb03ad.\n");
    }
    return info;
}
