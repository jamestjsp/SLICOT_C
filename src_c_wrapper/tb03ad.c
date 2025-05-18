/**
 * @file tb03ad.c
 * @brief C wrapper for SLICOT routine TB03AD.
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h> // Added for fprintf and stderr

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

// Helper function to transpose a 3D matrix block from Fortran column-major (read from f_arr_cm)
// to C row-major (written to c_arr_rm).
// d1, d2, d3 are the logical dimensions of the data block.
// ld1_f, ld2_f are the leading dimensions of the Fortran array f_arr_cm.
// ld1_c_data, ld2_c_data are the logical dimensions for the C array c_arr_rm (d1, d2 used for iteration, d3 is fastest varying).
static void transpose_3d_f2c_slice(const double* f_arr_cm, double* c_arr_rm,
                                   int d1, int d2, int d3, // Logical dimensions of the data block
                                   int ld1_f, int ld2_f,   // Leading dimensions of Fortran source
                                   int d1_c_buf, int d2_c_buf // Dimensions of C target buffer for indexing
                                  ) {
    if (!f_arr_cm || !c_arr_rm || d1 == 0 || d2 == 0 || d3 == 0) return;
    for (int k = 0; k < d3; ++k) {     // Loop over the third dimension (coefficients)
        for (int j = 0; j < d2; ++j) { // Loop over the second logical dimension
            for (int i = 0; i < d1; ++i) { // Loop over the first logical dimension
                // Fortran memory access: element (i,j,k) of the logical block
                // f_arr_cm is 0-indexed pointer to start of Fortran array
                // Fortran index: A(i+1, j+1, k+1) -> memory: val_ptr[i + j*ld1_f + k*ld1_f*ld2_f]
                double val = f_arr_cm[i + j * ld1_f + k * (size_t)ld1_f * ld2_f];
                
                // C row-major memory access: element (i,j,k) of the logical block
                // c_arr_rm is 0-indexed pointer to start of C array
                // C index: C[i][j][k] -> memory: c_ptr[i*d2_c_buf*d3 + j*d3 + k]
                c_arr_rm[i * d2_c_buf * d3 + j * d3 + k] = val;
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
    double* pcoeff_out, int ldpco1_c_in, int ldpco2_c_in,
    double* qcoeff_out, int ldqco1_c_in, int ldqco2_c_in,
    double* vcoeff_out, int ldvco1_c_in, int ldvco2_c_in,
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

    // Column-major copies for 2D matrices
    double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
    // Column-major buffers for 3D coefficient arrays (used only if row_major is true)
    double *pcoeff_f_buf = NULL, *qcoeff_f_buf = NULL, *vcoeff_f_buf = NULL;

    // Fortran leading dimensions
    int lda_f, ldb_f, ldc_f, ldd_f;
    int ldpco1_f, ldpco2_f, ldqco1_f, ldqco2_f, ldvco1_f, ldvco2_f;

    int porm, porp; // Actual dimensions of P(s) and Q(s) data blocks
    if (leri_f == 'L') {
        porm = p_param; porp = m_param;
    } else if (leri_f == 'R') {
        porm = m_param; porp = p_param;
    } else {
        info = -1; goto cleanup; 
    }

    // --- Parameter Validation ---
    if (equil_f != 'S' && equil_f != 'N') { info = -2; goto cleanup; }
    if (n_param < 0) { info = -3; goto cleanup; }
    if (m_param < 0) { info = -4; goto cleanup; }
    if (p_param < 0) { info = -5; goto cleanup; }
    if (nr_out == NULL) { info = -11; goto cleanup; }
    if (index_out == NULL && porm > 0) { info = -12; goto cleanup; }
    if (pcoeff_out == NULL && porm > 0) { info = -15; goto cleanup; }
    if (qcoeff_out == NULL && porm > 0 && porp > 0) { info = -18; goto cleanup; }
    if (vcoeff_out == NULL && porm > 0 && (n_param +1) > 0 ) { info = -21; goto cleanup; }


    int max_mp = MAX(1, MAX(m_param, p_param));

    // Fortran leading dimensions for 2D matrices A, B, C, D
    lda_f = MAX(1, n_param);
    ldb_f = MAX(1, n_param); 
    ldc_f = MAX(1, max_mp); 
    ldd_f = MAX(1, max_mp); 

    // Validate C leading dimensions for 2D matrices
    if (row_major) {
        if (n_param > 0 && lda_c < n_param) { info = -7; goto cleanup; }
        if (n_param > 0 && m_param > 0 && ldb_c < m_param) { info = -9; goto cleanup; }
        if (p_param > 0 && n_param > 0 && ldc_c < n_param) { info = -11; goto cleanup; }
        if (p_param > 0 && m_param > 0 && ldd_c < m_param) { info = -13; goto cleanup; }
    } else { // Column-major C
        if (n_param > 0 && lda_c < lda_f) { info = -7; goto cleanup; }
        if (n_param > 0 && ldb_c < ldb_f) { info = -9; goto cleanup; }
        if (n_param > 0 && ldc_c < ldc_f) { info = -11; goto cleanup; } // Check against Fortran's LDC
        if ((p_param > 0 || m_param > 0) && ldd_c < ldd_f) { info = -13; goto cleanup; }
    }
    
    // Fortran leading dimensions for 3D coefficient arrays
    // These are the dimensions Fortran's TB03AD routine expects for its array arguments.
    if (leri_f == 'L') { 
        ldpco1_f = MAX(1, p_param); ldpco2_f = MAX(1, p_param);
        ldqco1_f = MAX(1, p_param); ldqco2_f = MAX(1, m_param);
        ldvco1_f = MAX(1, p_param); ldvco2_f = MAX(1, n_param); 
    } else { // LERI = 'R'
        ldpco1_f = MAX(1, m_param); ldpco2_f = MAX(1, m_param);
        ldqco1_f = MAX(1, max_mp);  ldqco2_f = MAX(1, max_mp); // QCOEFF uses max_mp for LERI='R'
        ldvco1_f = MAX(1, m_param); ldvco2_f = MAX(1, n_param); 
    }

    // --- Workspace Allocation ---
    liwork = n_param + max_mp;
    liwork = MAX(1, liwork);
    iwork = (int*)malloc((size_t)liwork * sizeof(int));
    CHECK_ALLOC(iwork);

    int pm_ws = (leri_f == 'L') ? p_param : m_param;
    ldwork_val = MAX(1, MAX(n_param + MAX(n_param, MAX(3 * m_param, 3 * p_param)), pm_ws * (pm_ws + 2)));
    if (n_param == 0 && m_param == 0 && p_param == 0) ldwork_val = 1; 
    dwork = (double*)malloc((size_t)ldwork_val * sizeof(double));
    CHECK_ALLOC(dwork);

    // --- Memory for Column-Major Copies (if row_major) ---
    // Sizes for 2D matrix copies
    size_t a_nelems_f = (size_t)lda_f * n_param; if(n_param == 0) a_nelems_f = 0;
    size_t b_nelems_f = (size_t)ldb_f * max_mp;  if(n_param == 0) b_nelems_f = 0;
    size_t c_nelems_f = (size_t)ldc_f * n_param; if(p_param == 0 && n_param == 0) c_nelems_f = 0; else if (p_param == 0 && n_param > 0) c_nelems_f = 0;
    size_t d_nelems_f = (size_t)ldd_f * max_mp;  if(p_param == 0 && m_param == 0) d_nelems_f = 0;

    // Sizes for 3D coefficient array Fortran buffers (used if row_major)
    size_t pcoeff_f_nelems = (size_t)ldpco1_f * ldpco2_f * (n_param + 1); if(porm == 0) pcoeff_f_nelems = 0;
    size_t qcoeff_f_nelems = (size_t)ldqco1_f * ldqco2_f * (n_param + 1); if(porm == 0 && porp == 0 && !(leri_f=='R' && (m_param>0 || p_param>0))) qcoeff_f_nelems = 0; else if (porm == 0 || porp == 0) qcoeff_f_nelems = 0;
    size_t vcoeff_f_nelems = (size_t)ldvco1_f * ldvco2_f * (n_param + 1); if(porm == 0 || n_param == 0) vcoeff_f_nelems = 0;


    if (row_major) {
        if (a_nelems_f > 0 && a_io) { a_cm = (double*)malloc(a_nelems_f * sizeof(double)); CHECK_ALLOC(a_cm); }
        if (b_nelems_f > 0 && b_io) { b_cm = (double*)malloc(b_nelems_f * sizeof(double)); CHECK_ALLOC(b_cm); memset(b_cm, 0, b_nelems_f * sizeof(double));}
        if (c_nelems_f > 0 && c_io) { c_cm = (double*)malloc(c_nelems_f * sizeof(double)); CHECK_ALLOC(c_cm); memset(c_cm, 0, c_nelems_f * sizeof(double));}
        if (d_nelems_f > 0 && d_in) { d_cm = (double*)malloc(d_nelems_f * sizeof(double)); CHECK_ALLOC(d_cm); memset(d_cm, 0, d_nelems_f * sizeof(double));}
        
        if (pcoeff_f_nelems > 0 && pcoeff_out) { pcoeff_f_buf = (double*)malloc(pcoeff_f_nelems * sizeof(double)); CHECK_ALLOC(pcoeff_f_buf); }
        if (qcoeff_f_nelems > 0 && qcoeff_out) { qcoeff_f_buf = (double*)malloc(qcoeff_f_nelems * sizeof(double)); CHECK_ALLOC(qcoeff_f_buf); }
        if (vcoeff_f_nelems > 0 && vcoeff_out) { vcoeff_f_buf = (double*)malloc(vcoeff_f_nelems * sizeof(double)); CHECK_ALLOC(vcoeff_f_buf); }
    }

    // --- Prepare Fortran parameters and perform conversions ---
    double *a_ptr = a_io, *b_ptr = b_io, *c_ptr = c_io;
    const double *d_ptr = d_in;
    double *pcoeff_ptr = pcoeff_out, *qcoeff_ptr = qcoeff_out, *vcoeff_ptr = vcoeff_out;

    if (row_major) {
        if (a_cm && a_io && n_param > 0) { slicot_transpose_to_fortran_with_ld(a_io, a_cm, n_param, n_param, lda_c, lda_f, sizeof(double)); }
        a_ptr = a_cm ? a_cm : (n_param == 0 ? NULL : a_io);

        if (b_cm && b_io && n_param > 0 && m_param > 0) { 
            slicot_transpose_to_fortran_with_ld(b_io, b_cm, n_param, m_param, ldb_c, ldb_f, sizeof(double));
        }
        b_ptr = b_cm ? b_cm : (n_param == 0 || m_param == 0 ? NULL : b_io);

        if (c_cm && c_io && p_param > 0 && n_param > 0) { 
             slicot_transpose_to_fortran_with_ld(c_io, c_cm, p_param, n_param, ldc_c, ldc_f, sizeof(double));
        }
        c_ptr = c_cm ? c_cm : (p_param == 0 || n_param == 0 ? NULL : c_io);
        
        if (d_cm && d_in && p_param > 0 && m_param > 0) { 
            slicot_transpose_to_fortran_with_ld(d_in, d_cm, p_param, m_param, ldd_c, ldd_f, sizeof(double));
        }
        d_ptr = d_cm ? d_cm : (p_param == 0 || m_param == 0 ? NULL : d_in);

        pcoeff_ptr = pcoeff_f_buf ? pcoeff_f_buf : NULL;
        qcoeff_ptr = qcoeff_f_buf ? qcoeff_f_buf : NULL;
        vcoeff_ptr = vcoeff_f_buf ? vcoeff_f_buf : NULL;
    } else { // Column-major C: Pass pointers directly
        if (n_param == 0) a_ptr = NULL;
        if (n_param == 0 || m_param == 0) b_ptr = NULL;
        if (p_param == 0 || n_param == 0) c_ptr = NULL;
        if (p_param == 0 || m_param == 0) d_ptr = NULL;
        if (porm == 0) { pcoeff_ptr = NULL; qcoeff_ptr = NULL; vcoeff_ptr = NULL;}
        else if (porp == 0 && qcoeff_ptr != NULL && porm !=0 ) { /* QCOEFF is porm x 0 */ } // allow qcoeff_ptr if porm >0
        else if (porp == 0 && porm > 0) { /* qcoeff_ptr might be non-NULL but data is 0-width */ }

        // For coeff arrays, if row_major=0, pcoeff_out etc. are already in column-major.
        // The Fortran routine will use ldpco1_f, ldpco2_f etc.
        // The C caller must ensure their buffers (pcoeff_out etc.) are large enough
        // for Fortran's access pattern using ldpco1_f, ldpco2_f.
        // ldpco1_c_in, ldpco2_c_in are not used by this wrapper when row_major=0 for coeff arrays.
    }
    
    int nr_val;

    F77_FUNC(tb03ad, TB03AD)(&leri_f, &equil_f, &n_param, &m_param, &p_param,
                             a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                             &nr_val, index_out,
                             pcoeff_ptr, &ldpco1_f, &ldpco2_f,
                             qcoeff_ptr, &ldqco1_f, &ldqco2_f,
                             vcoeff_ptr, &ldvco1_f, &ldvco2_f, // Fortran expects LDVCO2 = N
                             &tol_param, iwork, dwork, &ldwork_val, &info,
                             1, 1);

    *nr_out = nr_val;
    int current_nr = nr_val;
    int kpcoef_actual = 0;
    if (info == 0 && porm > 0 && index_out) { // Check index_out
        for(int i=0; i<porm; ++i) {
            if (index_out[i] > kpcoef_actual) kpcoef_actual = index_out[i];
        }
        kpcoef_actual += 1;
    }
    if (porm == 0) kpcoef_actual = 0;
    kpcoef_actual = MAX(0, kpcoef_actual);


    if (row_major && info == 0) {
        if (a_cm && a_io && current_nr > 0) {
            slicot_transpose_to_c_with_ld(a_cm, a_io, current_nr, current_nr, lda_f, lda_c, sizeof(double));
        }
        if (b_cm && b_io && current_nr > 0 && m_param > 0) {
            slicot_transpose_to_c_with_ld(b_cm, b_io, current_nr, m_param, ldb_f, ldb_c, sizeof(double));
        }
        if (c_cm && c_io && p_param > 0 && current_nr > 0) {
            slicot_transpose_to_c_with_ld(c_cm, c_io, p_param, current_nr, ldc_f, ldc_c, sizeof(double));
        }

        // Transpose 3D coefficient arrays back
        // PCOEFF: data porm x porm x kpcoef_actual
        if (pcoeff_f_buf && pcoeff_out && porm > 0 && kpcoef_actual > 0) {
            transpose_3d_f2c_slice(pcoeff_f_buf, pcoeff_out, porm, porm, kpcoef_actual, ldpco1_f, ldpco2_f, ldpco1_c_in, ldpco2_c_in);
        }

        // QCOEFF: data (q_d1 x q_d2 x kpcoef_actual)
        int q_d1_data = (leri_f == 'L') ? p_param : p_param; // Q is PxM
        int q_d2_data = (leri_f == 'L') ? m_param : m_param;
        if (qcoeff_f_buf && qcoeff_out && q_d1_data > 0 && q_d2_data > 0 && kpcoef_actual > 0) {
             transpose_3d_f2c_slice(qcoeff_f_buf, qcoeff_out, q_d1_data, q_d2_data, kpcoef_actual, ldqco1_f, ldqco2_f, ldqco1_c_in, ldqco2_c_in);
        }
        
        // VCOEFF: data porm x current_nr x kpcoef_actual
        if (vcoeff_f_buf && vcoeff_out && porm > 0 && current_nr > 0 && kpcoef_actual > 0) {
            transpose_3d_f2c_slice(vcoeff_f_buf, vcoeff_out, porm, current_nr, kpcoef_actual, ldvco1_f, ldvco2_f, ldvco1_c_in, ldvco2_c_in);
        }
    }

cleanup:
    free(iwork);
    free(dwork);
    if (row_major) {
        free(a_cm); free(b_cm); free(c_cm); free(d_cm);
        free(pcoeff_f_buf); free(qcoeff_f_buf); free(vcoeff_f_buf);
    }
    
    if (info == SLICOT_MEMORY_ERROR && n_param >= 0) { 
       fprintf(stderr, "Error: Memory allocation failed in slicot_tb03ad.\n");
    }
    return info;
}
