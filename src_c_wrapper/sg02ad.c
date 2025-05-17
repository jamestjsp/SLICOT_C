/**
 * @file sg02ad.c
 * @brief C wrapper for SLICOT routine SG02AD.
 * @details Solves continuous- or discrete-time algebraic Riccati equations
 * for descriptor systems.
 * Workspace is allocated internally. Input/output matrix format is handled
 * via the row_major parameter.
 */

#include <stdlib.h> // For malloc, free
#include <string.h> // For memcpy, toupper
#include <ctype.h>  // For toupper
#include <math.h>   // For fmax (from C standard library)

#include "sg02ad.h"       // Public header for this wrapper
#include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions etc.
#include "slicot_f77.h"   // Provides F77_FUNC macro

/* External Fortran routine declaration */
extern void F77_FUNC(sg02ad, SG02AD)(
    const char* dico, const char* jobb, const char* fact, const char* uplo,
    const char* jobl, const char* scal, const char* sort, const char* acc,
    const int* n, const int* m, const int* p,
    const double* a, const int* lda,
    const double* e, const int* lde,
    const double* b, const int* ldb, /* If JOBB='G', B is G */
    const double* q, const int* ldq, /* If FACT='C'or'B', Q is C_factor */
    const double* r, const int* ldr, /* If FACT='D'or'B', R is D_factor */
    const double* l, const int* ldl,
    double* rcondu,
    double* x, const int* ldx,
    double* alfar, double* alfai, double* beta,
    double* s, const int* lds,
    double* t, const int* ldt,
    double* u, const int* ldu,
    const double* tol,
    int* iwork, double* dwork, const int* ldwork,
    int* bwork, /* Fortran LOGICAL */
    int* iwarn, int* info,
    size_t dico_len, size_t jobb_len, size_t fact_len, size_t uplo_len,
    size_t jobl_len, size_t scal_len, size_t sort_len, size_t acc_len
);

SLICOT_EXPORT
int slicot_sg02ad(
    char dico_param, char jobb_param, char fact_param, char uplo_param,
    char jobl_param, char scal_param, char sort_param, char acc_param,
    int n_param, int m_param, int p_param,
    double* a_io, int lda,   /* A is input */
    double* e_io, int lde,   /* E is input */
    double* b_io, int ldb,   /* B (or G) is input */
    double* q_io, int ldq,   /* Q (or C_factor) is input */
    double* r_io, int ldr,   /* R (or D_factor) is input */
    double* l_io, int ldl,   /* L is input */
    double* rcondu,       /* output */
    double* x_out, int ldx,  /* X is output */
    double* alfar_out, double* alfai_out, double* beta_out, /* output */
    double* s_out, int lds,  /* S is output */
    double* t_out, int ldt,  /* T is output */
    double* u_out, int ldu,  /* U is output */
    double tol_param,
    int* iwarn_out,          /* output */
    int row_major)
{
    // 1. Variable declarations
    int info = 0;
    int local_iwarn = 0;

    // Workspace arrays
    int *iwork = NULL;
    double *dwork = NULL;
    int *bwork_c = NULL; // C int array for Fortran LOGICAL

    // Workspace sizes
    int liwork = 0;
    int ldwork = 0;
    int lbwork = 0; // For BWORK (2*N elements)

    // Column-major copies for input matrices (if row_major)
    double *a_cm = NULL, *e_cm = NULL, *b_cm = NULL, *q_cm = NULL, *r_cm = NULL, *l_cm = NULL;
    // Column-major copies for output matrices (if row_major)
    double *x_cm = NULL, *s_cm = NULL, *t_cm = NULL, *u_cm = NULL;

    // Fortran-style leading dimensions
    int lda_f, lde_f, ldb_f, ldq_f, ldr_f, ldl_f;
    int ldx_f, lds_f, ldt_f, ldu_f;
    
    double tol_f = tol_param; // Pass tol by reference

    // Character parameters for Fortran
    char dico_f = toupper(dico_param);
    char jobb_f = toupper(jobb_param);
    char fact_f = toupper(fact_param);
    char uplo_f = toupper(uplo_param);
    char jobl_f = toupper(jobl_param);
    char scal_f = toupper(scal_param);
    char sort_f = toupper(sort_param);
    char acc_f  = toupper(acc_param);

    // 2. Input parameter validation
    if (strchr("CD", dico_f) == NULL) { info = -1; goto cleanup; }
    if (strchr("BG", jobb_f) == NULL) { info = -2; goto cleanup; }
    if (strchr("NCDB", fact_f) == NULL) { info = -3; goto cleanup; }
    if (strchr("UL", uplo_f) == NULL) { info = -4; goto cleanup; }
    if (strchr("ZN", jobl_f) == NULL) { info = -5; goto cleanup; }
    if (strchr("GN", scal_f) == NULL) { info = -6; goto cleanup; }
    if (strchr("SU", sort_f) == NULL) { info = -7; goto cleanup; }
    if (strchr("RN", acc_f) == NULL) { info = -8; goto cleanup; }

    if (n_param < 0) { info = -9; goto cleanup; }
    if (jobb_f == 'B' && m_param < 0) { info = -10; goto cleanup; }
    if (strchr("CDB", fact_f) != NULL && p_param < 0) { info = -11; goto cleanup; }


    // Validate pointers for required arrays
    if (a_io == NULL && n_param > 0) { info = -12; goto cleanup; }
    if (e_io == NULL && n_param > 0) { info = -14; goto cleanup; }
    if (b_io == NULL && n_param > 0 && ( (jobb_f == 'B' && m_param > 0) || (jobb_f == 'G') ) ) { info = -16; goto cleanup; }
    if (q_io == NULL && n_param > 0 && (strchr("NCDB", fact_f) != NULL) ) { info = -18; goto cleanup; } // Q or C_factor
    if (jobb_f == 'B' && r_io == NULL && m_param > 0 && (strchr("NCDB", fact_f) != NULL) ) { info = -20; goto cleanup; } // R or D_factor
    if (jobb_f == 'B' && jobl_f == 'N' && l_io == NULL && n_param > 0 && m_param > 0) { info = -22; goto cleanup; }
    
    if (rcondu == NULL) { info = -23; goto cleanup; }
    if (x_out == NULL && n_param > 0) { info = -24; goto cleanup; }
    if (alfar_out == NULL && n_param > 0) { info = -26; goto cleanup; } // ALFAR, ALFAI, BETA are 2*N
    if (alfai_out == NULL && n_param > 0) { info = -27; goto cleanup; }
    if (beta_out == NULL && n_param > 0) { info = -28; goto cleanup; }
    if (s_out == NULL && n_param > 0) { info = -29; goto cleanup; }
    if (t_out == NULL && n_param > 0) { info = -31; goto cleanup; }
    if (u_out == NULL && n_param > 0) { info = -33; goto cleanup; }
    if (iwarn_out == NULL) { info = -40; goto cleanup; } // Assuming IWARN is arg 40 based on typical position

    // Validate leading dimensions
    int min_lda_f = MAX(1, n_param);
    int min_lde_f = MAX(1, n_param);
    int min_ldb_f = MAX(1, n_param); // For B (NxM or G NxN)
    int min_ldq_f = (fact_f == 'N' || fact_f == 'D') ? MAX(1,n_param) : MAX(1,p_param);
    int min_ldr_f = (jobb_f == 'B' && (fact_f == 'N' || fact_f == 'C')) ? MAX(1,m_param) : ((jobb_f == 'B' && (fact_f == 'D' || fact_f == 'B')) ? MAX(1,p_param) : 1);
    int min_ldl_f = (jobb_f == 'B' && jobl_f == 'N') ? MAX(1,n_param) : 1;
    int min_ldx_f = MAX(1, n_param);
    int min_lds_f = (jobb_f == 'B') ? MAX(1, 2*n_param + m_param) : MAX(1, 2*n_param);
    int min_ldt_f = (jobb_f == 'B') ? MAX(1, 2*n_param + m_param) : MAX(1, 2*n_param);
    int min_ldu_f = MAX(1, 2*n_param);

    // Effective dimensions for row-major C arrays (number of columns)
    int lda_c_cols = n_param;
    int lde_c_cols = n_param;
    int ldb_c_cols = (jobb_f == 'B') ? m_param : n_param; // B is NxM, G is NxN
    int ldq_c_cols = (fact_f == 'N' || fact_f == 'D') ? n_param : n_param; // Q is NxN, C_factor is PxN
    int ldr_c_cols = (jobb_f == 'B' && (fact_f == 'N' || fact_f == 'C')) ? m_param : ((jobb_f == 'B' && (fact_f == 'D' || fact_f == 'B')) ? m_param : 0); // R is MxM, D_factor is PxM
    int ldl_c_cols = (jobb_f == 'B' && jobl_f == 'N') ? m_param : 0; // L is NxM
    int ldx_c_cols = n_param;
    int lds_c_cols = (jobb_f == 'B') ? (2*n_param + m_param) : (2*n_param);
    int ldt_c_cols = 2*n_param;
    int ldu_c_cols = 2*n_param;


    if (row_major) {
        if (n_param > 0 && lda < lda_c_cols) { info = -13; goto cleanup; }
        if (n_param > 0 && lde < lde_c_cols) { info = -15; goto cleanup; }
        if (n_param > 0 && ( (jobb_f == 'B' && m_param > 0) || jobb_f == 'G' ) && ldb < ldb_c_cols) { info = -17; goto cleanup; }
        if (n_param > 0 && ldq < ldq_c_cols ) { info = -19; goto cleanup; }
        if (jobb_f == 'B' && m_param > 0 && ldr < ldr_c_cols && ldr_c_cols > 0 ) { info = -21; goto cleanup; }
        if (jobb_f == 'B' && jobl_f == 'N' && n_param > 0 && m_param > 0 && ldl < ldl_c_cols ) { info = -22; goto cleanup; }
        if (n_param > 0 && ldx < ldx_c_cols) { info = -25; goto cleanup; }
        if (n_param > 0 && lds < lds_c_cols) { info = -30; goto cleanup; }
        if (n_param > 0 && ldt < ldt_c_cols) { info = -32; goto cleanup; }
        if (n_param > 0 && ldu < ldu_c_cols) { info = -34; goto cleanup; }
    } else { // Column-major C
        if (n_param > 0 && lda < min_lda_f) { info = -13; goto cleanup; }
        if (n_param > 0 && lde < min_lde_f) { info = -15; goto cleanup; }
        if (n_param > 0 && ( (jobb_f == 'B' && m_param > 0) || jobb_f == 'G' ) && ldb < min_ldb_f) { info = -17; goto cleanup; }
        if (n_param > 0 && ldq < min_ldq_f ) { info = -19; goto cleanup; }
        if (jobb_f == 'B' && m_param > 0 && ldr < min_ldr_f && min_ldr_f > 0) { info = -21; goto cleanup; }
        if (jobb_f == 'B' && jobl_f == 'N' && n_param > 0 && m_param > 0 && ldl < min_ldl_f) { info = -22; goto cleanup; }
        if (n_param > 0 && ldx < min_ldx_f) { info = -25; goto cleanup; }
        if (n_param > 0 && lds < min_lds_f) { info = -30; goto cleanup; }
        if (n_param > 0 && ldt < min_ldt_f) { info = -32; goto cleanup; }
        if (n_param > 0 && ldu < min_ldu_f) { info = -34; goto cleanup; }
    }
    if (info != 0) { goto cleanup; }

    // 3. Workspace Allocation
    if (jobb_f == 'B') {
        liwork = MAX(1, MAX(m_param, 2 * n_param));
    } else { // JOBB = 'G'
        liwork = MAX(1, 2 * n_param);
    }
    iwork = (int*)malloc((size_t)liwork * sizeof(int));
    CHECK_ALLOC(iwork);

    lbwork = 2 * n_param;
    if (lbwork > 0) {
        bwork_c = (int*)malloc((size_t)lbwork * sizeof(int));
        CHECK_ALLOC(bwork_c);
    } else {
        bwork_c = NULL; // Pass NULL if N=0
    }


    if (jobb_f == 'G') {
        ldwork = MAX(1, MAX(7 * (2 * n_param + 1) + 16, 16 * n_param));
    } else { // JOBB = 'B'
        ldwork = MAX(1, MAX(7 * (2 * n_param + 1) + 16, MAX(16 * n_param, MAX(2 * n_param + m_param, 3 * m_param))));
    }
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);

    // 4. Memory for column-major copies
    size_t a_nelems = (size_t)n_param * n_param; if(n_param == 0) a_nelems = 0;
    size_t e_nelems = (size_t)n_param * n_param; if(n_param == 0) e_nelems = 0;
    size_t b_nelems = (jobb_f == 'B') ? (size_t)n_param * m_param : (size_t)n_param * n_param; // B or G
    if(n_param == 0 || (jobb_f == 'B' && m_param == 0)) b_nelems = 0;
    if(jobb_f == 'G' && n_param == 0) b_nelems = 0;

    size_t q_nelems = (fact_f == 'N' || fact_f == 'D') ? (size_t)n_param*n_param : (size_t)p_param*n_param; // Q or C_factor
    if((fact_f == 'N' || fact_f == 'D') && n_param == 0) q_nelems = 0;
    if((fact_f == 'C' || fact_f == 'B') && (p_param == 0 || n_param == 0)) q_nelems = 0;
    
    size_t r_nelems = 0;
    if (jobb_f == 'B') {
        r_nelems = (fact_f == 'N' || fact_f == 'C') ? (size_t)m_param*m_param : (size_t)p_param*m_param; // R or D_factor
        if((fact_f == 'N' || fact_f == 'C') && m_param == 0) r_nelems = 0;
        if((fact_f == 'D' || fact_f == 'B') && (p_param == 0 || m_param == 0)) r_nelems = 0;
    }
    
    size_t l_nelems = 0;
    if (jobb_f == 'B' && jobl_f == 'N') {
        l_nelems = (size_t)n_param * m_param;
        if(n_param == 0 || m_param == 0) l_nelems = 0;
    }

    size_t x_nelems = (size_t)n_param * n_param; if(n_param == 0) x_nelems = 0;
    size_t s_nelems = (jobb_f == 'B') ? (size_t)(2*n_param+m_param)*(2*n_param+m_param) : (size_t)(2*n_param)*(2*n_param); // Approx, actual cols vary
    if(n_param == 0 && (jobb_f=='G' || (jobb_f=='B' && m_param == 0))) s_nelems = 0;
    size_t t_nelems = (jobb_f == 'B') ? (size_t)(2*n_param+m_param)*(2*n_param) : (size_t)(2*n_param)*(2*n_param);
    if(n_param == 0) t_nelems = 0;
    size_t u_nelems = (size_t)(2*n_param)*(2*n_param);
    if(n_param == 0) u_nelems = 0;


    if (row_major) {
        if (a_nelems > 0) { a_cm = (double*)malloc(a_nelems * sizeof(double)); CHECK_ALLOC(a_cm); }
        if (e_nelems > 0) { e_cm = (double*)malloc(e_nelems * sizeof(double)); CHECK_ALLOC(e_cm); }
        if (b_nelems > 0) { b_cm = (double*)malloc(b_nelems * sizeof(double)); CHECK_ALLOC(b_cm); }
        if (q_nelems > 0) { q_cm = (double*)malloc(q_nelems * sizeof(double)); CHECK_ALLOC(q_cm); }
        if (r_nelems > 0) { r_cm = (double*)malloc(r_nelems * sizeof(double)); CHECK_ALLOC(r_cm); }
        if (l_nelems > 0) { l_cm = (double*)malloc(l_nelems * sizeof(double)); CHECK_ALLOC(l_cm); }
        
        if (x_nelems > 0) { x_cm = (double*)malloc(x_nelems * sizeof(double)); CHECK_ALLOC(x_cm); }
        // S, T, U output buffers
        int s_cols_f = (jobb_f == 'B') ? (2*n_param + m_param) : (2*n_param);
        int t_cols_f = 2*n_param;
        int u_cols_f = 2*n_param;
        if (n_param > 0) { // Only allocate if N > 0, as S,T,U are 2N-dimensional
             s_cm = (double*)malloc((size_t)min_lds_f * s_cols_f * sizeof(double)); CHECK_ALLOC(s_cm);
             t_cm = (double*)malloc((size_t)min_ldt_f * t_cols_f * sizeof(double)); CHECK_ALLOC(t_cm);
             u_cm = (double*)malloc((size_t)min_ldu_f * u_cols_f * sizeof(double)); CHECK_ALLOC(u_cm);
        }
    }

    // 5. Prepare Fortran parameters and conversions
    double* a_ptr = a_io, *e_ptr = e_io, *b_ptr = b_io, *q_ptr = q_io, *r_ptr = r_io, *l_ptr = l_io;
    double* x_ptr = x_out, *s_ptr = s_out, *t_ptr = t_out, *u_ptr = u_out;

    lda_f = lda; lde_f = lde; ldb_f = ldb; ldq_f = ldq; ldr_f = ldr; ldl_f = ldl;
    ldx_f = ldx; lds_f = lds; ldt_f = ldt; ldu_f = ldu;

    if (row_major) {
        lda_f = MAX(1,n_param); lde_f = MAX(1,n_param); ldb_f = MAX(1,n_param);
        ldq_f = (fact_f == 'N' || fact_f == 'D') ? MAX(1,n_param) : MAX(1,p_param);
        ldr_f = (jobb_f == 'B' && (fact_f == 'N' || fact_f == 'C')) ? MAX(1,m_param) : ((jobb_f == 'B' && (fact_f == 'D' || fact_f == 'B')) ? MAX(1,p_param) : 1);
        ldl_f = (jobb_f == 'B' && jobl_f == 'N') ? MAX(1,n_param) : 1;
        ldx_f = MAX(1,n_param); 
        lds_f = min_lds_f; ldt_f = min_ldt_f; ldu_f = min_ldu_f;


        if (a_nelems > 0) { slicot_transpose_to_fortran_with_ld(a_io, a_cm, n_param, n_param, lda, lda_f, sizeof(double)); a_ptr = a_cm; } else {a_ptr = NULL;}
        if (e_nelems > 0) { slicot_transpose_to_fortran_with_ld(e_io, e_cm, n_param, n_param, lde, lde_f, sizeof(double)); e_ptr = e_cm; } else {e_ptr = NULL;}
        if (b_nelems > 0) { 
            int b_rows = n_param;
            int b_cols = (jobb_f == 'B') ? m_param : n_param;
            slicot_transpose_to_fortran_with_ld(b_io, b_cm, b_rows, b_cols, ldb, ldb_f, sizeof(double)); b_ptr = b_cm; 
        } else {b_ptr = NULL;}
        if (q_nelems > 0) {
            int q_rows = (fact_f == 'N' || fact_f == 'D') ? n_param : p_param;
            int q_cols = n_param;
            slicot_transpose_to_fortran_with_ld(q_io, q_cm, q_rows, q_cols, ldq, ldq_f, sizeof(double)); q_ptr = q_cm;
        } else {q_ptr = NULL;}
        if (r_nelems > 0 && jobb_f == 'B') {
            int r_rows = (fact_f == 'N' || fact_f == 'C') ? m_param : p_param;
            int r_cols = m_param;
            slicot_transpose_to_fortran_with_ld(r_io, r_cm, r_rows, r_cols, ldr, ldr_f, sizeof(double)); r_ptr = r_cm;
        } else if (jobb_f == 'B') {r_ptr = NULL;} // R not used if JOBB='G' or r_nelems is 0
        
        if (l_nelems > 0 && jobb_f == 'B' && jobl_f == 'N') {
            slicot_transpose_to_fortran_with_ld(l_io, l_cm, n_param, m_param, ldl, ldl_f, sizeof(double)); l_ptr = l_cm;
        } else if (jobb_f == 'B' && jobl_f == 'N') {l_ptr = NULL;}
        // L not used if JOBL='Z' or JOBB='G'

        if (x_nelems > 0) x_ptr = x_cm; else x_ptr = NULL;
        if (n_param > 0) { // S, T, U only if N > 0
             s_ptr = s_cm; t_ptr = t_cm; u_ptr = u_cm;
        } else {
             s_ptr = NULL; t_ptr = NULL; u_ptr = NULL;
        }

    } else { // Column-major C
        if(a_nelems == 0) a_ptr = NULL; if(e_nelems == 0) e_ptr = NULL; if(b_nelems == 0) b_ptr = NULL;
        if(q_nelems == 0) q_ptr = NULL; if(r_nelems == 0 && jobb_f == 'B') r_ptr = NULL; 
        if(l_nelems == 0 && jobb_f == 'B' && jobl_f == 'N') l_ptr = NULL;
        if(x_nelems == 0) x_ptr = NULL;
        if(n_param == 0) { s_ptr = NULL; t_ptr = NULL; u_ptr = NULL; }
    }
    
    int n_f_call = n_param, m_f_call = m_param, p_f_call = p_param;
    int ldwork_f_call = ldwork;

    // 7. Call Fortran
    F77_FUNC(sg02ad, SG02AD)(&dico_f, &jobb_f, &fact_f, &uplo_f, &jobl_f, &scal_f, &sort_f, &acc_f,
                             &n_f_call, &m_f_call, &p_f_call,
                             a_ptr, &lda_f, e_ptr, &lde_f, b_ptr, &ldb_f,
                             q_ptr, &ldq_f, r_ptr, &ldr_f, l_ptr, &ldl_f,
                             rcondu, x_ptr, &ldx_f,
                             alfar_out, alfai_out, beta_out,
                             s_ptr, &lds_f, t_ptr, &ldt_f, u_ptr, &ldu_f,
                             &tol_f, iwork, dwork, &ldwork_f_call, bwork_c,
                             &local_iwarn, &info,
                             1,1,1,1,1,1,1,1); // Lengths of char params

    if(iwarn_out != NULL) *iwarn_out = local_iwarn;

    // 8. Convert output matrices back if row_major
    if (row_major && info == 0) {
        if (x_nelems > 0 && x_cm != NULL && x_out != NULL) {
            slicot_transpose_to_c_with_ld(x_cm, x_out, n_param, n_param, ldx_f, ldx, sizeof(double));
        }
        if (n_param > 0) { // S, T, U only if N > 0
            int s_rows_f = 2*n_param; // Based on text "leading 2N-by-2N part"
            int s_cols_f = (jobb_f == 'B') ? (2*n_param + m_param) : (2*n_param);
            int t_rows_f = 2*n_param;
            int t_cols_f = 2*n_param;
            int u_rows_f = 2*n_param;
            int u_cols_f = 2*n_param;

            if (s_cm != NULL && s_out != NULL) slicot_transpose_to_c_with_ld(s_cm, s_out, s_rows_f, s_cols_f, lds_f, lds, sizeof(double));
            if (t_cm != NULL && t_out != NULL) slicot_transpose_to_c_with_ld(t_cm, t_out, t_rows_f, t_cols_f, ldt_f, ldt, sizeof(double));
            if (u_cm != NULL && u_out != NULL) slicot_transpose_to_c_with_ld(u_cm, u_out, u_rows_f, u_cols_f, ldu_f, ldu, sizeof(double));
        }
    }

cleanup:
    free(iwork);
    free(dwork);
    free(bwork_c);
    if(row_major){
        free(a_cm); free(e_cm); free(b_cm); free(q_cm); free(r_cm); free(l_cm);
        free(x_cm); free(s_cm); free(t_cm); free(u_cm);
    }
    return info;
}
