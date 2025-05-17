/**
 * @file sg03bd.c
 * @brief C wrapper for SLICOT routine SG03BD.
 * @details Solves for Cholesky factor of generalized stable continuous- or
 * discrete-time Lyapunov equations.
 * Matrices A, E, Q, Z are input/output. Matrix B is input (RHS term)
 * and output (Cholesky factor U).
 * Workspace (DWORK) is allocated internally.
 * Input/output matrix format is handled via the row_major parameter.
 */

#include <stdlib.h> // For malloc, free
#include <string.h> // For memcpy
#include <ctype.h>  // For toupper
#include <math.h>   // For fmax, fabs (from C standard library)

#include "sg03bd.h"       // Public header for this wrapper
#include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions etc.
#include "slicot_f77.h"   // Provides F77_FUNC macro

/* External Fortran routine declaration */
extern void F77_FUNC(sg03bd, SG03BD)(
    const char* dico, const char* fact, const char* trans,
    const int* n, const int* m,
    double* a, const int* lda, /* input/output */
    double* e, const int* lde, /* input/output */
    double* q, const int* ldq, /* input/output */
    double* z, const int* ldz, /* input/output */
    double* b_u, const int* ldb, /* input: B_rhs, output: U_cholesky (in B_U) */
    double* scale,             /* output */
    double* alphar, double* alfai, double* beta, /* output */
    double* dwork, const int* ldwork, /* workspace */
    int* info, /* output */
    size_t dico_len, size_t fact_len, size_t trans_len
);

SLICOT_EXPORT
int slicot_sg03bd(
    char dico_param, char fact_param, char trans_param,
    int n_param, int m_param,
    double* a_io, int lda,
    double* e_io, int lde,
    double* q_io, int ldq, 
    double* z_io, int ldz, 
    double* b_in_u_out, int ldb, 
    double* scale_out,
    double* alphar_out, double* alfai_out, double* beta_out,
    int row_major)
{
    // 1. Variable declarations
    int info = 0;
    double *dwork = NULL;
    int ldwork = 0;

    double *a_cm = NULL, *e_cm = NULL, *q_cm = NULL, *z_cm = NULL, *b_u_cm = NULL;

    int lda_f, lde_f, ldq_f, ldz_f, ldb_f;
    
    char dico_f = toupper(dico_param);
    char fact_f = toupper(fact_param);
    char trans_f= toupper(trans_param);

    // 2. Input parameter validation
    if (strchr("CD", dico_f) == NULL) { info = -1; goto cleanup; }
    if (strchr("NF", fact_f) == NULL) { info = -2; goto cleanup; }
    if (strchr("NT", trans_f) == NULL) { info = -3; goto cleanup; }
    if (n_param < 0) { info = -4; goto cleanup; }
    if (m_param < 0) { info = -5; goto cleanup; }

    if (a_io == NULL && n_param > 0) { info = -6; goto cleanup; }
    if (e_io == NULL && n_param > 0) { info = -8; goto cleanup; }
    if (q_io == NULL && n_param > 0) { info = -10; goto cleanup; }
    if (z_io == NULL && n_param > 0) { info = -12; goto cleanup; }
    // b_in_u_out: input B part needs M rows if TRANS='N', N rows if TRANS='T'. Output U is N x N.
    // Pointer must be valid if either input B part is non-empty or output U part is non-empty.
    if (b_in_u_out == NULL && ( (m_param > 0 && n_param > 0) || (n_param > 0) ) ) { info = -14; goto cleanup; }
    
    if (scale_out == NULL) { info = -16; goto cleanup; }
    if (alphar_out == NULL && n_param > 0) { info = -17; goto cleanup; } // ALPHAR, ALFAI, BETA are N-dim
    if (alfai_out == NULL && n_param > 0) { info = -18; goto cleanup; }
    if (beta_out == NULL && n_param > 0) { info = -19; goto cleanup; }

    int min_n_dim_f = MAX(1, n_param);
    int ldb_f_req; // Required Fortran LDB based on TRANS and output U
    if (trans_f == 'T') { // Input B is N x M, Output U is N x N
        ldb_f_req = MAX(1, n_param); 
    } else { // TRANS = 'N', Input B is M x N, Output U is N x N
        ldb_f_req = MAX(1, MAX(m_param, n_param));
    }

    if (row_major) {
        if (n_param > 0 && lda < n_param) { info = -7; goto cleanup; }
        if (n_param > 0 && lde < n_param) { info = -9; goto cleanup; }
        if (n_param > 0 && ldq < n_param) { info = -11; goto cleanup; }
        if (n_param > 0 && ldz < n_param) { info = -13; goto cleanup; }
        // For b_in_u_out (row-major), LDB is number of columns.
        // Input B part: if TRANS='T', B is NxM, LDB should be M. If TRANS='N', B is MxN, LDB should be N.
        // Output U part: U is NxN, LDB should be N.
        // So, LDB must be N for output. If input B requires different col count, user must handle.
        if (n_param > 0 && ldb < n_param) { info = -15; goto cleanup; } 
    } else { // Column-major C
        if (n_param > 0 && lda < min_n_dim_f) { info = -7; goto cleanup; }
        if (n_param > 0 && lde < min_n_dim_f) { info = -9; goto cleanup; }
        if (n_param > 0 && ldq < min_n_dim_f) { info = -11; goto cleanup; }
        if (n_param > 0 && ldz < min_n_dim_f) { info = -13; goto cleanup; }
        if ( ( (trans_f == 'T' && n_param > 0 && m_param > 0) || (trans_f == 'N' && m_param > 0 && n_param > 0) || (n_param > 0) ) && ldb < ldb_f_req ) { info = -15; goto cleanup; }
    }
    if (info != 0) { goto cleanup; }

    // 3. Workspace Allocation
    int ldwork_query_val = -1;
    double dwork_query_result[1];
    int temp_info_ws_query = 0; // Use a temporary info for workspace query
    F77_FUNC(sg03bd, SG03BD)(&dico_f, &fact_f, &trans_f, &n_param, &m_param,
                             NULL, &min_n_dim_f, NULL, &min_n_dim_f, NULL, &min_n_dim_f, NULL, &min_n_dim_f,
                             NULL, &ldb_f_req, 
                             scale_out, NULL, NULL, NULL,
                             dwork_query_result, &ldwork_query_val, &temp_info_ws_query,
                             1,1,1);
    
    if (temp_info_ws_query == 0 && ldwork_query_val == -1) { 
        ldwork = (int)dwork_query_result[0];
    } else { 
        if (fact_f == 'N') {
            ldwork = MAX(1, MAX(4 * n_param, 6 * n_param - 6));
        } else { // FACT = 'F'
            ldwork = MAX(1, MAX(2 * n_param, 6 * n_param - 6));
        }
    }
    ldwork = MAX(1, ldwork);
    dwork = (double*)malloc((size_t)ldwork * sizeof(double));
    CHECK_ALLOC(dwork);

    // 4. Memory for column-major copies
    size_t n_sq_nelems = (size_t)n_param * n_param; if(n_param == 0) n_sq_nelems = 0;
    
    // b_u_cm is for Fortran's B array, which is input B then output U.
    // Its Fortran dimensions are LDB_F x N1 (see HTML for N1).
    // LDB_F = ldb_f_req. N1 = (TRANS='T' ? MAX(M,N) : N).
    // For simplicity, allocate b_u_cm to hold an N x N matrix, as U is N x N.
    // The input B part will be copied into this.
    size_t b_u_cm_nelems = n_sq_nelems; 

    if (row_major) {
        if (n_sq_nelems > 0) {
            a_cm = (double*)malloc(n_sq_nelems * sizeof(double)); CHECK_ALLOC(a_cm);
            e_cm = (double*)malloc(n_sq_nelems * sizeof(double)); CHECK_ALLOC(e_cm);
            q_cm = (double*)malloc(n_sq_nelems * sizeof(double)); CHECK_ALLOC(q_cm);
            z_cm = (double*)malloc(n_sq_nelems * sizeof(double)); CHECK_ALLOC(z_cm);
            b_u_cm = (double*)malloc(b_u_cm_nelems * sizeof(double)); CHECK_ALLOC(b_u_cm);
        }
    }

    // 5. Prepare Fortran parameters and conversions
    double* a_ptr = a_io; double* e_ptr = e_io; double* q_ptr = q_io;
    double* z_ptr = z_io; double* b_u_ptr = b_in_u_out;

    lda_f=lda; lde_f=lde; ldq_f=ldq; ldz_f=ldz; ldb_f=ldb;

    if (row_major && n_param > 0) {
        lda_f=min_n_dim_f; lde_f=min_n_dim_f; ldq_f=min_n_dim_f; ldz_f=min_n_dim_f;
        ldb_f=ldb_f_req; // Fortran LDB

        slicot_transpose_to_fortran_with_ld(a_io,a_cm,n_param,n_param,lda,lda_f,sizeof(double)); a_ptr=a_cm;
        slicot_transpose_to_fortran_with_ld(e_io,e_cm,n_param,n_param,lde,lde_f,sizeof(double)); e_ptr=e_cm;
        
        if (fact_f == 'F') { 
            slicot_transpose_to_fortran_with_ld(q_io,q_cm,n_param,n_param,ldq,ldq_f,sizeof(double)); q_ptr=q_cm;
            slicot_transpose_to_fortran_with_ld(z_io,z_cm,n_param,n_param,ldz,ldz_f,sizeof(double)); z_ptr=z_cm;
        } else { 
            q_ptr=q_cm; z_ptr=z_cm; // Fortran writes to these
        }
        
        // Copy input B from b_in_u_out (MxN or NxM C array) to b_u_cm (Fortran LDB_F x N1_F)
        // The b_in_u_out C array has LDB = N_param (cols).
        b_u_ptr = b_u_cm; // Fortran will use b_u_cm
        if (m_param > 0) { // Only copy if B is non-empty
            if (trans_f == 'T') { // Input B is N x M (row-major)
                slicot_transpose_to_fortran_with_ld(b_in_u_out, b_u_cm, n_param, m_param, ldb, ldb_f, sizeof(double));
            } else { // Input B is M x N (row-major)
                slicot_transpose_to_fortran_with_ld(b_in_u_out, b_u_cm, m_param, n_param, ldb, ldb_f, sizeof(double));
            }
        } else if (n_param > 0) { // M=0, N>0. B is empty. U will be zero. b_u_cm is for output U.
             memset(b_u_cm, 0, b_u_cm_nelems * sizeof(double)); // Initialize for safety
        }


    } else { // Column-major C or N=0
        if (n_param == 0) {
            a_ptr=(a_io==NULL?NULL:a_io); lda_f=(a_io==NULL?1:lda);
            e_ptr=(e_io==NULL?NULL:e_io); lde_f=(e_io==NULL?1:lde);
            q_ptr=(q_io==NULL?NULL:q_io); ldq_f=(q_io==NULL?1:ldq);
            z_ptr=(z_io==NULL?NULL:z_io); ldz_f=(z_io==NULL?1:ldz);
            b_u_ptr=(b_in_u_out==NULL?NULL:b_in_u_out); ldb_f=(b_in_u_out==NULL?1:ldb);
        } else { // N > 0, column major C
            ldb_f = ldb; // Use user's LDB, already validated against ldb_f_req
        }
    }
    
    int n_f_call = n_param, m_f_call = m_param;
    int ldwork_f_call = ldwork;

    F77_FUNC(sg03bd, SG03BD)(&dico_f, &fact_f, &trans_f, &n_f_call, &m_f_call,
                             a_ptr, &lda_f, e_ptr, &lde_f, q_ptr, &ldq_f, z_ptr, &ldz_f,
                             b_u_ptr, &ldb_f, 
                             scale_out, alphar_out, alfai_out, beta_out,
                             dwork, &ldwork_f_call, &info,
                             1,1,1); 

    if (row_major && info == 0 && n_param > 0) {
        if (a_cm != NULL && a_io != NULL) slicot_transpose_to_c_with_ld(a_cm,a_io,n_param,n_param,lda_f,lda,sizeof(double));
        if (e_cm != NULL && e_io != NULL) slicot_transpose_to_c_with_ld(e_cm,e_io,n_param,n_param,lde_f,lde,sizeof(double));
        if (fact_f == 'N') { 
            if (q_cm != NULL && q_io != NULL) slicot_transpose_to_c_with_ld(q_cm,q_io,n_param,n_param,ldq_f,ldq,sizeof(double));
            if (z_cm != NULL && z_io != NULL) slicot_transpose_to_c_with_ld(z_cm,z_io,n_param,n_param,ldz_f,ldz,sizeof(double));
        }
        // Output U (Cholesky factor) is in b_u_ptr (which was b_u_cm). Transpose to b_in_u_out.
        if (b_u_cm != NULL && b_in_u_out != NULL) {
            slicot_transpose_to_c_with_ld(b_u_cm, b_in_u_out, n_param, n_param, ldb_f, ldb, sizeof(double));
        }
    }

cleanup:
    free(dwork);
    if(row_major){
        free(a_cm); free(e_cm); free(q_cm); free(z_cm); free(b_u_cm);
    }
    return info;
}
