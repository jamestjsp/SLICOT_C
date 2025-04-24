/**
 * @file ib01bd.c
 * @brief C wrapper implementation for SLICOT routine IB01BD.
 * @details Estimates system matrices (A, B, C, D), noise covariances (Q, Ry, S),
 * and Kalman gain (K) using processed data from IB01AD.
 * Workspace (IWORK, DWORK, BWORK) is allocated internally by this wrapper.
 * Input/output matrix format is handled via the row_major parameter, except for R.
 */

 #include <stdlib.h> // For malloc, free
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t
 #include <stdbool.h>// For bool type
 #include <math.h>   // For MAX/MIN if needed (often provided by slicot_utils.h)
 #include <stdio.h>  // For error logging (optional)
 
 #include "ib01bd.h"       // Public header for this wrapper
 #include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX/MIN, transpose functions etc.
 #include "slicot_f77.h"   // Provides F77_FUNC macro
 
 /* External Fortran routine declaration */
 extern void F77_FUNC(ib01bd, IB01BD)(
     const char* meth, const char* job, const char* jobck,
     const int* nobr, const int* n, const int* m, const int* l,
     const int* nsmpl, double* r, const int* ldr,
     double* a, const int* lda, double* c, const int* ldc,
     double* b, const int* ldb, double* d, const int* ldd,
     double* q, const int* ldq, double* ry, const int* ldry,
     double* s, const int* lds, double* k, const int* ldk,
     const double* tol, int* iwork, double* dwork, const int* ldwork,
     int* bwork, // Note: Fortran LOGICAL maps to int* for wrapper usage
     int* iwarn, int* info,
     int meth_len, int job_len, int jobck_len // Hidden string lengths
 );
 
 /* --- Internal Helper Functions for Workspace Calculation --- */
 
 // Calculate LIWORK based on documentation formulas
 static int calculate_liwork_internal(char meth, char job, char jobck, int n, int m, int l, int nobr) {
     int liw1 = 1;
     int liw2 = 0;
     int mnobr = m * nobr;
     int lnobr = l * nobr;
 
     if (meth != 'N') { // METH = 'M' or 'C'
         if (m == 0 || job == 'C') {
             liw1 = (jobck == 'N') ? n : MAX(n, mnobr + n);
         } else { // M > 0 and JOB <> 'C'
              liw1 = (jobck == 'N') ? MAX(lnobr, mnobr) : MAX(lnobr, mnobr + n);
         }
         // If METH is 'C', the N4SID part for B/D might increase requirement
         if (meth == 'C' && job != 'C') {
              liw1 = MAX(liw1, MAX(mnobr + n, m * (n + l)));
         }
     } else { // METH = 'N'
          liw1 = MAX(mnobr + n, m * (n + l));
     }
 
     if (jobck == 'K') {
         liw2 = n * n;
     }
 
     return MAX(1, MAX(liw1, liw2)); // Ensure minimum size 1
 }
 
 // Calculate LDWORK based on documentation formulas (complex!)
 // Corrected MAX usage.
 static int calculate_ldwork_internal(char meth, char job, char jobck, int n, int m, int l, int nobr) {
     int ldw1 = 1, ldw2 = 1, ldw3 = 1;
     int lnobr_l = l * nobr - l;
     int mnobr = m * nobr;
     int lnobr = l * nobr;
     int n_l = n + l; // Used in METH='N' part
     int m_nl = m * n_l; // Used in METH='N' part
 
     if (meth == 'M') {
         int term1_m = 2 * lnobr_l * n + 2 * n;
         int term2_m = lnobr_l * n + n * n + 7 * n;
         if (job == 'C' || (job == 'A' && m == 0)) {
             ldw1 = MAX(term1_m, term2_m);
         } else { // M > 0 and JOB = 'A', 'B', or 'D'
             int term3_m = lnobr_l * n + n + 6 * mnobr;
             int term4_m_sub = MAX(3 * lnobr + 1, m);
             int term4_m = lnobr_l * n + n + MAX(l + mnobr, lnobr + term4_m_sub);
             ldw1 = MAX(MAX(term1_m, term2_m), MAX(term3_m, term4_m));
         }
 
         if (jobck != 'N') {
             int aw = (m == 0 || job == 'C') ? (n + n * n) : 0;
             int term5_m_sub1 = MAX(5 * n, (2 * m + l) * nobr + l);
             int term5_m_sub2 = MAX(4 * (mnobr + n) + 1, mnobr + 2 * n + l);
             int term5_m = lnobr * n + MAX(lnobr_l * n + aw + 2 * n + term5_m_sub1, term5_m_sub2);
             ldw2 = term5_m;
         } else {
             ldw2 = 0;
         }
     } else if (meth == 'N') {
         int term1_n_sub = MAX(lnobr_l * n + 2 * n + (2 * m + l) * nobr + l,
                               2 * lnobr_l * n + n * n + 8 * n);
         int term1_n = lnobr * n + term1_n_sub;
         int term2_n = MAX(n + 4 * (mnobr + n) + 1, mnobr + 3 * n + l);
         ldw1 = MAX(term1_n, term2_n);
 
         if (m == 0 || job == 'C') {
             ldw2 = 0;
         } else { // M > 0 and JOB = 'A', 'B', or 'D'
             ldw2 = lnobr * n + mnobr * n_l * (m_nl + 1) + MAX(n_l * n_l, 4 * m_nl + 1);
         }
     } else { // METH = 'C'
         // LDW1 is max of (METH=M, JOB=C) and (METH=N) parts
         int ldw1_m_jobc = MAX(2 * lnobr_l * n + 2 * n, lnobr_l * n + n * n + 7 * n);
         int term1_n_sub_c = MAX(lnobr_l * n + 2 * n + (2 * m + l) * nobr + l,
                                 2 * lnobr_l * n + n * n + 8 * n);
         int term1_n_c = lnobr * n + term1_n_sub_c;
         int term2_n_c = MAX(n + 4 * (mnobr + n) + 1, mnobr + 3 * n + l);
         int ldw1_n_c = MAX(term1_n_c, term2_n_c);
         ldw1 = MAX(ldw1_m_jobc, ldw1_n_c);
 
         // LDW2 is from METH=N part (for B/D calculation)
         if (m == 0 || job == 'C') {
              ldw2 = 0;
         } else { // M > 0 and JOB = 'A', 'B', or 'D'
             ldw2 = lnobr * n + mnobr * n_l * (m_nl + 1) + MAX(n_l * n_l, 4 * m_nl + 1);
         }
     }
 
     if (jobck == 'K') {
         int term3_k_sub = MAX(3 * l, n * l);
         int term3_k = MAX(4 * n * n + 2 * n * l + l * l + term3_k_sub,
                           14 * n * n + 12 * n + 5);
         ldw3 = term3_k;
     } else {
         ldw3 = 0;
     }
 
     return MAX(1, MAX(MAX(ldw1, ldw2), ldw3)); // Ensure minimum size 1
 }
 
 
 /* --- C Wrapper Function Definition --- */
 SLICOT_EXPORT
 int slicot_ib01bd(char meth, char job, char jobck, int nobr, int n, int m, int l,
                   int nsmpl, double *r, int ldr,
                   double *a, int lda, double *c, int ldc,
                   double *b, int ldb, double *d, int ldd,
                   double *q, int ldq, double *ry, int ldry, double *s, int lds,
                   double *k, int ldk, double tol, int *iwarn, int row_major)
 {
     // 1. Variable declarations
     int info = 0;
     int iwarn_local = 0;
 
     // Internal workspace pointers
     int *iwork = NULL;
     double *dwork = NULL;
     int *bwork = NULL; // Use int for LOGICAL array
     int liwork = 0;
     int ldwork = 0;
     int lbwork = 0;
 
     // Pointers for column-major copies
     double *a_cm = NULL, *c_cm = NULL, *b_cm = NULL, *d_cm = NULL;
     double *q_cm = NULL, *ry_cm = NULL, *s_cm = NULL, *k_cm = NULL;
     // R is assumed to be column-major input
 
     // Fortran-compatible parameters
     char meth_upper = toupper(meth);
     char job_upper = toupper(job);
     char jobck_upper = toupper(jobck);
     const int meth_len = 1, job_len = 1, jobck_len = 1;
 
     // Pointers to pass to Fortran
     double *a_ptr, *c_ptr, *b_ptr, *d_ptr;
     double *q_ptr, *ry_ptr, *s_ptr, *k_ptr;
     int lda_f, ldc_f, ldb_f, ldd_f;
     int ldq_f, ldry_f, lds_f, ldk_f;
 
     // --- 2. Input Parameter Validation ---
     if (meth_upper != 'M' && meth_upper != 'N' && meth_upper != 'C') { info = -1; goto cleanup; }
     if (job_upper != 'A' && job_upper != 'C' && job_upper != 'B' && job_upper != 'D') { info = -2; goto cleanup; }
     if (jobck_upper != 'C' && jobck_upper != 'K' && jobck_upper != 'N') { info = -3; goto cleanup; }
     if (nobr <= 1) { info = -4; goto cleanup; }
     if (n <= 0 || n >= nobr) { info = -5; goto cleanup; }
     if (m < 0) { info = -6; goto cleanup; }
     if (l <= 0) { info = -7; goto cleanup; }
     if ((jobck_upper == 'C' || jobck_upper == 'K') && nsmpl < 2 * (m + l) * nobr) { info = -8; goto cleanup; }
     if (r == NULL) { info = -9; goto cleanup; } // R is required input
     if (ldr < 2 * (m + l) * nobr) { info = -10; goto cleanup; }
 
     // Determine if matrices are input or output based on flags
     bool is_A_input = (meth_upper == 'N' || meth_upper == 'C') && (job_upper == 'B' || job_upper == 'D');
     bool is_A_output = (job_upper == 'A' || job_upper == 'C');
     bool is_C_input = is_A_input; // Same condition for C input
     bool is_C_output = is_A_output; // Same condition for C output
     bool is_B_output = (m > 0) && (job_upper == 'A' || job_upper == 'B' || job_upper == 'D');
     bool is_D_output = (m > 0) && (job_upper == 'A' || job_upper == 'D');
     bool is_Q_output = (jobck_upper == 'C' || jobck_upper == 'K');
     bool is_RY_output = is_Q_output;
     bool is_S_output = is_Q_output;
     bool is_K_output = (jobck_upper == 'K');
 
     // Validate pointers for required arrays
     if ((is_A_input || is_A_output) && a == NULL) { info = -11; goto cleanup; }
     if ((is_C_input || is_C_output) && c == NULL) { info = -13; goto cleanup; }
     if (is_B_output && b == NULL) { info = -15; goto cleanup; }
     if (is_D_output && d == NULL) { info = -17; goto cleanup; }
     if (is_Q_output && q == NULL) { info = -19; goto cleanup; }
     if (is_RY_output && ry == NULL) { info = -21; goto cleanup; }
     if (is_S_output && s == NULL) { info = -23; goto cleanup; }
     if (is_K_output && k == NULL) { info = -25; goto cleanup; }
     // tol is scalar, iwarn is output pointer (can be NULL)
 
     // Validate leading dimensions
     int min_lda_f = MAX(1, n); int min_ldc_f = MAX(1, l); int min_ldb_f = MAX(1, n);
     int min_ldd_f = MAX(1, l); int min_ldq_f = MAX(1, n); int min_ldry_f = MAX(1, l);
     int min_lds_f = MAX(1, n); int min_ldk_f = MAX(1, n);
 
     if (row_major) {
         // C LD is columns
         if ((is_A_input || is_A_output) && lda < n) { info = -12; goto cleanup; }
         if ((is_C_input || is_C_output) && ldc < n) { info = -14; goto cleanup; }
         if (is_B_output && ldb < m) { info = -16; goto cleanup; }
         if (is_D_output && ldd < m) { info = -18; goto cleanup; }
         if (is_Q_output && ldq < n) { info = -20; goto cleanup; }
         if (is_RY_output && ldry < l) { info = -22; goto cleanup; }
         if (is_S_output && lds < l) { info = -24; goto cleanup; }
         if (is_K_output && ldk < l) { info = -26; goto cleanup; }
     } else {
         // C LD is rows
         if ((is_A_input || is_A_output) && lda < min_lda_f) { info = -12; goto cleanup; }
         if ((is_C_input || is_C_output) && ldc < min_ldc_f) { info = -14; goto cleanup; }
         if (is_B_output && ldb < min_ldb_f) { info = -16; goto cleanup; }
         if (is_D_output && ldd < min_ldd_f) { info = -18; goto cleanup; }
         if (is_Q_output && ldq < min_ldq_f) { info = -20; goto cleanup; }
         if (is_RY_output && ldry < min_ldry_f) { info = -22; goto cleanup; }
         if (is_S_output && lds < min_lds_f) { info = -24; goto cleanup; }
         if (is_K_output && ldk < min_ldk_f) { info = -26; goto cleanup; }
     }
 
     // --- 3. Internal Workspace Allocation ---
     liwork = calculate_liwork_internal(meth_upper, job_upper, jobck_upper, n, m, l, nobr);
     iwork = (int*)malloc((size_t)liwork * sizeof(int));
     CHECK_ALLOC(iwork);
 
     ldwork = calculate_ldwork_internal(meth_upper, job_upper, jobck_upper, n, m, l, nobr);
     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork);
 
     lbwork = (jobck_upper == 'K') ? 2 * n : 0;
     if (lbwork > 0) {
         bwork = (int*)malloc((size_t)lbwork * sizeof(int)); // Allocate int for LOGICAL
         CHECK_ALLOC(bwork);
     }
 
     // --- 4. Memory Allocation for Column-Major Copies ---
     size_t a_size = (size_t)n * n; if (n == 0) a_size = 0;
     size_t c_size = (size_t)l * n; if (l == 0 || n == 0) c_size = 0;
     size_t b_size = (size_t)n * m; if (n == 0 || m == 0) b_size = 0;
     size_t d_size = (size_t)l * m; if (l == 0 || m == 0) d_size = 0;
     size_t q_size = (size_t)n * n; if (n == 0) q_size = 0;
     size_t ry_size = (size_t)l * l; if (l == 0) ry_size = 0;
     size_t s_size = (size_t)n * l; if (n == 0 || l == 0) s_size = 0;
     size_t k_size = (size_t)n * l; if (n == 0 || l == 0) k_size = 0;
 
     if (row_major) {
         // Allocate only if needed for input or output and size > 0
         if ((is_A_input || is_A_output) && a_size > 0) { a_cm = (double*)malloc(a_size * sizeof(double)); CHECK_ALLOC(a_cm); }
         if ((is_C_input || is_C_output) && c_size > 0) { c_cm = (double*)malloc(c_size * sizeof(double)); CHECK_ALLOC(c_cm); }
         if (is_B_output && b_size > 0) { b_cm = (double*)malloc(b_size * sizeof(double)); CHECK_ALLOC(b_cm); }
         if (is_D_output && d_size > 0) { d_cm = (double*)malloc(d_size * sizeof(double)); CHECK_ALLOC(d_cm); }
         if (is_Q_output && q_size > 0) { q_cm = (double*)malloc(q_size * sizeof(double)); CHECK_ALLOC(q_cm); }
         if (is_RY_output && ry_size > 0) { ry_cm = (double*)malloc(ry_size * sizeof(double)); CHECK_ALLOC(ry_cm); }
         if (is_S_output && s_size > 0) { s_cm = (double*)malloc(s_size * sizeof(double)); CHECK_ALLOC(s_cm); }
         if (is_K_output && k_size > 0) { k_cm = (double*)malloc(k_size * sizeof(double)); CHECK_ALLOC(k_cm); }
     }
 
     // --- 5. Prepare Fortran Parameters and Conversions ---
     // Set default pointers and Fortran LDs
     a_ptr = a; lda_f = lda; c_ptr = c; ldc_f = ldc; b_ptr = b; ldb_f = ldb; d_ptr = d; ldd_f = ldd;
     q_ptr = q; ldq_f = ldq; ry_ptr = ry; ldry_f = ldry; s_ptr = s; lds_f = lds; k_ptr = k; ldk_f = ldk;
 
     if (row_major) {
         // Set Fortran LDs (rows)
         lda_f = MAX(1, n); ldc_f = MAX(1, l); ldb_f = MAX(1, n); ldd_f = MAX(1, l);
         ldq_f = MAX(1, n); ldry_f = MAX(1, l); lds_f = MAX(1, n); ldk_f = MAX(1, n);
 
         // Point to _cm buffers (or NULL if size is 0)
         // Point to copy if it's input OR output (Fortran needs col-major buffer)
         a_ptr = (a_size > 0 && (is_A_input || is_A_output)) ? a_cm : NULL;
         c_ptr = (c_size > 0 && (is_C_input || is_C_output)) ? c_cm : NULL;
         b_ptr = (b_size > 0 && is_B_output) ? b_cm : NULL;
         d_ptr = (d_size > 0 && is_D_output) ? d_cm : NULL;
         q_ptr = (q_size > 0 && is_Q_output) ? q_cm : NULL;
         ry_ptr = (ry_size > 0 && is_RY_output) ? ry_cm : NULL;
         s_ptr = (s_size > 0 && is_S_output) ? s_cm : NULL;
         k_ptr = (k_size > 0 && is_K_output) ? k_cm : NULL;
 
         // Transpose inputs from C row-major to _cm column-major
         if (is_A_input && a_size > 0 && a_cm != NULL) slicot_transpose_to_fortran(a, a_cm, n, n, sizeof(double));
         if (is_C_input && c_size > 0 && c_cm != NULL) slicot_transpose_to_fortran(c, c_cm, l, n, sizeof(double));
 
     } else {
         // Column-major: ensure NULL pointers if size is 0
         if (a_size == 0) a_ptr = NULL; if (c_size == 0) c_ptr = NULL; if (b_size == 0) b_ptr = NULL;
         if (d_size == 0) d_ptr = NULL; if (q_size == 0) q_ptr = NULL; if (ry_size == 0) ry_ptr = NULL;
         if (s_size == 0) s_ptr = NULL; if (k_size == 0) k_ptr = NULL;
     }
 
     // --- 6. Call Fortran Routine ---
     F77_FUNC(ib01bd, IB01BD)(
         &meth_upper, &job_upper, &jobck_upper, &nobr, &n, &m, &l, &nsmpl,
         r, &ldr, // R is always column-major input
         a_ptr, &lda_f, c_ptr, &ldc_f, b_ptr, &ldb_f, d_ptr, &ldd_f,
         q_ptr, &ldq_f, ry_ptr, &ldry_f, s_ptr, &lds_f, k_ptr, &ldk_f,
         &tol, iwork, dwork, &ldwork, bwork,
         &iwarn_local, &info,
         meth_len, job_len, jobck_len
     );
 
     // --- 7. Convert Results Back (if row_major and successful) ---
     if (row_major && info == 0) {
         // Transpose output matrices from _cm buffer back to original C array
         if (is_A_output && a_size > 0 && a_cm != NULL) slicot_transpose_to_c_with_ld(a_cm, a, n, n, lda_f, lda, sizeof(double));
         if (is_C_output && c_size > 0 && c_cm != NULL) slicot_transpose_to_c_with_ld(c_cm, c, l, n, ldc_f, ldc, sizeof(double));
         if (is_B_output && b_size > 0 && b_cm != NULL) slicot_transpose_to_c_with_ld(b_cm, b, n, m, ldb_f, ldb, sizeof(double));
         if (is_D_output && d_size > 0 && d_cm != NULL) slicot_transpose_to_c_with_ld(d_cm, d, l, m, ldd_f, ldd, sizeof(double));
         if (is_Q_output && q_size > 0 && q_cm != NULL) slicot_transpose_to_c_with_ld(q_cm, q, n, n, ldq_f, ldq, sizeof(double));
         if (is_RY_output && ry_size > 0 && ry_cm != NULL) slicot_transpose_to_c_with_ld(ry_cm, ry, l, l, ldry_f, ldry, sizeof(double));
         if (is_S_output && s_size > 0 && s_cm != NULL) slicot_transpose_to_c_with_ld(s_cm, s, n, l, lds_f, lds, sizeof(double));
         if (is_K_output && k_size > 0 && k_cm != NULL) slicot_transpose_to_c_with_ld(k_cm, k, n, l, ldk_f, ldk, sizeof(double));
     }
 
     // Update output warning flag if pointer provided
     if (iwarn != NULL) {
         *iwarn = iwarn_local;
     }
 
 cleanup:
     /* --- 8. Cleanup --- */
     free(iwork);
     free(dwork);
     free(bwork); // Safe to free NULL if lbwork was 0
     free(a_cm); free(c_cm); free(b_cm); free(d_cm);
     free(q_cm); free(ry_cm); free(s_cm); free(k_cm);
 
     if (info == SLICOT_MEMORY_ERROR) {
        fprintf(stderr, "Error: Memory allocation failed in slicot_ib01bd.\n");
     }
     return info;
 }
 