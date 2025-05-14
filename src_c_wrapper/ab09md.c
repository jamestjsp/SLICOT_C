/**
 * @file ab09md.c
 * @brief C wrapper implementation for SLICOT routine AB09MD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB09MD,
 * which computes a reduced order model (Ar,Br,Cr) for the ALPHA-stable
 * part of an original state-space representation (A,B,C) using
 * Balance & Truncate methods.
 */

 #include <stdlib.h> // For malloc, free
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t
 #include <math.h>   // For MAX/MIN if needed
 
 #include "ab09md.h"
 #include "slicot_utils.h" 
 #include "slicot_f77.h"   
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  */
 extern void F77_FUNC(ab09md, AB09MD)(
     const char* dico,      
     const char* job,       
     const char* equil,     
     const char* ordsel,    
     const int* n,          
     const int* m,          
     const int* p,          
     int* nr,               
     const double* alpha,   
     double* a,             
     const int* lda,        
     double* b,             
     const int* ldb,        
     double* c,             
     const int* ldc,        
     int* ns,               
     double* hsv,           
     const double* tol,     
     int* iwork,            
     double* dwork,         
     const int* ldwork,     
     int* iwarn,            
     int* info,             
     int dico_len,          
     int job_len,           
     int equil_len,         
     int ordsel_len         
 );
 
 
 SLICOT_EXPORT
 int slicot_ab09md(char dico_in, char job_in, char equil_in, char ordsel_in,
                   int n_in, int m_in, int p_in, int* nr_io, double alpha_in,
                   double* a_io, int lda_in,
                   double* b_io, int ldb_in,
                   double* c_io, int ldc_in,
                   int* ns_out, double* hsv_out, double tol_in, int* iwarn_out,
                   int row_major)
 {
     int info = 0;
     int* iwork_allocated_buffer = NULL;
     double* dwork_allocated_buffer = NULL;
     int liwork_actual_size = 0;
     int ldwork_for_query = -1;
     int ldwork_actual = 0;
     double dwork_val_from_query_arr[1]; 
 
     const int dico_len = 1, job_len = 1, equil_len = 1, ordsel_len = 1;
 
     char dico_upper = toupper(dico_in);
     char job_upper = toupper(job_in);
     char equil_upper = toupper(equil_in);
     char ordsel_upper = toupper(ordsel_in);
 
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL;
     double *a_ptr, *b_ptr, *c_ptr;
     double *hsv_ptr;
     int lda_f, ldb_f, ldc_f;
 
     /* --- Input Parameter Validation --- */
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (job_upper != 'B' && job_upper != 'N') { info = -2; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -3; goto cleanup; }
     if (ordsel_upper != 'F' && ordsel_upper != 'A') { info = -4; goto cleanup; }
     if (n_in < 0) { info = -5; goto cleanup; }
     if (m_in < 0) { info = -6; goto cleanup; }
     if (p_in < 0) { info = -7; goto cleanup; }
     if (nr_io == NULL) { info = -8; goto cleanup; }
     if (ordsel_upper == 'F' && (*nr_io < 0 || *nr_io > n_in) ) { info = -8; goto cleanup; }
     if (dico_upper == 'C' && alpha_in > 0.0) { info = -9; goto cleanup; }
     if (dico_upper == 'D' && (alpha_in < 0.0 || alpha_in > 1.0)) { info = -9; goto cleanup; }
 
     if (n_in > 0 && a_io == NULL) { info = -10; goto cleanup; }
     if (n_in > 0 && m_in > 0 && b_io == NULL) { info = -12; goto cleanup; }
     if (p_in > 0 && n_in > 0 && c_io == NULL) { info = -14; goto cleanup; }
     if (ns_out == NULL) { info = -16; goto cleanup; }
     if (n_in > 0 && hsv_out == NULL) { info = -17; goto cleanup; }
     if (iwarn_out == NULL) { info = -22; goto cleanup; }
 
     int min_lda_f_val = MAX(1, n_in);
     int min_ldb_f_val = MAX(1, n_in);
     int min_ldc_f_val = MAX(1, p_in);
 
     if (row_major) {
         if (n_in > 0 && lda_in < n_in) { info = -11; goto cleanup; } // C LDA (cols) >= N
         if (n_in > 0 && m_in > 0 && ldb_in < m_in) { info = -13; goto cleanup; } // C LDB (cols) >= M
         if (p_in > 0 && n_in > 0 && ldc_in < n_in) { info = -15; goto cleanup; } // C LDC (cols) >= N
         else if (p_in > 0 && n_in == 0 && ldc_in < 1) {info = -15; goto cleanup;} // C LDC must be >=1 if N=0, P>0
     } else { // Column-major
         if (n_in > 0 && lda_in < min_lda_f_val) { info = -11; goto cleanup; } // C LDA (rows) >= N
         if (n_in > 0 && ldb_in < min_ldb_f_val) { info = -13; goto cleanup; } // C LDB (rows) >= N
         if (p_in > 0 && ldc_in < min_ldc_f_val) { info = -15; goto cleanup; } // C LDC (rows) >= P
     }
 
     /* --- Workspace Allocation --- */
     if (job_upper == 'N') {
         if (n_in > 0) {
             liwork_actual_size = n_in;
             iwork_allocated_buffer = (int*)malloc((size_t)liwork_actual_size * sizeof(int));
             CHECK_ALLOC(iwork_allocated_buffer);
         } else { 
             liwork_actual_size = 0; 
             iwork_allocated_buffer = NULL;
         }
     } else { 
         liwork_actual_size = 0; 
         iwork_allocated_buffer = NULL;
     }
 
     double dummy_double_for_query = 0.0;
     int dummy_nr_for_query = (ordsel_upper == 'F' ? *nr_io : 0);
     int dummy_ns_for_query = 0;
 
     double* q_a_ptr = (n_in > 0) ? a_io : &dummy_double_for_query;
     int q_lda_f = (n_in > 0) ? (row_major ? MAX(1, n_in) : lda_in) : 1;
     double* q_b_ptr = (n_in > 0 && m_in > 0) ? b_io : &dummy_double_for_query;
     int q_ldb_f = (n_in > 0 && m_in > 0) ? (row_major ? MAX(1, n_in) : ldb_in) : 1;
     double* q_c_ptr = (p_in > 0 && n_in > 0) ? c_io : &dummy_double_for_query;
     int q_ldc_f = (p_in > 0 && n_in > 0) ? (row_major ? MAX(1, p_in) : ldc_in) : ( (p_in > 0) ? MAX(1,p_in) : 1 );
     double* q_hsv_ptr = (n_in > 0) ? hsv_out : &dummy_double_for_query; 
     
     int query_info_val = 0;
     F77_FUNC(ab09md, AB09MD)(&dico_upper, &job_upper, &equil_upper, &ordsel_upper,
                               &n_in, &m_in, &p_in, &dummy_nr_for_query, &alpha_in,
                               q_a_ptr, &q_lda_f, q_b_ptr, &q_ldb_f, q_c_ptr, &q_ldc_f,
                               &dummy_ns_for_query, q_hsv_ptr, &tol_in, 
                               iwork_allocated_buffer, 
                               dwork_val_from_query_arr, &ldwork_for_query, 
                               iwarn_out, &query_info_val,
                               dico_len, job_len, equil_len, ordsel_len);
 
     int min_ldwork_formula = 1;
     if (n_in > 0) {
          min_ldwork_formula = MAX(1, n_in * (2 * n_in + MAX(n_in, MAX(m_in, p_in)) + 5) + (n_in * (n_in + 1)) / 2);
     }
 
     if (query_info_val == 0 || query_info_val == -21) { // Parameter 21 is LDWORK for AB09MD
         ldwork_actual = (int)dwork_val_from_query_arr[0];
         ldwork_actual = MAX(ldwork_actual, min_ldwork_formula); 
     } else { 
         info = query_info_val; 
         goto cleanup; 
     }
     ldwork_actual = MAX(1, ldwork_actual); 
     
     dwork_allocated_buffer = (double*)malloc((size_t)ldwork_actual * sizeof(double));
     CHECK_ALLOC(dwork_allocated_buffer);
 
     /* --- Prepare Arrays for Computational Call --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         size_t a_size = (size_t)n_in * n_in;
         size_t b_size = (size_t)n_in * m_in;
         size_t c_size = (size_t)p_in * n_in;
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
 
         if (a_size > 0) slicot_transpose_to_fortran_with_ld(a_io, a_cm, n_in, n_in, lda_in, MAX(1,n_in), elem_size);
         if (b_size > 0) slicot_transpose_to_fortran_with_ld(b_io, b_cm, n_in, m_in, ldb_in, MAX(1,n_in), elem_size);
         if (c_size > 0) slicot_transpose_to_fortran_with_ld(c_io, c_cm, p_in, n_in, ldc_in, MAX(1,p_in), elem_size);
 
         lda_f = MAX(1, n_in); ldb_f = MAX(1, n_in); ldc_f = MAX(1, p_in);
 
         a_ptr = (n_in > 0) ? a_cm : NULL; 
         b_ptr = (n_in > 0 && m_in > 0) ? b_cm : NULL; 
         c_ptr = (p_in > 0 && n_in > 0) ? c_cm : NULL;
     } else {
         lda_f = lda_in; ldb_f = ldb_in; ldc_f = ldc_in;
         a_ptr = (n_in > 0) ? a_io : NULL;
         b_ptr = (n_in > 0 && m_in > 0) ? b_io : NULL;
         c_ptr = (p_in > 0 && n_in > 0) ? c_io : NULL;
     }
     hsv_ptr = (n_in > 0) ? hsv_out : NULL;
     
     // Reset info to 0 before computational call, unless CHECK_ALLOC failed.
     // The 'info' variable for the wrapper will be updated by the Fortran call.
     if (info != SLICOT_MEMORY_ERROR) info = 0; 
 
     F77_FUNC(ab09md, AB09MD)(&dico_upper, &job_upper, &equil_upper, &ordsel_upper,
                               &n_in, &m_in, &p_in, nr_io, &alpha_in,
                               a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f,
                               ns_out, hsv_ptr, &tol_in, 
                               iwork_allocated_buffer, dwork_allocated_buffer, &ldwork_actual,
                               iwarn_out, &info,
                               dico_len, job_len, equil_len, ordsel_len);
 
     if (row_major && info == 0) {
         int nr_val = *nr_io;
         if (nr_val >= 0) { 
             if (n_in > 0) { // Check if original dimension was > 0 before transposing
                 // A is overwritten with Ar (nr_val x nr_val)
                 if (a_cm != NULL && nr_val > 0) slicot_transpose_to_c_with_ld(a_cm, a_io, nr_val, nr_val, lda_f, lda_in, elem_size);
                 // B is overwritten with Br (nr_val x m_in)
                 if (m_in > 0 && b_cm != NULL && nr_val > 0) slicot_transpose_to_c_with_ld(b_cm, b_io, nr_val, m_in, ldb_f, ldb_in, elem_size);
                 // C is overwritten with Cr (p_in x nr_val)
                 if (p_in > 0 && c_cm != NULL && nr_val > 0) slicot_transpose_to_c_with_ld(c_cm, c_io, p_in, nr_val, ldc_f, ldc_in, elem_size);
             }
         }
     }
 
 cleanup:
     free(dwork_allocated_buffer);
     free(iwork_allocated_buffer); 
     free(a_cm);
     free(b_cm);
     free(c_cm);
 
     return info;
 }
