/**
 * @file ab09bd.c
 * @brief C wrapper implementation for SLICOT routine AB09BD
 *
 * This file provides a C wrapper implementation for the SLICOT routine AB09BD,
 * which computes a reduced order model (Ar,Br,Cr,Dr) for a stable
 * original state-space representation (A,B,C,D) using Singular
 * Perturbation Approximation (SPA) methods.
 * Workspace (IWORK, DWORK) is allocated internally by this wrapper.
 */

 #include <stdlib.h> // For malloc, free
 #include <ctype.h>  // For toupper
 #include <stddef.h> // For size_t
 #include <math.h>   // For MAX/MIN if needed
 
 // Include the header file for this wrapper
 #include "ab09bd.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX, MIN, transpose functions
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions
 
 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  */
 extern void F77_FUNC(ab09bd, AB09BD)(
     const char* dico,       // CHARACTER*1 DICO
     const char* job,        // CHARACTER*1 JOB
     const char* equil,      // CHARACTER*1 EQUIL
     const char* ordsel,     // CHARACTER*1 ORDSEL
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     int* nr,                // INTEGER NR (in/out)
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* c,              // DOUBLE PRECISION C(LDC,*) (in/out)
     const int* ldc,         // INTEGER LDC
     double* d,              // DOUBLE PRECISION D(LDD,*) (in/out)
     const int* ldd,         // INTEGER LDD
     double* hsv,            // DOUBLE PRECISION HSV(*) (output)
     const double* tol1,     // DOUBLE PRECISION TOL1
     const double* tol2,     // DOUBLE PRECISION TOL2
     int* iwork,             // INTEGER IWORK(*)
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* iwarn,             // INTEGER IWARN (output)
     int* info,              // INTEGER INFO (output)
     int dico_len,           // Hidden length
     int job_len,            // Hidden length
     int equil_len,          // Hidden length
     int ordsel_len          // Hidden length
 );
 
 
 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_ab09bd(char dico_in, char job_in, char equil_in, char ordsel_in,
                   int n_in, int m_in, int p_in, int* nr_io,
                   double* a_io, int lda_in,
                   double* b_io, int ldb_in,
                   double* c_io, int ldc_in,
                   double* d_io, int ldd_in,
                   double* hsv_out,
                   double tol1_in, double tol2_in,
                   int* iwarn_out,
                   int row_major)
 {
     /* Local variables */
     int info = 0; // This will hold the info from Fortran calls
     int ldwork_actual = 0;
     double* dwork_allocated_buffer = NULL;
     int* iwork_allocated_buffer = NULL;
     int liwork_actual_size = 0;
 
     const int dico_len = 1, job_len = 1, equil_len = 1, ordsel_len = 1;
 
     char dico_upper = toupper(dico_in);
     char job_upper = toupper(job_in);
     char equil_upper = toupper(equil_in);
     char ordsel_upper = toupper(ordsel_in);
 
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL;
     double *a_ptr, *b_ptr, *c_ptr, *d_ptr;
     double *hsv_ptr;
     int lda_f, ldb_f, ldc_f, ldd_f;
 
     /* --- Input Parameter Validation (as before) --- */
     if (dico_upper != 'C' && dico_upper != 'D') { info = -1; goto cleanup; }
     if (job_upper != 'B' && job_upper != 'N') { info = -2; goto cleanup; }
     if (equil_upper != 'S' && equil_upper != 'N') { info = -3; goto cleanup; }
     if (ordsel_upper != 'F' && ordsel_upper != 'A') { info = -4; goto cleanup; }
     if (n_in < 0) { info = -5; goto cleanup; }
     if (m_in < 0) { info = -6; goto cleanup; }
     if (p_in < 0) { info = -7; goto cleanup; }
     if (nr_io == NULL) { info = -8; goto cleanup; }
     if (ordsel_upper == 'F' && (*nr_io < 0 || *nr_io > n_in) ) { info = -8; goto cleanup; }
 
     if (n_in > 0 && a_io == NULL) { info = -9; goto cleanup; }
     // Parameter 10 is LDA
     if (n_in > 0 && m_in > 0 && b_io == NULL) { info = -11; goto cleanup; }
     // Parameter 12 is LDB
     if (p_in > 0 && n_in > 0 && c_io == NULL) { info = -13; goto cleanup; }
     // Parameter 14 is LDC
     if (p_in > 0 && m_in > 0 && d_io == NULL && !(n_in == 0 && p_in == 0 && m_in == 0) ) { info = -15; goto cleanup; }
     // Parameter 16 is LDD
     if (n_in > 0 && hsv_out == NULL) { info = -17; goto cleanup; }
     // TOL1 is param 18, TOL2 is param 19
     // IWORK is param 20
     // DWORK is param 21
     // LDWORK is param 22
     if (iwarn_out == NULL) { info = -23; goto cleanup; }
 
     int min_lda_f_val = MAX(1, n_in);
     int min_ldb_f_val = MAX(1, n_in);
     int min_ldc_f_val = MAX(1, p_in);
     int min_ldd_f_val = MAX(1, p_in);
 
     if (row_major) {
         if (n_in > 0 && lda_in < n_in) { info = -10; goto cleanup; }
         if (n_in > 0 && m_in > 0 && ldb_in < m_in) { info = -12; goto cleanup; }
         if (p_in > 0 && n_in > 0 && ldc_in < n_in) { info = -14; goto cleanup; }
         if (p_in > 0 && m_in > 0 && ldd_in < m_in) { info = -16; goto cleanup; }
     } else { 
         if (lda_in < min_lda_f_val) { info = -10; goto cleanup; }
         if (ldb_in < min_ldb_f_val) { info = -12; goto cleanup; }
         if (ldc_in < min_ldc_f_val) { info = -14; goto cleanup; }
         if (ldd_in < min_ldd_f_val) { info = -16; goto cleanup; }
     }
 
     /* --- Workspace Allocation --- */
     liwork_actual_size = MAX(1, 2 * n_in);
     iwork_allocated_buffer = (int*)malloc((size_t)liwork_actual_size * sizeof(int));
     CHECK_ALLOC(iwork_allocated_buffer);

     // Calculate LDWORK using the supplied formula
     if (n_in == 0) {
         ldwork_actual = 1;
     } else {
         ldwork_actual = MAX(1, n_in * (2 * n_in + MAX(n_in, MAX(m_in, p_in)) + 5) + (n_in * (n_in + 1)) / 2);
     }

     dwork_allocated_buffer = (double*)malloc((size_t)ldwork_actual * sizeof(double));
     CHECK_ALLOC(dwork_allocated_buffer);

     /* --- Prepare Arrays for Computational Call (info is 0 if we reached here without CHECK_ALLOC failing) --- */
     size_t elem_size = sizeof(double);
     // ... (rest of a_ptr, b_ptr, etc. setup as before) ...
     if (row_major) {
         size_t a_size = (size_t)n_in * n_in;
         size_t b_size = (size_t)n_in * m_in;
         size_t c_size = (size_t)p_in * n_in;
         size_t d_size = (size_t)p_in * m_in;
 
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
 
         if (a_size > 0) slicot_transpose_to_fortran_with_ld(a_io, a_cm, n_in, n_in, lda_in, MAX(1,n_in), elem_size);
         if (b_size > 0) slicot_transpose_to_fortran_with_ld(b_io, b_cm, n_in, m_in, ldb_in, MAX(1,n_in), elem_size);
         if (c_size > 0) slicot_transpose_to_fortran_with_ld(c_io, c_cm, p_in, n_in, ldc_in, MAX(1,p_in), elem_size);
         if (d_size > 0) slicot_transpose_to_fortran_with_ld(d_io, d_cm, p_in, m_in, ldd_in, MAX(1,p_in), elem_size);
 
         lda_f = MAX(1, n_in); ldb_f = MAX(1, n_in); ldc_f = MAX(1, p_in); ldd_f = MAX(1, p_in);
 
         a_ptr = (n_in > 0) ? a_cm : NULL; 
         b_ptr = (n_in > 0 && m_in > 0) ? b_cm : NULL; 
         c_ptr = (p_in > 0 && n_in > 0) ? c_cm : NULL;
         d_ptr = (p_in > 0 && m_in > 0) ? d_cm : NULL;
     } else {
         lda_f = lda_in; ldb_f = ldb_in; ldc_f = ldc_in; ldd_f = ldd_in;
         a_ptr = (n_in > 0) ? a_io : NULL;
         b_ptr = (n_in > 0 && m_in > 0) ? b_io : NULL;
         c_ptr = (p_in > 0 && n_in > 0) ? c_io : NULL;
         d_ptr = (p_in > 0 && m_in > 0) ? d_io : NULL;
     }
     hsv_ptr = (n_in > 0) ? hsv_out : NULL;
     
     // 'info' for the wrapper should be 0 before this computational call,
     // unless CHECK_ALLOC failed.
     // The 'info' variable passed to F77_FUNC will be set by the Fortran routine.
     F77_FUNC(ab09bd, AB09BD)(&dico_upper, &job_upper, &equil_upper, &ordsel_upper,
                               &n_in, &m_in, &p_in, nr_io, 
                               a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
                               hsv_ptr, &tol1_in, &tol2_in, 
                               iwork_allocated_buffer, dwork_allocated_buffer, &ldwork_actual,
                               iwarn_out, &info, // Fortran routine will set this 'info'
                               dico_len, job_len, equil_len, ordsel_len);
 
     if (row_major && info == 0) {
         int nr_val = *nr_io;
         if (nr_val >= 0) { 
             if (n_in > 0) { 
                 if (a_cm != NULL) slicot_transpose_to_c_with_ld(a_cm, a_io, nr_val, nr_val, lda_f, lda_in, elem_size);
                 if (m_in > 0 && b_cm != NULL) slicot_transpose_to_c_with_ld(b_cm, b_io, nr_val, m_in, ldb_f, ldb_in, elem_size);
                 if (p_in > 0 && c_cm != NULL) slicot_transpose_to_c_with_ld(c_cm, c_io, p_in, nr_val, ldc_f, ldc_in, elem_size);
             }
             if (p_in > 0 && m_in > 0 && d_cm != NULL) { // D is PxM
                  slicot_transpose_to_c_with_ld(d_cm, d_io, p_in, m_in, ldd_f, ldd_in, elem_size);
             }
         }
     }
 
 cleanup:
     free(dwork_allocated_buffer);
     free(iwork_allocated_buffer); 
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);
 
     return info; // Return the final info status
 }
