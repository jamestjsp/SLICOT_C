/**
 * @file sb02mt.c
 * @brief C wrapper implementation for SLICOT routine SB02MT
 *
 * This file provides a C wrapper implementation for the SLICOT routine SB02MT,
 * which converts linear-quadratic optimal control problems with
 * coupling weighting terms (matrix L) to standard problems by
 * computing modified matrices A, Q and optionally G = B*inv(R)*B'.
 */

 #include <stdlib.h>
 #include <ctype.h>
 #include <stddef.h> // For size_t
 #include <string.h> // For memcpy

 // Include the header file for this wrapper
 #include "sb02mt.h"
 // Include necessary SLICOT utility headers
 #include "slicot_utils.h" // Assumed to contain MAX, CHECK_ALLOC, SLICOT_MEMORY_ERROR, transpose routines, slicot_copy_symmetric_part, slicot_transpose_symmetric_to_fortran, slicot_transpose_symmetric_to_c
 #include "slicot_f77.h"   // For F77_FUNC macro and Fortran interface conventions

 /*
  * Declare the external Fortran routine using the F77_FUNC macro.
  * Arguments A, B, Q, R, L, IPIV, G, OUFACT can be modified depending on flags.
  */
 extern void F77_FUNC(sb02mt, SB02MT)(
     const char* jobg,       // CHARACTER*1 JOBG
     const char* jobl,       // CHARACTER*1 JOBL
     const char* fact,       // CHARACTER*1 FACT
     const char* uplo,       // CHARACTER*1 UPLO
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     double* a,              // DOUBLE PRECISION A(LDA,*) (in/out)
     const int* lda,         // INTEGER LDA
     double* b,              // DOUBLE PRECISION B(LDB,*) (in/out)
     const int* ldb,         // INTEGER LDB
     double* q,              // DOUBLE PRECISION Q(LDQ,*) (in/out)
     const int* ldq,         // INTEGER LDQ
     double* r,              // DOUBLE PRECISION R(LDR,*) (in/out)
     const int* ldr,         // INTEGER LDR
     double* l,              // DOUBLE PRECISION L(LDL,*) (in/out)
     const int* ldl,         // INTEGER LDL
     int* ipiv,              // INTEGER IPIV(*) (in/out if FACT='N')
     int* oufact,            // INTEGER OUFACT (output)
     double* g,              // DOUBLE PRECISION G(LDG,*) (output if JOBG='G')
     const int* ldg,         // INTEGER LDG
     int* iwork,             // INTEGER IWORK(*) (workspace if FACT='N')
     double* dwork,          // DOUBLE PRECISION DWORK(*)
     const int* ldwork,      // INTEGER LDWORK
     int* info,              // INTEGER INFO (output)
     int jobg_len,           // Hidden length
     int jobl_len,           // Hidden length
     int fact_len,           // Hidden length
     int uplo_len            // Hidden length
 );


 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_sb02mt(char jobg, char jobl, char fact, char uplo,
                   int n, int m,
                   double* a, int lda, double* b, int ldb,
                   double* q, int ldq, double* r, int ldr,
                   double* l, int ldl, int* ipiv, int* oufact,
                   double* g, int ldg, int row_major)
 {
     /* Local variables */
     int info = 0;
     int ldwork = -1; /* Use -1 for workspace query */
     double dwork_query;
     double* dwork = NULL;
     int* iwork = NULL; // Workspace for FACT='N'
     int iwork_size = 0;
     const int jobg_len = 1, jobl_len = 1, fact_len = 1, uplo_len = 1;

     char jobg_upper = toupper(jobg);
     char jobl_upper = toupper(jobl);
     char fact_upper = toupper(fact);
     char uplo_upper = toupper(uplo);

     /* Pointers for column-major copies if needed */
     double *a_cm = NULL, *b_cm = NULL, *q_cm = NULL;
     double *r_cm = NULL, *l_cm = NULL, *g_cm = NULL;
     double *a_ptr, *b_ptr, *q_ptr, *r_ptr, *l_ptr, *g_ptr;
     int lda_f, ldb_f, ldq_f, ldr_f, ldl_f, ldg_f;

     /* --- Input Parameter Validation --- */

     if (jobg_upper != 'G' && jobg_upper != 'N') { info = -1; goto cleanup; }
     if (jobl_upper != 'Z' && jobl_upper != 'N') { info = -2; goto cleanup; }
     if (fact_upper != 'N' && fact_upper != 'C' && fact_upper != 'U') { info = -3; goto cleanup; }
     if (uplo_upper != 'U' && uplo_upper != 'L') { info = -4; goto cleanup; }
     if (n < 0) { info = -5; goto cleanup; }
     if (m < 0) { info = -6; goto cleanup; }

     // Check leading dimensions based on storage order and flags
     int min_lda_f = (jobl_upper == 'N') ? MAX(1, n) : 1; // A not referenced if JOBL='Z'
     int min_ldb_f = MAX(1, n);
     int min_ldq_f = (jobl_upper == 'N') ? MAX(1, n) : 1; // Q not referenced if JOBL='Z'
     int min_ldr_f = MAX(1, m);
     int min_ldl_f = (jobl_upper == 'N') ? MAX(1, n) : 1; // L not referenced if JOBL='Z'
     int min_ldg_f = (jobg_upper == 'G') ? MAX(1, n) : 1; // G only computed if JOBG='G'

     if (row_major) {
         // For row-major C, LD is the number of columns
         int min_lda_rm_cols = (jobl_upper == 'N') ? n : 1;
         int min_ldb_rm_cols = m;
         int min_ldq_rm_cols = (jobl_upper == 'N') ? n : 1;
         int min_ldr_rm_cols = m;
         int min_ldl_rm_cols = (jobl_upper == 'N') ? m : 1;
         int min_ldg_rm_cols = (jobg_upper == 'G') ? n : 1;

         if (lda < min_lda_rm_cols) { info = -8; goto cleanup; }
         if (ldb < min_ldb_rm_cols) { info = -10; goto cleanup; }
         if (ldq < min_ldq_rm_cols) { info = -12; goto cleanup; }
         if (ldr < min_ldr_rm_cols) { info = -14; goto cleanup; }
         if (ldl < min_ldl_rm_cols) { info = -16; goto cleanup; }
         if (ldg < min_ldg_rm_cols) { info = -20; goto cleanup; }
     } else {
         // For column-major C, LD is the number of rows (Fortran style)
         if (lda < min_lda_f) { info = -8; goto cleanup; }
         if (ldb < min_ldb_f) { info = -10; goto cleanup; }
         if (ldq < min_ldq_f) { info = -12; goto cleanup; }
         if (ldr < min_ldr_f) { info = -14; goto cleanup; }
         if (ldl < min_ldl_f) { info = -16; goto cleanup; }
         if (ldg < min_ldg_f) { info = -20; goto cleanup; }
     }

     /* --- Prepare arrays for column-major format if using row-major --- */
     size_t elem_size = sizeof(double);
     if (row_major) {
         /* Allocate memory for column-major copies */
         size_t a_rows=n, a_cols=n, a_size = (size_t)a_rows*a_cols;
         size_t b_rows=n, b_cols=m, b_size = (size_t)b_rows*b_cols;
         size_t q_rows=n, q_cols=n, q_size = (size_t)q_rows*q_cols;
         size_t r_rows=m, r_cols=m, r_size = (size_t)r_rows*r_cols;
         size_t l_rows=n, l_cols=m, l_size = (size_t)l_rows*l_cols;
         size_t g_rows=n, g_cols=n, g_size = (size_t)g_rows*g_cols;

         // Allocate based on whether they are referenced and/or modified
         if (jobl_upper == 'N' && a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); } // B always in/out? Docs say modified if oufact=1.
         if (jobl_upper == 'N' && q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); }
         if (r_size > 0) { r_cm = (double*)malloc(r_size * elem_size); CHECK_ALLOC(r_cm); } // R always in/out? Docs say modified if oufact=1 or 2.
         if (jobl_upper == 'N' && l_size > 0) { l_cm = (double*)malloc(l_size * elem_size); CHECK_ALLOC(l_cm); } // L modified if oufact=1.
         if (jobg_upper == 'G' && g_size > 0) { g_cm = (double*)malloc(g_size * elem_size); CHECK_ALLOC(g_cm); } // G is output

         /* Transpose C inputs to Fortran copies */
         if (a_cm) slicot_transpose_to_fortran(a, a_cm, a_rows, a_cols, elem_size);
         if (b_cm) slicot_transpose_to_fortran(b, b_cm, b_rows, b_cols, elem_size);
         if (q_cm) slicot_transpose_symmetric_to_fortran(q, q_cm, n, uplo_upper, elem_size);
         if (r_cm) slicot_transpose_symmetric_to_fortran(r, r_cm, m, uplo_upper, elem_size); // R is symmetric
         if (l_cm) slicot_transpose_to_fortran(l, l_cm, l_rows, l_cols, elem_size);
         // G is output only, no input transpose needed.

         /* Fortran leading dimensions */
         lda_f = (a_rows > 0) ? a_rows : 1;
         ldb_f = (b_rows > 0) ? b_rows : 1;
         ldq_f = (q_rows > 0) ? q_rows : 1;
         ldr_f = (r_rows > 0) ? r_rows : 1;
         ldl_f = (l_rows > 0) ? l_rows : 1;
         ldg_f = (g_rows > 0) ? g_rows : 1;

         /* Set pointers */
         a_ptr = a_cm ? a_cm : a; // Use original if not referenced (JOBL='Z')
         b_ptr = b_cm;
         q_ptr = q_cm ? q_cm : q; // Use original if not referenced (JOBL='Z')
         r_ptr = r_cm;
         l_ptr = l_cm ? l_cm : l; // Use original if not referenced (JOBL='Z')
         g_ptr = g_cm ? g_cm : g; // Use original if not computed (JOBG='N')

     } else {
         /* Column-major case - use original arrays, potentially copy symmetric R, Q */
         lda_f = lda;
         ldb_f = ldb;
         ldq_f = ldq;
         ldr_f = ldr;
         ldl_f = ldl;
         ldg_f = ldg;
         a_ptr = a;
         b_ptr = b;
         q_ptr = q;
         r_ptr = r;
         l_ptr = l;
         g_ptr = g;

         // Need to copy symmetric inputs R, Q if they are modified and stored symmetrically
         // R is modified if oufact=1 or 2. Q is modified if JOBL='N'.
         size_t r_rows=m, r_cols=m, r_size = (size_t)r_rows*r_cols;
         size_t q_rows=n, q_cols=n, q_size = (size_t)q_rows*q_cols;
         // Copy R always, as it's potentially modified (stores factor).
         if (r_size > 0) { r_cm = (double*)malloc(r_size * elem_size); CHECK_ALLOC(r_cm); slicot_copy_symmetric_part(r, r_cm, m, uplo_upper, ldr, elem_size); r_ptr = r_cm; ldr_f = m;}
         // Copy Q only if JOBL='N' (it's modified).
         if (jobl_upper == 'N' && q_size > 0) { q_cm = (double*)malloc(q_size * elem_size); CHECK_ALLOC(q_cm); slicot_copy_symmetric_part(q, q_cm, n, uplo_upper, ldq, elem_size); q_ptr = q_cm; ldq_f = n;}
     }

     /* --- Workspace allocation --- */

     // Allocate IWORK (size M if FACT='N', else not referenced)
     if (fact_upper == 'N') {
         iwork_size = MAX(1, m); // Ensure size >= 1
         iwork = (int*)malloc((size_t)iwork_size * sizeof(int));
         CHECK_ALLOC(iwork);
     } else {
         iwork = NULL; // Not referenced
     }

     // Perform workspace query with small dwork array
     double dwork_temp[2]; // Small array to receive query result
     ldwork = -1; // Query mode
     F77_FUNC(sb02mt, SB02MT)(&jobg_upper, &jobl_upper, &fact_upper, &uplo_upper,
                              &n, &m, a_ptr, &lda_f, b_ptr, &ldb_f, q_ptr, &ldq_f,
                              r_ptr, &ldr_f, l_ptr, &ldl_f,
                              ipiv, oufact, g_ptr, &ldg_f,
                              iwork, dwork_temp, &ldwork,
                              &info,
                              jobg_len, jobl_len, fact_len, uplo_len);

     if (info < 0 && info != -23) { goto cleanup; } // Query failed due to invalid argument (allow INFO=-23 from query)
     info = 0; // Reset info after query

     // Get the required dwork size from query result
     ldwork = (int)dwork_temp[0]; // First element contains optimal ldwork
     // Check against minimum documented size
     int min_ldwork = 1;
     if (fact_upper == 'N' && (jobg_upper == 'G' || jobl_upper == 'N')) {
         min_ldwork = MAX(2, MAX(3 * m, n * m));
     } else if (fact_upper == 'N') { // JOBG='N' and JOBL='Z'
         min_ldwork = MAX(2, 3 * m);
     } else if (fact_upper == 'U' && (jobg_upper == 'G' || jobl_upper == 'N')) {
         min_ldwork = MAX(1, n * m);
     } // else min_ldwork = 1 (other FACT cases)

     ldwork = MAX(ldwork, min_ldwork);

     dwork = (double*)malloc((size_t)ldwork * sizeof(double));
     CHECK_ALLOC(dwork); // Sets info and jumps to cleanup on failure

     /* --- Call the computational routine --- */
     F77_FUNC(sb02mt, SB02MT)(&jobg_upper, &jobl_upper, &fact_upper, &uplo_upper,
                              &n, &m,
                              a_ptr, &lda_f, // A
                              b_ptr, &ldb_f, // B
                              q_ptr, &ldq_f, // Q
                              r_ptr, &ldr_f, // R
                              l_ptr, &ldl_f, // L
                              ipiv, oufact,
                              g_ptr, &ldg_f, // G
                              iwork, dwork, &ldwork, &info,
                              jobg_len, jobl_len, fact_len, uplo_len);

     /* --- Copy results back to row-major format if needed --- */
     if (row_major && (info == 0 || info == (m + 1))) {
         size_t a_rows=n, a_cols=n, a_size = (size_t)a_rows*a_cols;
         size_t b_rows=n, b_cols=m, b_size = (size_t)b_rows*b_cols;
         size_t q_rows=n, q_cols=n, q_size = (size_t)q_rows*q_cols;
         size_t r_rows=m, r_cols=m, r_size = (size_t)r_rows*r_cols;
         size_t l_rows=n, l_cols=m, l_size = (size_t)l_rows*l_cols;
         size_t g_rows=n, g_cols=n, g_size = (size_t)g_rows*g_cols;

         // Copy back based on modification flags (oufact) and job flags
         if (jobl_upper == 'N' && a_cm) slicot_transpose_to_c(a_cm, a, a_rows, a_cols, elem_size);
         if (*oufact == 1 && b_cm) slicot_transpose_to_c(b_cm, b, b_rows, b_cols, elem_size); // B modified if oufact=1
         if (jobl_upper == 'N' && q_cm) slicot_transpose_symmetric_to_c(q_cm, q, n, uplo_upper, elem_size);
         if ((*oufact == 1 || *oufact == 2) && r_cm) slicot_transpose_symmetric_to_c(r_cm, r, m, uplo_upper, elem_size); // R modified if oufact=1 or 2
         if (jobl_upper == 'N' && *oufact == 1 && l_cm) slicot_transpose_to_c(l_cm, l, l_rows, l_cols, elem_size); // L modified if oufact=1
         if (jobg_upper == 'G' && g_cm) slicot_transpose_symmetric_to_c(g_cm, g, n, uplo_upper, elem_size); // G is output
         // IPIV, OUFACT modified directly
     } else if (!row_major && (info == 0 || info == (m + 1))) {
         // Copy back results from temporary full copies (q_cm, r_cm) to original symmetric storage (q, r)
         size_t q_rows=n, q_cols=n, q_size = (size_t)q_rows*q_cols;
         size_t r_rows=m, r_cols=m, r_size = (size_t)r_rows*r_cols;
         if (jobl_upper == 'N' && q_cm) { // Check if q_cm was allocated
             slicot_copy_symmetric_part(q_cm, q, n, uplo_upper, ldq, elem_size);
         }
         if (r_cm) { // Check if r_cm was allocated
             slicot_copy_symmetric_part(r_cm, r, m, uplo_upper, ldr, elem_size);
         }
         // A, B, L, IPIV, OUFACT, G modified directly in original arrays.
     }

 cleanup:
     /* --- Cleanup --- */
     free(dwork);
     free(iwork); // Safe even if NULL
     free(a_cm);
     free(b_cm);
     free(q_cm); // Free Q copy if allocated
     free(r_cm); // Free R copy if allocated
     free(l_cm);
     free(g_cm);

     return info;
 }