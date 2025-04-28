/**
 * @file tf01md.c
 * @brief C wrapper implementation for SLICOT routine TF01MD
 * @details Computes the output response sequence of a linear time-invariant
 * discrete-time system given its state-space model (A,B,C,D), initial state,
 * and input sequence. Workspace (DWORK) is allocated internally.
 * Input/output matrix format is handled via the row_major parameter.
 */

 #include <stdlib.h> // For malloc, free
 #include <stddef.h> // For size_t
 #include <math.h>   // For isnan, isinf
 #include <stdio.h>  // For fprintf, stderr
 
 #include "tf01md.h"       // Public header for this wrapper
 #include "slicot_utils.h" // Provides CHECK_ALLOC, SLICOT_MEMORY_ERROR, MAX/MIN, transpose functions etc.
 #include "slicot_f77.h"   // Provides F77_FUNC macro for Fortran name mangling
 
 /* External Fortran routine declaration */
 extern void F77_FUNC(tf01md, TF01MD)(
     const int* n,           // INTEGER N
     const int* m,           // INTEGER M
     const int* p,           // INTEGER P
     const int* ny,          // INTEGER NY
     const double* a,        // DOUBLE PRECISION A(LDA,*) (in)
     const int* lda,         // INTEGER LDA
     const double* b,        // DOUBLE PRECISION B(LDB,*) (in)
     const int* ldb,         // INTEGER LDB
     const double* c,        // DOUBLE PRECISION C(LDC,*) (in)
     const int* ldc,         // INTEGER LDC
     const double* d,        // DOUBLE PRECISION D(LDD,*) (in)
     const int* ldd,         // INTEGER LDD
     const double* u,        // DOUBLE PRECISION U(LDU,*) (in)
     const int* ldu,         // INTEGER LDU
     double* x,              // DOUBLE PRECISION X(*) (in/out)
     double* y,              // DOUBLE PRECISION Y(LDY,*) (output)
     const int* ldy,         // INTEGER LDY
     double* dwork,          // DOUBLE PRECISION DWORK(*) (workspace)
     int* info               // INTEGER INFO (output)
 );
 
 /* C wrapper function definition */
 SLICOT_EXPORT
 int slicot_tf01md(int n, int m, int p, int ny,
                   const double* a, int lda, const double* b, int ldb,
                   const double* c, int ldc, const double* d, int ldd,
                   const double* u, int ldu,
                   double* x, double* y, int ldy, int row_major)
 {
     // 1. Variable declarations
     int info = 0;          // Return status code
     double *dwork = NULL;  // Internal double workspace pointer
     int ldwork = 0;        // Calculated size for dwork (Fortran requires pointer to size)
 
     // Pointers for column-major copies (if row_major is used)
     double *a_cm = NULL, *b_cm = NULL, *c_cm = NULL, *d_cm = NULL, *u_cm = NULL;
     double *y_cm = NULL; // Temporary buffer for output Y if row_major
 
     // Pointers to pass to Fortran
     const double *a_ptr, *b_ptr, *c_ptr, *d_ptr, *u_ptr;
     double *y_ptr; // Pointer to where Fortran should write Y output
 
     // Fortran-style leading dimensions (number of rows)
     int lda_f, ldb_f, ldc_f, ldd_f, ldu_f, ldy_f;
 
     // Size constants
     const size_t elem_size = sizeof(double);
 
     // 2. Input parameter validation (Check BEFORE allocating memory)
     if (n < 0) { info = -1; goto cleanup; }
     if (m < 0) { info = -2; goto cleanup; }
     if (p < 0) { info = -3; goto cleanup; }
     if (ny < 0) { info = -4; goto cleanup; }
 
     // Check pointers for required arrays (handle optional arrays carefully based on dimensions)
     if (a == NULL && n > 0) { info = -5; goto cleanup; }
     if (b == NULL && n > 0 && m > 0) { info = -7; goto cleanup; }
     if (c == NULL && p > 0 && n > 0) { info = -9; goto cleanup; }
     if (d == NULL && p > 0 && m > 0) { info = -11; goto cleanup; }
     if (u == NULL && m > 0 && ny > 0) { info = -13; goto cleanup; }
     if (x == NULL && n > 0) { info = -15; goto cleanup; } // X is required if N > 0
     if (y == NULL && p > 0 && ny > 0) { info = -16; goto cleanup; } // Y is required if P > 0 and NY > 0
 
     // Check leading dimensions based on row_major flag
     int min_lda_f = MAX(1, n); // Min Fortran LDA (rows)
     int min_ldb_f = MAX(1, n); // Min Fortran LDB (rows)
     int min_ldc_f = MAX(1, p); // Min Fortran LDC (rows)
     int min_ldd_f = MAX(1, p); // Min Fortran LDD (rows)
     int min_ldu_f = MAX(1, m); // Min Fortran LDU (rows)
     int min_ldy_f = MAX(1, p); // Min Fortran LDY (rows)
 
     if (row_major) {
         // C LDA is number of COLUMNS
         if (n > 0 && lda < n) { info = -6; goto cleanup; } // Check min columns for A
         if (m > 0 && ldb < m) { info = -8; goto cleanup; } // Check min columns for B
         if (n > 0 && ldc < n) { info = -10; goto cleanup; } // Check min columns for C
         if (m > 0 && ldd < m) { info = -12; goto cleanup; } // Check min columns for D
         if (ny > 0 && ldu < ny) { info = -14; goto cleanup; } // Check min columns for U
         if (ny > 0 && ldy < ny) { info = -17; goto cleanup; } // Check min columns for Y
     } else {
         // C LDA is number of ROWS (Fortran style)
         if (lda < min_lda_f) { info = -6; goto cleanup; } // Check min rows for A
         if (ldb < min_ldb_f) { info = -8; goto cleanup; } // Check min rows for B
         if (ldc < min_ldc_f) { info = -10; goto cleanup; } // Check min rows for C
         if (ldd < min_ldd_f) { info = -12; goto cleanup; } // Check min rows for D
         if (ldu < min_ldu_f) { info = -14; goto cleanup; } // Check min rows for U
         if (ldy < min_ldy_f) { info = -17; goto cleanup; } // Check min rows for Y
     }
 
     // Exit if any validation failed before allocating memory
     if (info != 0) { goto cleanup; }
 
     // 3. Internal Workspace Allocation
     // Calculate required size based on documentation: DWORK(N)
     ldwork = MAX(1, n); // Size is N, ensure at least 1
     dwork = (double*)malloc((size_t)ldwork * elem_size);
     // CHECK_ALLOC sets info = SLICOT_MEMORY_ERROR and jumps to cleanup on failure
     CHECK_ALLOC(dwork);
 
     // 4. Memory allocation for column-major copies (if row_major)
     // Calculate sizes, handling zero dimensions
     size_t a_rows = n; size_t a_cols = n; size_t a_size = (n > 0) ? (a_rows * a_cols) : 0;
     size_t b_rows = n; size_t b_cols = m; size_t b_size = (n > 0 && m > 0) ? (b_rows * b_cols) : 0;
     size_t c_rows = p; size_t c_cols = n; size_t c_size = (p > 0 && n > 0) ? (c_rows * c_cols) : 0;
     size_t d_rows = p; size_t d_cols = m; size_t d_size = (p > 0 && m > 0) ? (d_rows * d_cols) : 0;
     size_t u_rows = m; size_t u_cols = ny; size_t u_size = (m > 0 && ny > 0) ? (u_rows * u_cols) : 0;
     size_t y_rows = p; size_t y_cols = ny; size_t y_size = (p > 0 && ny > 0) ? (y_rows * y_cols) : 0;
 
     if (row_major) {
         // Allocate memory for column-major copies of inputs
         if (a_size > 0) { a_cm = (double*)malloc(a_size * elem_size); CHECK_ALLOC(a_cm); }
         if (b_size > 0) { b_cm = (double*)malloc(b_size * elem_size); CHECK_ALLOC(b_cm); }
         if (c_size > 0) { c_cm = (double*)malloc(c_size * elem_size); CHECK_ALLOC(c_cm); }
         if (d_size > 0) { d_cm = (double*)malloc(d_size * elem_size); CHECK_ALLOC(d_cm); }
         if (u_size > 0) { u_cm = (double*)malloc(u_size * elem_size); CHECK_ALLOC(u_cm); }
         // Allocate memory for column-major temporary output Y
         if (y_size > 0) { y_cm = (double*)malloc(y_size * elem_size); CHECK_ALLOC(y_cm); }
     }
 
     // 5. Prepare Fortran parameters and perform conversions
     if (row_major) {
         // Set Fortran LDA (number of ROWS)
         lda_f = min_lda_f; // MAX(1, n)
         ldb_f = min_ldb_f; // MAX(1, n)
         ldc_f = min_ldc_f; // MAX(1, p)
         ldd_f = min_ldd_f; // MAX(1, p)
         ldu_f = min_ldu_f; // MAX(1, m)
         ldy_f = min_ldy_f; // MAX(1, p)
 
         // Point to column-major copies if allocated, else NULL
         a_ptr = (a_size > 0) ? a_cm : NULL;
         b_ptr = (b_size > 0) ? b_cm : NULL;
         c_ptr = (c_size > 0) ? c_cm : NULL;
         d_ptr = (d_size > 0) ? d_cm : NULL;
         u_ptr = (u_size > 0) ? u_cm : NULL;
         y_ptr = (y_size > 0) ? y_cm : NULL; // Fortran writes to temp buffer
 
         // Convert input matrices from C row-major to Fortran column-major copies
         // Use C leading dimensions (cols) for transpose source (lda, ldb, ldc, ldd, ldu)
         // Use Fortran leading dimensions (rows) for transpose dest (lda_f, ldb_f, ldc_f, ldd_f, ldu_f)
         if (a_ptr) slicot_transpose_to_fortran_with_ld(a, a_cm, a_rows, a_cols, lda, lda_f, elem_size);
         if (b_ptr) slicot_transpose_to_fortran_with_ld(b, b_cm, b_rows, b_cols, ldb, ldb_f, elem_size);
         if (c_ptr) slicot_transpose_to_fortran_with_ld(c, c_cm, c_rows, c_cols, ldc, ldc_f, elem_size);
         if (d_ptr) slicot_transpose_to_fortran_with_ld(d, d_cm, d_rows, d_cols, ldd, ldd_f, elem_size);
         if (u_ptr) slicot_transpose_to_fortran_with_ld(u, u_cm, u_rows, u_cols, ldu, ldu_f, elem_size);
 
     } else {
         // Column-major C: Use original C pointers and LDs (which are rows)
         lda_f = lda;
         ldb_f = ldb;
         ldc_f = ldc;
         ldd_f = ldd;
         ldu_f = ldu;
         ldy_f = ldy;
 
         // Pass original pointers, ensuring NULL if size is 0
         a_ptr = (a_size > 0) ? a : NULL;
         b_ptr = (b_size > 0) ? b : NULL;
         c_ptr = (c_size > 0) ? c : NULL;
         d_ptr = (d_size > 0) ? d : NULL;
         u_ptr = (u_size > 0) ? u : NULL;
         y_ptr = (y_size > 0) ? y : NULL; // Fortran writes directly to user's Y
     }
 
     // X is 1D array (vector), no transposition needed, passed directly.
     // Pointer check for X already done. If n=0, x is NULL or ignored by Fortran.
 
     // 7. Call Fortran function
     F77_FUNC(tf01md, TF01MD)(
         &n, &m, &p, &ny,
         a_ptr, &lda_f, b_ptr, &ldb_f, c_ptr, &ldc_f, d_ptr, &ldd_f,
         u_ptr, &ldu_f,
         x, // Pass C pointer directly (in/out)
         y_ptr, // Pass pointer for Fortran output (y_cm or y)
         &ldy_f,
         dwork, // Pass internally allocated workspace
         &info
     );
 
     // 8. Convert results back to row-major (if needed)
     if (row_major && info == 0) {
         // X was modified in place (1D).
         // Y was written to y_cm, needs to be transposed back to user's y.
         if (y_ptr) { // Check if y_cm was allocated and used
              // Use Fortran LD (rows = ldy_f) for transpose source, C LD (cols = ldy) for dest
             slicot_transpose_to_c_with_ld(y_cm, y, y_rows, y_cols, ldy_f, ldy, elem_size);
         }
     }
     // If column-major, Y was written directly into user's y array.
 
 cleanup:
     /* --- Cleanup --- */
     // Free internally allocated workspace FIRST
     free(dwork); // Safe to call free(NULL)
 
     // Free temporary column-major arrays if allocated
     free(a_cm);
     free(b_cm);
     free(c_cm);
     free(d_cm);
     free(u_cm);
     free(y_cm); // Free the temporary output buffer
 
     // Check if info was set by CHECK_ALLOC during workspace/copy allocation
     if (info == SLICOT_MEMORY_ERROR) {
        fprintf(stderr, "Error: Memory allocation failed in slicot_tf01md.\n");
     }
     // Return the info code (either from validation, CHECK_ALLOC, or Fortran)
     return info;
 }
 