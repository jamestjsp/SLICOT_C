/**
 * @file slicot_utils.h
 * @brief Utility functions and definitions for SLICOT C wrappers
 *
 * This file contains utility function declarations, macro definitions,
 * and type definitions for the SLICOT C wrapper library, including array
 * transposition, complex number handling, and common macros.
 */

 #ifndef SLICOT_UTILS_H
 #define SLICOT_UTILS_H
 
 #include <stdlib.h> // For size_t, malloc (used in CHECK_ALLOC context)
 #include <stddef.h> // For size_t
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /* Complex number support */
 #ifdef __STDC_NO_COMPLEX__
 /* If complex numbers not supported in C99, define our own structure */
 typedef struct {
     double real;
     double imag;
 } slicot_complex_double;
 
 typedef struct {
     float real;
     float imag;
 } slicot_complex_float;
 #else
 /* Use C99 complex types if available */
 #include <complex.h>
 typedef double complex slicot_complex_double;
 typedef float complex slicot_complex_float;
 #endif
 
 /* Error Codes */
 /**
  * @brief Error code for memory allocation failure within SLICOT wrappers.
  * Chosen to be distinct from standard SLICOT INFO codes.
  */
 #define SLICOT_MEMORY_ERROR -1010
 
 /* Utility Macros */
 
 /**
  * @brief Macro to check the result of a memory allocation (e.g., malloc).
  *
  * If the pointer 'ptr' is NULL, it sets the 'info' variable (which must
  * be in scope) to SLICOT_MEMORY_ERROR and jumps to a 'cleanup' label
  * (which must exist in the calling function).
  *
  * Usage:
  * ptr = malloc(size);
  * CHECK_ALLOC(ptr);
  *
  * @param ptr The pointer returned by the allocation function.
  */
 #define CHECK_ALLOC(ptr)                                     \
     do {                                                     \
         if ((ptr) == NULL) {                                 \
             info = SLICOT_MEMORY_ERROR; /* Set error code */ \
             goto cleanup;           /* Jump to cleanup label */ \
         }                                                    \
     } while (0)
 
 /* MAX macro */
 #ifndef MAX
 #define MAX(a,b) (((a) > (b)) ? (a) : (b))
 #endif
 
 /* MIN macro */
 #ifndef MIN
 #define MIN(a,b) (((a) < (b)) ? (a) : (b))
 #endif
 
 /* Function Declarations for Transpose Utilities */
 
 /**
  * @brief Transpose a matrix from C (row-major) to Fortran (column-major) order.
  *
  * @param src Pointer to the source matrix (row-major).
  * @param dest Pointer to the destination matrix (column-major).
  * @param rows Number of rows in the matrix.
  * @param cols Number of columns in the matrix.
  * @param elem_size Size (in bytes) of a single matrix element.
  */
 void slicot_transpose_to_fortran(const void *src, void *dest, int rows, int cols, size_t elem_size);
 
 /**
  * @brief Transpose a matrix from Fortran (column-major) to C (row-major) order.
  *
  * @param src Pointer to the source matrix (column-major).
  * @param dest Pointer to the destination matrix (row-major).
  * @param rows Number of rows in the matrix.
  * @param cols Number of columns in the matrix.
  * @param elem_size Size (in bytes) of a single matrix element.
  */
 void slicot_transpose_to_c(const void *src, void *dest, int rows, int cols, size_t elem_size);
 
 /**
  * @brief In-place transpose of a square matrix (assumes row-major indexing for swap).
  *
  * @param matrix Pointer to the matrix to be transposed in-place.
  * @param rows Number of rows (must equal columns).
  * @param cols Number of columns (must equal rows).
  * @param elem_size Size (in bytes) of a single matrix element.
  * @return Returns 0 on success, -1 on error (if rows != cols or memory allocation fails).
  */
 int slicot_transpose_inplace(void *matrix, int rows, int cols, size_t elem_size);

/* Helper function to set a matrix to identity */
 static void set_identity(int n, double* mat, int ld, int row_major);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SLICOT_UTILS_H */
 