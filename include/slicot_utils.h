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
 
 #include "slicot_utils.h"  /* Complex number support */
 #ifdef __cplusplus
 /* Use C++ complex types when compiling with C++ */
 #include <complex>
 typedef std::complex<double> slicot_complex_double;
 typedef std::complex<float> slicot_complex_float;
 #elif defined(__STDC_NO_COMPLEX__)
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

 #ifdef __cplusplus
 extern "C" {
 #endif
 
 /* Macro to get the real part of slicot_complex_double */
 #ifdef __STDC_NO_COMPLEX__
 #define SLICOT_COMPLEX_REAL(z) ((z).real)
 #else
 #include <complex.h> // Ensure creal is available
 #define SLICOT_COMPLEX_REAL(z) (creal(z))
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
 
#ifndef SLICOT_C_WRAPPER_API // Prevent multiple definitions

  // Check if building a static library.
  // SLICOT_STATIC should be defined by the build system (e.g., CMake)
  // when BUILD_SHARED_LIBS is OFF.
  #ifdef SLICOT_STATIC
    #define SLICOT_C_WRAPPER_API

  // Building or using a shared library
  #else
    // Check if on Windows
    #ifdef _WIN32
      // Check if building the DLL.
      // SLICOT_C_WRAPPER_EXPORTS should be defined by the build system (e.g., CMake)
      // when building the shared library (BUILD_SHARED_LIBS is ON).
      #ifdef SLICOT_C_WRAPPER_EXPORTS
        #define SLICOT_C_WRAPPER_API __declspec(dllexport)
      // Using the DLL
      #else
        #define SLICOT_C_WRAPPER_API __declspec(dllimport)
      #endif // SLICOT_C_WRAPPER_EXPORTS

    // Non-Windows platforms (Linux, macOS, etc.)
    #else
      // Use GCC/Clang visibility attributes if available for better symbol handling
      #if defined(__GNUC__) && (__GNUC__ >= 4)
        #define SLICOT_C_WRAPPER_API __attribute__ ((visibility ("default")))
      #else
        #define SLICOT_C_WRAPPER_API // Default: empty definition
      #endif // __GNUC__

    #endif // _WIN32

  #endif // SLICOT_STATIC

#endif // SLICOT_C_WRAPPER_API



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
  SLICOT_C_WRAPPER_API
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
 SLICOT_C_WRAPPER_API
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
 SLICOT_C_WRAPPER_API
 int slicot_transpose_inplace(void *matrix, int rows, int cols, size_t elem_size);

/**
  * @brief Sets a square matrix to the identity matrix.
  *
  * This function sets the diagonal elements of the matrix to 1.0 and all
  * off-diagonal elements to 0.0. The matrix is assumed to be stored in either
  * row-major or column-major order, as specified by the 'row_major' parameter.
  *
  * @param n Order of the square matrix.
  * @param mat Pointer to the matrix (assumed to be allocated with sufficient size).
  * @param ld Leading dimension of the matrix (number of rows for column-major,
  *          number of columns for row-major).
  * @param row_major Integer flag indicating storage order:
  * = 0: Column-major (Fortran style).
  * = 1: Row-major (C style).
  */
  SLICOT_C_WRAPPER_API
 void set_identity(int n, double* mat, int ld, int row_major);
 
/**
 * @brief Copies the relevant triangle of a symmetric matrix to a full matrix.
 *
 * Used primarily for column-major wrappers where the Fortran routine might
 * access the unspecified triangle if the original storage is passed directly.
 * Copies the specified triangle (upper or lower) from src to dest,
 * filling the other triangle by symmetry. Assumes column-major storage
 * for both source and destination.
 *
 * @param src Pointer to the source symmetric matrix (column-major).
 * @param dest Pointer to the destination full matrix (column-major).
 * @param n Order of the square matrix.
 * @param uplo Specifies which triangle of src is stored ('U' or 'L').
 * @param ld Leading dimension of both src and dest.
 * @param elem_size Size (in bytes) of a single matrix element.
 */
  SLICOT_C_WRAPPER_API
void slicot_copy_symmetric_part(const void *src, void *dest, int n, char uplo, int ld, size_t elem_size);

/**
 * @brief Transpose a symmetric matrix from C (row-major, triangular storage)
 * to a full Fortran (column-major) matrix.
 *
 * @param src Pointer to the source matrix (row-major, stores upper or lower triangle).
 * @param dest Pointer to the destination matrix (column-major, will be filled fully).
 * @param n Order of the square matrix.
 * @param uplo Specifies which triangle of src is stored ('U' or 'L').
 * @param elem_size Size (in bytes) of a single matrix element.
 */
SLICOT_C_WRAPPER_API
void slicot_transpose_symmetric_to_fortran(const void *src, void *dest, int n, char uplo, size_t elem_size);

/**
 * @brief Transpose a symmetric matrix from Fortran (column-major, full storage)
 * to C (row-major, potentially triangular storage - copies only specified triangle).
 *
 * @param src Pointer to the source matrix (column-major, assumed fully populated).
 * @param dest Pointer to the destination matrix (row-major, will store only upper or lower triangle).
 * @param n Order of the square matrix.
 * @param uplo Specifies which triangle of src to copy to dest ('U' or 'L').
 * @param elem_size Size (in bytes) of a single matrix element.
 */
SLICOT_C_WRAPPER_API
void slicot_transpose_symmetric_to_c(const void *src, void *dest, int n, char uplo, size_t elem_size);


/**
 * @brief Transpose a matrix from Fortran (column-major) to C (row-major) order with custom leading dimensions.
 *
 * @param src Pointer to the source matrix (column-major).
 * @param dest Pointer to the destination matrix (row-major).
 * @param rows Number of rows to copy.
 * @param cols Number of columns to copy.
 * @param ld_src Leading dimension of the source matrix (number of rows).
 * @param ld_dest Leading dimension of the destination matrix (number of columns).
 * @param elem_size Size (in bytes) of a single matrix element.
 */
SLICOT_C_WRAPPER_API
void slicot_transpose_to_c_with_ld(const void *src, void *dest, int rows, int cols, 
                                  int ld_src, int ld_dest, size_t elem_size);
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* SLICOT_UTILS_H */
 