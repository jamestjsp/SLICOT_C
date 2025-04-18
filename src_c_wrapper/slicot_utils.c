/**
 * @file slicot_utils.c
 * @brief Implementation of utility functions for SLICOT C wrappers.
 *
 * Implements matrix transposition functions needed to interface C code
 * (typically using row-major order) with Fortran routines (using
 * column-major order).
 */

 #include "slicot_utils.h"
 #include <string.h> // For memcpy
 #include <stdlib.h> // For malloc, free
 
 /**
  * @brief Transpose a matrix from C (row-major) to Fortran (column-major) order.
  *
  * Copies elements from the source (row-major) matrix to the destination
  * (column-major) matrix, performing the necessary index transformation.
  *
  * @param src Pointer to the source matrix (row-major). Assumed size: rows * cols * elem_size.
  * @param dest Pointer to the destination matrix (column-major). Assumed size: rows * cols * elem_size.
  * @param rows Number of rows in the matrix.
  * @param cols Number of columns in the matrix.
  * @param elem_size Size (in bytes) of a single matrix element.
  */
 void slicot_transpose_to_fortran(const void *src, void *dest, int rows, int cols, size_t elem_size) {
     // Use char pointers for byte-level arithmetic
     const char *src_ptr = (const char *)src;
     char *dest_ptr = (char *)dest;
     int i, j; // Row and column indices
 
     // Check for invalid inputs (optional, but good practice)
     if (!src || !dest || rows <= 0 || cols <= 0 || elem_size <= 0) {
         return; // Or handle error appropriately
     }
 
     for (i = 0; i < rows; i++) {        // Iterate through rows
         for (j = 0; j < cols; j++) {    // Iterate through columns
             // Source index (row-major): element at [i][j] is at offset (i * cols + j)
             size_t src_offset = ((size_t)i * cols + j) * elem_size;
             // Destination index (column-major): element at [i][j] is at offset (i + j * rows)
             size_t dest_offset = (i + (size_t)j * rows) * elem_size;
 
             // Copy the element
             memcpy(dest_ptr + dest_offset, src_ptr + src_offset, elem_size);
         }
     }
 }
 
 /**
  * @brief Transpose a matrix from Fortran (column-major) to C (row-major) order.
  *
  * Copies elements from the source (column-major) matrix to the destination
  * (row-major) matrix, performing the necessary index transformation.
  *
  * @param src Pointer to the source matrix (column-major). Assumed size: rows * cols * elem_size.
  * @param dest Pointer to the destination matrix (row-major). Assumed size: rows * cols * elem_size.
  * @param rows Number of rows in the matrix.
  * @param cols Number of columns in the matrix.
  * @param elem_size Size (in bytes) of a single matrix element.
  */
 void slicot_transpose_to_c(const void *src, void *dest, int rows, int cols, size_t elem_size) {
     // Use char pointers for byte-level arithmetic
     const char *src_ptr = (const char *)src;
     char *dest_ptr = (char *)dest;
     int i, j; // Row and column indices
 
     // Check for invalid inputs (optional, but good practice)
     if (!src || !dest || rows <= 0 || cols <= 0 || elem_size <= 0) {
         return; // Or handle error appropriately
     }
 
     for (i = 0; i < rows; i++) {        // Iterate through rows
         for (j = 0; j < cols; j++) {    // Iterate through columns
             // Source index (column-major): element at [i][j] is at offset (i + j * rows)
             size_t src_offset = (i + (size_t)j * rows) * elem_size;
             // Destination index (row-major): element at [i][j] is at offset (i * cols + j)
             size_t dest_offset = ((size_t)i * cols + j) * elem_size;
 
             // Copy the element
             memcpy(dest_ptr + dest_offset, src_ptr + src_offset, elem_size);
         }
     }
 }
 
 /**
  * @brief In-place transpose of a square matrix (assumes row-major indexing for swap).
  *
  * Transposes the matrix by swapping elements across the main diagonal.
  * Requires temporary storage for one element.
  *
  * @param matrix Pointer to the matrix (assumed row-major layout).
  * @param rows Number of rows (must equal columns).
  * @param cols Number of columns (must equal rows).
  * @param elem_size Size (in bytes) of a single matrix element.
  * @return Returns 0 on success, -1 on error (if rows != cols or memory allocation fails).
  */
 int slicot_transpose_inplace(void *matrix, int rows, int cols, size_t elem_size) {
     // Use char pointer for byte-level arithmetic
     char *mat_ptr = (char *)matrix;
     int i, j;
     char *temp_buf = NULL; // Temporary buffer for one element
 
     /* In-place transpose only works for square matrices */
     if (rows != cols || rows <= 0 || elem_size <= 0 || !matrix) {
         return -1; // Invalid input or non-square matrix
     }
 
     /* Allocate temporary buffer for swapping */
     temp_buf = (char *)malloc(elem_size);
     if (!temp_buf) {
         return -1; // Memory allocation failed
     }
 
     /* Transpose the matrix in-place by swapping upper triangle with lower triangle */
     for (i = 0; i < rows; i++) {
         // Start j from i + 1 to only process the upper triangle (excluding diagonal)
         for (j = i + 1; j < cols; j++) {
             // Calculate offsets assuming row-major layout
             size_t offset_ij = ((size_t)i * cols + j) * elem_size;
             size_t offset_ji = ((size_t)j * cols + i) * elem_size; // Note: uses 'cols' here too for row-major
 
             /* Swap elements at [i][j] and [j][i] using the temporary buffer */
             memcpy(temp_buf,           mat_ptr + offset_ij, elem_size); // temp = mat[i][j]
             memcpy(mat_ptr + offset_ij, mat_ptr + offset_ji, elem_size); // mat[i][j] = mat[j][i]
             memcpy(mat_ptr + offset_ji, temp_buf,            elem_size); // mat[j][i] = temp
         }
     }
 
     /* Free the temporary buffer */
     free(temp_buf);
     return 0; // Success
 }
 
/* Helper function to set a matrix to identity */
 static void set_identity(int n, double* mat, int ld, int row_major) {
     if (!mat) return;
     memset(mat, 0, (size_t)ld * n * sizeof(double)); // Zero out the matrix first
     if (row_major) {
         // mat[i * ld + i] = 1.0;
         for (int i = 0; i < n; ++i) {
             if (i < ld) { // Check column index is within bounds
                mat[(size_t)i * ld + i] = 1.0;
             }
         }
     } else {
         // mat[i + i * ld] = 1.0;
         for (int i = 0; i < n; ++i) {
              if (i < n) { // Check row index is within bounds
                 mat[i + (size_t)i * ld] = 1.0;
              }
         }
     }
 }