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
#include <ctype.h> // Fot toupper
#include <stdio.h>

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
SLICOT_EXPORT
void slicot_transpose_to_fortran(const void *src, void *dest, int rows, int cols, size_t elem_size)
{
    // Use char pointers for byte-level arithmetic
    const char *src_ptr = (const char *)src;
    char *dest_ptr = (char *)dest;
    int i, j; // Row and column indices

    // Check for invalid inputs (optional, but good practice)
    if (!src || !dest || rows <= 0 || cols <= 0 || elem_size <= 0)
    {
        return; // Or handle error appropriately
    }

    for (i = 0; i < rows; i++)
    { // Iterate through rows
        for (j = 0; j < cols; j++)
        { // Iterate through columns
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
SLICOT_EXPORT
void slicot_transpose_to_c(const void *src, void *dest, int rows, int cols, size_t elem_size)
{
    // Use char pointers for byte-level arithmetic
    const char *src_ptr = (const char *)src;
    char *dest_ptr = (char *)dest;
    int i, j; // Row and column indices

    // Check for invalid inputs (optional, but good practice)
    if (!src || !dest || rows <= 0 || cols <= 0 || elem_size <= 0)
    {
        return; // Or handle error appropriately
    }

    for (i = 0; i < rows; i++)
    { // Iterate through rows
        for (j = 0; j < cols; j++)
        { // Iterate through columns
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
 * @brief Transpose a matrix from Fortran (column-major) to C (row-major) order with custom leading dimensions.
 *
 * Copies elements from the source (column-major) matrix to the destination
 * (row-major) matrix, performing the necessary index transformation and respecting
 * the provided leading dimensions.
 *
 * @param src Pointer to the source matrix (column-major).
 * @param dest Pointer to the destination matrix (row-major).
 * @param rows Number of rows to copy.
 * @param cols Number of columns to copy.
 * @param ld_src Leading dimension of the source matrix (number of rows).
 * @param ld_dest Leading dimension of the destination matrix (number of columns).
 * @param elem_size Size (in bytes) of a single matrix element.
 */
SLICOT_EXPORT
void slicot_transpose_to_c_with_ld(const void *src, void *dest, int rows, int cols, 
                                   int ld_src, int ld_dest, size_t elem_size)
{
    // Use char pointers for byte-level arithmetic
    const char *src_ptr = (const char *)src;
    char *dest_ptr = (char *)dest;
    int i, j; // Row and column indices

    // Check for invalid inputs
    if (!src || !dest || rows <= 0 || cols <= 0 || ld_src < rows || ld_dest < cols || elem_size <= 0)
    {
        return; // Or handle error appropriately
    }

    for (i = 0; i < rows; i++)
    { // Iterate through rows
        for (j = 0; j < cols; j++)
        { // Iterate through columns
            // Source index (column-major): element at [i][j] is at offset (i + j * ld_src)
            size_t src_offset = (i + (size_t)j * ld_src) * elem_size;
            // Destination index (row-major): element at [i][j] is at offset (i * ld_dest + j)
            size_t dest_offset = ((size_t)i * ld_dest + j) * elem_size;

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
SLICOT_EXPORT
int slicot_transpose_inplace(void *matrix, int rows, int cols, size_t elem_size)
{
    // Use char pointer for byte-level arithmetic
    char *mat_ptr = (char *)matrix;
    int i, j;
    char *temp_buf = NULL; // Temporary buffer for one element

    /* In-place transpose only works for square matrices */
    if (rows != cols || rows <= 0 || elem_size <= 0 || !matrix)
    {
        return -1; // Invalid input or non-square matrix
    }

    /* Allocate temporary buffer for swapping */
    temp_buf = (char *)malloc(elem_size);
    if (!temp_buf)
    {
        return -1; // Memory allocation failed
    }

    /* Transpose the matrix in-place by swapping upper triangle with lower triangle */
    for (i = 0; i < rows; i++)
    {
        // Start j from i + 1 to only process the upper triangle (excluding diagonal)
        for (j = i + 1; j < cols; j++)
        {
            // Calculate offsets assuming row-major layout
            size_t offset_ij = ((size_t)i * cols + j) * elem_size;
            size_t offset_ji = ((size_t)j * cols + i) * elem_size; // Note: uses 'cols' here too for row-major

            /* Swap elements at [i][j] and [j][i] using the temporary buffer */
            memcpy(temp_buf, mat_ptr + offset_ij, elem_size);            // temp = mat[i][j]
            memcpy(mat_ptr + offset_ij, mat_ptr + offset_ji, elem_size); // mat[i][j] = mat[j][i]
            memcpy(mat_ptr + offset_ji, temp_buf, elem_size);            // mat[j][i] = temp
        }
    }

    /* Free the temporary buffer */
    free(temp_buf);
    return 0; // Success
}

/* Helper function to set a matrix to identity */
SLICOT_EXPORT
void set_identity(int n, double *mat, int ld, int row_major)
{
    if (!mat)
        return;
    memset(mat, 0, (size_t)ld * n * sizeof(double)); // Zero out the matrix first
    if (row_major)
    {
        // mat[i * ld + i] = 1.0;
        for (int i = 0; i < n; ++i)
        {
            if (i < ld)
            { // Check column index is within bounds
                mat[(size_t)i * ld + i] = 1.0;
            }
        }
    }
    else
    {
        // mat[i + i * ld] = 1.0;
        for (int i = 0; i < n; ++i)
        {
            if (i < n)
            { // Check row index is within bounds
                mat[i + (size_t)i * ld] = 1.0;
            }
        }
    }
}


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
SLICOT_EXPORT
void slicot_copy_symmetric_part(const void *src, void *dest, int n, char uplo, int ld, size_t elem_size)
{
    const char *src_ptr = (const char *)src;
    char *dest_ptr = (char *)dest;
    int i, j;
    char uplo_upper = toupper(uplo);

    if (!src || !dest || n <= 0 || ld < n || elem_size <= 0) {
        return; // Invalid input
    }

    if (uplo_upper == 'U') {
        // Copy upper triangle and fill lower by symmetry
        for (j = 0; j < n; ++j) { // Column index
            for (i = 0; i <= j; ++i) { // Row index (up to diagonal)
                size_t src_offset = (i + (size_t)j * ld) * elem_size;
                size_t dest_offset_ij = (i + (size_t)j * ld) * elem_size;
                size_t dest_offset_ji = (j + (size_t)i * ld) * elem_size; // Symmetric position

                // Copy src[i][j] to dest[i][j]
                memcpy(dest_ptr + dest_offset_ij, src_ptr + src_offset, elem_size);
                // Copy src[i][j] to dest[j][i] (if not diagonal)
                if (i != j) {
                    memcpy(dest_ptr + dest_offset_ji, src_ptr + src_offset, elem_size);
                }
            }
        }
    } else { // Assume 'L'
        // Copy lower triangle and fill upper by symmetry
        for (j = 0; j < n; ++j) { // Column index
            for (i = j; i < n; ++i) { // Row index (from diagonal down)
                size_t src_offset = (i + (size_t)j * ld) * elem_size;
                size_t dest_offset_ij = (i + (size_t)j * ld) * elem_size;
                size_t dest_offset_ji = (j + (size_t)i * ld) * elem_size; // Symmetric position

                // Copy src[i][j] to dest[i][j]
                memcpy(dest_ptr + dest_offset_ij, src_ptr + src_offset, elem_size);
                // Copy src[i][j] to dest[j][i] (if not diagonal)
                if (i != j) {
                    memcpy(dest_ptr + dest_offset_ji, src_ptr + src_offset, elem_size);
                }
            }
        }
    }
}

/**
 * @brief Transpose a symmetric matrix from C (row-major, triangular storage)
 * to a full Fortran (column-major) matrix.
 *
 * Copies elements from the specified triangle (upper or lower) of the source
 * (row-major) matrix to the destination (column-major) matrix, filling
 * the other triangle by symmetry.
 *
 * @param src Pointer to the source matrix (row-major, stores upper or lower triangle).
 * @param dest Pointer to the destination matrix (column-major, will be filled fully).
 * @param n Order of the square matrix.
 * @param uplo Specifies which triangle of src is stored ('U' or 'L').
 * @param elem_size Size (in bytes) of a single matrix element.
 */
SLICOT_EXPORT
void slicot_transpose_symmetric_to_fortran(const void *src, void *dest, int n, char uplo, size_t elem_size)
{
    // Use char pointers for byte-level arithmetic
    const char *src_ptr = (const char *)src;
    char *dest_ptr = (char *)dest;
    int i, j; // Row and column indices
    char uplo_upper = toupper(uplo);

    // Check for invalid inputs (optional, but good practice)
    if (!src || !dest || n <= 0 || elem_size <= 0)
    {
        return; // Or handle error appropriately
    }

    if (uplo_upper == 'U') {
        // Source stores upper triangle (row-major)
        for (i = 0; i < n; i++) { // Iterate through rows
            for (j = i; j < n; j++) { // Iterate through columns (from diagonal right)
                // Source index (row-major): element at [i][j] is at offset (i * n + j)
                // Note: Assumes the source buffer has space for n*n elements, even if only triangle is used.
                // If source is packed, this offset calculation needs change. Assuming non-packed source.
                size_t src_offset = ((size_t)i * n + j) * elem_size;

                // Destination index (column-major): element at [i][j] is at offset (i + j * n)
                size_t dest_offset_ij = (i + (size_t)j * n) * elem_size;
                // Destination index for symmetric part [j][i]: offset (j + i * n)
                size_t dest_offset_ji = (j + (size_t)i * n) * elem_size;

                // Copy the element from src[i][j] to dest[i][j]
                memcpy(dest_ptr + dest_offset_ij, src_ptr + src_offset, elem_size);

                // Copy the element from src[i][j] to dest[j][i] (if not diagonal)
                if (i != j) {
                    memcpy(dest_ptr + dest_offset_ji, src_ptr + src_offset, elem_size);
                }
            }
        }
    } else { // Assume 'L'
        // Source stores lower triangle (row-major)
        for (i = 0; i < n; i++) { // Iterate through rows
            for (j = 0; j <= i; j++) { // Iterate through columns (up to diagonal)
                // Source index (row-major): element at [i][j] is at offset (i * n + j)
                size_t src_offset = ((size_t)i * n + j) * elem_size;

                // Destination index (column-major): element at [i][j] is at offset (i + j * n)
                size_t dest_offset_ij = (i + (size_t)j * n) * elem_size;
                // Destination index for symmetric part [j][i]: offset (j + i * n)
                size_t dest_offset_ji = (j + (size_t)i * n) * elem_size;

                // Copy the element from src[i][j] to dest[i][j]
                memcpy(dest_ptr + dest_offset_ij, src_ptr + src_offset, elem_size);

                // Copy the element from src[i][j] to dest[j][i] (if not diagonal)
                if (i != j) {
                    memcpy(dest_ptr + dest_offset_ji, src_ptr + src_offset, elem_size);
                }
            }
        }
    }
}

/**
 * @brief Transpose a symmetric matrix from Fortran (column-major, full storage)
 * to C (row-major, potentially triangular storage - copies only specified triangle).
 *
 * Copies elements from the specified triangle (upper or lower) of the source
 * (column-major, assumed full) matrix to the destination (row-major) matrix.
 * This is useful when the C code only needs to store one triangle of the result.
 *
 * @param src Pointer to the source matrix (column-major, assumed fully populated).
 * @param dest Pointer to the destination matrix (row-major, will store only upper or lower triangle).
 * @param n Order of the square matrix.
 * @param uplo Specifies which triangle of src to copy to dest ('U' or 'L').
 * @param elem_size Size (in bytes) of a single matrix element.
 */
SLICOT_EXPORT
void slicot_transpose_symmetric_to_c(const void *src, void *dest, int n, char uplo, size_t elem_size)
{
    // Use char pointers for byte-level arithmetic
    const char *src_ptr = (const char *)src;
    char *dest_ptr = (char *)dest;
    int i, j; // Row and column indices
    char uplo_upper = toupper(uplo);

    // Check for invalid inputs
    if (!src || !dest || n <= 0 || elem_size <= 0) {
        return;
    }

    if (uplo_upper == 'U') {
        // Copy upper triangle from Fortran (col-major) to C (row-major)
        for (i = 0; i < n; i++) { // Row index
            for (j = i; j < n; j++) { // Column index (from diagonal right)
                // Source index (col-major): element at [i][j] is at offset (i + j * n)
                size_t src_offset = (i + (size_t)j * n) * elem_size;
                // Destination index (row-major): element at [i][j] is at offset (i * n + j)
                size_t dest_offset = ((size_t)i * n + j) * elem_size;
                // Copy the element
                memcpy(dest_ptr + dest_offset, src_ptr + src_offset, elem_size);
            }
        }
    } else { // Assume 'L'
        // Copy lower triangle from Fortran (col-major) to C (row-major)
        for (i = 0; i < n; i++) { // Row index
            for (j = 0; j <= i; j++) { // Column index (up to diagonal)
                // Source index (col-major): element at [i][j] is at offset (i + j * n)
                size_t src_offset = (i + (size_t)j * n) * elem_size;
                // Destination index (row-major): element at [i][j] is at offset (i * n + j)
                size_t dest_offset = ((size_t)i * n + j) * elem_size;
                // Copy the element
                memcpy(dest_ptr + dest_offset, src_ptr + src_offset, elem_size);
            }
        }
    }
}

/**
 * @brief Prints a matrix of doubles to the standard output.
 *
 * This helper function is intended for debugging purposes. It takes a pointer to a matrix of doubles,
 * along with its dimensions, and prints the matrix in a readable format to the standard output (usually the console).
 * Each row of the matrix is printed on a new line, with elements separated by spaces or tabs for clarity.
 *
 * @param mat Pointer to the first element of the matrix (assumed to be stored in row-major order).
 * @param rows The number of rows in the matrix.
 * @param cols The number of columns in the matrix.
 *
 * @note This function is primarily intended for debugging and should not be used in production code
 * where performance or formatted output is critical.
 */
SLICOT_EXPORT
void printMatrixD(const char* name, const double* data, int rows, int cols, int ld, int rowMajor) {
    if (name)
        printf("%s (%dx%d, ld=%d, %s):\n", name, rows, cols, ld, rowMajor ? "RowMajor" : "ColMajor");
    for (int i = 0; i < rows; ++i) {
        printf("  [");
        for (int j = 0; j < cols; ++j) {
            double val = rowMajor ? data[i * ld + j] : data[i + j * ld];
            printf("%9.4f%s", val, (j == cols - 1 ? "" : ", "));
        }
        printf("]\n");
    }
}
