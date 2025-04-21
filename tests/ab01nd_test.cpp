#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits

// Include the specific header for the function being tested
#include "ab01nd.h" // Header for slicot_ab01nd

// Removed the custom ExpectNear helper function as it prevents using << for messages.
// We will use EXPECT_NEAR directly in the tests.


// =============================================================================
// Test Fixture and Test Case for Column-Major Input (ROW_MAJOR = 0)
// =============================================================================

// Define a fixture for the column-major test
class AB01NDTestColMajor : public ::testing::TestWithParam<char> {
    // No specific setup needed here for now
};

// Test case based on the example in the AB01ND documentation (Column-Major)
TEST_P(AB01NDTestColMajor, DocExampleColMajor) {
    // --- Input data from AB01ND documentation example ---
    char JOBZ = GetParam();  // Test all three values: 'N', 'I', 'F'
    int N = 3;        // Order of the system
    int M = 2;        // Number of inputs
    double TOL = 0.0; // Use default tolerance
    int ROW_MAJOR = 0; // Input arrays are column-major (Fortran style)
    double check_tol = 1e-4; // Tolerance for checking results

    // --- Leading Dimensions for Column-Major (Fortran style) ---
    int LDA = N;      // LDA >= MAX(1,N)
    int LDB = N;      // LDB >= MAX(1,N)
    int LDZ = N;      // LDZ >= MAX(1,N) if JOBZ != 'N', else >= 1

    // --- Input Matrix A (Column-major order as in Fortran example DATA section) ---
    std::vector<double> A = {
        -1.0, -2.0, -1.0, // Col 1
         0.0, -2.0,  0.0, // Col 2
         0.0, -2.0, -3.0  // Col 3
    };

    // --- Input Matrix B (Column-major order based on Fortran example DATA section) ---
    // Fortran DATA implies: B(1,1)=1.0, B(2,1)=0.0, B(3,1)=0.0, B(1,2)=0.0, B(2,2)=2.0, B(3,2)=1.0
     std::vector<double> B = {
         1.0,  0.0,  0.0, // Col 1
         0.0,  2.0,  1.0  // Col 2
    };

    // Output variables
    int NCONT;
    int INDCON;
    std::vector<int> NBLK(N);      // Block sizes
    std::vector<double> Z(N * N);  // Transformation matrix (LDZ*N = N*N)
    std::vector<double> TAU(N);    // Householder coefficients

    // --- Call the SLICOT C wrapper with ROW_MAJOR = 0 ---
    int INFO = slicot_ab01nd(JOBZ, N, M, A.data(), LDA, B.data(), LDB,
                            &NCONT, &INDCON, NBLK.data(), Z.data(), LDZ,
                            TAU.data(), TOL, ROW_MAJOR);

    // --- Check results ---
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab01nd returned error code: " << INFO;

    // --- Expected results from documentation RESULTS section ---
    int expected_NCONT = 2;
    int expected_INDCON = 1; // Doc says 1. Structure has 1 subdiag block A21.
    std::vector<int> expected_NBLK = { 2, 0, 0 };  // Only first element should be 2

    // Expected transformed A matrix (Acont, 2x2 - Column-major)
    double expected_A_cont_cm[4] = {
        -3.0000,  0.0000, // Col 1
         2.2361, -1.0000  // Col 2
    };

    // Expected transformed B matrix (Bcont, 2x2 - Column-major)
    double expected_B_cont_cm[4] = {
         0.0000,  1.0000, // Col 1
        -2.2361,  0.0000  // Col 2
    };

    // Expected transformation matrix Z (Column-major, 3x3)
    double expected_Z_values_cm[9] = {
         0.0000, -0.8944, -0.4472, // Col 1
         1.0000,  0.0000,  0.0000, // Col 2
         0.0000, -0.4472,  0.8944  // Col 3
    };

    // --- Verification of results ---
    EXPECT_EQ(NCONT, expected_NCONT) << "NCONT value doesn't match expected";
    EXPECT_EQ(INDCON, expected_INDCON) << "INDCON value doesn't match expected";
    EXPECT_EQ(NBLK[0], expected_NBLK[0]) << "First block size doesn't match expected";

    // --- Check controllable part of A (Column-major indexing) ---
    // Wrapper modifies A in place. Check top-left NCONT x NCONT block.
    for (int j = 0; j < NCONT; j++) { // col
        for (int i = 0; i < NCONT; i++) { // row
            // Access A[i, j] using column-major index: i + j * LDA
            // Use EXPECT_NEAR directly to allow message streaming with <<
            EXPECT_NEAR(A[i + j*LDA], expected_A_cont_cm[i + j*NCONT], check_tol)
                << "A[" << i << "," << j << "] doesn't match expected column-major result";
        }
    }

    // --- Check controllable part of B (Column-major indexing) ---
    // Wrapper modifies B in place. Check top NCONT x M block.
    for (int j = 0; j < M; j++) { // col
        for (int i = 0; i < NCONT; i++) { // row
            // Access B[i, j] using column-major index: i + j * LDB
            // Use EXPECT_NEAR directly to allow message streaming with <<
             EXPECT_NEAR(B[i + j*LDB], expected_B_cont_cm[i + j*M], check_tol) // Use M for expected indexing
                << "B[" << i << "," << j << "] doesn't match expected column-major result";
        }
    }

    // --- Verify transformation matrix Z (Column-major indexing) when applicable ---
    if (JOBZ == 'N') {
        // Z should not be computed or modified for 'N'
    } else if (JOBZ == 'I' || JOBZ == 'F') {
        // For 'I' or 'F', Z should contain the transformation in column-major format
        for (int j = 0; j < N; j++) { // col
            for (int i = 0; i < N; i++) { // row
                // Access Z[i, j] using column-major index: i + j * LDZ
                // Use EXPECT_NEAR directly to allow message streaming with <<
                EXPECT_NEAR(Z[i + j*LDZ], expected_Z_values_cm[i + j*N], check_tol)
                    << "Z[" << i << "," << j << "] doesn't match expected column-major result";
            }
        }
    }
}

// Instantiate the column-major test case with different JOBZ values
INSTANTIATE_TEST_SUITE_P(
    AB01NDJobzColMajor, // Renamed Suite P
    AB01NDTestColMajor,
    ::testing::Values('N', 'I', 'F') // Test all JOBZ options
);


// =============================================================================
// Test Fixture and Test Case for Row-Major Input (ROW_MAJOR = 1)
// =============================================================================

// Define a fixture for the row-major test
class AB01NDTestRowMajor : public ::testing::TestWithParam<char> {
     // No specific setup needed here for now
};


// Test case based on the example in the AB01ND documentation, using ROW_MAJOR = 1
TEST_P(AB01NDTestRowMajor, DocExampleRowMajor) {
    // --- Input data from AB01ND documentation example ---
    char JOBZ = GetParam();  // Test all three values: 'N', 'I', 'F'
    int N = 3;        // Order of the system
    int M = 2;        // Number of inputs
    double TOL = 0.0; // Use default tolerance
    int ROW_MAJOR = 1; // Input arrays are row-major
    double check_tol = 1e-4; // Tolerance for checking results

    // --- Leading Dimensions for Row-Major ---
    // LDA = number of columns in A for row-major input (N)
    // LDB = number of columns in B for row-major input (M)
    // LDZ = number of columns in Z for row-major input (N)
    int LDA = N;
    int LDB = M;
    int LDZ = N;

    // --- Input Matrix A (Row-major order) ---
    // Transposed from the Fortran column-major data
    std::vector<double> A = {
        -1.0,  0.0,  0.0, // Row 1
        -2.0, -2.0, -2.0, // Row 2
        -1.0,  0.0, -3.0  // Row 3
    };

    // --- Input Matrix B (Row-major order) ---
    // Transposed from the Fortran column-major data
    std::vector<double> B = {
         1.0,  0.0,  // Row 1
         0.0,  2.0,  // Row 2
         0.0,  1.0   // Row 3
    };


    // Output variables
    int NCONT;
    int INDCON;
    std::vector<int> NBLK(N);      // Block sizes
    std::vector<double> Z(N * N);  // Transformation matrix (N rows * LDZ cols = N*N)
    std::vector<double> TAU(N);    // Householder coefficients

    // --- Call the SLICOT C wrapper with ROW_MAJOR = 1 ---
    // The wrapper expects row-major A, B, Z and modifies them in place.
    int INFO = slicot_ab01nd(JOBZ, N, M, A.data(), LDA, B.data(), LDB,
                            &NCONT, &INDCON, NBLK.data(), Z.data(), LDZ,
                            TAU.data(), TOL, ROW_MAJOR);

    // --- Check results ---
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab01nd returned error code: " << INFO;

    // --- Expected results from documentation RESULTS section ---
    int expected_NCONT = 2;
    int expected_INDCON = 1;
    std::vector<int> expected_NBLK = { 2, 0, 0 };  // Only first element should be 2

    // Expected transformed A matrix (Acont, 2x2, Row-major)
    // Transposed from the expected column-major result
     double expected_A_cont_rm[4] = {
        -3.0000,  2.2361, // Row 1
         0.0000, -1.0000  // Row 2
    };

    // Expected transformed B matrix (Bcont, NCONT x M = 2x2, Row-major)
    // Transposed from the expected column-major result
    double expected_B_cont_rm[4] = {
         0.0000, -2.2361, // Row 1 (NCONT rows)
         1.0000,  0.0000  // Row 2 (NCONT rows)
    };

    // Expected transformation matrix Z (Row-major, 3x3)
    // Transposed from the expected column-major result
    double expected_Z_values_rm[9] = {
         0.0000,  1.0000,  0.0000, // Row 1
        -0.8944,  0.0000, -0.4472, // Row 2
        -0.4472,  0.0000,  0.8944  // Row 3
    };

    // --- Verification of results ---
    EXPECT_EQ(NCONT, expected_NCONT) << "NCONT value doesn't match expected";
    EXPECT_EQ(INDCON, expected_INDCON) << "INDCON value doesn't match expected";
    EXPECT_EQ(NBLK[0], expected_NBLK[0]) << "First block size doesn't match expected";

    // --- Check controllable part of A (Row-major indexing) ---
    // Wrapper modified A in place. Check top-left NCONT x NCONT block.
    for (int i = 0; i < NCONT; i++) { // row
        for (int j = 0; j < NCONT; j++) { // col
            // Access A[i, j] using row-major index: row * num_cols + col
            // LDA is the number of columns for row-major A (which is N)
            // Use EXPECT_NEAR directly to allow message streaming with <<
            EXPECT_NEAR(A[i * LDA + j], expected_A_cont_rm[i * NCONT + j], check_tol)
                << "A[" << i << "," << j << "] doesn't match expected row-major result";
        }
    }

    // --- Check controllable part of B (Row-major indexing) ---
    // Wrapper modified B in place. Check top NCONT x M block.
     for (int i = 0; i < NCONT; i++) { // row
        for (int j = 0; j < M; j++) { // col
            // Access B[i, j] using row-major index: row * num_cols + col
            // LDB is the number of columns for row-major B (which is M)
            // Use EXPECT_NEAR directly to allow message streaming with <<
            EXPECT_NEAR(B[i * LDB + j], expected_B_cont_rm[i * M + j], check_tol)
                << "B[" << i << "," << j << "] doesn't match expected row-major result";
        }
    }

    // --- Verify transformation matrix Z (Row-major indexing) when applicable ---
    if (JOBZ == 'N') {
        // Z should not be computed or modified for 'N'
    } else if (JOBZ == 'I' || JOBZ == 'F') {
        // For 'I' or 'F', Z should contain the transformation in row-major format
        for (int i = 0; i < N; i++) { // row
            for (int j = 0; j < N; j++) { // col
                // Access Z[i, j] using row-major index: row * num_cols + col
                // LDZ is the number of columns for row-major Z (which is N)
                // Use EXPECT_NEAR directly to allow message streaming with <<
                 EXPECT_NEAR(Z[i * LDZ + j], expected_Z_values_rm[i * N + j], check_tol)
                    << "Z[" << i << "," << j << "] doesn't match expected row-major result";
            }
        }
    }
}

// Instantiate the row-major test case with different JOBZ values
INSTANTIATE_TEST_SUITE_P(
    AB01NDJobzRowMajor, // Keep distinct name
    AB01NDTestRowMajor,
    ::testing::Values('N', 'I', 'F') // Test all JOBZ options
);
