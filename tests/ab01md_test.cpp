#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits

// Include the specific header for the function being tested
#include "ab01md.h" // Header for slicot_ab01md

// Removed helper functions - using EXPECT_NEAR directly allows message streaming <<

// =============================================================================
// Test Fixture and Test Case for Column-Major Input (ROW_MAJOR = 0)
// =============================================================================

// Define a fixture for the column-major test, parameterized by JOBZ
class AB01MDTestColMajor : public ::testing::TestWithParam<char> {
    // No specific setup needed here for now
};

// Test case based on the example in the AB01MD documentation (Column-Major)
TEST_P(AB01MDTestColMajor, DocExampleColMajor) {
    // --- Input data from AB01MD documentation example ---
    char JOBZ = GetParam(); // Parameterized JOBZ ('N', 'I', 'F')
    int N = 3;      // Order of the system
    double TOL = 0.0; // Use default tolerance
    int ROW_MAJOR = 0; // Input arrays are column-major (Fortran style)
    double check_tol = 1e-4; // Tolerance for checking results based on example precision

    // --- Leading Dimensions for Column-Major (Fortran style) ---
    int LDA = N;    // LDA >= MAX(1,N)
    int LDZ = N;    // LDZ >= MAX(1,N) if JOBZ != 'N', else >= 1

    // --- Input Matrix A (Column-major order as in Fortran example) ---
    std::vector<double> A = {
        1.0,  4.0,  0.0,  // Col 1
        2.0, -1.0,  0.0,  // Col 2
        0.0,  0.0,  1.0   // Col 3
    };

    // --- Input Vector B ---
    std::vector<double> B = { 1.0, 0.0, 1.0 };

    // Output variables
    int NCONT;
    std::vector<double> Z(N * N); // Allocate space for Z (LDZ*N = N*N)
    std::vector<double> TAU(N);   // Allocate space for TAU

    // --- Call the SLICOT C wrapper with ROW_MAJOR = 0 ---
    int INFO = slicot_ab01md(JOBZ, N, A.data(), LDA, B.data(), &NCONT, Z.data(), LDZ, TAU.data(), TOL, ROW_MAJOR);

    // --- Check results ---
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab01md returned error code: " << INFO;

    // --- Expected results from documentation ---
    int expected_NCONT = 3;
    // Expected A (controllable part, 3x3, column-major)
    double expected_A_cont_cm[9] = {
        1.0000,  2.8284,  0.0000, // Col 1
        1.4142, -1.0000,  1.4142, // Col 2
        0.0000,  2.8284,  1.0000  // Col 3
    };
    // Expected B (controllable part, first NCONT elements)
    double expected_B_cont[3] = { -1.4142, 0.0000, 0.0000 };
    // Expected Z (column-major, 3x3)
    double expected_Z_cm[9] = {
       -0.7071,  0.0000, -0.7071, // Col 1
        0.0000, -1.0000,  0.0000, // Col 2
       -0.7071,  0.0000,  0.7071  // Col 3
    };

    // --- Verification ---
    EXPECT_EQ(NCONT, expected_NCONT) << "NCONT mismatch";

    // Verify the controllable part of A (NCONT x NCONT) using column-major indexing
    if (NCONT > 0) {
        for (int j = 0; j < NCONT; ++j) { // col
            for (int i = 0; i < NCONT; ++i) { // row
                // Access A[i, j] using column-major index: i + j * LDA
                EXPECT_NEAR(A[i + j * LDA], expected_A_cont_cm[i + j * NCONT], check_tol)
                    << "A[" << i << "," << j << "] mismatch (Column-Major)";
            }
        }
    }

    // Verify the controllable part of B (first NCONT elements)
    for (int i = 0; i < NCONT; ++i) {
         EXPECT_NEAR(B[i], expected_B_cont[i], check_tol)
             << "B[" << i << "] mismatch (Column-Major)";
    }
    // Verify remaining elements of B are zero (optional, but good practice)
    for (int i = NCONT; i < N; ++i) {
        EXPECT_NEAR(B[i], 0.0, check_tol)
             << "B[" << i << "] should be zero (Column-Major)";
    }


    // Verify Z matrix using column-major indexing
    if (JOBZ == 'I' || JOBZ == 'F') {
        for (int j = 0; j < N; ++j) { // col
            for (int i = 0; i < N; ++i) { // row
                 // Access Z[i, j] using column-major index: i + j * LDZ
                 EXPECT_NEAR(Z[i + j * LDZ], expected_Z_cm[i + j * N], check_tol)
                     << "Z[" << i << "," << j << "] mismatch (Column-Major)";
            }
        }
    }
}

// Instantiate the column-major test case with different JOBZ values
INSTANTIATE_TEST_SUITE_P(
    AB01MDJobzColMajor,
    AB01MDTestColMajor,
    ::testing::Values('N', 'I', 'F') // Test all JOBZ options
);


// =============================================================================
// Test Fixture and Test Case for Row-Major Input (ROW_MAJOR = 1)
// =============================================================================

// Define a fixture for the row-major test, parameterized by JOBZ
class AB01MDTestRowMajor : public ::testing::TestWithParam<char> {
     // No specific setup needed here for now
};


// Test case based on the example in the AB01MD documentation, using ROW_MAJOR = 1
TEST_P(AB01MDTestRowMajor, DocExampleRowMajor) {
    // --- Input data from AB01MD documentation example ---
    char JOBZ = GetParam(); // Parameterized JOBZ ('N', 'I', 'F')
    int N = 3;      // Order of the system
    double TOL = 0.0; // Use default tolerance
    int ROW_MAJOR = 1; // Input arrays are row-major
    double check_tol = 1e-4; // Tolerance for checking results

    // --- Leading Dimensions for Row-Major ---
    // LDA = number of columns in A for row-major input (N)
    // LDZ = number of columns in Z for row-major input (N)
    int LDA = N;
    int LDZ = N;

    // --- Input Matrix A (Row-major order) ---
    // Transposed from the Fortran column-major data
    std::vector<double> A = {
        1.0,  2.0,  0.0,  // Row 1
        4.0, -1.0,  0.0,  // Row 2
        0.0,  0.0,  1.0   // Row 3
    };

    // --- Input Vector B (same as column-major) ---
    std::vector<double> B = { 1.0, 0.0, 1.0 };

    // Output variables
    int NCONT;
    std::vector<double> Z(N * N); // Allocate space for Z (N rows * LDZ cols = N*N)
    std::vector<double> TAU(N);   // Allocate space for TAU

    // --- Call the SLICOT C wrapper with ROW_MAJOR = 1 ---
    // The wrapper expects row-major A, Z and modifies A, B, Z in place.
    int INFO = slicot_ab01md(JOBZ, N, A.data(), LDA, B.data(), &NCONT, Z.data(), LDZ, TAU.data(), TOL, ROW_MAJOR);

    // --- Check results ---
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab01md returned error code: " << INFO;

    // --- Expected results from documentation (transposed to Row-Major) ---
    int expected_NCONT = 3;
    // Expected A (controllable part, 3x3, row-major)
    double expected_A_cont_rm[9] = {
        1.0000,  1.4142,  0.0000, // Row 1
        2.8284, -1.0000,  2.8284, // Row 2
        0.0000,  1.4142,  1.0000  // Row 3
    };
    // Expected B (controllable part, first NCONT elements - same as column-major)
    double expected_B_cont[3] = { -1.4142, 0.0000, 0.0000 };
    // Expected Z (row-major, 3x3)
    double expected_Z_rm[9] = {
       -0.7071,  0.0000, -0.7071, // Row 1
        0.0000, -1.0000,  0.0000, // Row 2
       -0.7071,  0.0000,  0.7071  // Row 3
    };

    // --- Verification ---
    EXPECT_EQ(NCONT, expected_NCONT) << "NCONT mismatch";

    // Verify the controllable part of A (NCONT x NCONT) using row-major indexing
    if (NCONT > 0) {
        for (int i = 0; i < NCONT; ++i) { // row
            for (int j = 0; j < NCONT; ++j) { // col
                 // Access A[i, j] using row-major index: i * LDA + j
                 // LDA is number of columns for row-major A (which is N)
                 EXPECT_NEAR(A[i * LDA + j], expected_A_cont_rm[i * NCONT + j], check_tol)
                     << "A[" << i << "," << j << "] mismatch (Row-Major)";
            }
        }
    }

    // Verify the controllable part of B (first NCONT elements)
    for (int i = 0; i < NCONT; ++i) {
         EXPECT_NEAR(B[i], expected_B_cont[i], check_tol)
             << "B[" << i << "] mismatch (Row-Major)";
    }
     // Verify remaining elements of B are zero (optional, but good practice)
    for (int i = NCONT; i < N; ++i) {
        EXPECT_NEAR(B[i], 0.0, check_tol)
             << "B[" << i << "] should be zero (Row-Major)";
    }

    // Verify Z matrix using row-major indexing
    if (JOBZ == 'I' || JOBZ == 'F') {
        for (int i = 0; i < N; ++i) { // row
            for (int j = 0; j < N; ++j) { // col
                 // Access Z[i, j] using row-major index: i * LDZ + j
                 // LDZ is number of columns for row-major Z (which is N)
                 EXPECT_NEAR(Z[i * LDZ + j], expected_Z_rm[i * N + j], check_tol)
                     << "Z[" << i << "," << j << "] mismatch (Row-Major)";
            }
        }
    }
}

// Instantiate the row-major test case with different JOBZ values
INSTANTIATE_TEST_SUITE_P(
    AB01MDJobzRowMajor,
    AB01MDTestRowMajor,
    ::testing::Values('N', 'I', 'F') // Test all JOBZ options
);
