#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits
#include <algorithm> // Required for std::min, std::max
#include <string> // Required for std::to_string
#include <iostream> // For printing matrix values during debugging (optional)
#include <iomanip>  // For std::setprecision, std::fixed

// Include the specific header for the function being tested
#include "ab13bd.h" // Header for slicot_ab13bd
#include "slicot_utils.h" // For slicot_transpose_to_c if needed by test helpers

// =============================================================================
// Test Fixture and Test Case for Column-Major Input (ROW_MAJOR = 0)
// =============================================================================

class AB13BDTestColMajor : public ::testing::Test {
protected:
    // Define common variables for tests
    int N = 7, M = 2, P = 3;
    char DICO = 'C', JOBN = 'L'; // L2-norm for Continuous system
    double TOL = 1.0e-10;
    int ROW_MAJOR = 0;
    // Use a tolerance appropriate for the expected result's precision
    double check_tol = 1e-5;

    // Expected results from documentation
    double expected_norm = 7.93948; // L2-norm from example
    int expected_IWARN = 0; // Implied
    int expected_INFO = 0;

     // Input Matrices (Column-major order, N=7, M=2, P=3)
    std::vector<double> A_in = {
        -0.04165, -5.2100,  0.0000,  0.5450,  0.0000,  0.0000,  0.0000, // Col 1
         0.0000, -12.500,  3.3300,  0.0000,  0.0000,  0.0000,  0.0000, // Col 2
         4.9200,  0.0000, -3.3300,  0.0000,  0.0000,  0.0000,  0.0000, // Col 3
         0.4920,  0.0000,  0.0000,  0.0000, -0.4920,  0.0000,  0.0000, // Col 4
         0.0000,  0.0000,  0.0000,  0.0545,  0.004165, 0.5210,  0.0000, // Col 5
         0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -12.500,  3.3300, // Col 6
         0.0000,  0.0000,  0.0000,  0.0000,  4.9200,  0.0000, -3.3300  // Col 7
    };
    std::vector<double> B_in = {
         0.0000, 12.500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, // Col 1
         0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 12.500,  0.0000  // Col 2
    };
     std::vector<double> C_in = {
        1.0000, 0.0000, 0.0000, // Col 1
        0.0000, 0.0000, 0.0000, // Col 2
        0.0000, 0.0000, 0.0000, // Col 3
        0.0000, 1.0000, 0.0000, // Col 4
        0.0000, 0.0000, 1.0000, // Col 5
        0.0000, 0.0000, 0.0000, // Col 6
        0.0000, 0.0000, 0.0000  // Col 7
    };
     std::vector<double> D_in = {
        0.0, 0.0, 0.0, // Col 1
        0.0, 0.0, 0.0  // Col 2
    }; // Zero matrix
};

TEST_F(AB13BDTestColMajor, DocExampleColMajor) {

    // --- Leading Dimensions for Column-Major (Fortran style) ---
    int LDA = N; // >= N
    int LDB = N; // >= N
    int LDC = P; // >= P
    int LDD = P; // >= P

    // Create copies of input matrices as they are overwritten
    std::vector<double> A = A_in;
    std::vector<double> B = B_in;
    std::vector<double> C = C_in;
    std::vector<double> D = D_in;

    ASSERT_EQ(A.size(), LDA * N);
    ASSERT_EQ(B.size(), LDB * M);
    ASSERT_EQ(C.size(), LDC * N);
    ASSERT_EQ(D.size(), LDD * M);

    // Output variables
    int nq_out;
    int iwarn_out;
    int info_out;

    // --- Call the SLICOT C wrapper ---
    double norm_value = slicot_ab13bd(DICO, JOBN, N, M, P,
                                      A.data(), LDA, B.data(), LDB, C.data(), LDC, D.data(), LDD,
                                      &nq_out, TOL, &iwarn_out, &info_out,
                                      ROW_MAJOR);

    // --- Check results ---
    ASSERT_EQ(info_out, expected_INFO) << "SLICOT routine ab13bd returned error code: " << info_out;

    // Check scalar outputs
    EXPECT_EQ(iwarn_out, expected_IWARN) << "Warning flag IWARN mismatch";
    EXPECT_GE(nq_out, 0) << "Output order NQ should be non-negative"; // Basic check

    // Check the computed norm
    EXPECT_NEAR(norm_value, expected_norm, check_tol) << "Computed norm mismatch";

    // Cannot check modified A, B, C, D as expected values are not in the doc example.
}

// =============================================================================
// Test Fixture and Test Case for Row-Major Input (ROW_MAJOR = 1)
// =============================================================================

// Inherit expected values and parameters from ColMajor fixture
class AB13BDTestRowMajor : public AB13BDTestColMajor {
public:
    AB13BDTestRowMajor() {
        ROW_MAJOR = 1; // Override ROW_MAJOR for this fixture
    }
};


TEST_F(AB13BDTestRowMajor, DocExampleRowMajor) {

    // --- Leading Dimensions for Row-Major (C style) ---
    int LDA = N; // Cols for A (NxN)
    int LDB = M; // Cols for B (NxM)
    int LDC = N; // Cols for C (PxN)
    int LDD = M; // Cols for D (PxM)

    // --- Input Matrices (Row-major order, N=7, M=2, P=3) ---
    // Transpose the column-major inputs defined in the fixture
    std::vector<double> A(N * LDA);
    std::vector<double> B(N * LDB);
    std::vector<double> C(P * LDC);
    std::vector<double> D(P * LDD);

    slicot_transpose_to_c(A_in.data(), A.data(), N, N, sizeof(double));
    slicot_transpose_to_c(B_in.data(), B.data(), N, M, sizeof(double));
    slicot_transpose_to_c(C_in.data(), C.data(), P, N, sizeof(double));
    slicot_transpose_to_c(D_in.data(), D.data(), P, M, sizeof(double));

    ASSERT_EQ(A.size(), N * LDA);
    ASSERT_EQ(B.size(), N * LDB);
    ASSERT_EQ(C.size(), P * LDC);
    ASSERT_EQ(D.size(), P * LDD);

    // Output variables
    int nq_out;
    int iwarn_out;
    int info_out;

    // --- Call the SLICOT C wrapper ---
     double norm_value = slicot_ab13bd(DICO, JOBN, N, M, P,
                                      A.data(), LDA, B.data(), LDB, C.data(), LDC, D.data(), LDD,
                                      &nq_out, TOL, &iwarn_out, &info_out,
                                      ROW_MAJOR);

    // --- Check results ---
    ASSERT_EQ(info_out, expected_INFO) << "SLICOT routine ab13bd returned error code: " << info_out;

    // Check scalar outputs
    EXPECT_EQ(iwarn_out, expected_IWARN) << "Warning flag IWARN mismatch";
    EXPECT_GE(nq_out, 0) << "Output order NQ should be non-negative"; // Basic check

    // Check the computed norm
    EXPECT_NEAR(norm_value, expected_norm, check_tol) << "Computed norm mismatch";

    // Cannot check modified A, B, C, D as expected values are not in the doc example.
    // If needed for debugging, could print the modified matrices A, B, C, D here.
}
