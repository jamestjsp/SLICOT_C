#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits

// Include the specific header for the function being tested
#include "ab01nd.h" // Header for slicot_ab01nd

// Define a fixture for common setup
class AB01NDTest : public ::testing::TestWithParam<char> {
protected:
    // Helper function for comparing doubles with tolerance
    void ExpectNear(double val1, double val2, double tol = std::numeric_limits<double>::epsilon() * 100) {
        EXPECT_NEAR(val1, val2, tol);
    }

    // Helper function to compare vectors of doubles
    void ExpectVectorNear(const std::vector<double>& vec1, const std::vector<double>& vec2, double tol = std::numeric_limits<double>::epsilon() * 100) {
        ASSERT_EQ(vec1.size(), vec2.size());
        for (size_t i = 0; i < vec1.size(); ++i) {
            EXPECT_NEAR(vec1[i], vec2[i], tol) << "Mismatch at index " << i;
        }
    }
};

// Test case based on the example in the AB01ND documentation
TEST_P(AB01NDTest, DocExample) {
    // --- Input data from AB01ND documentation example ---
    char JOBZ = GetParam();  // Test all three values: 'N', 'I', 'F'
    int N = 3;        // Order of the system
    int M = 2;        // Number of inputs
    double TOL = 0.0; // Use default tolerance
    int LDA = N;      // Leading dimension of A
    int LDB = N;      // Leading dimension of B
    int LDZ = N;      // Leading dimension of Z
    int ROW_MAJOR = 0; // Input arrays are column-major (Fortran style)

    // Matrix A (Column-major order as in Fortran example DATA section)
    std::vector<double> A = {
        -1.0, -2.0, -1.0, // Col 1
         0.0, -2.0,  0.0, // Col 2
         0.0, -2.0, -3.0  // Col 3
    };

    // Matrix B (Column-major order as in Fortran example DATA section)
    std::vector<double> B = {
         1.0,  0.0, // Col 1
         0.0,  2.0, // Col 2
         0.0,  1.0  // Col 3
    };

    // Create a copy of the original data to compare with documentation results later
    std::vector<double> A_orig = A;
    std::vector<double> B_orig = B;

    // Output variables
    int NCONT;
    int INDCON;
    std::vector<int> NBLK(N);      // Block sizes
    std::vector<double> Z(N * N);  // Transformation matrix
    std::vector<double> TAU(N);    // Householder coefficients

    // --- Call the SLICOT C wrapper --- 
    int INFO = slicot_ab01nd(JOBZ, N, M, A.data(), LDA, B.data(), LDB, 
                            &NCONT, &INDCON, NBLK.data(), Z.data(), LDZ, 
                            TAU.data(), TOL, ROW_MAJOR);

    // --- Check results --- 
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab01nd returned error code: " << INFO;

    // --- Expected results from documentation RESULTS section ---
    int expected_NCONT = 2;
    int expected_INDCON = 1;
    std::vector<int> expected_NBLK = { 2, 0, 0 };  // Only first element should be 2
    
    // Expected transformed A matrix (2x2 - controllable part)
    double expected_A_cont[4] = {
        -3.0000, 2.2361,
         0.0000, -1.0000
    };
    
    // Expected transformed B matrix (2x2)
    double expected_B_cont[4] = {
         0.0000, -2.2361,
         1.0000, 0.0000
    };

    // Expected transformation matrix Z (column-major, 3x3)
    double expected_Z_values[9] = {
         0.0000, -0.8944, -0.4472,
         1.0000,  0.0000,  0.0000,
         0.0000, -0.4472,  0.8944
    };
    
    // --- Verification of results ---
    EXPECT_EQ(NCONT, expected_NCONT) << "NCONT value doesn't match expected";
    EXPECT_EQ(INDCON, expected_INDCON) << "INDCON value doesn't match expected";
    EXPECT_EQ(NBLK[0], expected_NBLK[0]) << "First block size doesn't match expected";

    // Check controllable part of A
    for (int i = 0; i < expected_NCONT; i++) {
        for (int j = 0; j < expected_NCONT; j++) {
            EXPECT_NEAR(A[i + j*LDA], expected_A_cont[i + j*expected_NCONT], 1e-4)
                << "A[" << i << "," << j << "] doesn't match expected";
        }
    }

    // Check controllable part of B
    for (int i = 0; i < expected_NCONT; i++) {
        for (int j = 0; j < M; j++) {
            EXPECT_NEAR(B[i + j*LDB], expected_B_cont[i + j*expected_NCONT], 1e-4)
                << "B[" << i << "," << j << "] doesn't match expected";
        }
    }

    // Verify transformation matrix Z when applicable
    if (JOBZ == 'N') {
        // Z should not be computed for 'N'
    } else if (JOBZ == 'I' || JOBZ == 'F') {
        // For 'I' or 'F', Z should contain the transformation
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                EXPECT_NEAR(Z[i + j*LDZ], expected_Z_values[i + j*N], 1e-4)
                    << "Z[" << i << "," << j << "] doesn't match expected";
            }
        }
    }
}

// Run the test with different JOBZ values
INSTANTIATE_TEST_CASE_P(
    AB01NDJobz,
    AB01NDTest,
    ::testing::Values('N', 'I', 'F')
);
