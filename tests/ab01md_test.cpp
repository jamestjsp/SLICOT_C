#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits

// Include the specific header for the function being tested
#include "ab01md.h" // Assuming this is the correct header for slicot_ab01md

// Define a fixture for common setup, if needed
class AB01MDTest : public ::testing::Test {
protected:
    // You can add setup/teardown logic here if tests share data/state
    void SetUp() override {
        // Example setup
    }

    void TearDown() override {
        // Example teardown
    }

    // Helper function for comparing doubles with tolerance
    // Uses a default tolerance based on machine epsilon
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

// Test case based on the example in the AB01MD documentation
TEST_F(AB01MDTest, DocExample) {
    // --- Input data from AB01MD documentation example ---
    char JOBZ = 'I'; // Compute Z explicitly, initialize Z to identity
    int N = 3;      // Order of the system
    double TOL = 0.0; // Use default tolerance
    int LDA = N;    // Leading dimension of A (using N for simplicity in column-major)
    int LDZ = N;    // Leading dimension of Z (using N for simplicity in column-major)
    int ROW_MAJOR = 0; // Input arrays are column-major (Fortran style)

    // Matrix A (Column-major order as in Fortran example)
    std::vector<double> A = {
        1.0,  4.0,  0.0,  // Col 1
        2.0, -1.0,  0.0,  // Col 2
        0.0,  0.0,  1.0   // Col 3
    };

    // Matrix B
    std::vector<double> B = { 1.0, 0.0, 1.0 };

    // Output variables
    int NCONT;
    std::vector<double> Z(N * N); // Allocate space for Z
    std::vector<double> TAU(N);   // Allocate space for TAU

    // --- Call the SLICOT C wrapper --- 
    // The wrapper handles workspace allocation internally.
    int INFO = slicot_ab01md(JOBZ, N, A.data(), LDA, B.data(), &NCONT, Z.data(), LDZ, TAU.data(), TOL, ROW_MAJOR);

    // --- Check results --- 
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab01md returned error code: " << INFO;

    // --- Expected results from documentation ---
    int expected_NCONT = 3;
    // Expected A (controllable part, column-major)
    std::vector<double> expected_A_cont = {
        1.0000,  2.8284,  0.0000, // Col 1
        1.4142, -1.0000,  1.4142, // Col 2
        0.0000,  2.8284,  1.0000  // Col 3
    };
    // Expected B (controllable part)
    std::vector<double> expected_B_cont = { -1.4142, 0.0000, 0.0000 };
    // Expected Z (column-major)
    std::vector<double> expected_Z = {
       -0.7071,  0.0000, -0.7071, // Col 1
        0.0000, -1.0000,  0.0000, // Col 2
       -0.7071,  0.0000,  0.7071  // Col 3
    };

    // --- Verification ---
    EXPECT_EQ(NCONT, expected_NCONT);

    // Verify the controllable part of A (NCONT x NCONT)
    // Extract the NCONT x NCONT submatrix from the potentially modified A
    std::vector<double> A_cont_result;
    if (NCONT > 0) {
        A_cont_result.resize(NCONT * NCONT);
        for (int j = 0; j < NCONT; ++j) {
            for (int i = 0; i < NCONT; ++i) {
                A_cont_result[i + j * NCONT] = A[i + j * LDA]; // Column-major access
            }
        }
        ExpectVectorNear(A_cont_result, expected_A_cont, 1e-4); // Use tolerance from example output
    }

    // Verify the controllable part of B (first NCONT elements)
    std::vector<double> B_cont_result(B.begin(), B.begin() + NCONT);
    ExpectVectorNear(B_cont_result, expected_B_cont, 1e-4); // Use tolerance from example output

    // Verify Z
    ExpectVectorNear(Z, expected_Z, 1e-4); // Use tolerance from example output
}

// Add more tests for different parameters (e.g., JOBZ='N', JOBZ='F', different N, TOL) and edge cases
