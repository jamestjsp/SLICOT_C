#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits

// Include the specific header for the function being tested
#include "ab01od.h" // Header for slicot_ab01od

// Define a fixture for common setup
class AB01ODTest : public ::testing::Test {
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

// Test case based on the example in the AB01OD documentation (forward stage only)
TEST_F(AB01ODTest, DocExampleForwardStage) {
    // --- Input data from AB01OD documentation example ---
    char STAGES = 'F';  // Perform the forward stage only
    char JOBU = 'N';    // Do not form U
    char JOBV = 'N';    // Do not form V
    int N = 5;          // Order of the system
    int M = 2;          // Number of inputs
    double TOL = 0.0;   // Use default tolerance
    int LDA = N;        // Leading dimension of A (using N for simplicity in column-major)
    int LDB = N;        // Leading dimension of B (using N for simplicity in column-major)
    int LDU = 1;        // Leading dimension of U (not used, JOBU='N')
    int LDV = 1;        // Leading dimension of V (not used, JOBV='N')
    int ROW_MAJOR = 0;  // Input arrays are column-major (Fortran style)

    // Matrix A (Column-major order as in Fortran example)
    std::vector<double> A = {
        17.0, 23.0, 4.0, 10.0, 11.0,  // Col 1
        24.0, 5.0, 6.0, 12.0, 18.0,   // Col 2
        1.0, 7.0, 13.0, 19.0, 25.0,   // Col 3
        8.0, 14.0, 20.0, 21.0, 2.0,   // Col 4
        15.0, 16.0, 22.0, 3.0, 9.0    // Col 5
    };

    // Matrix B (Column-major order as in Fortran example)
    std::vector<double> B = {
        -1.0, 4.0, -9.0, 16.0, -25.0,  // Col 1
        -4.0, 9.0, -16.0, 25.0, -36.0  // Col 2
    };

    // Create dummy arrays for U and V (not used with JOBU='N' and JOBV='N')
    std::vector<double> U(1);
    std::vector<double> V(1);

    // Output variables
    int NCONT;
    int INDCON;
    std::vector<int> KSTAIR(N);

    // --- Call the SLICOT C wrapper --- 
    int INFO = slicot_ab01od(STAGES, JOBU, JOBV, N, M,
                            A.data(), LDA, B.data(), LDB,
                            U.data(), LDU, V.data(), LDV,
                            &NCONT, &INDCON, KSTAIR.data(),
                            TOL, ROW_MAJOR);

    // --- Check results --- 
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab01od returned error code: " << INFO;

    // Expected results from documentation RESULTS section
    int expected_NCONT = 5;  // All states are controllable in this case 
    int expected_INDCON = 3; // Number of stairs in staircase form
    std::vector<int> expected_KSTAIR = { 2, 2, 1 };  // Dimensions of stairs

    // Expected transformed A matrix (column-major)
    // Allow a bit larger tolerance since example results are rounded
    const double tol = 1e-4;
    std::vector<double> expected_A = {
        12.8848, 4.4741, 14.4576, 0.0000, 0.0000,  // Col 1
        3.2345, -12.5544, 7.6855, 1.4805, 0.0000,  // Col 2
        11.8211, 5.3509, 23.1452, 27.4668, -30.4822,  // Col 3
        3.3758, 5.9403, 26.3872, 22.6564, 0.6745,   // Col 4
        -0.8982, 1.4360, -29.9557, -0.0072, 18.8680  // Col 5
    };

    // Expected transformed B matrix (column-major)
    std::vector<double> expected_B = {
        31.1199, 3.2480, 0.0000, 0.0000, 0.0000,  // Col 1
        47.6865, 0.0000, 0.0000, 0.0000, 0.0000   // Col 2
    };

    // --- Verification ---
    EXPECT_EQ(NCONT, expected_NCONT);
    EXPECT_EQ(INDCON, expected_INDCON);

    // Check the stair dimensions
    for (int i = 0; i < expected_INDCON; i++) {
        EXPECT_EQ(KSTAIR[i], expected_KSTAIR[i]) << "Mismatch in KSTAIR at index " << i;
    }

    // Check the transformed A and B matrices
    ExpectVectorNear(A, expected_A, tol);
    ExpectVectorNear(B, expected_B, tol);
}

// Test case for both forward and backward stages
TEST_F(AB01ODTest, AllStages) {
    char STAGES = 'A';  // Perform both forward and backward stages
    char JOBU = 'I';    // Form transformation matrix U
    char JOBV = 'I';    // Form transformation matrix V
    int N = 5;
    int M = 2;
    double TOL = 0.0;
    int LDA = N;
    int LDB = N;
    int LDU = N;
    int LDV = M;
    int ROW_MAJOR = 0;  // Column-major

    // Matrix A (same as DocExampleForwardStage)
    std::vector<double> A = {
        17.0, 23.0, 4.0, 10.0, 11.0,  // Col 1
        24.0, 5.0, 6.0, 12.0, 18.0,   // Col 2
        1.0, 7.0, 13.0, 19.0, 25.0,   // Col 3
        8.0, 14.0, 20.0, 21.0, 2.0,   // Col 4
        15.0, 16.0, 22.0, 3.0, 9.0    // Col 5
    };

    // Matrix B (same as DocExampleForwardStage)
    std::vector<double> B = {
        -1.0, 4.0, -9.0, 16.0, -25.0,  // Col 1
        -4.0, 9.0, -16.0, 25.0, -36.0  // Col 2
    };

    // Matrices U and V will be initialized inside the function
    std::vector<double> U(N * N);
    std::vector<double> V(M * M);

    // Output variables
    int NCONT;
    int INDCON;
    std::vector<int> KSTAIR(N);

    // Call slicot_ab01od
    int INFO = slicot_ab01od(STAGES, JOBU, JOBV, N, M,
                            A.data(), LDA, B.data(), LDB,
                            U.data(), LDU, V.data(), LDV,
                            &NCONT, &INDCON, KSTAIR.data(),
                            TOL, ROW_MAJOR);

    // Check SLICOT execution
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab01od returned error code: " << INFO;

    // Verify controllability index is 3
    EXPECT_EQ(INDCON, 3);
    
    // Verify dimensions of the stairs
    std::vector<int> expected_KSTAIR = { 2, 2, 1 };
    for (int i = 0; i < INDCON; i++) {
        EXPECT_EQ(KSTAIR[i], expected_KSTAIR[i]) << "Mismatch in KSTAIR at index " << i;
    }

    // Verify that B has zeros in all rows except the first two
    for (int j = 0; j < M; j++) {
        for (int i = 2; i < N; i++) {
            EXPECT_NEAR(B[i + j*N], 0.0, 1e-10) << "B(" << i << "," << j << ") is not zero";
        }
    }
    
    // Verify that A is block upper Hessenberg with zeroes below the first block subdiagonal
    for (int j = 0; j < N; j++) {
        for (int i = j+3; i < N; i++) {  // Check elements below the first block subdiagonal
            EXPECT_NEAR(A[i + j*N], 0.0, 1e-10) << "A(" << i << "," << j << ") is not zero";
        }
    }
    
    // Verify that U and V are orthogonal (U'*U and V'*V should be identity matrices)
    // This is a crude but effective check
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double dot_product = 0.0;
            for (int k = 0; k < N; k++) {
                dot_product += U[k + i*N] * U[k + j*N];
            }
            if (i == j) {
                EXPECT_NEAR(dot_product, 1.0, 1e-10) << "U'*U(" << i << "," << j << ") is not 1";
            } else {
                EXPECT_NEAR(dot_product, 0.0, 1e-10) << "U'*U(" << i << "," << j << ") is not 0";
            }
        }
    }
    
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double dot_product = 0.0;
            for (int k = 0; k < M; k++) {
                dot_product += V[k + i*M] * V[k + j*M];
            }
            if (i == j) {
                EXPECT_NEAR(dot_product, 1.0, 1e-10) << "V'*V(" << i << "," << j << ") is not 1";
            } else {
                EXPECT_NEAR(dot_product, 0.0, 1e-10) << "V'*V(" << i << "," << j << ") is not 0";
            }
        }
    }
}

// Test case for backward stage only
TEST_F(AB01ODTest, BackwardStage) {
    // First perform forward stage to get inputs for backward stage
    char STAGES_F = 'F';
    char JOBU_F = 'I';
    char JOBV_F = 'N';
    int N = 5;
    int M = 2;
    double TOL = 0.0;
    int LDA = N;
    int LDB = N;
    int LDU = N;
    int LDV = 1;
    int ROW_MAJOR = 0;  // Column-major

    // Initial A and B matrices
    std::vector<double> A = {
        17.0, 23.0, 4.0, 10.0, 11.0,
        24.0, 5.0, 6.0, 12.0, 18.0,
        1.0, 7.0, 13.0, 19.0, 25.0,
        8.0, 14.0, 20.0, 21.0, 2.0,
        15.0, 16.0, 22.0, 3.0, 9.0
    };

    std::vector<double> B = {
        -1.0, 4.0, -9.0, 16.0, -25.0,
        -4.0, 9.0, -16.0, 25.0, -36.0
    };

    std::vector<double> U(N * N);
    std::vector<double> V_dummy(1);

    // Output variables for forward stage
    int NCONT_F;
    int INDCON_F;
    std::vector<int> KSTAIR_F(N);

    // Call forward stage
    int INFO_F = slicot_ab01od(STAGES_F, JOBU_F, JOBV_F, N, M,
                              A.data(), LDA, B.data(), LDB,
                              U.data(), LDU, V_dummy.data(), LDV,
                              &NCONT_F, &INDCON_F, KSTAIR_F.data(),
                              TOL, ROW_MAJOR);

    ASSERT_EQ(INFO_F, 0) << "Forward stage failed with error code: " << INFO_F;

    // Now perform backward stage
    char STAGES_B = 'B';
    char JOBU_B = 'I';
    char JOBV_B = 'I';
    int LDV_B = M;
    std::vector<double> V(M * M);

    // Output variables for backward stage
    int NCONT_B = NCONT_F;  // Pass forward results
    int INDCON_B = INDCON_F;
    std::vector<int> KSTAIR_B = KSTAIR_F;

    // Call backward stage
    int INFO_B = slicot_ab01od(STAGES_B, JOBU_B, JOBV_B, N, M,
                              A.data(), LDA, B.data(), LDB,
                              U.data(), LDU, V.data(), LDV_B,
                              &NCONT_B, &INDCON_B, KSTAIR_B.data(),
                              TOL, ROW_MAJOR);

    ASSERT_EQ(INFO_B, 0) << "Backward stage failed with error code: " << INFO_B;

    // Verify that NCONT and INDCON remain unchanged
    EXPECT_EQ(NCONT_B, NCONT_F);
    EXPECT_EQ(INDCON_B, INDCON_F);
    
    // Verify that KSTAIR remains unchanged
    for (int i = 0; i < INDCON_B; i++) {
        EXPECT_EQ(KSTAIR_B[i], KSTAIR_F[i]) << "KSTAIR changed at index " << i;
    }

    // Verify that B matrix is in upper triangular form for the first block
    // First block should be upper triangular (first KSTAIR[0] rows)
    for (int i = 0; i < KSTAIR_B[0]; i++) {
        for (int j = 0; j < M; j++) {
            if (i > j) {
                EXPECT_NEAR(B[i + j*N], 0.0, 1e-10) << "B(" << i << "," << j << ") is not zero in upper triangular part";
            }
        }
    }
}

// Test with row major storage
TEST_F(AB01ODTest, RowMajorStorage) {
    char STAGES = 'F';
    char JOBU = 'N';
    char JOBV = 'N';
    int N = 5;
    int M = 2;
    double TOL = 0.0;
    int LDA = N;  // For row-major, this is the number of columns
    int LDB = M;  // For row-major, this is the number of columns
    int LDU = 1;
    int LDV = 1;
    int ROW_MAJOR = 1;  // Row-major

    // Matrix A (Row-major order)
    std::vector<double> A = {
        17.0, 24.0, 1.0, 8.0, 15.0,
        23.0, 5.0, 7.0, 14.0, 16.0,
        4.0, 6.0, 13.0, 20.0, 22.0,
        10.0, 12.0, 19.0, 21.0, 3.0,
        11.0, 18.0, 25.0, 2.0, 9.0
    };

    // Matrix B (Row-major order)
    std::vector<double> B = {
        -1.0, -4.0,
        4.0, 9.0,
        -9.0, -16.0,
        16.0, 25.0,
        -25.0, -36.0
    };

    std::vector<double> U(1);
    std::vector<double> V(1);

    int NCONT;
    int INDCON;
    std::vector<int> KSTAIR(N);

    int INFO = slicot_ab01od(STAGES, JOBU, JOBV, N, M,
                            A.data(), LDA, B.data(), LDB,
                            U.data(), LDU, V.data(), LDV,
                            &NCONT, &INDCON, KSTAIR.data(),
                            TOL, ROW_MAJOR);

    ASSERT_EQ(INFO, 0) << "SLICOT routine ab01od failed with row-major storage, error code: " << INFO;

    // Expected results (from DocExampleForwardStage converted to row-major)
    int expected_INDCON = 3;
    std::vector<int> expected_KSTAIR = { 2, 2, 1 };

    // Expected transformed A matrix (row-major)
    const double tol = 1e-4;
    std::vector<double> expected_A = {
        12.8848, 3.2345, 11.8211, 3.3758, -0.8982,
        4.4741, -12.5544, 5.3509, 5.9403, 1.4360,
        14.4576, 7.6855, 23.1452, 26.3872, -29.9557,
        0.0000, 1.4805, 27.4668, 22.6564, -0.0072,
        0.0000, 0.0000, -30.4822, 0.6745, 18.8680
    };

    // Expected transformed B matrix (row-major)
    std::vector<double> expected_B = {
        31.1199, 47.6865,
        3.2480, 0.0000,
        0.0000, 0.0000,
        0.0000, 0.0000,
        0.0000, 0.0000
    };

    EXPECT_EQ(INDCON, expected_INDCON);
    for (int i = 0; i < expected_INDCON; i++) {
        EXPECT_EQ(KSTAIR[i], expected_KSTAIR[i]) << "Mismatch in KSTAIR at index " << i;
    }

    ExpectVectorNear(A, expected_A, tol);
    ExpectVectorNear(B, expected_B, tol);
}

// Test edge case with zero matrices
TEST_F(AB01ODTest, ZeroMatrices) {
    char STAGES = 'F';
    char JOBU = 'N';
    char JOBV = 'N';
    int N = 3;
    int M = 2;
    double TOL = 0.0;
    int LDA = N;
    int LDB = N;
    int LDU = 1;
    int LDV = 1;
    int ROW_MAJOR = 0;

    // Zero matrices
    std::vector<double> A(N * N, 0.0);
    std::vector<double> B(N * M, 0.0);
    std::vector<double> U(1);
    std::vector<double> V(1);

    int NCONT;
    int INDCON;
    std::vector<int> KSTAIR(N);

    int INFO = slicot_ab01od(STAGES, JOBU, JOBV, N, M,
                            A.data(), LDA, B.data(), LDB,
                            U.data(), LDU, V.data(), LDV,
                            &NCONT, &INDCON, KSTAIR.data(),
                            TOL, ROW_MAJOR);

    ASSERT_EQ(INFO, 0) << "SLICOT routine ab01od failed with zero matrices, error code: " << INFO;

    // For zero matrices, expect NCONT=0
    EXPECT_EQ(NCONT, 0);
    // There should be no stairs
    EXPECT_EQ(INDCON, 0);
    
    // A and B should remain zero
    for (size_t i = 0; i < A.size(); i++) {
        EXPECT_NEAR(A[i], 0.0, 1e-10);
    }
    for (size_t i = 0; i < B.size(); i++) {
        EXPECT_NEAR(B[i], 0.0, 1e-10);
    }
}
