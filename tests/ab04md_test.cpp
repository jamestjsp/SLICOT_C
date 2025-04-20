#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits

// Include the specific header for the function being tested
#include "ab04md.h" // Header for slicot_ab04md

// Define a fixture for AB04MD tests
class AB04MDTest : public ::testing::Test {
protected:
    // Helper function for comparing doubles with tolerance
    void ExpectNear(double val1, double val2, double tol = std::numeric_limits<double>::epsilon() * 100) {
        EXPECT_NEAR(val1, val2, tol);
    }

    // Helper function to compare matrices represented as vectors
    void ExpectMatricesNear(const std::vector<double>& mat1, const std::vector<double>& mat2, 
                           int rows, int cols, int ld1, int ld2, bool is_row_major,
                           double tol = 1e-4) {
        ASSERT_EQ(mat1.size(), ld1 * (is_row_major ? rows : cols));
        ASSERT_EQ(mat2.size(), ld2 * (is_row_major ? rows : cols));
        
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                int idx1 = is_row_major ? i * ld1 + j : j * ld1 + i;
                int idx2 = is_row_major ? i * ld2 + j : j * ld2 + i;
                EXPECT_NEAR(mat1[idx1], mat2[idx2], tol) << "Matrices differ at (" << i << "," << j << ")";
            }
        }
    }
    
    // Test data from the SLICOT documentation example
    // These are in column-major format (Fortran style)
    
    // Continuous-time system matrices
    const std::vector<double> Ac = {
        1.0, 0.5,  // Col 1
        0.5, 1.0   // Col 2
    };
    const std::vector<double> Bc = {
        0.0, -1.0, // Col 1
        1.0, 0.0   // Col 2
    };
    const std::vector<double> Cc = {
        -1.0, 0.0, // Col 1
        0.0, 1.0   // Col 2
    };
    const std::vector<double> Dc = {
        1.0, 0.0,  // Col 1
        0.0, -1.0  // Col 2
    };
    
    // Discrete-time system matrices
    const std::vector<double> Ad = {
        -1.0, -4.0, // Col 1
        -4.0, -1.0  // Col 2
    };
    const std::vector<double> Bd = {
        2.8284, 0.0,     // Col 1
        0.0, -2.8284     // Col 2
    };
    const std::vector<double> Cd = {
        0.0, 2.8284,     // Col 1
        -2.8284, 0.0     // Col 2
    };
    const std::vector<double> Dd = {
        -1.0, 0.0,       // Col 1
        0.0, -3.0        // Col 2
    };
    
    // System dimensions
    const int N = 2; // System order
    const int M = 2; // Number of inputs
    const int P = 2; // Number of outputs
    
    // Parameters for bilinear transformation
    const double ALPHA = 1.0;
    const double BETA = 1.0;
};

// Test continuous to discrete transformation
TEST_F(AB04MDTest, ContinuousToDiscrete) {
    // Copy the continuous time system matrices
    std::vector<double> A = Ac;
    std::vector<double> B = Bc;
    std::vector<double> C = Cc;
    std::vector<double> D = Dc;

    int LDA = N;
    int LDB = N;
    int LDC = P;
    int LDD = P;
    int ROW_MAJOR = 0; // Using column-major (Fortran-style) storage
    
    // Call the SLICOT wrapper function to convert continuous to discrete
    int INFO = slicot_ab04md('C', N, M, P, ALPHA, BETA, 
                             A.data(), LDA, B.data(), LDB, 
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);
                             
    // Check results
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab04md returned error code: " << INFO;
    
    // Verify output matches expected discrete-time matrices
    ExpectMatricesNear(A, Ad, N, N, LDA, N, false, 1e-4);
    ExpectMatricesNear(B, Bd, N, M, LDB, N, false, 1e-4);
    ExpectMatricesNear(C, Cd, P, N, LDC, P, false, 1e-4);
    ExpectMatricesNear(D, Dd, P, M, LDD, P, false, 1e-4);
}

// Test discrete to continuous transformation
TEST_F(AB04MDTest, DiscreteToContiuous) {
    // Copy the discrete time system matrices
    std::vector<double> A = Ad;
    std::vector<double> B = Bd;
    std::vector<double> C = Cd;
    std::vector<double> D = Dd;

    int LDA = N;
    int LDB = N;
    int LDC = P;
    int LDD = P;
    int ROW_MAJOR = 0; // Using column-major (Fortran-style) storage
    
    // Call the SLICOT wrapper function to convert discrete to continuous
    int INFO = slicot_ab04md('D', N, M, P, ALPHA, BETA, 
                             A.data(), LDA, B.data(), LDB, 
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);
                             
    // Check results
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab04md returned error code: " << INFO;
    
    // Verify output matches expected continuous-time matrices
    ExpectMatricesNear(A, Ac, N, N, LDA, N, false, 1e-4);
    ExpectMatricesNear(B, Bc, N, M, LDB, N, false, 1e-4);
    ExpectMatricesNear(C, Cc, P, N, LDC, P, false, 1e-4);
    ExpectMatricesNear(D, Dc, P, M, LDD, P, false, 1e-4);
}

// Test continuous-discrete-continuous cycle
TEST_F(AB04MDTest, ContinuousDiscreteContinuous) {
    // Start with continuous-time system
    std::vector<double> A = Ac;
    std::vector<double> B = Bc;
    std::vector<double> C = Cc;
    std::vector<double> D = Dc;
    
    // Save the original matrices
    std::vector<double> A_orig = A;
    std::vector<double> B_orig = B;
    std::vector<double> C_orig = C;
    std::vector<double> D_orig = D;

    int LDA = N;
    int LDB = N;
    int LDC = P;
    int LDD = P;
    int ROW_MAJOR = 0; // Using column-major (Fortran-style) storage
    
    // Convert to discrete-time
    int INFO = slicot_ab04md('C', N, M, P, ALPHA, BETA, 
                             A.data(), LDA, B.data(), LDB, 
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "First conversion failed with code: " << INFO;
    
    // Convert back to continuous-time
    INFO = slicot_ab04md('D', N, M, P, ALPHA, BETA, 
                         A.data(), LDA, B.data(), LDB, 
                         C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "Second conversion failed with code: " << INFO;
    
    // Verify we get the original matrices back
    ExpectMatricesNear(A, A_orig, N, N, LDA, LDA, false, 1e-4);
    ExpectMatricesNear(B, B_orig, N, M, LDB, LDB, false, 1e-4);
    ExpectMatricesNear(C, C_orig, P, N, LDC, LDC, false, 1e-4);
    ExpectMatricesNear(D, D_orig, P, M, LDD, LDD, false, 1e-4);
}

// Test discrete-continuous-discrete cycle
TEST_F(AB04MDTest, DiscreteContiuousDiscrete) {
    // Start with discrete-time system
    std::vector<double> A = Ad;
    std::vector<double> B = Bd;
    std::vector<double> C = Cd;
    std::vector<double> D = Dd;
    
    // Save the original matrices
    std::vector<double> A_orig = A;
    std::vector<double> B_orig = B;
    std::vector<double> C_orig = C;
    std::vector<double> D_orig = D;

    int LDA = N;
    int LDB = N;
    int LDC = P;
    int LDD = P;
    int ROW_MAJOR = 0; // Using column-major (Fortran-style) storage
    
    // Convert to continuous-time
    int INFO = slicot_ab04md('D', N, M, P, ALPHA, BETA, 
                             A.data(), LDA, B.data(), LDB, 
                             C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "First conversion failed with code: " << INFO;
    
    // Convert back to discrete-time
    INFO = slicot_ab04md('C', N, M, P, ALPHA, BETA, 
                         A.data(), LDA, B.data(), LDB, 
                         C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "Second conversion failed with code: " << INFO;
    
    // Verify we get the original matrices back
    ExpectMatricesNear(A, A_orig, N, N, LDA, LDA, false, 1e-4);
    ExpectMatricesNear(B, B_orig, N, M, LDB, LDB, false, 1e-4);
    ExpectMatricesNear(C, C_orig, P, N, LDC, LDC, false, 1e-4);
    ExpectMatricesNear(D, D_orig, P, M, LDD, LDD, false, 1e-4);
}

// Test row-major storage option
TEST_F(AB04MDTest, RowMajorStorage) {
    // Convert the test data to row-major format
    std::vector<double> A_rm(N * N);
    std::vector<double> B_rm(N * M);
    std::vector<double> C_rm(P * N);
    std::vector<double> D_rm(P * M);
    
    // Copy from column-major to row-major
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A_rm[i*N + j] = Ac[j*N + i]; // Transpose
        }
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            B_rm[i*M + j] = Bc[j*N + i]; // Transpose
        }
    }
    for (int i = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j) {
            C_rm[i*N + j] = Cc[j*P + i]; // Transpose
        }
    }
    for (int i = 0; i < P; ++i) {
        for (int j = 0; j < M; ++j) {
            D_rm[i*M + j] = Dc[j*P + i]; // Transpose
        }
    }
    
    // Save copies for comparison after conversion
    std::vector<double> A_rm_orig = A_rm;
    std::vector<double> B_rm_orig = B_rm;
    std::vector<double> C_rm_orig = C_rm;
    std::vector<double> D_rm_orig = D_rm;
    
    // Row-major leading dimensions
    int LDA = N; // Number of columns in each row
    int LDB = M;
    int LDC = N;
    int LDD = M;
    int ROW_MAJOR = 1; // Using row-major storage
    
    // Convert to discrete-time and back
    int INFO = slicot_ab04md('C', N, M, P, ALPHA, BETA, 
                             A_rm.data(), LDA, B_rm.data(), LDB, 
                             C_rm.data(), LDC, D_rm.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "First conversion failed with code: " << INFO;
    
    // Convert back to continuous-time
    INFO = slicot_ab04md('D', N, M, P, ALPHA, BETA, 
                        A_rm.data(), LDA, B_rm.data(), LDB, 
                        C_rm.data(), LDC, D_rm.data(), LDD, ROW_MAJOR);
    ASSERT_EQ(INFO, 0) << "Second conversion failed with code: " << INFO;
    
    // Verify we get the original matrices back
    ExpectMatricesNear(A_rm, A_rm_orig, N, N, LDA, LDA, true, 1e-4);
    ExpectMatricesNear(B_rm, B_rm_orig, N, M, LDB, LDB, true, 1e-4);
    ExpectMatricesNear(C_rm, C_rm_orig, P, N, LDC, LDC, true, 1e-4);
    ExpectMatricesNear(D_rm, D_rm_orig, P, M, LDD, LDD, true, 1e-4);
}
