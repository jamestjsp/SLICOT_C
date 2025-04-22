#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits
#include <iomanip> // For std::setw in printMatrixD if needed

// Include the header for the function being tested
#include "ab05od.h" // Header for slicot_ab05od

// Define a fixture for AB05OD tests in column-major format (original Fortran ordering)
class AB05ODTestColMajor : public ::testing::Test {
protected:
    // Define common variables and expected values from documentation
    int N1 = 3;      // Order of first system
    int M1 = 2;      // Inputs to first system
    int P1 = 2;      // Outputs from first system
    int N2 = 3;      // Order of second system
    int M2 = 2;      // Inputs to second system
    double ALPHA = 1.0; // Scalar multiplier for second system
    int ROW_MAJOR = 0;  // Column-major
    double check_tol = 1e-4; // Tolerance for result verification
    
    // Expected dimensions for resulting system
    int N_expected = 6; // N1 + N2
    int M_expected = 4; // M1 + M2
    int expected_info = 0;
    
    // Input matrices in column-major order (as they appear in SLICOT docs)
    // Note: In the documentation, matrices are presented row by row for readability
    // but stored column by column in Fortran
    
    // First system matrices (A1, B1, C1, D1)
    std::vector<double> A1_in = {
        // Column 1
        1.0, 0.0, 1.0,
        // Column 2
        0.0, -1.0, 1.0,
        // Column 3
        -1.0, 1.0, 2.0
    };
    std::vector<double> B1_in = {
        // Column 1
        1.0, 1.0, 0.0,
        // Column 2
        2.0, 0.0, 1.0
    };
    std::vector<double> C1_in = {
        // Column 1
        3.0, 0.0,
        // Column 2
        -2.0, 1.0,
        // Column 3
        1.0, 0.0
    };
    std::vector<double> D1_in = {
        // Column 1
        1.0, 0.0,
        // Column 2
        0.0, 1.0
    };
    
    // Second system matrices (A2, B2, C2, D2)
    std::vector<double> A2_in = {
        // Column 1
        -3.0, 1.0, 0.0,
        // Column 2
        0.0, 0.0, -1.0,
        // Column 3
        0.0, 1.0, 2.0
    };
    std::vector<double> B2_in = {
        // Column 1
        0.0, 1.0, 0.0,
        // Column 2
        -1.0, 0.0, 2.0
    };
    std::vector<double> C2_in = {
        // Column 1
        1.0, 1.0,
        // Column 2
        1.0, 1.0,
        // Column 3
        0.0, -1.0
    };
    std::vector<double> D2_in = {
        // Column 1
        1.0, 0.0,
        // Column 2
        1.0, 1.0
    };
    
    // Expected output matrices from documentation in column-major order
    std::vector<double> A_expected = {
        // Column 1
        1.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000,
        // Column 2
        0.0000, -1.0000, 1.0000, 0.0000, 0.0000, 0.0000,
        // Column 3
        -1.0000, 1.0000, 2.0000, 0.0000, 0.0000, 0.0000,
        // Column 4
        0.0000, 0.0000, 0.0000, -3.0000, 1.0000, 0.0000,
        // Column 5
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, -1.0000,
        // Column 6
        0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 2.0000
    };
    std::vector<double> B_expected = {
        // Column 1
        1.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000,
        // Column 2
        2.0000, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000,
        // Column 3
        0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.0000,
        // Column 4
        0.0000, 0.0000, 0.0000, -1.0000, 0.0000, 2.0000
    };
    std::vector<double> C_expected = {
        // Column 1
        3.0000, 0.0000,
        // Column 2
        -2.0000, 1.0000,
        // Column 3
        1.0000, 0.0000,
        // Column 4
        1.0000, 1.0000,
        // Column 5
        1.0000, 1.0000,
        // Column 6
        0.0000, -1.0000
    };
    std::vector<double> D_expected = {
        // Column 1
        1.0000, 0.0000,
        // Column 2
        0.0000, 1.0000,
        // Column 3
        1.0000, 0.0000,
        // Column 4
        1.0000, 1.0000
    };
};

// Row-major test fixture (inherit from column-major and convert data)
class AB05ODTestRowMajor : public AB05ODTestColMajor {
public:
    AB05ODTestRowMajor() {
        ROW_MAJOR = 1;  // Override for row-major tests
        
        // Convert input matrices from column-major to row-major
        A1_rm.resize(N1 * N1);
        B1_rm.resize(N1 * M1);
        C1_rm.resize(P1 * N1);
        D1_rm.resize(P1 * M1);
        A2_rm.resize(N2 * N2);
        B2_rm.resize(N2 * M2);
        C2_rm.resize(P1 * N2);
        D2_rm.resize(P1 * M2);
        
        // Convert expected output matrices to row-major
        A_expected_rm.resize(N_expected * N_expected);
        B_expected_rm.resize(N_expected * M_expected);
        C_expected_rm.resize(P1 * N_expected);
        D_expected_rm.resize(P1 * M_expected);
        
        // Helper function to transpose data (column-major to row-major)
        auto transpose = [](const std::vector<double>& colMajor, std::vector<double>& rowMajor, int rows, int cols) {
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    rowMajor[i * cols + j] = colMajor[i + j * rows];
                }
            }
        };
        
        // Convert all input matrices
        transpose(A1_in, A1_rm, N1, N1);
        transpose(B1_in, B1_rm, N1, M1);
        transpose(C1_in, C1_rm, P1, N1);
        transpose(D1_in, D1_rm, P1, M1);
        transpose(A2_in, A2_rm, N2, N2);
        transpose(B2_in, B2_rm, N2, M2);
        transpose(C2_in, C2_rm, P1, N2);
        transpose(D2_in, D2_rm, P1, M2);
        
        // Convert expected outputs
        transpose(A_expected, A_expected_rm, N_expected, N_expected);
        transpose(B_expected, B_expected_rm, N_expected, M_expected);
        transpose(C_expected, C_expected_rm, P1, N_expected);
        transpose(D_expected, D_expected_rm, P1, M_expected);
    }
    
    std::vector<double> A1_rm, B1_rm, C1_rm, D1_rm;
    std::vector<double> A2_rm, B2_rm, C2_rm, D2_rm;
    std::vector<double> A_expected_rm, B_expected_rm, C_expected_rm, D_expected_rm;
};


// Test case for column-major format (Fortran ordering)
TEST_F(AB05ODTestColMajor, DocExample) {
    // Define leading dimensions for column-major format
    int LDA1 = N1;
    int LDB1 = N1;
    int LDC1 = P1;
    int LDD1 = P1;
    int LDA2 = N2;
    int LDB2 = N2;
    int LDC2 = P1;
    int LDD2 = P1;
    
    // Output matrices leading dimensions
    int LDA = N_expected;
    int LDB = N_expected;
    int LDC = P1;
    int LDD = P1;
    
    // Copy input matrices (they may be modified by the function)
    std::vector<double> A1 = A1_in;
    std::vector<double> B1 = B1_in;
    std::vector<double> C1 = C1_in;
    std::vector<double> D1 = D1_in;
    std::vector<double> A2 = A2_in;
    std::vector<double> B2 = B2_in;
    std::vector<double> C2 = C2_in;
    std::vector<double> D2 = D2_in;
    
    // Output matrices
    int N_out, M_out;
    std::vector<double> A(N_expected * N_expected);
    std::vector<double> B(N_expected * M_expected);
    std::vector<double> C(P1 * N_expected);
    std::vector<double> D(P1 * M_expected);
    
    // Call the C wrapper
    char OVER = 'N';  // No overlap mode
    int info = slicot_ab05od(OVER, 
                             N1, M1, P1, N2, M2, ALPHA,
                             A1.data(), LDA1, B1.data(), LDB1,
                             C1.data(), LDC1, D1.data(), LDD1,
                             A2.data(), LDA2, B2.data(), LDB2,
                             C2.data(), LDC2, D2.data(), LDD2,
                             &N_out, &M_out,
                             A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD,
                             ROW_MAJOR);
    
    // Verify return code and output dimensions
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(N_out, N_expected);
    ASSERT_EQ(M_out, M_expected);
    
    // Verify result matrices (column-major format)
    for (int j = 0; j < N_expected; ++j) {  // Column
        for (int i = 0; i < N_expected; ++i) {  // Row
            EXPECT_NEAR(A[i + j*LDA], A_expected[i + j*N_expected], check_tol)
                << "A[" << i << "," << j << "] mismatch (Column-Major)";
        }
    }
    
    for (int j = 0; j < M_expected; ++j) {  // Column
        for (int i = 0; i < N_expected; ++i) {  // Row
            EXPECT_NEAR(B[i + j*LDB], B_expected[i + j*N_expected], check_tol)
                << "B[" << i << "," << j << "] mismatch (Column-Major)";
        }
    }
    
    for (int j = 0; j < N_expected; ++j) {  // Column
        for (int i = 0; i < P1; ++i) {  // Row
            EXPECT_NEAR(C[i + j*LDC], C_expected[i + j*P1], check_tol)
                << "C[" << i << "," << j << "] mismatch (Column-Major)";
        }
    }
    
    for (int j = 0; j < M_expected; ++j) {  // Column
        for (int i = 0; i < P1; ++i) {  // Row
            EXPECT_NEAR(D[i + j*LDD], D_expected[i + j*P1], check_tol)
                << "D[" << i << "," << j << "] mismatch (Column-Major)";
        }
    }
}

// Test case for row-major format (C ordering)
TEST_F(AB05ODTestRowMajor, DocExample) {
    // Define leading dimensions for row-major format
    int LDA1 = N1;
    int LDB1 = M1;
    int LDC1 = N1;
    int LDD1 = M1;
    int LDA2 = N2;
    int LDB2 = M2;
    int LDC2 = N2;
    int LDD2 = M2;
    
    // Output matrices leading dimensions
    int LDA = N_expected;
    int LDB = M_expected;
    int LDC = N_expected;
    int LDD = M_expected;
    
    // Copy input matrices
    std::vector<double> A1 = A1_rm;
    std::vector<double> B1 = B1_rm;
    std::vector<double> C1 = C1_rm;
    std::vector<double> D1 = D1_rm;
    std::vector<double> A2 = A2_rm;
    std::vector<double> B2 = B2_rm;
    std::vector<double> C2 = C2_rm;
    std::vector<double> D2 = D2_rm;
    
    // Output matrices
    int N_out, M_out;
    std::vector<double> A(N_expected * LDA);
    std::vector<double> B(N_expected * LDB);
    std::vector<double> C(P1 * LDC);
    std::vector<double> D(P1 * LDD);
    
    // Call the C wrapper
    char OVER = 'N';  // No overlap mode
    int info = slicot_ab05od(OVER, 
                             N1, M1, P1, N2, M2, ALPHA,
                             A1.data(), LDA1, B1.data(), LDB1,
                             C1.data(), LDC1, D1.data(), LDD1,
                             A2.data(), LDA2, B2.data(), LDB2,
                             C2.data(), LDC2, D2.data(), LDD2,
                             &N_out, &M_out,
                             A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD,
                             ROW_MAJOR);
    
    // Verify return code and output dimensions
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(N_out, N_expected);
    ASSERT_EQ(M_out, M_expected);
    
    // Verify result matrices (row-major format)
    for (int i = 0; i < N_expected; ++i) {  // Row
        for (int j = 0; j < N_expected; ++j) {  // Column
            EXPECT_NEAR(A[i*LDA + j], A_expected_rm[i*N_expected + j], check_tol)
                << "A[" << i << "," << j << "] mismatch (Row-Major)";
        }
    }
    
    for (int i = 0; i < N_expected; ++i) {  // Row
        for (int j = 0; j < M_expected; ++j) {  // Column
            EXPECT_NEAR(B[i*LDB + j], B_expected_rm[i*M_expected + j], check_tol)
                << "B[" << i << "," << j << "] mismatch (Row-Major)";
        }
    }
    
    for (int i = 0; i < P1; ++i) {  // Row
        for (int j = 0; j < N_expected; ++j) {  // Column
            EXPECT_NEAR(C[i*LDC + j], C_expected_rm[i*N_expected + j], check_tol)
                << "C[" << i << "," << j << "] mismatch (Row-Major)";
        }
    }
    
    for (int i = 0; i < P1; ++i) {  // Row
        for (int j = 0; j < M_expected; ++j) {  // Column
            EXPECT_NEAR(D[i*LDD + j], D_expected_rm[i*M_expected + j], check_tol)
                << "D[" << i << "," << j << "] mismatch (Row-Major)";
        }
    }
}

// Test input validation
TEST_F(AB05ODTestColMajor, InputValidation) {
    // Define valid leading dimensions
    int LDA1 = MAX(1, N1);
    int LDB1 = MAX(1, N1);
    int LDC1 = MAX(1, P1);
    int LDD1 = MAX(1, P1);
    int LDA2 = MAX(1, N2);
    int LDB2 = MAX(1, N2);
    int LDC2 = MAX(1, P1);
    int LDD2 = MAX(1, P1);
    int LDA = MAX(1, N1 + N2);
    int LDB = MAX(1, N1 + N2);
    int LDC = MAX(1, P1);
    int LDD = MAX(1, P1);
    
    // Output matrices
    int N_out, M_out;
    std::vector<double> A((N1 + N2) * LDA);
    std::vector<double> B((N1 + N2) * LDB);
    std::vector<double> C(P1 * LDC);
    std::vector<double> D(P1 * LDD);
    
    // Test negative dimensions
    int info = slicot_ab05od('N', -1, M1, P1, N2, M2, ALPHA,
                           A1_in.data(), LDA1, B1_in.data(), LDB1,
                           C1_in.data(), LDC1, D1_in.data(), LDD1,
                           A2_in.data(), LDA2, B2_in.data(), LDB2,
                           C2_in.data(), LDC2, D2_in.data(), LDD2,
                           &N_out, &M_out, A.data(), LDA, B.data(), LDB,
                           C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    EXPECT_NE(info, 0) << "Expected error for negative N1";
    
    // Test invalid OVER value (should be handled gracefully in the wrapper)
    info = slicot_ab05od('X', N1, M1, P1, N2, M2, ALPHA,
                       A1_in.data(), LDA1, B1_in.data(), LDB1,
                       C1_in.data(), LDC1, D1_in.data(), LDD1,
                       A2_in.data(), LDA2, B2_in.data(), LDB2,
                       C2_in.data(), LDC2, D2_in.data(), LDD2,
                       &N_out, &M_out, A.data(), LDA, B.data(), LDB,
                       C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    EXPECT_NE(info, 0) << "Expected error for invalid OVER value";
    
    // Test insufficient leading dimensions
    info = slicot_ab05od('N', N1, M1, P1, N2, M2, ALPHA,
                       A1_in.data(), 0, B1_in.data(), LDB1,
                       C1_in.data(), LDC1, D1_in.data(), LDD1,
                       A2_in.data(), LDA2, B2_in.data(), LDB2,
                       C2_in.data(), LDC2, D2_in.data(), LDD2,
                       &N_out, &M_out, A.data(), LDA, B.data(), LDB,
                       C.data(), LDC, D.data(), LDD, ROW_MAJOR);
    EXPECT_NE(info, 0) << "Expected error for insufficient LDA1";
}

// Test 'O' option which gets converted to 'N' in the wrapper
TEST_F(AB05ODTestColMajor, OverlapOptionHandling) {
    // Define leading dimensions
    int LDA1 = N1;
    int LDB1 = N1;
    int LDC1 = P1;
    int LDD1 = P1;
    int LDA2 = N2;
    int LDB2 = N2;
    int LDC2 = P1;
    int LDD2 = P1;
    
    // Output matrices leading dimensions
    int LDA = N_expected;
    int LDB = N_expected;
    int LDC = P1;
    int LDD = P1;
    
    // Copy input matrices
    std::vector<double> A1 = A1_in;
    std::vector<double> B1 = B1_in;
    std::vector<double> C1 = C1_in;
    std::vector<double> D1 = D1_in;
    std::vector<double> A2 = A2_in;
    std::vector<double> B2 = B2_in;
    std::vector<double> C2 = C2_in;
    std::vector<double> D2 = D2_in;
    
    // Output matrices
    int N_out, M_out;
    std::vector<double> A(N_expected * N_expected);
    std::vector<double> B(N_expected * M_expected);
    std::vector<double> C(P1 * N_expected);
    std::vector<double> D(P1 * M_expected);
    
    // Call the wrapper with OVER='O' - should be handled as 'N'
    char OVER = 'O';
    int info = slicot_ab05od(OVER, 
                             N1, M1, P1, N2, M2, ALPHA,
                             A1.data(), LDA1, B1.data(), LDB1,
                             C1.data(), LDC1, D1.data(), LDD1,
                             A2.data(), LDA2, B2.data(), LDB2,
                             C2.data(), LDC2, D2.data(), LDD2,
                             &N_out, &M_out,
                             A.data(), LDA, B.data(), LDB,
                             C.data(), LDC, D.data(), LDD,
                             ROW_MAJOR);
    
    // Should still produce correct results
    ASSERT_EQ(info, expected_info);
    ASSERT_EQ(N_out, N_expected);
    ASSERT_EQ(M_out, M_expected);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
