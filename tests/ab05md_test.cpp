#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <limits> // Required for std::numeric_limits

// Include the specific header for the function being tested
#include "ab05md.h" // Header for slicot_ab05md

// Define a fixture for AB05MD tests
class AB05MDTest : public ::testing::Test {
    // No specific setup needed for these tests
};

// =============================================================================
// Test Case: Documentation Example (UPLO='L', OVER='N') - ROW-MAJOR
// =============================================================================
TEST_F(AB05MDTest, DocExample_RowMajor) {
    // --- Input Parameters ---
    char UPLO = 'L'; // Lower block diagonal form
    char OVER = 'N'; // No array overlapping
    int N1 = 3;      // Order of system 1
    int M1 = 2;      // Inputs to system 1
    int P1 = 2;      // Outputs from system 1 / Inputs to system 2
    int N2 = 3;      // Order of system 2
    int P2 = 2;      // Outputs from system 2
    int ROW_MAJOR = 1; // Test Row-Major storage
    double check_tol = 1e-4; // Tolerance for checking results

    // --- Leading Dimensions (Row-Major: Number of Columns) ---
    int LDA1 = N1;
    int LDB1 = M1;
    int LDC1 = N1;
    int LDD1 = M1;
    int LDA2 = N2;
    int LDB2 = P1;
    int LDC2 = N2;
    int LDD2 = P1;
    int N_expected = N1 + N2; // Expected resulting system order (6)
    int LDA = N_expected;     // LD for resulting A
    int LDB = M1;             // LD for resulting B
    int LDC = N_expected;     // LD for resulting C
    int LDD = M1;             // LD for resulting D

    // --- Input Matrices (Row-Major format, from AB05MD.html example DATA) ---
    std::vector<double> A1 = {
        1.0, 0.0, -1.0,  // Row 1
        0.0, -1.0, 1.0,  // Row 2
        1.0, 1.0, 2.0    // Row 3
    };
    std::vector<double> B1 = {
        1.0, 1.0,  // Row 1
        2.0, 0.0,  // Row 2
        1.0, 0.0   // Row 3
    };
    std::vector<double> C1 = {
        3.0, -2.0, 1.0,  // Row 1
        0.0, 1.0, 0.0    // Row 2
    };
    std::vector<double> D1 = {
        1.0, 0.0,  // Row 1
        0.0, 1.0   // Row 2
    };
    std::vector<double> A2 = {
        -3.0, 0.0, 0.0,  // Row 1
        1.0, 0.0, 1.0,   // Row 2
        0.0, -1.0, 2.0   // Row 3
    };
    std::vector<double> B2 = {
        0.0, -1.0,  // Row 1
        1.0, 0.0,   // Row 2
        0.0, 2.0    // Row 3
    };
    std::vector<double> C2 = {
        1.0, 1.0, 0.0,  // Row 1
        1.0, 1.0, -1.0  // Row 2
    };
    std::vector<double> D2 = {
        1.0, 1.0,  // Row 1
        0.0, 1.0   // Row 2
    };

    // --- Output Variables ---
    int N_result;
    std::vector<double> A(N_expected * LDA); // 6x6
    std::vector<double> B(N_expected * LDB); // 6x2
    std::vector<double> C(P2 * LDC);         // 2x6
    std::vector<double> D(P2 * LDD);         // 2x2
    std::vector<double> DWORK(1); // Dummy workspace for OVER='N'

    // --- Call the C wrapper function ---
    // Assuming the C wrapper slicot_ab05md handles the DWORK/LDWORK arguments internally
    // or they are not needed when OVER='N'. Modify if your wrapper requires them explicitly.
    int INFO = slicot_ab05md(UPLO, OVER,
                            N1, M1, P1, N2, P2,
                            A1.data(), LDA1, B1.data(), LDB1,
                            C1.data(), LDC1, D1.data(), LDD1,
                            A2.data(), LDA2, B2.data(), LDB2,
                            C2.data(), LDC2, D2.data(), LDD2,
                            &N_result,
                            A.data(), LDA, B.data(), LDB,
                            C.data(), LDC, D.data(), LDD,
                            ROW_MAJOR);

    // --- Check results ---
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab05md returned error code: " << INFO;
    ASSERT_EQ(N_result, N_expected) << "Output dimension N does not match expected value";

    // --- Expected Output Matrices (Row-Major format) ---
    // Using calculated values based on formula.
     std::vector<double> A_expected = {
        1.0000, 0.0000, -1.0000, 0.0000, 0.0000, 0.0000,  // Row 1
        0.0000, -1.0000, 1.0000, 0.0000, 0.0000, 0.0000,  // Row 2
        1.0000, 1.0000, 2.0000, 0.0000, 0.0000, 0.0000,   // Row 3
        0.0000, -1.0000, 0.0000, -3.0000, 0.0000, 0.0000, // Row 4 (B2*C1)
        3.0000, -2.0000, 1.0000, 1.0000, 0.0000, 1.0000,  // Row 5 (B2*C1 | A2)
        0.0000, 2.0000, 0.0000, 0.0000, -1.0000, 2.0000   // Row 6 (B2*C1 | A2)
    };
     std::vector<double> B_expected = {
        1.0000, 1.0000,  // Row 1 (B1)
        2.0000, 0.0000,  // Row 2 (B1)
        1.0000, 0.0000,  // Row 3 (B1)
        0.0000, -1.0000, // Row 4 (B2*D1) - Calculated
        1.0000, 0.0000,  // Row 5 (B2*D1) - Calculated
        0.0000, 2.0000   // Row 6 (B2*D1) - Calculated
    };
    std::vector<double> C_expected = {
        3.0000, -1.0000, 1.0000, 1.0000, 1.0000, 0.0000,  // Row 1 (D2*C1 | C2)
        0.0000, 1.0000, 0.0000, 1.0000, 1.0000, -1.0000   // Row 2 (D2*C1 | C2)
    };
    std::vector<double> D_expected = {
        1.0000, 1.0000,  // Row 1 (D2*D1)
        0.0000, 1.0000   // Row 2 (D2*D1)
    };


    // --- Verification (Row-Major Indexing) ---
    // Verify A
    for (int i = 0; i < N_expected; ++i) { // row
        for (int j = 0; j < N_expected; ++j) { // col
            EXPECT_NEAR(A[i * LDA + j], A_expected[i * N_expected + j], check_tol)
                << "A[" << i << "," << j << "] mismatch (Row-Major)";
        }
    }
    // Verify B
    for (int i = 0; i < N_expected; ++i) { // row
        for (int j = 0; j < M1; ++j) { // col
            EXPECT_NEAR(B[i * LDB + j], B_expected[i * M1 + j], check_tol)
                << "B[" << i << "," << j << "] mismatch (Row-Major)";
        }
    }
    // Verify C
    for (int i = 0; i < P2; ++i) { // row
        for (int j = 0; j < N_expected; ++j) { // col
            EXPECT_NEAR(C[i * LDC + j], C_expected[i * N_expected + j], check_tol)
                << "C[" << i << "," << j << "] mismatch (Row-Major)";
        }
    }
    // Verify D
    for (int i = 0; i < P2; ++i) { // row
        for (int j = 0; j < M1; ++j) { // col
            EXPECT_NEAR(D[i * LDD + j], D_expected[i * M1 + j], check_tol)
                << "D[" << i << "," << j << "] mismatch (Row-Major)";
        }
    }
}


// =============================================================================
// Test Case: Documentation Example (UPLO='L', OVER='N') - COLUMN-MAJOR
// =============================================================================
TEST_F(AB05MDTest, DocExample_ColMajor) {
    // --- Input Parameters ---
    char UPLO = 'L'; // Lower block diagonal form
    char OVER = 'N'; // No array overlapping
    int N1 = 3;      // Order of system 1
    int M1 = 2;      // Inputs to system 1
    int P1 = 2;      // Outputs from system 1 / Inputs to system 2
    int N2 = 3;      // Order of system 2
    int P2 = 2;      // Outputs from system 2
    int ROW_MAJOR = 0; // Test Column-Major storage
    double check_tol = 1e-4; // Tolerance for checking results

    // --- Leading Dimensions (Column-Major: Number of Rows) ---
    int LDA1 = N1;
    int LDB1 = N1;
    int LDC1 = P1; // >= MAX(1,P1)
    int LDD1 = P1; // >= MAX(1,P1)
    int LDA2 = N2;
    int LDB2 = N2;
    int LDC2 = P2; // >= MAX(1,P2)
    int LDD2 = P2; // >= MAX(1,P2)
    int N_expected = N1 + N2; // Expected resulting system order (6)
    int LDA = N_expected;     // LD for resulting A
    int LDB = N_expected;     // LD for resulting B
    int LDC = P2;             // LD for resulting C // >= MAX(1,P2)
    int LDD = P2;             // LD for resulting D // >= MAX(1,P2)

    // --- Input Matrices (Column-Major format, transposed from example DATA) ---
    std::vector<double> A1 = {
        1.0, 0.0, 1.0,   // Col 1
        0.0, -1.0, 1.0,  // Col 2
       -1.0, 1.0, 2.0    // Col 3
    };
    std::vector<double> B1 = {
        1.0, 2.0, 1.0,   // Col 1
        1.0, 0.0, 0.0    // Col 2
    };
    std::vector<double> C1 = {
        3.0, 0.0,   // Col 1
       -2.0, 1.0,   // Col 2
        1.0, 0.0    // Col 3
    };
    std::vector<double> D1 = {
        1.0, 0.0,   // Col 1
        0.0, 1.0    // Col 2
    };
    std::vector<double> A2 = {
       -3.0, 1.0, 0.0,   // Col 1
        0.0, 0.0, -1.0,  // Col 2
        0.0, 1.0, 2.0    // Col 3
    };
    std::vector<double> B2 = {
        0.0, 1.0, 0.0,   // Col 1
       -1.0, 0.0, 2.0    // Col 2
    };
    std::vector<double> C2 = {
        1.0, 1.0,   // Col 1
        1.0, 1.0,   // Col 2
        0.0, -1.0   // Col 3
    };
    std::vector<double> D2 = {
        1.0, 0.0,   // Col 1
        1.0, 1.0    // Col 2
    };

    // --- Output Variables ---
    int N_result;
    std::vector<double> A(LDA * N_expected); // 6x6
    std::vector<double> B(LDB * M1);         // 6x2
    std::vector<double> C(LDC * N_expected); // 2x6
    std::vector<double> D(LDD * M1);         // 2x2
    std::vector<double> DWORK(1); // Dummy workspace for OVER='N'


    // --- Call the C wrapper function ---
    // Assuming the C wrapper slicot_ab05md handles the DWORK/LDWORK arguments internally
    // or they are not needed when OVER='N'. Modify if your wrapper requires them explicitly.
    int INFO = slicot_ab05md(UPLO, OVER,
                            N1, M1, P1, N2, P2,
                            A1.data(), LDA1, B1.data(), LDB1,
                            C1.data(), LDC1, D1.data(), LDD1,
                            A2.data(), LDA2, B2.data(), LDB2,
                            C2.data(), LDC2, D2.data(), LDD2,
                            &N_result,
                            A.data(), LDA, B.data(), LDB,
                            C.data(), LDC, D.data(), LDD,
                            ROW_MAJOR);

    // --- Check results ---
    ASSERT_EQ(INFO, 0) << "SLICOT routine ab05md returned error code: " << INFO;
    ASSERT_EQ(N_result, N_expected) << "Output dimension N does not match expected value";

    // --- Expected Output Matrices (Column-Major format, calculated from formulas) ---
     std::vector<double> A_expected = {
        // A1          | B2*C1 block
        1.0,  0.0,  1.0,  0.0,  3.0,  0.0, // Col 1
        0.0, -1.0,  1.0, -1.0, -2.0,  2.0, // Col 2
       -1.0,  1.0,  2.0,  0.0,  1.0,  0.0, // Col 3
        // 0           | A2 block
        0.0,  0.0,  0.0, -3.0,  1.0,  0.0, // Col 4
        0.0,  0.0,  0.0,  0.0,  0.0, -1.0, // Col 5
        0.0,  0.0,  0.0,  0.0,  1.0,  2.0  // Col 6
    };
     std::vector<double> B_expected = {
        // B1 | B2*D1
        1.0,  2.0,  1.0,  0.0,  1.0,  0.0, // Col 1
        1.0,  0.0,  0.0, -1.0,  0.0,  2.0  // Col 2
    };
     std::vector<double> C_expected = {
        // D2*C1 | C2
        3.0,  0.0,   // Col 1
       -1.0,  1.0,   // Col 2
        1.0,  0.0,   // Col 3
        1.0,  1.0,   // Col 4
        1.0,  1.0,   // Col 5
        0.0, -1.0    // Col 6
    };
    std::vector<double> D_expected = {
        // D2*D1
        1.0, 0.0,   // Col 1
        1.0, 1.0    // Col 2
    };

    // --- Verification (Column-Major Indexing) ---
    // Verify A
    for (int j = 0; j < N_expected; ++j) { // col
        for (int i = 0; i < N_expected; ++i) { // row
            EXPECT_NEAR(A[i + j * LDA], A_expected[i + j * N_expected], check_tol)
                << "A[" << i << "," << j << "] mismatch (Col-Major)";
        }
    }
    // Verify B
    for (int j = 0; j < M1; ++j) { // col
        for (int i = 0; i < N_expected; ++i) { // row
            // Access B_expected using correct stride N_expected (number of rows)
            EXPECT_NEAR(B[i + j * LDB], B_expected[i + j * N_expected], check_tol)
                << "B[" << i << "," << j << "] mismatch (Col-Major)";
        }
    }
    // Verify C
    for (int j = 0; j < N_expected; ++j) { // col
        for (int i = 0; i < P2; ++i) { // row
            // Access C_expected using correct stride P2 (number of rows)
            EXPECT_NEAR(C[i + j * LDC], C_expected[i + j * P2], check_tol)
                << "C[" << i << "," << j << "] mismatch (Col-Major)";
        }
    }
    // Verify D
    for (int j = 0; j < M1; ++j) { // col
        for (int i = 0; i < P2; ++i) { // row
             // Access D_expected using correct stride P2 (number of rows)
            EXPECT_NEAR(D[i + j * LDD], D_expected[i + j * P2], check_tol)
                << "D[" << i << "," << j << "] mismatch (Col-Major)";
        }
    }
}

// TODO: Add tests for UPLO='U' and OVER='O' for Column-Major if needed,
//       being mindful of the potential discrepancies noted for Row-Major.

