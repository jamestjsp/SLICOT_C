#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <cstring>   // For memcpy
#include <iomanip>   // For std::fixed, std::setprecision
#include <iostream>  // For debug prints (optional)

#include "sb04qd.h"
#include "slicot_utils.h" // For transpose functions

// --- Test Fixture ---
class SB04QDTest : public ::testing::Test {
protected:
    int N_in, M_in;

    std::vector<double> A_io;
    std::vector<double> B_io;
    std::vector<double> C_io; // Input C, Output X
    std::vector<double> Z_out;

    int LDA_in, LDB_in, LDC_in, LDZ_in;

    // Expected results from SB04QD.html example (stored column-major)
    std::vector<double> X_expected_cm;
    std::vector<double> Z_expected_cm;
    // A and B are also modified by the routine as per documentation:
    // A -> upper Hessenberg H + U factorization data
    // B -> Schur factor S of B'
    // We will primarily check X and Z for simplicity, as A and B are intermediate/transformed.

    int INFO_expected = 0;
    double check_tol = 1e-4; // Based on example precision (F8.4)

    void InitializeDimensions(int n_val, int m_val, bool row_major_layout) {
        N_in = n_val;
        M_in = m_val;

        if (row_major_layout) {
            LDA_in = std::max(1, N_in); // A is N x N, RM LDA is cols of A.
            LDB_in = std::max(1, M_in); // B is M x M, RM LDB is cols of B.
            LDC_in = std::max(1, M_in); // C (X) is N x M, RM LDC is cols of C.
            LDZ_in = std::max(1, M_in); // Z is M x M, RM LDZ is cols of Z.

            A_io.resize((size_t)N_in * LDA_in); 
            B_io.resize((size_t)M_in * LDB_in); 
            C_io.resize((size_t)N_in * LDC_in); 
            Z_out.resize((size_t)M_in * LDZ_in);
        } else { // Column-major
            LDA_in = std::max(1, N_in); // CM LDA is rows of A.
            LDB_in = std::max(1, M_in); // CM LDB is rows of B.
            LDC_in = std::max(1, N_in); // CM LDC is rows of C.
            LDZ_in = std::max(1, M_in); // CM LDZ is rows of Z.

            A_io.resize((size_t)LDA_in * N_in); 
            B_io.resize((size_t)LDB_in * M_in); 
            C_io.resize((size_t)LDC_in * M_in);
            Z_out.resize((size_t)LDZ_in * M_in);
        }
    }

    void SetUpExampleData(bool row_major_layout) {
        N_in = 3; M_in = 3;
        InitializeDimensions(N_in, M_in, row_major_layout);

        // Example data (column-major format from HTML)
        std::vector<double> a_example_cm = {
            1.0, 6.0, 9.0, // col 1
            2.0, 7.0, 2.0, // col 2
            3.0, 8.0, 3.0  // col 3
        };
        std::vector<double> b_example_cm = {
            7.0, 2.0, 3.0, // col 1
            2.0, 1.0, 4.0, // col 2
            3.0, 2.0, 1.0  // col 3
        };
        std::vector<double> c_example_cm = { // This is input C
            271.0, 923.0, 578.0, // col 1
            135.0, 494.0, 383.0, // col 2
            147.0, 482.0, 287.0  // col 3
        };

        X_expected_cm = { // This is C on output (solution X)
             2.0, 4.0, 5.0, // col 1 of X
             3.0, 7.0, 3.0, // col 2 of X
             6.0, 1.0, 2.0  // col 3 of X
        };
        Z_expected_cm = {
             0.8337,  0.3881,  0.3928, // col 1 of Z
             0.5204, -0.7900, -0.3241, // col 2 of Z
            -0.1845, -0.4746,  0.8606  // col 3 of Z
        };
        INFO_expected = 0;

        if (row_major_layout) {
            // Transpose A, B, C from column-major example data to row-major A_io, B_io, C_io
            slicot_transpose_to_c_with_ld(a_example_cm.data(), A_io.data(), N_in, N_in, N_in, LDA_in, sizeof(double));
            slicot_transpose_to_c_with_ld(b_example_cm.data(), B_io.data(), M_in, M_in, M_in, LDB_in, sizeof(double));
            slicot_transpose_to_c_with_ld(c_example_cm.data(), C_io.data(), N_in, M_in, N_in, LDC_in, sizeof(double));
        } else { // Column-major
            memcpy(A_io.data(), a_example_cm.data(), a_example_cm.size() * sizeof(double));
            memcpy(B_io.data(), b_example_cm.data(), b_example_cm.size() * sizeof(double));
            memcpy(C_io.data(), c_example_cm.data(), c_example_cm.size() * sizeof(double));
        }
    }

    void CheckSolution(const std::vector<double>& computed_x_in_c, const std::vector<double>& computed_z, bool is_row_major_computed) {
        // Check X (solution in C_io)
        if (N_in > 0 && M_in > 0) {
            if (is_row_major_computed) {
                std::vector<double> x_expected_rm(N_in * M_in);
                slicot_transpose_to_c_with_ld(X_expected_cm.data(), x_expected_rm.data(), N_in, M_in, N_in, M_in, sizeof(double));
                for (int i = 0; i < N_in; ++i) {
                    for (int j = 0; j < M_in; ++j) {
                        EXPECT_NEAR(computed_x_in_c[i * LDC_in + j], x_expected_rm[i * M_in + j], check_tol)
                            << "X_rm(" << i << "," << j << ") mismatch";
                    }
                }
            } else { // Column-major computed
                for (int j = 0; j < M_in; ++j) { // col
                    for (int i = 0; i < N_in; ++i) { // row
                        EXPECT_NEAR(computed_x_in_c[j * LDC_in + i], X_expected_cm[j * N_in + i], check_tol)
                            << "X_cm(" << i << "," << j << ") mismatch";
                    }
                }
            }
        }

        // Check Z
        if (M_in > 0) {
            if (is_row_major_computed) {
                std::vector<double> z_expected_rm(M_in * M_in);
                slicot_transpose_to_c_with_ld(Z_expected_cm.data(), z_expected_rm.data(), M_in, M_in, M_in, M_in, sizeof(double));
                for (int i = 0; i < M_in; ++i) {
                    for (int j = 0; j < M_in; ++j) {
                        EXPECT_NEAR(computed_z[i * LDZ_in + j], z_expected_rm[i * M_in + j], check_tol)
                            << "Z_rm(" << i << "," << j << ") mismatch";
                    }
                }
            } else { // Column-major computed
                for (int j = 0; j < M_in; ++j) { // col
                    for (int i = 0; i < M_in; ++i) { // row
                        EXPECT_NEAR(computed_z[j * LDZ_in + i], Z_expected_cm[j * M_in + i], check_tol)
                            << "Z_cm(" << i << "," << j << ") mismatch";
                    }
                }
            }
        }
    }
};

// Test: Documentation Example (Column-Major)
TEST_F(SB04QDTest, DocExample_ColMajor) {
    SetUpExampleData(false); // false for column-major

    int info = slicot_sb04qd(N_in, M_in,
                             A_io.data(), LDA_in, B_io.data(), LDB_in,
                             C_io.data(), LDC_in, Z_out.data(), LDZ_in,
                             0 /* col-major */);

    ASSERT_EQ(info, INFO_expected);
    CheckSolution(C_io, Z_out, false);
}

// Test: Documentation Example (Row-Major)
TEST_F(SB04QDTest, DocExample_RowMajor) {
    SetUpExampleData(true); // true for row-major

    int info = slicot_sb04qd(N_in, M_in,
                             A_io.data(), LDA_in, B_io.data(), LDB_in,
                             C_io.data(), LDC_in, Z_out.data(), LDZ_in,
                             1 /* row-major */);
    
    ASSERT_EQ(info, INFO_expected);
    CheckSolution(C_io, Z_out, true);
}

// Test: Zero Dimensions (N=0)
TEST_F(SB04QDTest, ZeroDimensions_N0) {
    N_in = 0; M_in = 2;
    InitializeDimensions(N_in, M_in, false); // col-major

    // A, C are N x N and N x M, so effectively empty or not referenced if N=0.
    // B, Z are M x M.
    std::vector<double> b_dummy_cm = {1,0,0,1}; // M x M identity for B
    if (M_in > 0) memcpy(B_io.data(), b_dummy_cm.data(), b_dummy_cm.size()*sizeof(double));

    int info = slicot_sb04qd(N_in, M_in,
                             nullptr, LDA_in, B_io.data(), LDB_in,
                             nullptr, LDC_in, Z_out.data(), LDZ_in,
                             0 /* col-major */);
    EXPECT_EQ(info, 0);
    // If N=0, X is 0xM; C_io (output) is not well-defined for checking specific values.
    // Z should be computed if M > 0. B is modified.
    // Primary check is that INFO=0 for valid zero-dimension inputs.
}

// Test: Zero Dimensions (M=0)
TEST_F(SB04QDTest, ZeroDimensions_M0) {
    N_in = 2; M_in = 0;
    InitializeDimensions(N_in, M_in, false); // col-major

    // B, C, Z are M x M or N x M, so effectively empty or not referenced if M=0.
    // A is N x N.
    std::vector<double> a_dummy_cm = {1,0,0,1}; // N x N identity for A
    if (N_in > 0) memcpy(A_io.data(), a_dummy_cm.data(), a_dummy_cm.size()*sizeof(double));

    int info = slicot_sb04qd(N_in, M_in,
                             A_io.data(), LDA_in, nullptr, LDB_in,
                             nullptr, LDC_in, nullptr, LDZ_in,
                             0 /* col-major */);
    EXPECT_EQ(info, 0);
    // If M=0, X is Nx0, Z is 0x0. C_io, Z_out (output) are not well-defined for checking.
}

// Test: Zero Dimensions (N=0, M=0)
TEST_F(SB04QDTest, ZeroDimensions_N0_M0) {
    N_in = 0; M_in = 0;
    InitializeDimensions(N_in, M_in, false); // col-major

    int info = slicot_sb04qd(N_in, M_in,
                             nullptr, LDA_in, nullptr, LDB_in,
                             nullptr, LDC_in, nullptr, LDZ_in,
                             0 /* col-major */);
    EXPECT_EQ(info, 0);
}

// Test: Parameter Validation (Selected)
TEST_F(SB04QDTest, ParameterValidation) {
    SetUpExampleData(false); // Initialize with valid N, M etc. for non-NULL pointers

    // Test invalid N
    int info = slicot_sb04qd(-1, M_in, A_io.data(), LDA_in, B_io.data(), LDB_in, C_io.data(), LDC_in, Z_out.data(), LDZ_in, 0);
    EXPECT_EQ(info, -1);

    // Test invalid M
    info = slicot_sb04qd(N_in, -1, A_io.data(), LDA_in, B_io.data(), LDB_in, C_io.data(), LDC_in, Z_out.data(), LDZ_in, 0);
    EXPECT_EQ(info, -2);
    
    // Test invalid LDA (assuming N_in > 0 from SetUpExampleData)
    if (N_in > 0) {
        info = slicot_sb04qd(N_in, M_in, A_io.data(), 0, B_io.data(), LDB_in, C_io.data(), LDC_in, Z_out.data(), LDZ_in, 0);
        EXPECT_EQ(info, -4); // LDA is arg 4
    }
    
    // Test invalid LDB (assuming M_in > 0 from SetUpExampleData)
    if (M_in > 0) {
        info = slicot_sb04qd(N_in, M_in, A_io.data(), LDA_in, B_io.data(), 0, C_io.data(), LDC_in, Z_out.data(), LDZ_in, 0);
        EXPECT_EQ(info, -6); // LDB is arg 6
    }
    
    // Test NULL A when N > 0
    if (N_in > 0) {
        info = slicot_sb04qd(N_in, M_in, nullptr, LDA_in, B_io.data(), LDB_in, C_io.data(), LDC_in, Z_out.data(), LDZ_in, 0);
        // The wrapper should catch this (or the Fortran routine will signal for arg 3).
        // The wrapper assigns error -3 based on Fortran argument position for A pointer.
        EXPECT_EQ(info, -3); 
    }
}

