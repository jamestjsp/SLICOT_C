#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <cstring>   // For memcpy
#include <iomanip>   // For std::fixed, std::setprecision
#include <iostream>  // For debug prints (optional)

#include "sb04md.h"
#include "slicot_utils.h" // For transpose functions

// --- Test Fixture ---
class SB04MDTest : public ::testing::Test {
protected:
    int N_in, M_in;

    std::vector<double> A_io;
    std::vector<double> B_io;
    std::vector<double> C_io; // Input C, Output X
    std::vector<double> Z_out;

    int LDA_in, LDB_in, LDC_in, LDZ_in;

    // Expected results from SB04MD.html example (stored column-major)
    std::vector<double> X_expected_cm;
    std::vector<double> Z_expected_cm;
    // A and B are also modified by the routine, but the example doesn't show their final state.
    // We will primarily check X and Z.

    int INFO_expected = 0;
    double check_tol = 1e-4; // Based on example precision

    void InitializeDimensions(int n_val, int m_val, bool row_major_layout) {
        N_in = n_val;
        M_in = m_val;

        if (row_major_layout) {
            LDA_in = std::max(1, N_in); // A is N x N, C is N x M. RM LDA is cols of A.
            LDB_in = std::max(1, M_in); // B is M x M, RM LDB is cols of B.
            LDC_in = std::max(1, M_in); // RM LDC is cols of C.
            LDZ_in = std::max(1, M_in); // Z is M x M, RM LDZ is cols of Z.

            A_io.resize((size_t)N_in * LDA_in); // N_in rows, LDA_in (N_in) cols
            B_io.resize((size_t)M_in * LDB_in); // M_in rows, LDB_in (M_in) cols
            C_io.resize((size_t)N_in * LDC_in); // N_in rows, LDC_in (M_in) cols
            Z_out.resize((size_t)M_in * LDZ_in); // M_in rows, LDZ_in (M_in) cols
        } else { // Column-major
            LDA_in = std::max(1, N_in); // CM LDA is rows of A.
            LDB_in = std::max(1, M_in); // CM LDB is rows of B.
            LDC_in = std::max(1, N_in); // CM LDC is rows of C.
            LDZ_in = std::max(1, M_in); // CM LDZ is rows of Z.

            A_io.resize((size_t)LDA_in * N_in); // LDA_in (N_in) rows, N_in cols
            B_io.resize((size_t)LDB_in * M_in); // LDB_in (M_in) rows, M_in cols
            C_io.resize((size_t)LDC_in * M_in); // LDC_in (N_in) rows, M_in cols
            Z_out.resize((size_t)LDZ_in * M_in); // LDZ_in (M_in) rows, M_in cols
        }
    }

    void SetUpExampleData(bool row_major_layout) {
        N_in = 3; M_in = 2;
        InitializeDimensions(N_in, M_in, row_major_layout);

        // Example data (column-major format from HTML interpretation)
        std::vector<double> a_example_cm = {
            2.0, 0.0, 6.0, // col 1
            1.0, 2.0, 1.0, // col 2
            3.0, 1.0, 2.0  // col 3
        };
        std::vector<double> b_example_cm = {
            2.0, 1.0, // col 1
            1.0, 6.0  // col 2
        };
        std::vector<double> c_example_cm = {
            2.0, 1.0, 0.0, // col 1
            1.0, 4.0, 5.0  // col 2
        };

        X_expected_cm = {
            -2.7685, -1.0531,  4.5257, // col 1 of X
             0.5498,  0.6865, -0.4389  // col 2 of X
        };
        Z_expected_cm = {
            -0.9732,  0.2298, // col 1 of Z
            -0.2298, -0.9732  // col 2 of Z
        };
        INFO_expected = 0;

        if (row_major_layout) {
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
        if (N_in == 0 || M_in == 0) { // If N or M is 0, X is N x M, so it's empty or not computed. Z is M x M.
            if (M_in == 0) return; // Z is also not computed if M=0
        }

        // Check X (solution in C_io)
        if (N_in > 0 && M_in > 0) {
            if (is_row_major_computed) {
                std::vector<double> x_expected_rm(N_in * M_in);
                slicot_transpose_to_c_with_ld(X_expected_cm.data(), x_expected_rm.data(), N_in, M_in, N_in, LDC_in, sizeof(double));
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
                slicot_transpose_to_c_with_ld(Z_expected_cm.data(), z_expected_rm.data(), M_in, M_in, M_in, LDZ_in, sizeof(double));
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
TEST_F(SB04MDTest, DocExample_ColMajor) {
    SetUpExampleData(false); // false for column-major

    int info = slicot_sb04md(N_in, M_in,
                             A_io.data(), LDA_in, B_io.data(), LDB_in,
                             C_io.data(), LDC_in, Z_out.data(), LDZ_in,
                             0 /* col-major */);

    ASSERT_EQ(info, INFO_expected);
    CheckSolution(C_io, Z_out, false);
}

// Test: Documentation Example (Row-Major)
TEST_F(SB04MDTest, DocExample_RowMajor) {
    SetUpExampleData(true); // true for row-major

    int info = slicot_sb04md(N_in, M_in,
                             A_io.data(), LDA_in, B_io.data(), LDB_in,
                             C_io.data(), LDC_in, Z_out.data(), LDZ_in,
                             1 /* row-major */);
    
    ASSERT_EQ(info, INFO_expected);
    CheckSolution(C_io, Z_out, true);
}

// Test: Zero Dimensions (N=0)
TEST_F(SB04MDTest, ZeroDimensions_N0) {
    N_in = 0; M_in = 2;
    InitializeDimensions(N_in, M_in, false); // col-major

    // A, C are N x N and N x M, so effectively empty or not referenced if N=0.
    // B, Z are M x M.
    std::vector<double> b_dummy_cm = {1,0,0,1}; // M x M identity
    if (M_in > 0) memcpy(B_io.data(), b_dummy_cm.data(), b_dummy_cm.size()*sizeof(double));


    int info = slicot_sb04md(N_in, M_in,
                             nullptr, LDA_in, B_io.data(), LDB_in,
                             nullptr, LDC_in, Z_out.data(), LDZ_in,
                             0 /* col-major */);
    EXPECT_EQ(info, 0);
    // If N=0, X is 0xM, so C_io content is not defined as solution.
    // Z should be computed if M > 0. B is modified.
    // The example doesn't specify Z for N=0, but it should be an M x M orthogonal matrix.
    // For simplicity, just check info.
}

// Test: Zero Dimensions (M=0)
TEST_F(SB04MDTest, ZeroDimensions_M0) {
    N_in = 2; M_in = 0;
    InitializeDimensions(N_in, M_in, false); // col-major

    // B, C, Z are M x M or N x M, so effectively empty or not referenced if M=0.
    // A is N x N.
    std::vector<double> a_dummy_cm = {1,0,0,1}; // N x N identity
    if (N_in > 0) memcpy(A_io.data(), a_dummy_cm.data(), a_dummy_cm.size()*sizeof(double));

    int info = slicot_sb04md(N_in, M_in,
                             A_io.data(), LDA_in, nullptr, LDB_in,
                             nullptr, LDC_in, nullptr, LDZ_in,
                             0 /* col-major */);
    EXPECT_EQ(info, 0);
    // If M=0, X is Nx0, Z is 0x0. C_io, Z_out content not defined as solution.
}

// Test: Zero Dimensions (N=0, M=0)
TEST_F(SB04MDTest, ZeroDimensions_N0_M0) {
    N_in = 0; M_in = 0;
    InitializeDimensions(N_in, M_in, false); // col-major

    int info = slicot_sb04md(N_in, M_in,
                             nullptr, LDA_in, nullptr, LDB_in,
                             nullptr, LDC_in, nullptr, LDZ_in,
                             0 /* col-major */);
    EXPECT_EQ(info, 0);
}

// Test: Parameter Validation (Selected)
TEST_F(SB04MDTest, ParameterValidation) {
    SetUpExampleData(false); // Initialize with valid N, M etc. for non-NULL pointers

    // Test invalid N
    int info = slicot_sb04md(-1, M_in, A_io.data(), LDA_in, B_io.data(), LDB_in, C_io.data(), LDC_in, Z_out.data(), LDZ_in, 0);
    EXPECT_EQ(info, -1);

    // Test invalid M
    info = slicot_sb04md(N_in, -1, A_io.data(), LDA_in, B_io.data(), LDB_in, C_io.data(), LDC_in, Z_out.data(), LDZ_in, 0);
    EXPECT_EQ(info, -2);
    
    // Test invalid LDA (assuming N_in > 0 from SetUpExampleData)
    if (N_in > 0) {
        info = slicot_sb04md(N_in, M_in, A_io.data(), 0, B_io.data(), LDB_in, C_io.data(), LDC_in, Z_out.data(), LDZ_in, 0);
        EXPECT_EQ(info, -4); // LDA is arg 4
    }
    
    // Test invalid LDB (assuming M_in > 0 from SetUpExampleData)
    if (M_in > 0) {
        info = slicot_sb04md(N_in, M_in, A_io.data(), LDA_in, B_io.data(), 0, C_io.data(), LDC_in, Z_out.data(), LDZ_in, 0);
        EXPECT_EQ(info, -6); // LDB is arg 6
    }
    // Note: SLICOT routines usually don't check for NULL pointers if dimensions are zero.
    // The wrapper itself should handle NULLs appropriately if dimensions allow.
    // The Fortran routine expects valid pointers if dimensions > 0.
    // The C wrapper's internal NULL checks are not explicitly tested here, but covered by zero-dim tests.
}
