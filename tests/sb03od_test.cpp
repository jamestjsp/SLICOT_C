#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <cstring>   // For memcpy

#include "sb03od.h"
#include "slicot_utils.h" // For transpose functions

// Helper to print matrix (col-major)
void print_matrix_sb03od_cm(const char* name, const double* mat, int rows, int cols, int ld) {
    std::cout << name << " (Col-Major, " << rows << "x" << cols << ", LD=" << ld << "):" << std::endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(4) << mat[j * ld + i] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Helper to print matrix (row-major)
void print_matrix_sb03od_rm(const char* name, const double* mat, int rows, int cols, int ld) {
    std::cout << name << " (Row-Major, " << rows << "x" << cols << ", LD=" << ld << "):" << std::endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(4) << mat[i * ld + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


// --- Test Fixture ---
class SB03ODTest : public ::testing::Test {
protected:
    char DICO_in = 'C';
    char FACT_in = 'N';
    char TRANS_in = 'N';
    int N_in = 4;
    int M_in = 5;

    std::vector<double> A_io;
    std::vector<double> Q_io;
    std::vector<double> B_io; // Input B, Output U

    double SCALE_out = 0.0;
    std::vector<double> WR_out;
    std::vector<double> WI_out;

    int LDA_in, LDQ_in, LDB_in;

    std::vector<double> U_expected_cm; // Expected U (col-major)
    double SCALE_expected = 1.0;
    int INFO_expected = 0;

    double check_tol = 1e-4; // Based on example precision

    void InitializeDimensions(int n, int m, bool row_major_layout) {
        N_in = n;
        M_in = m;

        WR_out.resize(N_in > 0 ? N_in : 1); // Min size 1 for safety if N=0
        WI_out.resize(N_in > 0 ? N_in : 1);

        if (row_major_layout) {
            LDA_in = std::max(1, N_in); 
            LDQ_in = std::max(1, N_in); 
            
            // For B_io (row-major):
            // LDB_in is the number of columns. For output U (N x N), LDB_in must be N_in.
            // The input B can be M_in x N_in (TRANS='N') or N_in x M_in (TRANS='T').
            // B_io must accommodate the larger of input rows or N_in (for U's rows),
            // and LDB_in (N_in) columns.
            LDB_in = std::max(1, N_in); // Number of columns for B_io will be N_in.

            int b_io_rows;
            if (TRANS_in == 'N') {
                b_io_rows = std::max(M_in, N_in); // Max of input B rows (M_in) or output U rows (N_in)
            } else { // TRANS_in == 'T'
                b_io_rows = N_in; // Input B rows (N_in), output U rows (N_in)
            }
            B_io.resize((size_t)b_io_rows * LDB_in);

        } else { // Column-major
            LDA_in = std::max(1, N_in); // rows for A (NxN)
            LDQ_in = std::max(1, N_in); // rows for Q (NxN)
            // LDB for B_io (col-major):
            // Fortran LDB: MAX(1,N,M) if TRANS='N'; MAX(1,N) if TRANS='T'.
            if (TRANS_in == 'N') {
                LDB_in = std::max(1, std::max(N_in, M_in));
            } else { // TRANS_in == 'T'
                LDB_in = std::max(1, N_in);
            }
        }
        
        A_io.resize((size_t)LDA_in * N_in); // A is N x N
        Q_io.resize((size_t)LDQ_in * N_in); // Q is N x N

        // B_io size depends on layout and TRANS for input, but must hold N x N for output U.
        if (row_major_layout) {
            // B_io (row-major) for input:
            // if TRANS='N', M_in rows, N_in cols. LDB_in (cols) = N_in.
            // if TRANS='T', N_in rows, M_in cols. LDB_in (cols) = M_in.
            // B_io (row-major) for output U: N_in rows, N_in cols. LDB_in (cols) = N_in.
            // B_io must be large enough for input and for the N_in x N_in output U.
            // LDB_in passed to slicot_sb03od should be the number of columns of B_io
            // as it's prepared for output U, i.e., N_in.
            // The input data will be copied into this structure.

            int b_input_actual_rows = (TRANS_in == 'N') ? M_in : N_in;
            int b_input_actual_cols = (TRANS_in == 'N') ? N_in : M_in;
            
            // LDB_in for the B_io array that will receive U (NxN)
            LDB_in = std::max(1, N_in); 
            // B_io must be at least N_in x N_in for U.
            // If input is larger (e.g. M_in > N_in for TRANS='N'), B_io needs to be M_in x N_in.
            B_io.resize((size_t)std::max(b_input_actual_rows, N_in) * LDB_in);

        } else { // Column-major
            // B_io (col-major) for input:
            // if TRANS='N', LDB_in x N cols.
            // if TRANS='T', LDB_in x M cols.
            // B_io (col-major) for output U: LDB_in x N cols (U is NxN part).
            int b_cm_cols = (TRANS_in == 'N') ? N_in : M_in;
            if (TRANS_in == 'T' && N_in > M_in) b_cm_cols = N_in; // Doc: B(LDB,max(M,N)) if TRANS='T'
                                                                 // And output U is NxN.
            b_cm_cols = N_in; // To be safe for output U.
            B_io.resize((size_t)LDB_in * b_cm_cols);
        }
    }

    // Setup based on SB03OD.html example
    void SetUpExampleData(bool row_major_layout) {
        DICO_in = 'C'; FACT_in = 'N'; TRANS_in = 'N';
        N_in = 4; M_in = 5;
        InitializeDimensions(N_in, M_in, row_major_layout);

        std::vector<double> a_cm_data = {
            -1.0, -1.0,  2.0,  2.0,
            37.0,-10.0, -4.0,  2.0,
           -12.0,  0.0,  7.0,  7.0,
           -12.0,  4.0, -6.0, -9.0
        };
        // B input (TRANS='N', M x N = 5 x 4, col-major)
        std::vector<double> b_input_cm_data = {
            1.0,  0.0, -1.0,  1.0, -1.0, // col 1
            2.5,  1.0, -2.5,  2.5, -2.5, // col 2
            1.0,  0.0, -1.0,  4.0, -4.0, // col 3
            3.5,  1.0, -1.5, -5.5,  3.5  // col 4
        };

        U_expected_cm = { // Col-major
            1.0, 0.0, 0.0, 0.0,
            3.0, 1.0, 0.0, 0.0,
            2.0,-1.0, 1.0, 0.0,
           -1.0, 1.0,-2.0, 1.0
        };
        SCALE_expected = 1.0;
        INFO_expected = 0;

        if (row_major_layout) {
            slicot_transpose_to_c_with_ld(a_cm_data.data(), A_io.data(), N_in, N_in, N_in, LDA_in, sizeof(double));
            // Q_io is output if FACT='N'
            
            // Prepare B_io for row-major input (M_in x N_in if TRANS='N')
            // LDB_in is N_in (cols of B_io, which will hold U)
            // b_input_cm_data is M_in x N_in (col-major), ld_cm = M_in
            // temp_b_rm_input is M_in x N_in (row-major), ld_rm = N_in
            std::vector<double> temp_b_rm_input(M_in * N_in);
            // Manual copy from column-major b_input_cm_data to row-major temp_b_rm_input
            // to ensure temp_b_rm_input[i][j] (math) = b_input_cm_data[i][j] (math)
            int ld_cm_b_input = M_in; // Leading dimension of b_input_cm_data
            int ld_rm_b_temp = N_in;  // Leading dimension of temp_b_rm_input
            for (int i = 0; i < M_in; ++i) { // iterate mathematical rows
                for (int j = 0; j < N_in; ++j) { // iterate mathematical columns
                    temp_b_rm_input[i * ld_rm_b_temp + j] = b_input_cm_data[i + j * ld_cm_b_input];
                }
            }

            for(int i=0; i<M_in; ++i) { 
                for(int j=0; j<N_in; ++j) {
                    if (i < std::max(M_in, N_in) && j < LDB_in) { // Check bounds of B_io
                         B_io[i * LDB_in + j] = temp_b_rm_input[i * N_in + j];
                    }
                }
            }
        } else { // Column-major
            memcpy(A_io.data(), a_cm_data.data(), a_cm_data.size() * sizeof(double));
            // Q_io is output
            // B_io (col-major, LDB_in rows, N_in cols for input)
            memcpy(B_io.data(), b_input_cm_data.data(), b_input_cm_data.size() * sizeof(double));
        }
    }

    void CheckSolutionU(const std::vector<double>& u_computed, bool is_row_major_computed) {
        if (N_in == 0) return;
        
        int ld_computed = LDB_in; // For row-major, LDB_in is N_in (cols of U)
                                  // For col-major, LDB_in is rows of B_io (>=N_in)

        if (is_row_major_computed) {
            std::vector<double> u_expected_rm(N_in * N_in);
            slicot_transpose_to_c_with_ld(U_expected_cm.data(), u_expected_rm.data(), N_in, N_in, N_in, N_in, sizeof(double));
            for (int i = 0; i < N_in; ++i) {
                for (int j = 0; j < N_in; ++j) {
                    EXPECT_NEAR(u_computed[i * ld_computed + j], u_expected_rm[i * N_in + j], check_tol)
                        << "U_rm(" << i << "," << j << ") mismatch";
                }
            }
        } else { // Column-major computed (U is upper triangular in B_io)
            for (int j = 0; j < N_in; ++j) { // col
                for (int i = 0; i <= j; ++i) { // row (upper triangle including diagonal)
                     EXPECT_NEAR(u_computed[j * ld_computed + i], U_expected_cm[j * N_in + i], check_tol)
                        << "U_cm(" << i << "," << j << ") mismatch";
                }
                 for (int i = j + 1; i < N_in; ++i) { // row (lower triangle, should be zero)
                     EXPECT_NEAR(u_computed[j * ld_computed + i], 0.0, check_tol)
                        << "U_cm(" << i << "," << j << ") non-zero in lower triangle";
                }
            }
        }
    }
};

// Test: Basic Continuous-Time System (Col-Major) - From SB03OD.html example
TEST_F(SB03ODTest, ContinuousLyapunov_ColMajor_FactN_TransN) {
    SetUpExampleData(false); // false for column-major

    int info = slicot_sb03od(DICO_in, FACT_in, TRANS_in, N_in, M_in,
                             A_io.data(), LDA_in, Q_io.data(), LDQ_in,
                             B_io.data(), LDB_in, &SCALE_out,
                             WR_out.data(), WI_out.data(), 0 /* col-major */);

    ASSERT_EQ(info, INFO_expected);
    EXPECT_NEAR(SCALE_out, SCALE_expected, check_tol);
    CheckSolutionU(B_io, false); // B_io now contains U
}

// Test: Basic Continuous-Time System (Row-Major)
TEST_F(SB03ODTest, ContinuousLyapunov_RowMajor_FactN_TransN) {
    SetUpExampleData(true); // true for row-major

    int info = slicot_sb03od(DICO_in, FACT_in, TRANS_in, N_in, M_in,
                             A_io.data(), LDA_in, Q_io.data(), LDQ_in,
                             B_io.data(), LDB_in, &SCALE_out,
                             WR_out.data(), WI_out.data(), 1 /* row-major */);
    
    ASSERT_EQ(info, INFO_expected);
    EXPECT_NEAR(SCALE_out, SCALE_expected, check_tol);
    CheckSolutionU(B_io, true); // B_io now contains U
}


// Test: Zero Dimensions (N=0)
TEST_F(SB03ODTest, ZeroDimensions_N0) {
    DICO_in = 'C'; FACT_in = 'N'; TRANS_in = 'N';
    N_in = 0; M_in = 2; // M can be > 0
    InitializeDimensions(N_in, M_in, false); // col-major. This sets LDB_in correctly.

    double* null_a = nullptr; double* null_q = nullptr; double* null_b = nullptr;
    double* null_wr = nullptr; double* null_wi = nullptr;
    // LDA_in, LDQ_in, LDB_in are set by InitializeDimensions.
    // For N=0, M=2, TRANS='N', LDB_in should be MAX(1,M_in) = 2.
    // LDA_in, LDQ_in should be 1.

    int info = slicot_sb03od(DICO_in, FACT_in, TRANS_in, N_in, M_in,
                             null_a, LDA_in, null_q, LDQ_in,
                             null_b, LDB_in, &SCALE_out,
                             null_wr, null_wi, 0 /* col-major */);

    EXPECT_EQ(info, 0);
    EXPECT_EQ(SCALE_out, 1.0); // Fortran SB03OD sets SCALE=1.0 if N=0
}

// Test: Zero Dimensions (M=0)
TEST_F(SB03ODTest, ZeroDimensions_M0) {
    DICO_in = 'C'; FACT_in = 'N'; TRANS_in = 'N';
    N_in = 2; M_in = 0;
    InitializeDimensions(N_in, M_in, false); // col-major
    A_io = {1,0,0,1}; // Dummy A
    // B can be NULL if M=0. Fortran SB03OD sets U=0 if M=0, N>0.
    
    int info = slicot_sb03od(DICO_in, FACT_in, TRANS_in, N_in, M_in,
                             A_io.data(), LDA_in, Q_io.data(), LDQ_in,
                             B_io.data(), LDB_in, &SCALE_out, // B_io will contain U
                             WR_out.data(), WI_out.data(), 0 /* col-major */);

    EXPECT_EQ(info, 0);
    EXPECT_EQ(SCALE_out, 1.0); // Fortran SB03OD sets SCALE=1.0 if M=0
    // Check if U (in B_io) is zero
    for(int j=0; j<N_in; ++j) {
        for(int i=0; i<N_in; ++i) {
            EXPECT_NEAR(B_io[j*LDB_in + i], 0.0, 1e-9);
        }
    }
}


// Test: Parameter Validation (Selected)
TEST_F(SB03ODTest, ParameterValidation) {
    SetUpExampleData(false); // Initialize with valid N, M etc.

    // Test invalid DICO
    int info = slicot_sb03od('X', FACT_in, TRANS_in, N_in, M_in, A_io.data(), LDA_in, Q_io.data(), LDQ_in, B_io.data(), LDB_in, &SCALE_out, WR_out.data(), WI_out.data(), 0);
    EXPECT_EQ(info, -1);

    // Test invalid N
    info = slicot_sb03od(DICO_in, FACT_in, TRANS_in, -1, M_in, A_io.data(), LDA_in, Q_io.data(), LDQ_in, B_io.data(), LDB_in, &SCALE_out, WR_out.data(), WI_out.data(), 0);
    EXPECT_EQ(info, -4);
    
    // Test invalid LDA
    if (N_in > 0) {
        info = slicot_sb03od(DICO_in, FACT_in, TRANS_in, N_in, M_in, A_io.data(), 0, Q_io.data(), LDQ_in, B_io.data(), LDB_in, &SCALE_out, WR_out.data(), WI_out.data(), 0);
        EXPECT_EQ(info, -7);
    }
    
    // Test NULL A when N > 0
    if (N_in > 0) {
        info = slicot_sb03od(DICO_in, FACT_in, TRANS_in, N_in, M_in, nullptr, LDA_in, Q_io.data(), LDQ_in, B_io.data(), LDB_in, &SCALE_out, WR_out.data(), WI_out.data(), 0);
        EXPECT_EQ(info, -6);
    }
     // Test NULL SCALE
    info = slicot_sb03od(DICO_in, FACT_in, TRANS_in, N_in, M_in, A_io.data(), LDA_in, Q_io.data(), LDQ_in, B_io.data(), LDB_in, nullptr, WR_out.data(), WI_out.data(), 0);
    EXPECT_EQ(info, -12);

}

// TODO: Add more tests:
// - Discrete-time
// - TRANS='T'
// - FACT='F' (requires providing Schur A and Q)
// - Cases leading to INFO = 1, 2, 3 etc.

