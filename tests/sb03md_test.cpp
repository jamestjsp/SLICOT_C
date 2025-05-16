#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <cstring>   // For memcpy

#include "sb03md.h"
#include "slicot_utils.h" // For slicot_transpose_to_c_with_ld etc.

// Helper function to print matrices (for debugging)
void print_matrix_col_major(const char* name, const double* mat, int rows, int cols, int ld) {
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

void print_matrix_row_major(const char* name, const double* mat, int rows, int cols, int ld) {
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
class SB03MDTest : public ::testing::Test {
protected:
    // Default parameters
    char DICO_in = 'D';    // Discrete-time from example
    char JOB_in = 'X';     // Compute solution X only
    char FACT_in = 'N';    // Compute Schur factorization
    char TRANA_in = 'N';   // op(A) = A
    int N_in = 3;

    // Input matrices (will be initialized in SetUp or specific tests)
    std::vector<double> A_io;
    std::vector<double> U_io; // Used if FACT='F' or output if FACT='N'
    std::vector<double> C_io; // Input C, output X

    // Output scalars & vectors
    double SCALE_out = 0.0;
    double SEP_out = 0.0;
    double FERR_out = 0.0;
    std::vector<double> WR_out;
    std::vector<double> WI_out;

    // Leading dimensions
    int LDA_in, LDU_in, LDC_in;

    // Expected results (for specific test cases)
    std::vector<double> X_expected;
    double SCALE_expected = 1.0;
    int INFO_expected = 0;

    // Test verification tolerance
    double check_tol = 1e-4; // Based on example output precision

    void InitializeDimensions(int n, bool row_major_layout) {
        N_in = n;
        if (row_major_layout) {
            LDA_in = std::max(1, N_in); // cols
            LDU_in = std::max(1, N_in); // cols
            LDC_in = (JOB_in == 'S') ? 1 : std::max(1, N_in); // cols
        } else { // Column-major
            LDA_in = std::max(1, N_in); // rows
            LDU_in = std::max(1, N_in); // rows
            LDC_in = (JOB_in == 'S') ? 1 : std::max(1, N_in); // rows
        }

        A_io.resize((size_t)LDA_in * N_in); // Actual size depends on layout
        U_io.resize((size_t)LDU_in * N_in);
        if (JOB_in != 'S') {
            C_io.resize((size_t)LDC_in * N_in);
        } else {
            C_io.clear(); // Not referenced
        }
        WR_out.resize(N_in);
        WI_out.resize(N_in);
    }

    // Setup based on SB03MD.html example
    void SetUpExampleData(bool row_major_layout) {
        DICO_in = 'D'; JOB_in = 'X'; FACT_in = 'N'; TRANA_in = 'N';
        N_in = 3;
        InitializeDimensions(N_in, row_major_layout);

        std::vector<double> a_cm_data = {3.0, 1.0, 0.0, 1.0, 3.0, 0.0, 1.0, 0.0, 3.0}; // A col-major
        std::vector<double> c_cm_data = {25.0, 24.0, 15.0, 24.0, 32.0, 8.0, 15.0, 8.0, 40.0}; // C col-major (symmetric)
        X_expected = {2.0, 1.0, 1.0, 1.0, 3.0, 0.0, 1.0, 0.0, 4.0}; // X col-major (symmetric)
        SCALE_expected = 1.0;
        INFO_expected = 0;

        if (row_major_layout) {
            slicot_transpose_to_c_with_ld(a_cm_data.data(), A_io.data(), N_in, N_in, N_in, LDA_in, sizeof(double));
            if (JOB_in != 'S') {
                // For row-major input C, provide upper triangle (or full)
                // The wrapper will use slicot_transpose_symmetric_to_fortran with UPLO='U'
                // So, ensure A_io has at least the upper triangle correctly for that.
                // Here, c_cm_data is full, so direct transpose is fine.
                slicot_transpose_to_c_with_ld(c_cm_data.data(), C_io.data(), N_in, N_in, N_in, LDC_in, sizeof(double));
            }
        } else {
            memcpy(A_io.data(), a_cm_data.data(), a_cm_data.size() * sizeof(double));
            if (JOB_in != 'S') {
                // For col-major input C, provide upper triangle (or full)
                // The wrapper will use slicot_copy_symmetric_part with UPLO='U'
                memcpy(C_io.data(), c_cm_data.data(), c_cm_data.size() * sizeof(double));
            }
        }
        // U_io is output if FACT='N', so no need to initialize its content here.
        // WR_out, WI_out are outputs.
    }

    void CheckSolution(const std::vector<double>& x_computed, const std::vector<double>& x_expected_cm, bool is_row_major_computed) {
        if (N_in == 0) return;
        ASSERT_EQ(x_computed.size(), (size_t)LDC_in * N_in); // LDC_in is cols for RM, rows for CM

        if (is_row_major_computed) {
            std::vector<double> x_expected_rm(N_in * N_in);
            slicot_transpose_to_c_with_ld(x_expected_cm.data(), x_expected_rm.data(), N_in, N_in, N_in, N_in, sizeof(double));
            for (int i = 0; i < N_in; ++i) {
                for (int j = 0; j < N_in; ++j) {
                    EXPECT_NEAR(x_computed[i * LDC_in + j], x_expected_rm[i * N_in + j], check_tol)
                        << "X_rm(" << i << "," << j << ") mismatch";
                }
            }
        } else { // Column-major computed (wrapper fills upper triangle)
            for (int j = 0; j < N_in; ++j) { // col
                for (int i = 0; i <= j; ++i) { // row (upper triangle including diagonal)
                     EXPECT_NEAR(x_computed[j * LDC_in + i], x_expected_cm[j * N_in + i], check_tol)
                        << "X_cm(" << i << "," << j << ") mismatch";
                }
            }
        }
    }
};

// Test: Basic Discrete-Time System (Column-Major) - From SB03MD.html example
TEST_F(SB03MDTest, DiscreteLyapunov_ColMajor_FactN) {
    SetUpExampleData(false); // false for column-major

    int info = slicot_sb03md(DICO_in, JOB_in, FACT_in, TRANA_in, N_in,
                             A_io.data(), LDA_in, U_io.data(), LDU_in,
                             C_io.data(), LDC_in, &SCALE_out, &SEP_out, &FERR_out,
                             WR_out.data(), WI_out.data(), 0 /* col-major */);

    ASSERT_EQ(info, INFO_expected);
    EXPECT_NEAR(SCALE_out, SCALE_expected, check_tol);
    CheckSolution(C_io, X_expected, false); // C_io now contains X
}

// Test: Basic Discrete-Time System (Row-Major)
TEST_F(SB03MDTest, DiscreteLyapunov_RowMajor_FactN) {
    SetUpExampleData(true); // true for row-major

    int info = slicot_sb03md(DICO_in, JOB_in, FACT_in, TRANA_in, N_in,
                             A_io.data(), LDA_in, U_io.data(), LDU_in,
                             C_io.data(), LDC_in, &SCALE_out, &SEP_out, &FERR_out,
                             WR_out.data(), WI_out.data(), 1 /* row-major */);

    ASSERT_EQ(info, INFO_expected);
    EXPECT_NEAR(SCALE_out, SCALE_expected, check_tol);
    CheckSolution(C_io, X_expected, true); // C_io now contains X
}

// Test: Continuous-Time System (Column-Major)
TEST_F(SB03MDTest, ContinuousLyapunov_ColMajor_FactN) {
    DICO_in = 'C'; JOB_in = 'X'; FACT_in = 'N'; TRANA_in = 'N';
    N_in = 2;
    InitializeDimensions(N_in, false); // col-major

    // A = [0, 1; -2, -3], C_in = [0, 0; 0, -1] (symmetric)
    // Solves A'X + XA = C_in
    // Solution X = [1/3, 0; 0, 1/6]
    A_io = {0.0, -2.0, 1.0, -3.0}; // A column-major: A(0,0)=0, A(1,0)=-2, A(0,1)=1, A(1,1)=-3
    C_io = {0.0, 0.0, 0.0, -1.0};  // C_in column-major (symmetric): C(0,0)=0, C(1,0)=0, C(0,1)=0, C(1,1)=-1
                                   // Wrapper expects upper triangle for C input.
    X_expected = {1.0/3.0, 0.0, 0.0, 1.0/6.0}; // X_expected column-major
    
    SCALE_expected = 1.0; // Expect scale to be 1 for this well-behaved problem
    INFO_expected = 0;

    int info = slicot_sb03md(DICO_in, JOB_in, FACT_in, TRANA_in, N_in,
                             A_io.data(), LDA_in, U_io.data(), LDU_in,
                             C_io.data(), LDC_in, &SCALE_out, &SEP_out, &FERR_out,
                             WR_out.data(), WI_out.data(), 0 /* col-major */);
    
    ASSERT_EQ(info, INFO_expected);
    EXPECT_NEAR(SCALE_out, SCALE_expected, check_tol);
    CheckSolution(C_io, X_expected, false);
}


// Test: Zero Dimensions (N=0)
TEST_F(SB03MDTest, ZeroDimensions) {
    DICO_in = 'C'; JOB_in = 'X'; FACT_in = 'N'; TRANA_in = 'N';
    N_in = 0;
    InitializeDimensions(N_in, false); // col-major

    // For N=0, A, U, C, WR, WI can be NULL or point to dummy.
    // SCALE, SEP, FERR must be valid pointers.
    // Fortran SB03MD sets SCALE=1.0, SEP=0.0 (if JOB='S'/'B') for N=0.
    
    double dummy_val = 0.0;
    double* null_a = nullptr; // Or &dummy_val if routine requires non-NULL for LD >= 1
    double* null_u = nullptr;
    double* null_c = nullptr;
    double* null_wr = nullptr;
    double* null_wi = nullptr;
    
    // If N=0, LDA, LDU, LDC must be >= 1.
    LDA_in = 1; LDU_in = 1; LDC_in = 1;


    int info = slicot_sb03md(DICO_in, JOB_in, FACT_in, TRANA_in, N_in,
                             null_a, LDA_in, null_u, LDU_in,
                             null_c, LDC_in, &SCALE_out, &SEP_out, &FERR_out,
                             null_wr, null_wi, 0 /* col-major */);

    EXPECT_EQ(info, 0);
    EXPECT_EQ(SCALE_out, 1.0);
    // SEP and FERR are not referenced if JOB='X'
}

// Test: Parameter Validation (Selected)
TEST_F(SB03MDTest, ParameterValidation) {
    SetUpExampleData(false); // Initialize with valid N, M, P etc.

    // Test invalid DICO
    int info = slicot_sb03md('Z', JOB_in, FACT_in, TRANA_in, N_in, A_io.data(), LDA_in, U_io.data(), LDU_in, C_io.data(), LDC_in, &SCALE_out, &SEP_out, &FERR_out, WR_out.data(), WI_out.data(), 0);
    EXPECT_EQ(info, -1);

    // Test invalid N
    info = slicot_sb03md(DICO_in, JOB_in, FACT_in, TRANA_in, -1, A_io.data(), LDA_in, U_io.data(), LDU_in, C_io.data(), LDC_in, &SCALE_out, &SEP_out, &FERR_out, WR_out.data(), WI_out.data(), 0);
    EXPECT_EQ(info, -5);
    
    // Test invalid LDA
    if (N_in > 0) {
        info = slicot_sb03md(DICO_in, JOB_in, FACT_in, TRANA_in, N_in, A_io.data(), 0, U_io.data(), LDU_in, C_io.data(), LDC_in, &SCALE_out, &SEP_out, &FERR_out, WR_out.data(), WI_out.data(), 0);
        EXPECT_EQ(info, -7); // LDA is arg 7
    }
    
    // Test NULL A when N > 0
    if (N_in > 0) {
       // The C wrapper should have checks for this.
       // A is arg 6.
        info = slicot_sb03md(DICO_in, JOB_in, FACT_in, TRANA_in, N_in, nullptr, LDA_in, U_io.data(), LDU_in, C_io.data(), LDC_in, &SCALE_out, &SEP_out, &FERR_out, WR_out.data(), WI_out.data(), 0);
        EXPECT_EQ(info, -6); // Expect error for NULL A (arg 6)
    }
}

// TODO: Add more tests:
// - JOB='S', JOB='B'
// - FACT='F'
// - TRANA='T', TRANA='C'
// - Cases leading to INFO = N+1 (singular/nearly singular)
// - Test with actual Schur inputs for FACT='F'

