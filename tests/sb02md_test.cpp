#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max

#include "sb02md.h"
#include "slicot_utils.h" // For transpose functions

// --- Column-Major Test Fixture ---
class SB02MDTestColMajor : public ::testing::Test {
protected:
    // Parameters from SB02MD.html example
    char DICO_in = 'C';
    char HINV_in = 'D'; // Not used for DICO='C'
    char UPLO_in = 'U';
    char SCAL_in = 'N';
    char SORT_in = 'S';
    int N_in = 2;

    // Input matrices (column-major format)
    std::vector<double> A_in;
    std::vector<double> G_in;
    std::vector<double> Q_in_out; // Q on input, X (solution) on output

    // Output scalars
    double RCOND_out = -1.0;

    // Output arrays (will be sized in SetUp)
    std::vector<double> WR_out;
    std::vector<double> WI_out;
    std::vector<double> S_out;
    std::vector<double> U_out;

    // Expected results
    std::vector<double> X_expected_colmajor; // Solution X
    double RCOND_expected = 0.30901699437494745; // Updated from 0.31
    int INFO_expected = 0;

    // Leading dimensions
    int LDA_in = 0;
    int LDG_in = 0;
    int LDQ_in_out = 0;
    int LDS_out = 0;
    int LDU_out = 0;

    // Verification tolerance
    double check_tol = 1e-4; // As per example output precision

    void SetUp() override {
        // Define leading dimensions for Column Major (Fortran-style: LD is number of rows)
        LDA_in = std::max(1, N_in);
        LDG_in = std::max(1, N_in);
        LDQ_in_out = std::max(1, N_in);
        LDS_out = std::max(1, 2 * N_in);
        LDU_out = std::max(1, 2 * N_in);

        // Initialize input matrices (column-major format from SB02MD.html example)
        // A = [[0,1],[0,0]]
        A_in = {0.0, 0.0, 1.0, 0.0};
        // G = [[0,0],[0,1]] (symmetric, UPLO='U')
        G_in = {0.0, 0.0, 0.0, 1.0}; // Full matrix provided for simplicity
        // Q = [[1,0],[0,2]] (symmetric, UPLO='U')
        Q_in_out = {1.0, 0.0, 0.0, 2.0}; // Full matrix provided for simplicity

        // Resize output arrays
        WR_out.resize(2 * (size_t)N_in);
        WI_out.resize(2 * (size_t)N_in);
        S_out.resize((size_t)LDS_out * (2 * N_in));
        U_out.resize((size_t)LDU_out * (2 * N_in));

        // Expected solution X = [[2,1],[1,2]] (column-major)
        X_expected_colmajor = {2.0, 1.0, 1.0, 2.0};
    }
};

// --- Row-Major Test Fixture ---
class SB02MDTestRowMajor : public SB02MDTestColMajor {
protected:
    std::vector<double> A_rm_in;
    std::vector<double> G_rm_in;
    std::vector<double> Q_rm_in_out;
    std::vector<double> S_rm_out;
    std::vector<double> U_rm_out;
    std::vector<double> X_expected_rowmajor;

    void SetUp() override {
        SB02MDTestColMajor::SetUp(); // Base class SetUp called first

        // Adjust leading dimensions for Row Major (C-style: LD is number of columns)
        LDA_in = std::max(1, N_in);     // N cols for A (NxN)
        LDG_in = std::max(1, N_in);     // N cols for G (NxN)
        LDQ_in_out = std::max(1, N_in); // N cols for Q (NxN)
        LDS_out = std::max(1, 2 * N_in); // 2N cols for S (2Nx2N)
        LDU_out = std::max(1, 2 * N_in); // 2N cols for U (2Nx2N)

        // Resize row-major matrices
        A_rm_in.resize((size_t)N_in * LDA_in);
        G_rm_in.resize((size_t)N_in * LDG_in);
        Q_rm_in_out.resize((size_t)N_in * LDQ_in_out);
        S_rm_out.resize((size_t)(2 * N_in) * LDS_out);
        U_rm_out.resize((size_t)(2 * N_in) * LDU_out);
        X_expected_rowmajor.resize(X_expected_colmajor.size());

        // Transpose column-major inputs to row-major
        if (!A_in.empty()) slicot_transpose_to_c_with_ld(A_in.data(), A_rm_in.data(), N_in, N_in, std::max(1,N_in), LDA_in, sizeof(double));
        // For G and Q, which are symmetric and UPLO is used, provide full matrices for simplicity in test setup.
        // The wrapper will handle them correctly.
        if (!G_in.empty()) slicot_transpose_to_c_with_ld(G_in.data(), G_rm_in.data(), N_in, N_in, std::max(1,N_in), LDG_in, sizeof(double));
        if (!Q_in_out.empty()) slicot_transpose_to_c_with_ld(Q_in_out.data(), Q_rm_in_out.data(), N_in, N_in, std::max(1,N_in), LDQ_in_out, sizeof(double));
        
        // Transpose expected X to row-major
        if (!X_expected_colmajor.empty()) slicot_transpose_to_c_with_ld(X_expected_colmajor.data(), X_expected_rowmajor.data(), N_in, N_in, std::max(1,N_in), N_in, sizeof(double));
    }
};

// --- Test Cases ---

TEST_F(SB02MDTestColMajor, DocExample) {
    // Q_in_out is used for both input Q and output X
    std::vector<double> Q_test = Q_in_out; // Make a copy for the test

    int info_out = slicot_sb02md(
        DICO_in, HINV_in, UPLO_in, SCAL_in, SORT_in, N_in,
        A_in.data(), LDA_in, G_in.data(), LDG_in, Q_test.data(), LDQ_in_out,
        &RCOND_out, WR_out.data(), WI_out.data(),
        S_out.data(), LDS_out, U_out.data(), LDU_out,
        0 // row_major = false
    );

    ASSERT_EQ(info_out, INFO_expected);
    EXPECT_NEAR(RCOND_out, RCOND_expected, check_tol);
    for (size_t i = 0; i < X_expected_colmajor.size(); ++i) {
        // If UPLO_in is 'U', solution X is in upper triangle of Q_test.
        // If UPLO_in is 'L', solution X is in lower triangle of Q_test.
        // The wrapper ensures Q_test (X) is fully populated if col-major and UPLO is used for input.
        // For this example, X_expected_colmajor is full.
        EXPECT_NEAR(Q_test[i], X_expected_colmajor[i], check_tol) << "X mismatch at index " << i;
    }
}

TEST_F(SB02MDTestRowMajor, DocExample) {
    // Q_rm_in_out is used for both input Q and output X
    std::vector<double> Q_rm_test = Q_rm_in_out; // Make a copy

    int info_out = slicot_sb02md(
        DICO_in, HINV_in, UPLO_in, SCAL_in, SORT_in, N_in,
        A_rm_in.data(), LDA_in, G_rm_in.data(), LDG_in, Q_rm_test.data(), LDQ_in_out,
        &RCOND_out, WR_out.data(), WI_out.data(), // WR/WI are 1D, no RM/CM issue
        S_rm_out.data(), LDS_out, U_rm_out.data(), LDU_out,
        1 // row_major = true
    );

    ASSERT_EQ(info_out, INFO_expected);
    EXPECT_NEAR(RCOND_out, RCOND_expected, check_tol);
    for (size_t i = 0; i < X_expected_rowmajor.size(); ++i) {
        // The wrapper ensures Q_rm_test (X) is correctly populated in row-major.
        EXPECT_NEAR(Q_rm_test[i], X_expected_rowmajor[i], check_tol) << "X_rm mismatch at index " << i;
    }
}

TEST_F(SB02MDTestColMajor, ZeroDimensions) {
    int n_zero = 0;
    int lda_z = 1, ldg_z = 1, ldq_z = 1, lds_z = 1, ldu_z = 1;
    double rcond_z;
    
    // For N=0, WR, WI, S, U are not meaningfully sized beyond 1 for WR/WI, 2*N for S/U.
    // The wrapper should handle NULLs or minimal allocations internally.
    // We pass NULLs for arrays that would be zero-sized.
    // Workspace arrays are handled internally.

    int info_out = slicot_sb02md(
        DICO_in, HINV_in, UPLO_in, SCAL_in, SORT_in, n_zero,
        nullptr, lda_z, nullptr, ldg_z, nullptr, ldq_z,
        &rcond_z, nullptr, nullptr, // WR, WI
        nullptr, lds_z, nullptr, ldu_z, // S, U
        0 // row_major = false
    );
    EXPECT_EQ(info_out, 0); // Expect success for N=0
}

TEST_F(SB02MDTestColMajor, ParameterValidation) {
    int n_val = N_in; // Use N_in for valid dimensions where needed
    std::vector<double> a_val( (size_t)LDA_in * n_val );
    std::vector<double> g_val( (size_t)LDG_in * n_val );
    std::vector<double> q_val( (size_t)LDQ_in_out * n_val );
    std::vector<double> s_val( (size_t)LDS_out * (2*n_val) );
    std::vector<double> u_val( (size_t)LDU_out * (2*n_val) );
    std::vector<double> wr_val( 2 * (size_t)n_val );
    std::vector<double> wi_val( 2 * (size_t)n_val );
    double rcond_val;

    // Test invalid DICO
    EXPECT_EQ(slicot_sb02md('X', HINV_in, UPLO_in, SCAL_in, SORT_in, n_val, a_val.data(), LDA_in, g_val.data(), LDG_in, q_val.data(), LDQ_in_out, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), LDS_out, u_val.data(), LDU_out, 0), -1);
    // Test invalid HINV for DICO='D'
    EXPECT_EQ(slicot_sb02md('D', 'X', UPLO_in, SCAL_in, SORT_in, n_val, a_val.data(), LDA_in, g_val.data(), LDG_in, q_val.data(), LDQ_in_out, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), LDS_out, u_val.data(), LDU_out, 0), -2);
    // Test invalid UPLO
    EXPECT_EQ(slicot_sb02md(DICO_in, HINV_in, 'X', SCAL_in, SORT_in, n_val, a_val.data(), LDA_in, g_val.data(), LDG_in, q_val.data(), LDQ_in_out, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), LDS_out, u_val.data(), LDU_out, 0), -3);
    // Test invalid SCAL
    EXPECT_EQ(slicot_sb02md(DICO_in, HINV_in, UPLO_in, 'X', SORT_in, n_val, a_val.data(), LDA_in, g_val.data(), LDG_in, q_val.data(), LDQ_in_out, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), LDS_out, u_val.data(), LDU_out, 0), -4);
    // Test invalid SORT
    EXPECT_EQ(slicot_sb02md(DICO_in, HINV_in, UPLO_in, SCAL_in, 'X', n_val, a_val.data(), LDA_in, g_val.data(), LDG_in, q_val.data(), LDQ_in_out, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), LDS_out, u_val.data(), LDU_out, 0), -5);
    // Test invalid N
    EXPECT_EQ(slicot_sb02md(DICO_in, HINV_in, UPLO_in, SCAL_in, SORT_in, -1, a_val.data(), LDA_in, g_val.data(), LDG_in, q_val.data(), LDQ_in_out, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), LDS_out, u_val.data(), LDU_out, 0), -6);

    // Test invalid LDA (col-major)
    if (n_val > 0) EXPECT_EQ(slicot_sb02md(DICO_in, HINV_in, UPLO_in, SCAL_in, SORT_in, n_val, a_val.data(), 0, g_val.data(), LDG_in, q_val.data(), LDQ_in_out, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), LDS_out, u_val.data(), LDU_out, 0), -8);
    // Test invalid LDG (col-major)
    if (n_val > 0) EXPECT_EQ(slicot_sb02md(DICO_in, HINV_in, UPLO_in, SCAL_in, SORT_in, n_val, a_val.data(), LDA_in, g_val.data(), 0, q_val.data(), LDQ_in_out, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), LDS_out, u_val.data(), LDU_out, 0), -10);
    // Test invalid LDQ (col-major)
    if (n_val > 0) EXPECT_EQ(slicot_sb02md(DICO_in, HINV_in, UPLO_in, SCAL_in, SORT_in, n_val, a_val.data(), LDA_in, g_val.data(), LDG_in, q_val.data(), 0, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), LDS_out, u_val.data(), LDU_out, 0), -12);
    // Test invalid LDS (col-major)
    if (n_val > 0) EXPECT_EQ(slicot_sb02md(DICO_in, HINV_in, UPLO_in, SCAL_in, SORT_in, n_val, a_val.data(), LDA_in, g_val.data(), LDG_in, q_val.data(), LDQ_in_out, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), 0, u_val.data(), LDU_out, 0), -17);
    // Test invalid LDU (col-major)
    if (n_val > 0) EXPECT_EQ(slicot_sb02md(DICO_in, HINV_in, UPLO_in, SCAL_in, SORT_in, n_val, a_val.data(), LDA_in, g_val.data(), LDG_in, q_val.data(), LDQ_in_out, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), LDS_out, u_val.data(), 0, 0), -19);

    // Test NULL pointers (assuming wrapper adds these checks)
    // These are not currently in sb02md.c, so they might pass if not checked, or crash.
    // Example: NULL A
    // if (n_val > 0) EXPECT_EQ(slicot_sb02md(DICO_in, HINV_in, UPLO_in, SCAL_in, SORT_in, n_val, nullptr, LDA_in, g_val.data(), LDG_in, q_val.data(), LDQ_in_out, &rcond_val, wr_val.data(), wi_val.data(), s_val.data(), LDS_out, u_val.data(), LDU_out, 0), -7); // Assuming -7 for A
}

