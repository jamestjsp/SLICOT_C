#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max

#include "sb02od.h"
#include "slicot_utils.h"

// --- Test Fixture ---
class SB02ODTest : public ::testing::Test {
protected:
    // Default parameters for a simple continuous-time case
    char DICO_in = 'C';
    char JOBB_in = 'B';
    char FACT_in = 'N';
    char UPLO_in = 'U';
    char JOBL_in = 'Z';
    char SORT_in = 'S';
    int N_in = 2;
    int M_in = 1;
    int P_in = 1; // Not strictly used if FACT='N' but good to define
    double TOL_in = 0.0; // Use default tolerance

    // Input matrices (will be initialized in SetUp or specific tests)
    std::vector<double> A_io;
    std::vector<double> B_io;
    std::vector<double> Q_io;
    std::vector<double> R_io;
    std::vector<double> L_io;

    // Output matrices & scalars
    double RCOND_out = 0.0;
    std::vector<double> X_out;
    std::vector<double> ALFAR_out;
    std::vector<double> ALFAI_out;
    std::vector<double> BETA_out;
    std::vector<double> S_out;
    std::vector<double> T_out;
    std::vector<double> U_out;

    // Leading dimensions (will be set based on N, M, P, JOBB, DICO)
    int LDA_in, LDB_in, LDQ_in, LDR_in, LDL_in;
    int LDX_out, LDS_out, LDT_out, LDU_out;

    // Expected results (for specific test cases)
    std::vector<double> X_expected;
    double RCOND_expected_min = 0.0; // Minimum expected RCOND
    double RCOND_expected_val_n0 = 1.0; // Expected RCOND for N=0
    int INFO_expected = 0;

    // Test verification tolerance
    double check_tol = 1e-5;


    void InitializeDimensions(int n, int m, int p, char jobb, char dico, bool row_major_layout) {
        N_in = n;
        M_in = m;
        P_in = p;
        JOBB_in = jobb;
        DICO_in = dico;

        int s_t_dim = (JOBB_in == 'B') ? (2 * N_in + M_in) : (2 * N_in);
        int t_rows_fortran = (JOBB_in == 'B') ? (2*N_in+M_in) : ((DICO_in == 'D') ? (2*N_in) : 1);
        int t_cols_fortran = (DICO_in == 'C' && JOBB_in == 'G') ? 1 : (2*N_in);


        if (row_major_layout) {
            LDA_in = std::max(1, N_in);
            LDB_in = std::max(1, (JOBB_in == 'B' ? M_in : N_in));
            LDQ_in = std::max(1, N_in); // Assuming Q(N,N) or C(P,N) -> N cols for Q/C
            LDR_in = std::max(1, (JOBB_in == 'B' ? M_in : 1)); // R(M,M) or D(P,M) -> M cols for R/D
            LDL_in = std::max(1, (JOBB_in == 'B' && JOBL_in == 'N' ? M_in : 1));
            LDX_out = std::max(1, N_in);
            LDS_out = std::max(1, s_t_dim); // S is s_t_dim x s_t_dim
            LDT_out = std::max(1, t_cols_fortran); // T is t_rows_fortran x t_cols_fortran
            LDU_out = std::max(1, 2 * N_in); // U is 2N x 2N
        } else { // Column-major
            LDA_in = std::max(1, N_in);
            LDB_in = std::max(1, N_in); // B(N,M) or G(N,N)
            LDQ_in = std::max(1, (FACT_in == 'N' || FACT_in == 'D' ? N_in : P_in));
            LDR_in = std::max(1, (JOBB_in == 'B' ? (FACT_in == 'N' || FACT_in == 'C' ? M_in : P_in) : 1));
            LDL_in = std::max(1, (JOBB_in == 'B' && JOBL_in == 'N' ? N_in : 1));
            LDX_out = std::max(1, N_in);
            LDS_out = std::max(1, s_t_dim);
            LDT_out = std::max(1, t_rows_fortran);
            LDU_out = std::max(1, 2 * N_in);
        }

        A_io.resize((size_t)LDA_in * (row_major_layout ? N_in : N_in)); // N rows for CM, N cols for RM
        if (JOBB_in == 'B') B_io.resize((size_t)LDB_in * (row_major_layout ? N_in : M_in));
        else B_io.resize((size_t)LDB_in * (row_major_layout ? N_in : N_in)); // G
        
        if (FACT_in == 'N' || FACT_in == 'D') Q_io.resize((size_t)LDQ_in * (row_major_layout ? N_in : N_in)); // Q
        else Q_io.resize((size_t)LDQ_in * (row_major_layout ? P_in : N_in)); // C

        if (JOBB_in == 'B') {
            if (FACT_in == 'N' || FACT_in == 'C') R_io.resize((size_t)LDR_in * (row_major_layout ? M_in : M_in)); // R
            else R_io.resize((size_t)LDR_in * (row_major_layout ? P_in : M_in)); // D
            if (JOBL_in == 'N') L_io.resize((size_t)LDL_in * (row_major_layout ? N_in : M_in));
            else L_io.clear();
        } else {
            R_io.clear();
            L_io.clear();
        }
        
        X_out.resize((size_t)LDX_out * N_in);
        ALFAR_out.resize(2 * std::max(1,N_in));
        ALFAI_out.resize(2 * std::max(1,N_in));
        BETA_out.resize(2 * std::max(1,N_in));
        S_out.resize((size_t)LDS_out * s_t_dim); // S is s_t_dim x s_t_dim
        if (!(DICO_in == 'C' && JOBB_in == 'G')) {
             T_out.resize((size_t)LDT_out * (row_major_layout ? t_rows_fortran : t_cols_fortran));
        } else {
            T_out.clear();
        }
        U_out.resize((size_t)LDU_out * (2 * N_in)); // U is 2N x 2N
    }

    void SetUpSimpleContinuousExample() {
        // Example from SB02OD.html (slightly modified for N=2, M=1)
        // A = [0, 1; 0, 0], B = [0; 1], Q = [1, 0; 0, 1], R = [1], L = [0;0]
        // Expected X for this system (approx): [sqrt(3), 1; 1, sqrt(3)]
        DICO_in = 'C'; JOBB_in = 'B'; FACT_in = 'N'; UPLO_in = 'U'; JOBL_in = 'Z'; SORT_in = 'S';
        InitializeDimensions(2, 1, 1, JOBB_in, DICO_in, false); // col-major for setup

        A_io = {0.0, 0.0, 1.0, 0.0}; // Fortran: A(1,1), A(2,1), A(1,2), A(2,2)
        B_io = {0.0, 1.0};           // Fortran: B(1,1), B(2,1)
        Q_io = {1.0, 0.0, 0.0, 1.0}; // Fortran: Q(1,1), Q(2,1), Q(1,2), Q(2,2) (upper)
                                     // For UPLO='U', Q(2,1) is not strictly needed by Fortran but good for consistency
        R_io = {1.0};                // Fortran: R(1,1) (upper)
        L_io.assign((size_t)LDL_in * M_in, 0.0); // Zero L

        // X_expected is full, column-major
        X_expected = {1.73205081, 1.0, 1.0, 1.73205081}; // sqrt(3), 1, 1, sqrt(3)
        RCOND_expected_min = 1e-2; // Expect reasonably well-conditioned
        INFO_expected = 0;
    }
};

// Test: Basic Continuous-Time System (Column-Major)
TEST_F(SB02ODTest, ContinuousARE_ColMajor_Basic) {
    SetUpSimpleContinuousExample(); // Sets up for col-major

    int info = slicot_sb02od(DICO_in, JOBB_in, FACT_in, UPLO_in, JOBL_in, SORT_in,
                             N_in, M_in, P_in,
                             A_io.data(), LDA_in, B_io.data(), LDB_in,
                             Q_io.data(), LDQ_in, R_io.data(), LDR_in,
                             (JOBL_in == 'N' ? L_io.data() : nullptr), LDL_in,
                             &RCOND_out, X_out.data(), LDX_out,
                             ALFAR_out.data(), ALFAI_out.data(), BETA_out.data(),
                             S_out.data(), LDS_out, T_out.data(), LDT_out,
                             U_out.data(), LDU_out, TOL_in, 0 /* col-major */);

    ASSERT_EQ(info, INFO_expected);
    EXPECT_GT(RCOND_out, RCOND_expected_min);

    // Check X (symmetric, check upper triangle)
    // X_expected is in column major format for this setup
    EXPECT_NEAR(X_out[0], X_expected[0], check_tol); // X(0,0)
    EXPECT_NEAR(X_out[1], X_expected[1], check_tol); // X(1,0)
    // X_out[2] is X(0,1)
    EXPECT_NEAR(X_out[LDX_out * 1 + 0], X_expected[2], check_tol); // X(0,1)
    EXPECT_NEAR(X_out[LDX_out * 1 + 1], X_expected[3], check_tol); // X(1,1)

    // Verify symmetry of X_out if UPLO_in was 'U'
    if (UPLO_in == 'U' && N_in > 1) {
         EXPECT_NEAR(X_out[LDX_out * 0 + 1], X_out[LDX_out*1 + 0], 1e-9); // X(1,0) vs X(0,1)
    }
}

// Test: Basic Continuous-Time System (Row-Major)
TEST_F(SB02ODTest, ContinuousARE_RowMajor_Basic) {
    SetUpSimpleContinuousExample(); // Sets up col-major data initially

    // Transpose inputs to Row-Major
    std::vector<double> A_rm(A_io.size());
    std::vector<double> B_rm(B_io.size());
    std::vector<double> Q_rm(Q_io.size());
    std::vector<double> R_rm(R_io.size());
    std::vector<double> L_rm;
    if (JOBL_in == 'N') L_rm.resize(L_io.size());

    slicot_transpose_to_c_with_ld(A_io.data(), A_rm.data(), N_in, N_in, LDA_in, N_in, sizeof(double));
    slicot_transpose_to_c_with_ld(B_io.data(), B_rm.data(), N_in, M_in, LDB_in, M_in, sizeof(double));
    // Q is symmetric, for UPLO='U', only upper is significant in input Q_io.
    // For row-major, we'd pass the full Q_rm.
    // The wrapper handles taking the UPLO part for Fortran.
    // So, create full Q_rm from Q_io (which is upper CM)
    std::vector<double> Q_full_cm = Q_io; // Q_io is already full for this example
    // if (N_in > 1 && UPLO_in == 'U') Q_full_cm[0*LDQ_in + 1] = Q_full_cm[1*LDQ_in+0]; // Q(1,0) = Q(0,1) if Q_io was only upper
    // else if (N_in > 1 && UPLO_in == 'L') Q_full_cm[1*LDQ_in+0] = Q_full_cm[0*LDQ_in+1]; // Q(0,1) = Q(1,0) if Q_io was only lower
    slicot_transpose_to_c_with_ld(Q_full_cm.data(), Q_rm.data(), N_in, N_in, std::max(1,N_in), N_in, sizeof(double));

    // R is symmetric
    std::vector<double> R_full_cm = R_io; // R_io is already full for this example
    // if M_in > 1 ... (similar for R if not scalar)
    slicot_transpose_to_c_with_ld(R_full_cm.data(), R_rm.data(), M_in, M_in, std::max(1,M_in), M_in, sizeof(double));


    if (JOBL_in == 'N') {
      slicot_transpose_to_c_with_ld(L_io.data(), L_rm.data(), N_in, M_in, LDL_in, M_in, sizeof(double));
    }

    // Update leading dimensions for row-major call
    InitializeDimensions(N_in, M_in, P_in, JOBB_in, DICO_in, true); // true for row-major

    int info = slicot_sb02od(DICO_in, JOBB_in, FACT_in, UPLO_in, JOBL_in, SORT_in,
                             N_in, M_in, P_in,
                             A_rm.data(), LDA_in, B_rm.data(), LDB_in,
                             Q_rm.data(), LDQ_in, R_rm.data(), LDR_in,
                             (JOBL_in == 'N' ? L_rm.data() : nullptr), LDL_in,
                             &RCOND_out, X_out.data(), LDX_out,
                             ALFAR_out.data(), ALFAI_out.data(), BETA_out.data(),
                             S_out.data(), LDS_out, 
                             (DICO_in == 'C' && JOBB_in == 'G' ? nullptr : T_out.data()), LDT_out,
                             U_out.data(), LDU_out, TOL_in, 1 /* row-major */);

    ASSERT_EQ(info, INFO_expected);
    EXPECT_GT(RCOND_out, RCOND_expected_min);

    // X_expected is in column-major. Transpose it to row-major for comparison.
    std::vector<double> X_expected_rm(X_expected.size());
    slicot_transpose_to_c_with_ld(X_expected.data(), X_expected_rm.data(), N_in, N_in, std::max(1,N_in), std::max(1,N_in), sizeof(double));

    // Check X_out (row-major) - should be full symmetric
    for (int i = 0; i < N_in; ++i) {
        for (int j = 0; j < N_in; ++j) {
            EXPECT_NEAR(X_out[i * LDX_out + j], X_expected_rm[i * N_in + j], check_tol)
                << "X_rm mismatch at (" << i << "," << j << ")";
        }
    }
     // Verify symmetry of X_out explicitly
    if (N_in > 1) {
         EXPECT_NEAR(X_out[1 * LDX_out + 0], X_out[0 * LDX_out + 1], check_tol); // X_rm(1,0) vs X_rm(0,1)
    }
}


// Test: Zero Dimensions
TEST_F(SB02ODTest, ZeroDimensions) {
    InitializeDimensions(0, 0, 0, 'B', 'C', false); // N=0, M=0, P=0
    double dummy = 0.0; // For arrays that cannot be nullptr

    // For N=0, many arrays can be nullptr.
    // RCOND, ALFAR, ALFAI, BETA, S, T, U need valid pointers even if not used.
    // X also needs a valid pointer.
    // The wrapper should handle NULLs for A,B,Q,R,L if N=0/M=0.
    // However, the Fortran routine might not like NULL for X, S, T, U, ALFAR etc.
    // Let's provide dummy pointers for outputs.

    ALFAR_out.assign(2, 0.0); // Min size for 2*N
    ALFAI_out.assign(2, 0.0);
    BETA_out.assign(2, 0.0);
    X_out.assign(1,0.0); LDX_out = 1;
    S_out.assign(1,0.0); LDS_out = 1; // S_T_DIM will be 0 if M=0,N=0. Max(1,0) = 1.
    T_out.assign(1,0.0); LDT_out = 1; // T_ROWS_F will be 1 if N=0,M=0,DICO=C,JOBB=B. Max(1,1)=1.
    U_out.assign(1,0.0); LDU_out = 1; // 2*N = 0. Max(1,0)=1.


    int info = slicot_sb02od(DICO_in, JOBB_in, FACT_in, UPLO_in, JOBL_in, SORT_in,
                             N_in, M_in, P_in,
                             nullptr, 1, nullptr, 1, // A, B
                             nullptr, 1, nullptr, 1, // Q, R
                             nullptr, 1,              // L
                             &RCOND_out, X_out.data(), LDX_out,
                             ALFAR_out.data(), ALFAI_out.data(), BETA_out.data(),
                             S_out.data(), LDS_out, T_out.data(), LDT_out,
                             U_out.data(), LDU_out, TOL_in, 0 /* col-major */);
     EXPECT_EQ(info, 0);
     EXPECT_EQ(RCOND_out, RCOND_expected_val_n0); 
}
 
 // Test: Parameter Validation (Selected)
TEST_F(SB02ODTest, ParameterValidation) {
    SetUpSimpleContinuousExample(); // Initialize with valid N, M, P etc.

    // Test invalid DICO
    int info = slicot_sb02od('X', JOBB_in, FACT_in, UPLO_in, JOBL_in, SORT_in, N_in, M_in, P_in, A_io.data(), LDA_in, B_io.data(), LDB_in, Q_io.data(), LDQ_in, R_io.data(), LDR_in, L_io.data(), LDL_in, &RCOND_out, X_out.data(), LDX_out, ALFAR_out.data(), ALFAI_out.data(), BETA_out.data(), S_out.data(), LDS_out, T_out.data(), LDT_out, U_out.data(), LDU_out, TOL_in, 0);
    EXPECT_EQ(info, -1);

    // Test invalid N
    info = slicot_sb02od(DICO_in, JOBB_in, FACT_in, UPLO_in, JOBL_in, SORT_in, -1, M_in, P_in, A_io.data(), LDA_in, B_io.data(), LDB_in, Q_io.data(), LDQ_in, R_io.data(), LDR_in, L_io.data(), LDL_in, &RCOND_out, X_out.data(), LDX_out, ALFAR_out.data(), ALFAI_out.data(), BETA_out.data(), S_out.data(), LDS_out, T_out.data(), LDT_out, U_out.data(), LDU_out, TOL_in, 0);
    EXPECT_EQ(info, -7);
    
    // Test invalid LDA
    if (N_in > 0) {
        info = slicot_sb02od(DICO_in, JOBB_in, FACT_in, UPLO_in, JOBL_in, SORT_in, N_in, M_in, P_in, A_io.data(), 0, B_io.data(), LDB_in, Q_io.data(), LDQ_in, R_io.data(), LDR_in, L_io.data(), LDL_in, &RCOND_out, X_out.data(), LDX_out, ALFAR_out.data(), ALFAI_out.data(), BETA_out.data(), S_out.data(), LDS_out, T_out.data(), LDT_out, U_out.data(), LDU_out, TOL_in, 0);
        EXPECT_EQ(info, -11);
    }

    // Test NULL A when N > 0
    if (N_in > 0) {
        info = slicot_sb02od(DICO_in, JOBB_in, FACT_in, UPLO_in, JOBL_in, SORT_in, N_in, M_in, P_in, nullptr, LDA_in, B_io.data(), LDB_in, Q_io.data(), LDQ_in, R_io.data(), LDR_in, L_io.data(), LDL_in, &RCOND_out, X_out.data(), LDX_out, ALFAR_out.data(), ALFAI_out.data(), BETA_out.data(), S_out.data(), LDS_out, T_out.data(), LDT_out, U_out.data(), LDU_out, TOL_in, 0);
        EXPECT_EQ(info, -10); // Expect error for NULL A (arg 10)
     }
 }
 
// TODO: Add more tests:
// - JOBB='G'
// - FACT='C', 'D', 'B'
// - UPLO='L'
// - JOBL='N'
// - SORT='U'
// - DICO='D'
// - Cases leading to INFO > 0 (e.g., singular pencil, QZ fail, reorder fail, etc. - might be hard to construct)

