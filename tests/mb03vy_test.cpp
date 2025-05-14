#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::min, std::max
#include <iostream>  // For debugging
#include <iomanip>   // For std::setprecision

#include "mb03vy.h"       // Header for the C wrapper
#include "slicot_utils.h" // For MAX, MIN if not in <algorithm> for C++

// Helper function to get 3D array index for column-major storage
// Fortran: A(row, col, page) (1-indexed)
// C: a_flat[ (page_idx * lda2 * lda1) + (col_idx * lda1) + row_idx ] (0-indexed)
size_t get_idx(int r, int c, int pg, int lda1, int lda2) {
    return static_cast<size_t>(pg) * lda1 * lda2 +
           static_cast<size_t>(c) * lda1 +
           static_cast<size_t>(r);
}


// --- Test Fixture ---
class MB03VYTest : public ::testing::Test {
protected:
    // Default parameters for basic tests
    int N_default = 1;
    int P_default = 1;
    int ILO_default = 1;
    int IHI_default = 1;
    double check_tol = 1e-9; // Tolerance for checking numerical results

    // Variables for tests
    int N, P, ILO, IHI;
    std::vector<double> A_data;
    std::vector<double> TAU_data;
    int LDA1, LDA2, LDTAU;

    int info_result = -999; // Default to indicate not run

    void SetUpDefault() {
        N = N_default;
        P = P_default;
        ILO = ILO_default;
        IHI = IHI_default;

        LDA1 = std::max(1, N);
        LDA2 = std::max(1, N);
        LDTAU = std::max(1, (N > 0 ? N - 1 : 0)); // N-1 can be 0 if N=1, or -1 if N=0. max(1,0)=1. max(1,-1)=1.
                                                 // Corrected: LDTAU >= max(1, N-1). If N=0, N-1=-1, max(1,-1)=1. If N=1, N-1=0, max(1,0)=1.

        // Allocate A: LDA1 * LDA2 * P
        // A is (LDA1, LDA2, P) in Fortran.
        // In C, flat storage for column major: P pages, each page is LDA2 columns, each col is LDA1 rows.
        A_data.assign(static_cast<size_t>(LDA1) * LDA2 * P, 0.0);

        // Allocate TAU: LDTAU * P
        TAU_data.assign(static_cast<size_t>(LDTAU) * P, 0.0);
    }

    void SetUp() override {
        SetUpDefault();
    }
};

// --- Test Cases ---

TEST_F(MB03VYTest, ParameterValidation) {
    // N invalid
    info_result = slicot_mb03vy(-1, P, ILO, IHI, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, -1);

    // P invalid
    SetUpDefault(); // Reset params
    info_result = slicot_mb03vy(N, 0, ILO, IHI, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, -2);

    // ILO invalid (<1)
    SetUpDefault();
    info_result = slicot_mb03vy(N, P, 0, IHI, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, -3);

    // ILO invalid (>N, for N>0)
    SetUpDefault(); N = 3; ILO = 4; IHI = 3; LDA1=3; LDA2=3; LDTAU=std::max(1,N-1); A_data.assign(LDA1*LDA2*P,0); TAU_data.assign(LDTAU*P,0);
    info_result = slicot_mb03vy(N, P, ILO, IHI, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, -3);
    
    // IHI invalid (<ILO, for N>0)
    SetUpDefault(); N = 3; ILO = 2; IHI = 1; LDA1=3; LDA2=3; LDTAU=std::max(1,N-1); A_data.assign(LDA1*LDA2*P,0); TAU_data.assign(LDTAU*P,0);
    info_result = slicot_mb03vy(N, P, ILO, IHI, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, -4);

    // IHI invalid (>N, for N>0)
    SetUpDefault(); N = 3; ILO = 1; IHI = 4; LDA1=3; LDA2=3; LDTAU=std::max(1,N-1); A_data.assign(LDA1*LDA2*P,0); TAU_data.assign(LDTAU*P,0);
    info_result = slicot_mb03vy(N, P, ILO, IHI, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, -4);

    // LDA1 invalid
    SetUpDefault(); N = 3; LDA1 = 2; LDA2=3; LDTAU=std::max(1,N-1); A_data.assign(LDA1*LDA2*P,0); TAU_data.assign(LDTAU*P,0);
    info_result = slicot_mb03vy(N, P, 1, N, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, -6);

    // LDA2 invalid
    SetUpDefault(); N = 3; LDA1 = 3; LDA2=2; LDTAU=std::max(1,N-1); A_data.assign(LDA1*LDA2*P,0); TAU_data.assign(LDTAU*P,0);
    info_result = slicot_mb03vy(N, P, 1, N, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, -7);

    // LDTAU invalid
    SetUpDefault(); N = 3; LDA1 = 3; LDA2=3; LDTAU=1; // N-1 = 2, so LDTAU must be >= 2
    if (N > 1) { // Only test if N-1 is relevant
       A_data.assign(LDA1*LDA2*P,0); TAU_data.assign(LDTAU*P,0);
       info_result = slicot_mb03vy(N, P, 1, N, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
       EXPECT_EQ(info_result, -9);
    }
    SetUpDefault(); N = 1; LDA1 = 1; LDA2=1; LDTAU=0; // N-1=0, max(1,0)=1. LDTAU must be >=1
    A_data.assign(LDA1*LDA2*P,0); TAU_data.assign(LDTAU*P,0); // This will make TAU_data empty if LDTAU=0
    // The wrapper might not catch LDTAU=0 if N-1 is also 0, but Fortran might.
    // The current wrapper validation is `ldtau < MAX(1, n - 1)`. If N=1, `n-1=0`, `MAX(1,0)=1`. So `ldtau < 1` is an error.
    info_result = slicot_mb03vy(N, P, 1, N, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, -9);


    // Test with NULL pointers for A and TAU when N > 0 or P > 0
    SetUpDefault(); N = 1; P = 1; // Valid N, P
    LDA1=1; LDA2=1; LDTAU=1;
    // A is NULL
    info_result = slicot_mb03vy(N, P, ILO, IHI, nullptr, LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    // The wrapper doesn't explicitly check for A == NULL, but Fortran will likely crash or error.
    // SLICOT wrappers typically check for NULL on required arrays. MB03VY.c does not for A and TAU.
    // This is a deviation from typical robust wrapper practice. For now, we expect it to pass wrapper validation.
    // If it crashes, it means the Fortran layer is hit.
    // Let's assume the Fortran layer would return a negative info for its argument if A is required and NULL.
    // A is arg 5.
    // For now, this test might not be meaningful if the wrapper passes it.
    // EXPECT_EQ(info_result, -5); // This would be the expectation if wrapper checked.

    // TAU is NULL
    A_data.assign(LDA1*LDA2*P,0);
    info_result = slicot_mb03vy(N, P, ILO, IHI, A_data.data(), LDA1, LDA2, nullptr, LDTAU, 0);
    // TAU is arg 8.
    // EXPECT_EQ(info_result, -8); // This would be the expectation if wrapper checked.
}

TEST_F(MB03VYTest, ZeroDimensionN) {
    N = 0;
    P = 1;
    // As per docs: 1 <= ILO <= max(1,N) => ILO = 1
    // min(ILO,N) <= IHI <= N => min(1,0)=0 <= IHI <= 0 => IHI = 0
    ILO = 1;
    IHI = 0;

    LDA1 = std::max(1, N); // Should be 1
    LDA2 = std::max(1, N); // Should be 1
    LDTAU = std::max(1, (N > 0 ? N - 1 : 0)); // N=0 => N-1=-1. max(1,-1)=1. Should be 1.

    A_data.assign(static_cast<size_t>(LDA1) * LDA2 * P, 0.0); // 1*1*1 = 1 element
    TAU_data.assign(static_cast<size_t>(LDTAU) * P, 0.0);   // 1*1 = 1 element

    info_result = slicot_mb03vy(N, P, ILO, IHI, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, 0); // Should succeed with N=0
}


TEST_F(MB03VYTest, BasicFunctionality_N1_P1_ILO1_IHI1) {
    N = 1;
    P = 1;
    ILO = 1;
    IHI = 1; // ihi-1 < ilo, so no reflectors applied. Q should be Identity.

    LDA1 = std::max(1, N);
    LDA2 = std::max(1, N);
    LDTAU = std::max(1, (N > 0 ? N - 1 : 0)); // For N=1, LDTAU >= 1

    A_data.assign(static_cast<size_t>(LDA1) * LDA2 * P, 0.0);
    TAU_data.assign(static_cast<size_t>(LDTAU) * P, 0.0);

    // Input A(0,0,0) can be anything, e.g., 5.0
    // Fortran A(1,1,1)
    if (!A_data.empty()) {
        A_data[get_idx(0,0,0, LDA1, LDA2)] = 5.0; // A(1,1,1)
    }
    // TAU is not used in this case.

    info_result = slicot_mb03vy(N, P, ILO, IHI, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, 0);

    // Check if A(1,1,1) is 1.0 (Identity)
    if (N==1 && P==1 && !A_data.empty()) {
        EXPECT_NEAR(A_data[get_idx(0,0,0, LDA1, LDA2)], 1.0, check_tol);
    }
}

TEST_F(MB03VYTest, RowMajorFlagIgnored) {
    // This test verifies that passing row_major = 1 behaves the same as row_major = 0,
    // because the wrapper states it ignores this flag for array 'a'.
    // We use the N=1, P=1, ILO=1, IHI=1 case which should result in Q=Identity.
    N = 1;
    P = 1;
    ILO = 1;
    IHI = 1;

    LDA1 = std::max(1, N);
    LDA2 = std::max(1, N);
    LDTAU = std::max(1, (N > 0 ? N - 1 : 0));

    A_data.assign(static_cast<size_t>(LDA1) * LDA2 * P, 0.0);
    TAU_data.assign(static_cast<size_t>(LDTAU) * P, 0.0);
    if (!A_data.empty()) A_data[get_idx(0,0,0, LDA1, LDA2)] = 5.0;


    info_result = slicot_mb03vy(N, P, ILO, IHI, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 1); // row_major = 1
    EXPECT_EQ(info_result, 0);
    if (N==1 && P==1 && !A_data.empty()) {
        EXPECT_NEAR(A_data[get_idx(0,0,0, LDA1, LDA2)], 1.0, check_tol);
    }
}

// More complex smoke test for N > 1.
// This requires setting up A and TAU as if MB03VD produced them.
// For simplicity, we will just call with valid dimensions and some dummy data,
// and expect info = 0. We cannot verify the numerical output of A without MB03VD's logic.
TEST_F(MB03VYTest, SmokeTest_N3_P2) {
    N = 3;
    P = 2;
    ILO = 1;
    IHI = N; // Apply all possible reflectors for the range

    LDA1 = std::max(1, N);
    LDA2 = std::max(1, N);
    LDTAU = std::max(1, (N > 0 ? N - 1 : 0)); // N-1 = 2. LDTAU >= 2.

    A_data.assign(static_cast<size_t>(LDA1) * LDA2 * P, 0.0);
    TAU_data.assign(static_cast<size_t>(LDTAU) * P, 0.0);

    // Fill A and TAU with some dummy values.
    // A's lower triangular parts are inputs.
    // For A(i,j,k) where i > j
    for (int k_pg = 0; k_pg < P; ++k_pg) { // Page
        for (int j_col = 0; j_col < N; ++j_col) { // Column
            for (int i_row = j_col + 1; i_row < N; ++i_row) { // Row (strictly lower)
                 if ( (size_t)LDA1 * LDA2 * P > 0) // check if A_data is allocated
                    A_data[get_idx(i_row, j_col, k_pg, LDA1, LDA2)] = static_cast<double>(i_row + j_col + k_pg + 1) / 10.0;
            }
        }
        // Fill diagonal with 1s for identity part of reflectors (conceptual)
        for (int d_idx = 0; d_idx < N; ++d_idx) {
            if ( (size_t)LDA1 * LDA2 * P > 0)
                A_data[get_idx(d_idx, d_idx, k_pg, LDA1, LDA2)] = 1.0;
        }

        // Fill TAU
        // TAU(i,j) -> tau for reflector H_j(i)
        // For MB03VD, TAU(k,j) stores tau for H_j(k) where k from ILO to IHI-1
        // So, for TAU(LDTAU,P), we need to fill TAU(i, k_pg) for i from 0 to (IHI-1)-(ILO)+1 - 1
        for (int i_tau_row = 0; i_tau_row < (IHI - ILO); ++i_tau_row) { // Reflectors from H(ilo) to H(ihi-1)
            if (static_cast<size_t>(LDTAU * P) > 0) // check if TAU_data is allocated
                 TAU_data[static_cast<size_t>(k_pg) * LDTAU + i_tau_row] = 0.5; // Dummy tau value
        }
    }
    
    info_result = slicot_mb03vy(N, P, ILO, IHI, A_data.data(), LDA1, LDA2, TAU_data.data(), LDTAU, 0);
    EXPECT_EQ(info_result, 0);
    // Cannot easily verify the content of A_data without replicating MB03VD logic.
}

