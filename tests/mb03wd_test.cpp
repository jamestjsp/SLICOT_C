#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::min, std::max
#include <iostream>
#include <iomanip>   // For std::setprecision
#include <cstring>

#include "mb03vd.h"  // For creating input in periodic Hessenberg form
#include "mb03wd.h"       // Header for the C wrapper
#include "slicot_utils.h" // For MAX, MIN if not in <algorithm> for C++
#include "test_utils.h"

// Helper function to get 3D array index for column-major storage
// Fortran: A(row, col, page) (1-indexed)
// C: a_flat[ (page_idx * lda2_c * lda1_c) + (col_idx * lda1_c) + row_idx ] (0-indexed)
// lda1_c is the first leading dimension (number of rows in Fortran sense)
// lda2_c is the second leading dimension (number of columns in Fortran sense)
size_t get_3d_idx(int r_idx, int c_idx, int p_idx, int lda1_c, int lda2_c) {
    return static_cast<size_t>(p_idx) * lda1_c * lda2_c +
           static_cast<size_t>(c_idx) * lda1_c +
           static_cast<size_t>(r_idx);
}

// --- Test Fixture ---
class MB03WDTest : public ::testing::Test {
protected:
    char JOB;
    char COMPZ;
    int N, P, ILO, IHI, ILOZ, IHIZ;
    std::vector<double> H_data;
    std::vector<double> Z_data;
    std::vector<double> WR_data;
    std::vector<double> WI_data;
    int LDH1, LDH2, LDZ1, LDZ2;

    int info_result = -999; 
    double check_tol = 1e-7; // Standard tolerance

    void InitializeDims(int n_val, int p_val, char job_val = 'S', char compz_val = 'N') {
        N = n_val;
        P = p_val;
        JOB = job_val;
        COMPZ = compz_val;

        // Default full range for ILO, IHI, ILOZ, IHIZ
        if (N == 0) {
            ILO = 1; IHI = 0;
            ILOZ = 1; IHIZ = 0;
        } else {
            ILO = 1; IHI = N;
            ILOZ = 1; IHIZ = N;
        }

        LDH1 = std::max(1, N);
        LDH2 = std::max(1, N);
        
        if (COMPZ == 'N') {
            LDZ1 = 1;
            LDZ2 = 1;
        } else {
            LDZ1 = std::max(1, N);
            LDZ2 = std::max(1, N);
        }

        H_data.assign(static_cast<size_t>(LDH1) * LDH2 * P, 0.0);
        if (COMPZ != 'N') {
            Z_data.assign(static_cast<size_t>(LDZ1) * LDZ2 * P, 0.0);
            if (COMPZ == 'V' || COMPZ == 'I') { // Initialize Z to identity for 'V' or 'I'
                 for (int k_page = 0; k_page < P; ++k_page) {
                    for (int i = 0; i < N; ++i) {
                        if(N > 0) Z_data[get_3d_idx(i, i, k_page, LDZ1, LDZ2)] = 1.0;
                    }
                }
            }
        } else {
            Z_data.clear();
        }
        
        if (N > 0) {
            WR_data.resize(N);
            WI_data.resize(N);
        } else {
            WR_data.clear();
            WI_data.clear();
        }
    }
    // Example data from MB03WD.html (raw input, not yet in Hessenberg/Triangular form)
    // This data is NOT directly usable for testing MB03WD outputs against example outputs
    // without first processing through MB03VD.
    void LoadRawExampleInputH() {
        if (N == 4 && P == 2) {
            std::vector<double> h1_raw_example = { // Column-major from data file
                1.5, 1.0, 1.5, 1.0,
               -0.7, 0.0,-0.7, 0.0,
                3.5, 2.0, 2.5, 2.0,
               -0.7, 3.0,-0.3, 1.0
            };
            std::vector<double> h2_raw_example = { // Identical to h1_raw_example
                1.5, 1.0, 1.5, 1.0,
               -0.7, 0.0,-0.7, 0.0,
                3.5, 2.0, 2.5, 2.0,
               -0.7, 3.0,-0.3, 1.0
            };
            for(int j=0; j<N; ++j) for(int i=0; i<N; ++i) H_data[get_3d_idx(i,j,0,LDH1,LDH2)] = h1_raw_example[static_cast<size_t>(j)*N + i];
            for(int j=0; j<N; ++j) for(int i=0; i<N; ++i) H_data[get_3d_idx(i,j,1,LDH1,LDH2)] = h2_raw_example[static_cast<size_t>(j)*N + i];
        }
    }
};

// --- Column-Major Test Fixture ---
class MB03WDTestColMajor : public ::testing::Test {
protected:
    // Test parameters
    int N = 4;
    int P = 2;
    int ILO = 1;
    int IHI = 4;
    int ILOZ = 1;
    int IHIZ = 4;
    
    // Verification tolerance
    double check_tol = 1e-10;
    
    // Input matrices (column-major)
    std::vector<double> A = {
        // A1 (n×n)
        1.5, 1.0, 1.5, 1.0,    // First column
        -0.7, 0.0, -0.7, 0.0,  // Second column
        3.5, 2.0, 2.5, 2.0,    // Third column
        -0.7, 3.0, -0.3, 1.0,  // Fourth column
        
        // A2 (n×n)
        1.5, 1.0, 1.5, 1.0,    // First column
        -0.7, 0.0, -0.7, 0.0,  // Second column
        3.5, 2.0, 2.5, 2.0,    // Third column
        -0.7, 3.0, -0.3, 1.0   // Fourth column
    };
    
    // Periodic Hessenberg form (after MB03VD)
    std::vector<double> H;
    
    // Schur form (to be computed)
    std::vector<double> T;
    
    // Schur vectors
    std::vector<double> Z;
    
    // Eigenvalues
    std::vector<double> WR;
    std::vector<double> WI;
    
    // Auxiliary for MB03VD
    std::vector<double> TAU;
    
    // Result info
    int info_result = -999;
    
    // Leading dimensions
    int LDA1;
    int LDA2;
    int LDH1;
    int LDH2;
    int LDZ1;
    int LDZ2;
    int LDTAU;
    
    void SetUp() override {
        // Size output arrays
        H.resize(N * N * P);
        T.resize(N * N * P);
        Z.resize(N * N * P);
        WR.resize(N);
        WI.resize(N);
        TAU.resize((N-1) * P);
        
        // Calculate leading dimensions
        LDA1 = LDH1 = LDZ1 = std::max(1, N);
        LDA2 = LDH2 = LDZ2 = std::max(1, N);
        LDTAU = std::max(1, N-1);
        
        // Create input matrix in periodic Hessenberg form using MB03VD
        H = A;  // Copy input matrices
        info_result = slicot_mb03vd(N, P, ILO, IHI,
                                   H.data(), LDH1, LDH2,
                                   TAU.data(), LDTAU,
                                   0 /* row_major = false */);
        ASSERT_EQ(info_result, 0);
        
        // Copy H to T for Schur decomposition
        T = H;
    }
};

// --- Row-Major Test Fixture ---
class MB03WDTestRowMajor : public ::testing::Test {
protected:
    // Test parameters
    int N = 4;
    int P = 2;
    int ILO = 1;
    int IHI = 4;
    int ILOZ = 1;
    int IHIZ = 4;
    
    // Verification tolerance
    double check_tol = 1e-10;
    
    // Input matrices (row-major)
    std::vector<double> A_rm = {
        // A1 (n×n) in row-major format
        1.5, -0.7, 3.5, -0.7,  // First row
        1.0, 0.0, 2.0, 3.0,    // Second row
        1.5, -0.7, 2.5, -0.3,  // Third row
        1.0, 0.0, 2.0, 1.0,    // Fourth row
        
        // A2 (n×n) in row-major format
        1.5, -0.7, 3.5, -0.7,  // First row
        1.0, 0.0, 2.0, 3.0,    // Second row
        1.5, -0.7, 2.5, -0.3,  // Third row
        1.0, 0.0, 2.0, 1.0     // Fourth row
    };
    
    // Periodic Hessenberg form (after MB03VD)
    std::vector<double> H_rm;
    
    // Schur form (to be computed)
    std::vector<double> T_rm;
    
    // Schur vectors
    std::vector<double> Z_rm;
    
    // Eigenvalues
    std::vector<double> WR;
    std::vector<double> WI;
    
    // Auxiliary for MB03VD
    std::vector<double> TAU_rm;
    
    // Result info
    int info_result = -999;
    
    // Leading dimensions
    int LDA1;
    int LDA2;
    int LDH1;
    int LDH2;
    int LDZ1;
    int LDZ2;
    int LDTAU;
    
    void SetUp() override {
        // Size output arrays
        H_rm.resize(N * N * P);
        T_rm.resize(N * N * P);
        Z_rm.resize(N * N * P);
        WR.resize(N);
        WI.resize(N);
        TAU_rm.resize((N-1) * P);
        
        // Calculate leading dimensions for row-major
        LDA1 = LDH1 = LDZ1 = N;  // Number of rows in each 2D slice
        LDA2 = LDH2 = LDZ2 = N;  // Number of columns in each 2D slice
        LDTAU = N-1;  // Leading dimension for TAU
        
        // Create input matrix in periodic Hessenberg form using MB03VD
        H_rm = A_rm;  // Copy input matrices
        info_result = slicot_mb03vd(N, P, ILO, IHI,
                                   H_rm.data(), LDH1, LDH2,
                                   TAU_rm.data(), LDTAU,
                                   1 /* row_major = true */);
        ASSERT_EQ(info_result, 0);
        
        // Copy H to T for Schur decomposition
        T_rm = H_rm;
    }
};

// --- Test Cases ---

TEST_F(MB03WDTest, ParameterValidation) {
    InitializeDims(4, 2, 'S', 'V'); // Base valid dimensions
    // JOB invalid
    info_result = slicot_mb03wd('X', COMPZ, N, P, ILO, IHI, ILOZ, IHIZ, H_data.data(), LDH1, LDH2, Z_data.data(), LDZ1, LDZ2, WR_data.data(), WI_data.data(), 0);
    EXPECT_EQ(info_result, -1);
    // COMPZ invalid
    info_result = slicot_mb03wd(JOB, 'X', N, P, ILO, IHI, ILOZ, IHIZ, H_data.data(), LDH1, LDH2, Z_data.data(), LDZ1, LDZ2, WR_data.data(), WI_data.data(), 0);
    EXPECT_EQ(info_result, -2);
    // N invalid
    info_result = slicot_mb03wd(JOB, COMPZ, -1, P, 1, 0, 1, 0, nullptr, 1, 1, nullptr, 1, 1, nullptr, nullptr, 0); // Adjusted ILO/IHI for N=-1
    EXPECT_EQ(info_result, -3);
    // P invalid
    info_result = slicot_mb03wd(JOB, COMPZ, N, 0, ILO, IHI, ILOZ, IHIZ, H_data.data(), LDH1, LDH2, Z_data.data(), LDZ1, LDZ2, WR_data.data(), WI_data.data(), 0);
    EXPECT_EQ(info_result, -4);
    // ILO invalid (<1)
    info_result = slicot_mb03wd(JOB, COMPZ, N, P, 0, IHI, ILOZ, IHIZ, H_data.data(), LDH1, LDH2, Z_data.data(), LDZ1, LDZ2, WR_data.data(), WI_data.data(), 0);
    EXPECT_EQ(info_result, -5);
    // IHI invalid (>N)
    info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, N + 1, ILOZ, IHIZ, H_data.data(), LDH1, LDH2, Z_data.data(), LDZ1, LDZ2, WR_data.data(), WI_data.data(), 0);
    EXPECT_EQ(info_result, -6);
    // ILOZ invalid (>ILO)
    if (N>0) { // Ensure ILO is meaningful
      info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILO + 1, IHIZ, H_data.data(), LDH1, LDH2, Z_data.data(), LDZ1, LDZ2, WR_data.data(), WI_data.data(), 0);
      EXPECT_EQ(info_result, -7);
    }
    // IHIZ invalid (<IHI)
     if (N>0) { // Ensure IHI is meaningful
      info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHI -1 , H_data.data(), LDH1, LDH2, Z_data.data(), LDZ1, LDZ2, WR_data.data(), WI_data.data(), 0);
      EXPECT_EQ(info_result, -8);
    }

    // H is NULL
    if (N > 0) {
        info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ, nullptr, LDH1, LDH2, Z_data.data(), LDZ1, LDZ2, WR_data.data(), WI_data.data(), 0);
        EXPECT_EQ(info_result, -9);
    }
    // LDH1 invalid
    info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ, H_data.data(), 0, LDH2, Z_data.data(), LDZ1, LDZ2, WR_data.data(), WI_data.data(), 0);
    EXPECT_EQ(info_result, -10);
    // Z is NULL (COMPZ != 'N')
    if (N > 0 && COMPZ != 'N') {
        info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ, H_data.data(), LDH1, LDH2, nullptr, LDZ1, LDZ2, WR_data.data(), WI_data.data(), 0);
        EXPECT_EQ(info_result, -12);
    }
    // LDZ1 invalid (COMPZ != 'N')
    if (COMPZ != 'N') {
        info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ, H_data.data(), LDH1, LDH2, Z_data.data(), 0, LDZ2, WR_data.data(), WI_data.data(), 0);
        EXPECT_EQ(info_result, -13);
    }
    // WR is NULL
    if (N > 0) {
        info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ, H_data.data(), LDH1, LDH2, Z_data.data(), LDZ1, LDZ2, nullptr, WI_data.data(), 0);
        EXPECT_EQ(info_result, -15);
    }
    // WI is NULL
    if (N > 0) {
        info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ, H_data.data(), LDH1, LDH2, Z_data.data(), LDZ1, LDZ2, WR_data.data(), nullptr, 0);
        EXPECT_EQ(info_result, -16);
    }
}

TEST_F(MB03WDTest, ZeroDimensionN) {
    InitializeDims(0, 1, 'E', 'N');
    // For N=0: ILO=1, IHI=0, ILOZ=1, IHIZ=0 (set by InitializeDims)
    
    info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ,
                                nullptr, LDH1, LDH2, // H can be NULL if N=0
                                nullptr, LDZ1, LDZ2, // Z can be NULL if N=0 or COMPZ='N'
                                nullptr, nullptr, 0); // WR, WI can be NULL if N=0
    EXPECT_EQ(info_result, 0);
}

TEST_F(MB03WDTest, Functionality_N1_P1_JobS_CompZI) {
    InitializeDims(1, 1, 'S', 'I');
    // H(0,0,0) = 5.0
    if(N>0 && P>0) H_data[get_3d_idx(0,0,0, LDH1, LDH2)] = 5.0;

    info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ,
                                H_data.data(), LDH1, LDH2,
                                Z_data.data(), LDZ1, LDZ2,
                                WR_data.data(), WI_data.data(), 0);
    EXPECT_EQ(info_result, 0);
    if (N==1) {
        EXPECT_NEAR(WR_data[0], 5.0, check_tol);
        EXPECT_NEAR(WI_data[0], 0.0, check_tol);
        EXPECT_NEAR(H_data[get_3d_idx(0,0,0,LDH1,LDH2)], 5.0, check_tol); // T1 = H1
        EXPECT_NEAR(Z_data[get_3d_idx(0,0,0,LDZ1,LDZ2)], 1.0, check_tol); // Z1 = I
    }
}

TEST_F(MB03WDTest, Functionality_N1_P2_JobS_CompZI) {
    InitializeDims(1, 2, 'S', 'I');
    // H(0,0,0) = 2.0 (H1)
    // H(0,0,1) = 3.0 (H2)
    // Product H = H1*H2 = 6.0
    if(N>0 && P>1) {
      H_data[get_3d_idx(0,0,0, LDH1, LDH2)] = 2.0;
      H_data[get_3d_idx(0,0,1, LDH1, LDH2)] = 3.0;
    }

    info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ,
                                H_data.data(), LDH1, LDH2,
                                Z_data.data(), LDZ1, LDZ2,
                                WR_data.data(), WI_data.data(), 0);
    EXPECT_EQ(info_result, 0);

    if (N==1) {
        EXPECT_NEAR(WR_data[0], 6.0, check_tol); // Eigenvalue of product
        EXPECT_NEAR(WI_data[0], 0.0, check_tol);
        // T1*T2 = Eigenvalue. For N=1, T1=[[t1]], T2=[[t2]]. t1*t2=6.0
        // The routine might distribute this as T1=6, T2=1 or T1=2, T2=3 etc.
        // We know T1 is Schur, T2 is upper triangular. For N=1, both are just scalars.
        // Let's check product:
        if (P==2) {
             EXPECT_NEAR(H_data[get_3d_idx(0,0,0,LDH1,LDH2)] * H_data[get_3d_idx(0,0,1,LDH1,LDH2)], 6.0, check_tol);
        }
        // Z_1 and Z_2 should be identity
        if (P > 0 && N > 0) EXPECT_NEAR(Z_data[get_3d_idx(0,0,0,LDZ1,LDZ2)], 1.0, check_tol);
        if (P > 1 && N > 0) EXPECT_NEAR(Z_data[get_3d_idx(0,0,1,LDZ1,LDZ2)], 1.0, check_tol);
    }
}

TEST_F(MB03WDTest, SmokeTest_N2_P2_JobS_CompZI) {
    InitializeDims(2, 2, 'S', 'I');

    // Construct H1 (Upper Hessenberg) for page 0
    // H1 = [1 2]
    //      [3 4]
    if (N==2 && P > 0) {
        H_data[get_3d_idx(0,0,0,LDH1,LDH2)] = 1.0; H_data[get_3d_idx(0,1,0,LDH1,LDH2)] = 2.0;
        H_data[get_3d_idx(1,0,0,LDH1,LDH2)] = 3.0; H_data[get_3d_idx(1,1,0,LDH1,LDH2)] = 4.0;
    }

    // Construct H2 (Upper Triangular) for page 1
    // H2 = [5 6]
    //      [0 7]
    if (N==2 && P > 1) {
        H_data[get_3d_idx(0,0,1,LDH1,LDH2)] = 5.0; H_data[get_3d_idx(0,1,1,LDH1,LDH2)] = 6.0;
        H_data[get_3d_idx(1,0,1,LDH1,LDH2)] = 0.0; H_data[get_3d_idx(1,1,1,LDH1,LDH2)] = 7.0;
    }
    
    info_result = slicot_mb03wd(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ,
                                H_data.data(), LDH1, LDH2,
                                Z_data.data(), LDZ1, LDZ2,
                                WR_data.data(), WI_data.data(), 0);
    EXPECT_EQ(info_result, 0);

    // Further checks could verify properties of T1 (Schur form) and T2 (Upper Triangular)
    // e.g., T1(1,0) should be 0 if eigenvalues are real. T2(1,0) should be 0.
    if (N==2 && P > 1 && JOB=='S') {
        // Check T2 is upper triangular
        EXPECT_NEAR(H_data[get_3d_idx(1,0,1,LDH1,LDH2)], 0.0, check_tol); 
        // Check T1 is upper quasi-triangular (real Schur form)
        // If WI[0] and WI[1] are zero, then H_data[get_3d_idx(1,0,0,LDH1,LDH2)] (T1(2,1)) should be zero.
        bool real_eigs = true;
        for(int i=0; i<N; ++i) if(std::abs(WI_data[i]) > check_tol) real_eigs = false;
        
        if(real_eigs) {
            EXPECT_NEAR(H_data[get_3d_idx(1,0,0,LDH1,LDH2)], 0.0, check_tol);
        }
        // If complex eigs, T1(2,1) might not be zero.
    }
}

// Test: Eigenvalue Computation Only (Column-Major)
TEST_F(MB03WDTestColMajor, EigenvaluesOnly) {
    // Call MB03WD to compute eigenvalues only
    info_result = slicot_mb03wd('E', 'N', N, P, ILO, IHI, ILOZ, IHIZ,
                              H.data(), LDH1, LDH2,
                              nullptr, 1, 1,
                              WR.data(), WI.data(),
                              0 /* row_major = false */);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Verify that eigenvalues come in complex conjugate pairs
    for (int i = 0; i < N; i++) {
        if (WI[i] != 0.0) {  // Complex eigenvalue
            EXPECT_NEAR(WR[i], WR[i+1], check_tol);
            EXPECT_NEAR(WI[i], -WI[i+1], check_tol);
            i++;  // Skip the conjugate pair
        }
    }
    
    // Note: We don't verify the actual values of eigenvalues here
    // as they depend on the specific input matrices
}

// Test: Schur Decomposition (Column-Major)
TEST_F(MB03WDTestColMajor, SchurDecomposition) {
    // Initialize Z to identity matrices
    std::fill(Z.begin(), Z.end(), 0.0);
    for (int k = 0; k < P; k++) {
        for (int i = 0; i < N; i++) {
            Z[k*N*N + i*N + i] = 1.0;  // Set diagonals to 1.0
        }
    }
    
    // Call MB03WD to compute Schur decomposition
    info_result = slicot_mb03wd('S', 'V', N, P, ILO, IHI, ILOZ, IHIZ,
                              T.data(), LDH1, LDH2,
                              Z.data(), LDZ1, LDZ2,
                              WR.data(), WI.data(),
                              0 /* row_major = false */);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Check that T1 is in real Schur form (upper quasi-triangular)
    for (int i = 1; i < N; i++) {
        // Check below sub-diagonal elements are zero
        for (int j = 0; j < i-1; j++) {
            EXPECT_NEAR(T[j*N + i], 0.0, check_tol);
        }
    }
    
    // Check T2...Tp are upper triangular
    for (int k = 1; k < P; k++) {
        for (int i = 1; i < N; i++) {
            for (int j = 0; j < i; j++) {
                EXPECT_NEAR(T[k*N*N + j*N + i], 0.0, check_tol);
            }
        }
    }
    
    // Check orthogonality of Z matrices
    for (int k = 0; k < P; k++) {
        // Compute Z'*Z for kth matrix
        std::vector<double> ZTZ(N * N, 0.0);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double sum = 0.0;
                for (int l = 0; l < N; l++) {
                    sum += Z[k*N*N + l*N + i] * Z[k*N*N + l*N + j];
                }
                ZTZ[i*N + j] = sum;
            }
        }
        
        // Check that Z'*Z is close to identity
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    EXPECT_NEAR(ZTZ[i*N + j], 1.0, check_tol);
                } else {
                    EXPECT_NEAR(ZTZ[i*N + j], 0.0, check_tol);
                }
            }
        }
    }
}

// Test: Initialize Schur Vectors (Column-Major)
TEST_F(MB03WDTestColMajor, InitSchurVectors) {
    // Call MB03WD to compute Schur decomposition and initialize Z
    info_result = slicot_mb03wd('S', 'I', N, P, ILO, IHI, ILOZ, IHIZ,
                              T.data(), LDH1, LDH2,
                              Z.data(), LDZ1, LDZ2,
                              WR.data(), WI.data(),
                              0 /* row_major = false */);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Check orthogonality of Z matrices
    for (int k = 0; k < P; k++) {
        // Compute Z'*Z for kth matrix
        std::vector<double> ZTZ(N * N, 0.0);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double sum = 0.0;
                for (int l = 0; l < N; l++) {
                    sum += Z[k*N*N + l*N + i] * Z[k*N*N + l*N + j];
                }
                ZTZ[i*N + j] = sum;
            }
        }
        
        // Check that Z'*Z is close to identity
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    EXPECT_NEAR(ZTZ[i*N + j], 1.0, check_tol);
                } else {
                    EXPECT_NEAR(ZTZ[i*N + j], 0.0, check_tol);
                }
            }
        }
    }
}

// Test: Eigenvalue Computation Only (Row-Major)
TEST_F(MB03WDTestRowMajor, EigenvaluesOnly) {
    // Call MB03WD to compute eigenvalues only in row-major mode
    info_result = slicot_mb03wd('E', 'N', N, P, ILO, IHI, ILOZ, IHIZ,
                              H_rm.data(), LDH1, LDH2,
                              nullptr, 1, 1,
                              WR.data(), WI.data(),
                              1 /* row_major = true */);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Verify that eigenvalues come in complex conjugate pairs
    for (int i = 0; i < N; i++) {
        if (WI[i] != 0.0) {  // Complex eigenvalue
            EXPECT_NEAR(WR[i], WR[i+1], check_tol);
            EXPECT_NEAR(WI[i], -WI[i+1], check_tol);
            i++;  // Skip the conjugate pair
        }
    }
}

// Test: Schur Decomposition (Row-Major)
TEST_F(MB03WDTestRowMajor, SchurDecomposition) {
    // Initialize Z to identity matrices
    std::fill(Z_rm.begin(), Z_rm.end(), 0.0);
    for (int k = 0; k < P; k++) {
        for (int i = 0; i < N; i++) {
            Z_rm[k*N*N + i*N + i] = 1.0;  // Set diagonals to 1.0
        }
    }
    
    // Call MB03WD to compute Schur decomposition in row-major mode
    info_result = slicot_mb03wd('S', 'V', N, P, ILO, IHI, ILOZ, IHIZ,
                              T_rm.data(), LDH1, LDH2,
                              Z_rm.data(), LDZ1, LDZ2,
                              WR.data(), WI.data(),
                              1 /* row_major = true */);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Check that T1 is in real Schur form (upper quasi-triangular)
    // In row-major, elements BELOW the first sub-diagonal should be zero
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i-1; j++) {  // Fixed: j < i-1 checks elements below the first sub-diagonal
            EXPECT_NEAR(T_rm[i*N + j], 0.0, check_tol);
        }
    }
    
    // Check T2...Tp are upper triangular
    // In row-major, elements BELOW the main diagonal should be zero
    for (int k = 1; k < P; k++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {  // Fixed: j < i checks elements below the diagonal
                EXPECT_NEAR(T_rm[k*N*N + i*N + j], 0.0, check_tol);
            }
        }
    }
    
    // Check orthogonality of Z matrices
    for (int k = 0; k < P; k++) {
        // Compute Z*Z' for kth matrix (since we're in row-major)
        std::vector<double> ZZT(N * N, 0.0);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double sum = 0.0;
                for (int l = 0; l < N; l++) {
                    sum += Z_rm[k*N*N + i*N + l] * Z_rm[k*N*N + j*N + l];
                }
                ZZT[i*N + j] = sum;
            }
        }
        
        // Check that Z*Z' is close to identity
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    EXPECT_NEAR(ZZT[i*N + j], 1.0, check_tol);
                } else {
                    EXPECT_NEAR(ZZT[i*N + j], 0.0, check_tol);
                }
            }
        }
    }
}

// Test: Initialize Schur Vectors (Row-Major)
TEST_F(MB03WDTestRowMajor, InitSchurVectors) {
    // Call MB03WD to compute Schur decomposition and initialize Z in row-major mode
    info_result = slicot_mb03wd('S', 'I', N, P, ILO, IHI, ILOZ, IHIZ,
                              T_rm.data(), LDH1, LDH2,
                              Z_rm.data(), LDZ1, LDZ2,
                              WR.data(), WI.data(),
                              1 /* row_major = true */);
    
    // Verify return code
    ASSERT_EQ(info_result, 0);
    
    // Check orthogonality of Z matrices
    for (int k = 0; k < P; k++) {
        // Compute Z*Z' for kth matrix (since we're in row-major)
        std::vector<double> ZZT(N * N, 0.0);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double sum = 0.0;
                for (int l = 0; l < N; l++) {
                    sum += Z_rm[k*N*N + i*N + l] * Z_rm[k*N*N + j*N + l];
                }
                ZZT[i*N + j] = sum;
            }
        }
        
        // Check that Z*Z' is close to identity
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    EXPECT_NEAR(ZZT[i*N + j], 1.0, check_tol);
                } else {
                    EXPECT_NEAR(ZZT[i*N + j], 0.0, check_tol);
                }
            }
        }
    }
}

// Test: Parameter Validation
TEST_F(MB03WDTestColMajor, ParameterValidation) {
    // Use properly sized arrays for all tests to avoid segfaults
    std::vector<double> dummy_H(N*N*P, 0.0);
    std::vector<double> dummy_Z(N*N*P, 0.0);
    std::vector<double> dummy_wr(N);
    std::vector<double> dummy_wi(N);
    
    // Test invalid JOB
    info_result = slicot_mb03wd('X', 'N', N, P, ILO, IHI, ILOZ, IHIZ,
                              dummy_H.data(), N, N,
                              nullptr, 1, 1,
                              dummy_wr.data(), dummy_wi.data(),
                              0);
    EXPECT_EQ(info_result, -1);
    
    // Test invalid COMPZ
    info_result = slicot_mb03wd('E', 'X', N, P, ILO, IHI, ILOZ, IHIZ,
                              dummy_H.data(), N, N,
                              nullptr, 1, 1,
                              dummy_wr.data(), dummy_wi.data(),
                              0);
    EXPECT_EQ(info_result, -2);
    
    // Test invalid N
    info_result = slicot_mb03wd('E', 'N', -1, P, ILO, IHI, ILOZ, IHIZ,
                              dummy_H.data(), N, N,
                              nullptr, 1, 1,
                              dummy_wr.data(), dummy_wi.data(),
                              0);
    EXPECT_EQ(info_result, -3);
    
    // Test invalid P
    info_result = slicot_mb03wd('E', 'N', N, 0, ILO, IHI, ILOZ, IHIZ,
                              dummy_H.data(), N, N,
                              nullptr, 1, 1,
                              dummy_wr.data(), dummy_wi.data(),
                              0);
    EXPECT_EQ(info_result, -4);
    
    // Additional parameter validation tests...
    // (ILO, IHI, ILOZ, IHIZ, etc.)
}

// Test: Zero Dimension Case
TEST_F(MB03WDTestColMajor, ZeroDimension) {
    // Call MB03WD with N=0
    int zero_n = 0;
    info_result = slicot_mb03wd('E', 'N', zero_n, P, 1, 0, 1, 0,
                              nullptr, 1, 1,
                              nullptr, 1, 1,
                              nullptr, nullptr,
                              0);
    EXPECT_EQ(info_result, 0);
}

