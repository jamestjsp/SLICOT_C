#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <iostream>  // For debugging

#include "sb10jd.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions

// --- Test Fixture Base ---
class SB10JDTest : public ::testing::Test {
protected:
    // From Python example test_sb10jd
    int N_param = 6;
    int M_param = 1;
    int NP_param = 6; // C is 6x6, D is 7x1 in example, but routine takes NP for C rows and D rows.
                      // The python D is zeros((7,1)), so NP should be 7 if D is NPxM.
                      // However, C is eye(6), implying NP=6 for C.
                      // Let's assume NP refers to the rows of C and D.
                      // The Fortran interface implies C is NPxN and D is NPxM.
                      // If C is 6xN and D is 7xM, this is inconsistent.
                      // The Python example has D=zeros((7,1)) but C=eye(6).
                      // Let's use NP_param = 6 for C and adjust D to be 6x1 for consistency with C.
                      // If D needs to be 7x1, then NP_param should be 7, and C would need to be 7xN.
                      // The expected C_r is 6xNSYS, D_r is 6xM. This confirms NP=6.

    double check_tol = 1e-5; // As per Python assert_allclose

    // Input data vectors (will be initialized in SetUp)
    std::vector<double> A_data, B_data, C_data, D_data, E_data;
    
    // Output data (will be modified in place)
    std::vector<double> A_out, B_out, C_out, D_out; // Will copy from _data and pass to func
    int NSYS_out;

    // Expected results (from Python example)
    std::vector<double> A_expected, B_expected, C_expected, D_expected;
    int NSYS_expected; // From example, A_r is 1x1, so NSYS_expected = 1.
    int expected_info = 0;

    // Leading dimensions
    int LDA, LDB, LDC, LDD, LDE;


    void SetUpBase() {
        // Input data from Python example test_sb10jd
        A_data = {
             0,  0,  0, -1,  1,  0,
             0, 32,  0,  0, -1,  1,
             0,  0,  1,  0,  0,  0,
             0,  0,  0,  1,  0,  0,
            -1,  1,  0,  0,  0,  0,
             0, -1,  1,  0,  0,  0
        };
        E_data = {
             0,  0,  0,  0,  0,  0,
             0,  0,  0,  0,  0,  0,
             0,  0,  0,-10,  0, 10,
             0,  0,  0,  0,  0,  0,
             0,  0,  0,  0,  0,  0,
             0,  0,  0,  0,  0,  0
        };
        B_data = {
            -7.1,
             0.0,
             0.0,
             0.0,
             0.0,
             0.0
        };
        // C = eye(6) -> NP_param = 6
        C_data.assign(NP_param * N_param, 0.0);
        for(int i=0; i < std::min(NP_param, N_param); ++i) C_data[i * N_param + i] = 1.0;

        // D = zeros((NP_param, M_param)) -> 6x1
        D_data.assign(NP_param * M_param, 0.0);
        // The python example has D = zeros((7,1)) but C_exp is 6xNSYS.
        // To match C_exp, NP must be 6. So D should be 6x1.
        // The D_exp from python is 6x1.

        // Expected results
        NSYS_expected = 1; // Aexp is 1x1
        A_expected = {-0.003125};
        B_expected = {0.059000};
        C_expected = { // 6x1 (NP_param x NSYS_expected)
            -1.17519e-02,
            -1.17519e-02,
            -1.17519e-02,
             0.0,
             0.0,
             3.76060e-01
        };
        D_expected = { // 6x1 (NP_param x M_param)
             2.21875e-01,
             2.21875e-01,
             2.21875e-01,
             0.0,
             7.10000e+00, // Corrected from 7.100000+00 to 7.1
             0.0
        };

        // Prepare output vectors by copying input data, as they are modified in-place
        A_out = A_data;
        B_out = B_data;
        C_out = C_data;
        D_out = D_data; 
        // E_data is also input/output but its final content is not useful.
        // We still need to pass a modifiable copy.
    }
};

class SB10JDTestColMajor : public SB10JDTest {
protected:
    std::vector<double> A_cm, B_cm, C_cm, D_cm, E_cm; 
    void SetUp() override {
        SetUpBase();
        // Fortran-style Leading Dimensions (number of rows)
        LDA = std::max(1, N_param);    LDB = std::max(1, N_param);
        LDC = std::max(1, NP_param);   LDD = std::max(1, NP_param);
        LDE = std::max(1, N_param);

        // Transpose row-major data to column-major for Fortran call
        if (!A_data.empty()) { A_cm.resize(A_data.size()); slicot_transpose_to_fortran_with_ld(A_data.data(), A_cm.data(), N_param, N_param, N_param, LDA, sizeof(double));}
        if (!B_data.empty()) { B_cm.resize(B_data.size()); slicot_transpose_to_fortran_with_ld(B_data.data(), B_cm.data(), N_param, M_param, M_param, LDB, sizeof(double));}
        if (!C_data.empty()) { C_cm.resize(C_data.size()); slicot_transpose_to_fortran_with_ld(C_data.data(), C_cm.data(), NP_param, N_param, N_param, LDC, sizeof(double));}
        if (!D_data.empty()) { D_cm.resize(D_data.size()); slicot_transpose_to_fortran_with_ld(D_data.data(), D_cm.data(), NP_param, M_param, M_param, LDD, sizeof(double));}
        if (!E_data.empty()) { E_cm.resize(E_data.size()); slicot_transpose_to_fortran_with_ld(E_data.data(), E_cm.data(), N_param, N_param, N_param, LDE, sizeof(double));}
        
        // Output vectors will use these column-major copies
        A_out = A_cm; B_out = B_cm; C_out = C_cm; D_out = D_cm;
    }
};

class SB10JDTestRowMajor : public SB10JDTest {
protected:
    void SetUp() override {
        SetUpBase();
        // C Row-Major Leading Dimensions (number of columns)
        LDA = std::max(1, N_param);    LDB = std::max(1, M_param);
        LDC = std::max(1, N_param);    LDD = std::max(1, M_param);
        LDE = std::max(1, N_param);
        // A_out, B_out etc are already row-major from SetUpBase copies.
    }
};

TEST_F(SB10JDTestColMajor, PythonExampleColMajor) {
    std::vector<double> e_copy = E_cm; // E is modified by the routine

    int info_result = slicot_sb10jd(N_param, M_param, NP_param,
                                   A_out.data(), LDA,
                                   B_out.data(), LDB,
                                   C_out.data(), LDC,
                                   D_out.data(), LDD,
                                   e_copy.data(), LDE, // Pass copy of E
                                   &NSYS_out,
                                   0 /* row_major = false */);
    ASSERT_EQ(info_result, expected_info);
    ASSERT_EQ(NSYS_out, NSYS_expected);

    // A_out, B_out, C_out, D_out are now in column-major (as returned by Fortran)
    // Expected results are currently row-major. Convert them to column-major for comparison.
    std::vector<double> A_exp_cm(A_expected.size()), B_exp_cm(B_expected.size());
    std::vector<double> C_exp_cm(C_expected.size()), D_exp_cm(D_expected.size());

    if(NSYS_expected > 0) {
      slicot_transpose_to_fortran_with_ld(A_expected.data(), A_exp_cm.data(), NSYS_expected, NSYS_expected, NSYS_expected, NSYS_expected, sizeof(double));
      for(size_t i=0; i < (size_t)NSYS_expected * NSYS_expected; ++i) EXPECT_NEAR(A_out[i], A_exp_cm[i], check_tol) << "A_out["<<i<<"]";
    }
    if(NSYS_expected > 0 && M_param > 0) {
      slicot_transpose_to_fortran_with_ld(B_expected.data(), B_exp_cm.data(), NSYS_expected, M_param, M_param, NSYS_expected, sizeof(double));
      for(size_t i=0; i < (size_t)NSYS_expected * M_param; ++i) EXPECT_NEAR(B_out[i], B_exp_cm[i], check_tol) << "B_out["<<i<<"]";
    }
    if(NP_param > 0 && NSYS_expected > 0) {
      slicot_transpose_to_fortran_with_ld(C_expected.data(), C_exp_cm.data(), NP_param, NSYS_expected, NSYS_expected, NP_param, sizeof(double));
      for(size_t i=0; i < (size_t)NP_param * NSYS_expected; ++i) EXPECT_NEAR(C_out[i], C_exp_cm[i], check_tol) << "C_out["<<i<<"]";
    }
    if(NP_param > 0 && M_param > 0) {
      slicot_transpose_to_fortran_with_ld(D_expected.data(), D_exp_cm.data(), NP_param, M_param, M_param, NP_param, sizeof(double));
      for(size_t i=0; i < (size_t)NP_param * M_param; ++i) EXPECT_NEAR(D_out[i], D_exp_cm[i], check_tol) << "D_out["<<i<<"]";
    }
}

TEST_F(SB10JDTestRowMajor, PythonExampleRowMajor) {
    std::vector<double> e_copy = E_data; // E is modified

    int info_result = slicot_sb10jd(N_param, M_param, NP_param,
                                   A_out.data(), LDA,
                                   B_out.data(), LDB,
                                   C_out.data(), LDC,
                                   D_out.data(), LDD,
                                   e_copy.data(), LDE,
                                   &NSYS_out,
                                   1 /* row_major = true */);
    ASSERT_EQ(info_result, expected_info);
    ASSERT_EQ(NSYS_out, NSYS_expected);

    // A_out, B_out, C_out, D_out are now in row-major.
    // Compare directly with expected (which are row-major).
    // Need to consider that only the first NSYSxNSYS (for A), NSYSxM (for B), NPxNSYS (for C)
    // elements of A_out, B_out, C_out are meaningful. D_out is NPxM.

    if(NSYS_expected > 0) {
      for(int i=0; i < NSYS_expected; ++i) {
        for(int j=0; j < NSYS_expected; ++j) {
          EXPECT_NEAR(A_out[i*LDA + j], A_expected[i*NSYS_expected+j], check_tol) << "A_out["<<i<<"]["<<j<<"]";
        }
      }
    }
     if(NSYS_expected > 0 && M_param > 0) {
      for(int i=0; i < NSYS_expected; ++i) {
        for(int j=0; j < M_param; ++j) {
          EXPECT_NEAR(B_out[i*LDB + j], B_expected[i*M_param+j], check_tol) << "B_out["<<i<<"]["<<j<<"]";
        }
      }
    }
    if(NP_param > 0 && NSYS_expected > 0) {
      for(int i=0; i < NP_param; ++i) {
        for(int j=0; j < NSYS_expected; ++j) {
          EXPECT_NEAR(C_out[i*LDC + j], C_expected[i*NSYS_expected+j], check_tol) << "C_out["<<i<<"]["<<j<<"]";
        }
      }
    }
    if(NP_param > 0 && M_param > 0) {
       for(int i=0; i < NP_param; ++i) {
        for(int j=0; j < M_param; ++j) {
          EXPECT_NEAR(D_out[i*LDD + j], D_expected[i*M_param+j], check_tol) << "D_out["<<i<<"]["<<j<<"]";
        }
      }
    }
}

TEST_F(SB10JDTestColMajor, ParameterValidation) {
    int nsys_val;
    std::vector<double> e_dummy_cm;
    if (N_param > 0) e_dummy_cm.resize(N_param*N_param);

    int info_result = slicot_sb10jd(-1, M_param, NP_param, nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1,&nsys_val,0);
    EXPECT_EQ(info_result, -1); // Invalid N

    info_result = slicot_sb10jd(N_param, M_param, NP_param, A_cm.data(),0,B_cm.data(),1,C_cm.data(),1,D_cm.data(),1,e_dummy_cm.data(),1,&nsys_val,0);
    if (N_param > 0) EXPECT_EQ(info_result, -5); // Invalid LDA
    else EXPECT_EQ(info_result, 0); // Or other if N=0 allows LDA=0
}

TEST_F(SB10JDTestColMajor, ZeroDimensionN) {
    int n_zero = 0;
    int m_zd = 1, np_zd = 1; 
    int nsys_zd_out;

    // For N=0, A,B,C,E are effectively empty. D is NPxM.
    // Ad, Bd, Cd become 0x0, 0xM, NPx0. NSYS should be 0.
    // Dd = D.
    std::vector<double> d_zd_data = {1.23}; // 1x1 D matrix
    std::vector<double> d_zd_out = d_zd_data; // D is modified in place to Dd

    int ldd_zd = std::max(1, np_zd);

    int info_result = slicot_sb10jd(n_zero, m_zd, np_zd,
                                   nullptr, 1,       // A
                                   nullptr, 1,       // B
                                   nullptr, std::max(1,np_zd), // C
                                   d_zd_out.data(), ldd_zd, // D
                                   nullptr, 1,       // E
                                   &nsys_zd_out,
                                   0 /* row_major = false */);
    EXPECT_EQ(info_result, 0);
    EXPECT_EQ(nsys_zd_out, 0); // Expect system order to become 0
    if (!d_zd_out.empty() && !d_zd_data.empty()) {
        EXPECT_NEAR(d_zd_out[0], d_zd_data[0], check_tol); // Dd should be D if N=0
    }
}

