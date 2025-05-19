#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // For std::max, std::min
#include <numeric>   // For std::accumulate
#include <cctype>    // For toupper C++ style
#include <iostream>  // For debug prints

#include "tg01fd.h"       // Wrapper header
#include "slicot_utils.h" // For transpose functions

// --- Base Test Fixture for TG01FD ---
class Tg01fdBaseTest : public ::testing::Test {
protected:
    char COMPQ_val, COMPZ_val, JOBA_val;
    int L_val, N_val, M_val, P_val;
    double TOL_val;

    // Input matrices (initially in Fortran col-major flat for convenience from test data)
    std::vector<double> A_init_f, E_init_f, B_init_f, C_init_f;
    std::vector<double> Q_init_f, Z_init_f; // For COMPQ/Z = 'U'

    // Actual C-layout matrices passed to the wrapper (will be modified)
    std::vector<double> A_c, E_c, B_c, C_c, Q_c, Z_c;
    
    // Expected outputs (Fortran column-major for matrices for comparison)
    std::vector<double> A_expected_f, E_expected_f, B_expected_f, C_expected_f;
    std::vector<double> Q_expected_f, Z_expected_f;
    int RANKE_expected, RNKA22_expected;
    int INFO_expected;

    // Actual outputs from wrapper
    int RANKE_out, RNKA22_out;
    int INFO_out;

    int LDA_c, LDE_c, LDB_c, LDC_c, LDQ_c, LDZ_c;
    double check_tol = 1e-7; 

    virtual void SetUpTestData() = 0;

    // Helper to prepare a C-layout matrix from a Fortran-layout flat vector
    void PrepareMatrixForCall(const std::vector<double>& source_f_flat, 
                              std::vector<double>& target_c_layout, 
                              int rows, int cols, bool is_row_major_dest) {
        if (rows == 0 || cols == 0) {
            target_c_layout.clear();
            return;
        }
        target_c_layout.resize((size_t)rows * cols);
        if (source_f_flat.empty() && (rows*cols > 0) ) { // Source is empty but target shouldn't be
             // This might happen if an optional input matrix for 'U' mode is not provided
             // but space still needs to be passed. Fill with a known pattern or zeros.
            std::fill(target_c_layout.begin(), target_c_layout.end(), 0.0); // Example: fill with zeros
            return;
        }
        if (source_f_flat.size() != (size_t)rows * cols) {
            // This indicates an issue with test data setup if source_f_flat is not empty
            // For now, if sizes mismatch but target needs to be non-empty, fill with 0.
            if (!target_c_layout.empty()) std::fill(target_c_layout.begin(), target_c_layout.end(), 0.0);
            return;
        }

        if (is_row_major_dest) {
            slicot_transpose_to_c_with_ld(source_f_flat.data(), target_c_layout.data(), 
                                          rows, cols, rows, cols, sizeof(double));
        } else { // Destination is column-major C
            target_c_layout = source_f_flat;
        }
    }


    void SetUpBase(bool is_row_major) {
        SetUpTestData(); 

        // Prepare C-layout matrices that will be passed to the wrapper
        PrepareMatrixForCall(A_init_f, A_c, L_val, N_val, is_row_major);
        PrepareMatrixForCall(E_init_f, E_c, L_val, N_val, is_row_major);
        if (M_val > 0) PrepareMatrixForCall(B_init_f, B_c, L_val, M_val, is_row_major); else B_c.clear();
        if (P_val > 0) PrepareMatrixForCall(C_init_f, C_c, P_val, N_val, is_row_major); else C_c.clear();
        
        if (COMPQ_val == 'U' && L_val > 0) PrepareMatrixForCall(Q_init_f, Q_c, L_val, L_val, is_row_major);
        else if (COMPQ_val == 'I' && L_val > 0) Q_c.assign((size_t)L_val * L_val, 0.0); // Allocate for output
        else Q_c.clear();

        if (COMPZ_val == 'U' && N_val > 0) PrepareMatrixForCall(Z_init_f, Z_c, N_val, N_val, is_row_major);
        else if (COMPZ_val == 'I' && N_val > 0) Z_c.assign((size_t)N_val * N_val, 0.0); // Allocate for output
        else Z_c.clear();
        
        // Calculate C leading dimensions
        if (is_row_major) {
            LDA_c = std::max(1, N_val); LDE_c = std::max(1, N_val);
            LDB_c = (M_val > 0) ? std::max(1, M_val) : 1;
            LDC_c = (P_val > 0) ? std::max(1, N_val) : 1;
            LDQ_c = (L_val > 0) ? std::max(1, L_val) : 1; 
            LDZ_c = (N_val > 0) ? std::max(1, N_val) : 1; 
        } else { // Column-major C
            LDA_c = std::max(1, L_val); LDE_c = std::max(1, L_val);
            LDB_c = (M_val > 0) ? std::max(1, L_val) : 1;
            LDC_c = (P_val > 0) ? std::max(1, P_val) : 1;
            LDQ_c = (L_val > 0) ? std::max(1, L_val) : 1;
            LDZ_c = (N_val > 0) ? std::max(1, N_val) : 1;
        }
    }

    void VerifyMatrix(const std::string& name, const std::vector<double>& actual_c_data, 
                      const std::vector<double>& expected_f_data, 
                      int rows, int cols, int ld_c_val, bool is_row_major_input) {
        if (rows == 0 || cols == 0) {
            EXPECT_TRUE(expected_f_data.empty()) << name << " expected to be empty.";
            return;
        }
        ASSERT_EQ(expected_f_data.size(), (size_t)rows * cols) << name << " expected size mismatch.";
        
        // actual_c_data should have been prepared by SetUpBase or the test itself
        if (actual_c_data.size() != (size_t)rows * cols) {
             FAIL() << name << " actual_c_data size mismatch. Expected: " << (size_t)rows*cols << " Got: " << actual_c_data.size();
        }


        std::vector<double> actual_f(expected_f_data.size());
        if (is_row_major_input) {
            slicot_transpose_to_fortran_with_ld(actual_c_data.data(), actual_f.data(), 
                                                rows, cols, ld_c_val, std::max(1, rows), sizeof(double));
        } else {
            if (ld_c_val == rows && !actual_c_data.empty()) { actual_f = actual_c_data; } 
            else if (!actual_c_data.empty()) { 
                 for (int j = 0; j < cols; ++j) {
                    for (int i = 0; i < rows; ++i) {
                        actual_f[i + j * rows] = actual_c_data[i + j * ld_c_val];
                    }
                }
            } else { // actual_c_data is empty but rows/cols > 0
                 FAIL() << name << " actual_c_data is empty but rows/cols > 0";
            }
        }
        for (size_t i = 0; i < expected_f_data.size(); ++i) {
            EXPECT_NEAR(actual_f[i], expected_f_data[i], check_tol)
                << name << " matrix mismatch at flat index " << i 
                << ", actual_f: " << actual_f[i] << ", expected_f: " << expected_f_data[i];
        }
    }
    
    void VerifyOutputs(bool is_row_major_input) {
        ASSERT_EQ(INFO_out, INFO_expected);
        if (INFO_expected != 0) return;

        VerifyMatrix("A", A_c, A_expected_f, L_val, N_val, LDA_c, is_row_major_input);
        VerifyMatrix("E", E_c, E_expected_f, L_val, N_val, LDE_c, is_row_major_input);
        
        if (M_val > 0) VerifyMatrix("B", B_c, B_expected_f, L_val, M_val, LDB_c, is_row_major_input);
        else EXPECT_TRUE(B_expected_f.empty());
        
        if (P_val > 0) VerifyMatrix("C", C_c, C_expected_f, P_val, N_val, LDC_c, is_row_major_input);
        else EXPECT_TRUE(C_expected_f.empty());

        if (COMPQ_val != 'N' && L_val > 0) VerifyMatrix("Q", Q_c, Q_expected_f, L_val, L_val, LDQ_c, is_row_major_input);
        else if (COMPQ_val != 'N') EXPECT_TRUE(Q_expected_f.empty());


        if (COMPZ_val != 'N' && N_val > 0) VerifyMatrix("Z", Z_c, Z_expected_f, N_val, N_val, LDZ_c, is_row_major_input);
        else if (COMPZ_val != 'N') EXPECT_TRUE(Z_expected_f.empty());

        EXPECT_EQ(RANKE_out, RANKE_expected);
        if (JOBA_val != 'N') {
            EXPECT_EQ(RNKA22_out, RNKA22_expected);
        }
    }
};

// --- Test Fixture for TG01FD Documentation Example (from Python test1_tg01fd) ---
class Tg01fdDocExampleTest : public Tg01fdBaseTest {
protected:
    void SetUpTestData() override {
        L_val = 4; N_val = 4; M_val = 2; P_val = 2;
        COMPQ_val = 'I'; COMPZ_val = 'I'; JOBA_val = 'T'; 
        TOL_val = 0.0;
        INFO_expected = 0;

        // Input data (Fortran column-major layout)
        A_init_f = {-1.0, 0.0, 1.0, 0.0,  0.0, 0.0, 1.0, 0.0,  0.0, 1.0, 0.0, 0.0,  3.0, 2.0, 4.0, 0.0};
        E_init_f = {1.0, 0.0, 3.0, 0.0,  2.0, 1.0, 9.0, 0.0,  0.0, 0.0, 6.0, 2.0,  0.0, 1.0, 3.0, 0.0};
        B_init_f = {1.0, 0.0, 0.0, 1.0,  0.0, 0.0, 1.0, 1.0};
        C_init_f = {-1.0, 0.0, 0.0, 1.0, 1.0, -1.0, 0.0, 1.0};
        Q_init_f.clear(); Z_init_f.clear(); // Not used for COMPQ/Z = 'I'

        // Expected outputs (Fortran column-major) from Python test1_tg01fd
        A_expected_f = {2.02781052, -0.09804588, 0.27131089, 0.06900656, 
                        0.10783277, 0.25437761, 0.77603837, -0.56694671,
                        3.90616686, 1.60529591, -0.36920735, -2.19740106,
                       -2.15710472, -0.12692683, -0.48533567, 0.3086067};
        E_expected_f = {10.15874008, 0.0, 0.0, 0.0, 
                        5.82296975, -2.468405, 0.0, 0.0,
                        1.30205562, -0.18960188, 1.03378058, 0.0,
                        0.0, 0.0, 0.0, 0.0};
        B_expected_f = {-0.21566555, 0.30148458, 0.75952691, 1.13389342,
                        -0.97049496, 0.95156071, 0.09906873, 0.37796447};
        C_expected_f = {0.36514837, -1.09544512, -1.00000000,  1.00000000, 
                        -0.44721360, -0.89442719, -0.81649658,  2.22044605e-16};
        Q_expected_f = {-0.21566555, -0.10783277, -0.97049496,  0.0,
                        -0.50875523, -0.25437761,  0.1413209,   0.81023981,
                         0.61092382, -0.77603837, -0.04953436,  0.14860309,
                         0.56694671,  0.56694671, -0.18898224,  0.56694671};
        Z_expected_f = {-0.36514837, -0.91287093,  6.19714937e-17, -0.18257419,
                        -1.35772740e-16, 0.0, -1.0, -6.78863700e-17,
                         0.44721360,  0.0,  0.0, -0.89442719,
                         0.81649658, -0.40824829, -1.38572473e-16,  0.40824829};

        RANKE_expected = 3;
        RNKA22_expected = 1;
    }
};

TEST_F(Tg01fdDocExampleTest, ColMajor_I_I_T) {
    SetUpBase(false); // Column-major C setup
    
    INFO_out = slicot_tg01fd(COMPQ_val, COMPZ_val, JOBA_val, 
                             L_val, N_val, M_val, P_val,
                             A_c.data(), LDA_c, E_c.data(), LDE_c,
                             B_c.data(), LDB_c, C_c.data(), LDC_c,
                             Q_c.data(), LDQ_c, Z_c.data(), LDZ_c,
                             &RANKE_out, &RNKA22_out, TOL_val, 0);
    VerifyOutputs(false);
}

TEST_F(Tg01fdDocExampleTest, RowMajor_I_I_T) {
    SetUpBase(true); // Row-major C setup

    INFO_out = slicot_tg01fd(COMPQ_val, COMPZ_val, JOBA_val, 
                             L_val, N_val, M_val, P_val,
                             A_c.data(), LDA_c, E_c.data(), LDE_c,
                             B_c.data(), LDB_c, C_c.data(), LDC_c,
                             Q_c.data(), LDQ_c, Z_c.data(), LDZ_c,
                             &RANKE_out, &RNKA22_out, TOL_val, 1);
    VerifyOutputs(true);
}


// --- Test Fixture for Parameter Validation ---
class Tg01fdParamValidationTest : public Tg01fdBaseTest {
protected:
    void SetUpTestData() override {
        L_val = 1; N_val = 1; M_val = 0; P_val = 0; TOL_val = 0.0;
        COMPQ_val = 'N'; COMPZ_val = 'N'; JOBA_val = 'N'; 
        A_init_f = {1.0}; E_init_f = {1.0}; 
        B_init_f.clear(); C_init_f.clear(); Q_init_f.clear(); Z_init_f.clear();
        A_expected_f.clear(); E_expected_f.clear(); B_expected_f.clear(); C_expected_f.clear();
        Q_expected_f.clear(); Z_expected_f.clear();
        RANKE_expected = 0; RNKA22_expected = 0; 
    }
};

TEST_F(Tg01fdParamValidationTest, InvalidCOMPQ) {
    INFO_expected = -1; SetUpBase(false); COMPQ_val = 'X';
    INFO_out = slicot_tg01fd(COMPQ_val,COMPZ_val,JOBA_val,L_val,N_val,M_val,P_val,A_c.data(),LDA_c,E_c.data(),LDE_c,B_c.data(),LDB_c,C_c.data(),LDC_c,Q_c.data(),LDQ_c,Z_c.data(),LDZ_c,&RANKE_out,&RNKA22_out,TOL_val,0);
    ASSERT_EQ(INFO_out, INFO_expected);
}
// ... Add more validation tests ...


// --- Test Fixture for Zero Dimensions ---
class Tg01fdZeroDimTest : public Tg01fdBaseTest {
protected:
    void SetUpTestData() override { 
        COMPQ_val = 'N'; COMPZ_val = 'N'; JOBA_val = 'N'; 
        TOL_val = 0.0; INFO_expected = 0;
        RANKE_expected = 0; RNKA22_expected = 0;
        A_init_f.clear(); E_init_f.clear(); B_init_f.clear(); C_init_f.clear();
        Q_init_f.clear(); Z_init_f.clear();
        A_expected_f.clear(); E_expected_f.clear(); B_expected_f.clear(); C_expected_f.clear();
        Q_expected_f.clear(); Z_expected_f.clear();
    }
};

TEST_F(Tg01fdZeroDimTest, AllZeros) {
    L_val = 0; N_val = 0; M_val = 0; P_val = 0;
    SetUpBase(false); 
    INFO_out = slicot_tg01fd(COMPQ_val,COMPZ_val,JOBA_val,L_val,N_val,M_val,P_val,nullptr,LDA_c,nullptr,LDE_c,nullptr,LDB_c,nullptr,LDC_c,nullptr,LDQ_c,nullptr,LDZ_c,&RANKE_out,&RNKA22_out,TOL_val,0);
    VerifyOutputs(false); 
}

