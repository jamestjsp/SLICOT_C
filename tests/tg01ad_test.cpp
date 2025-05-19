#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // For std::max
#include <numeric>   // For std::accumulate
#include <cctype>    // For toupper C++ style

#include "tg01ad.h"       // Wrapper header
#include "slicot_utils.h" // For transpose functions

// --- Base Test Fixture for TG01AD ---
class Tg01adBaseTest : public ::testing::Test {
protected:
    // Parameters for slicot_tg01ad
    char JOB_val;
    int L_val, N_val, M_val, P_val;
    double THRESH_val;

    std::vector<double> A_c_layout;
    std::vector<double> E_c_layout;
    std::vector<double> B_c_layout;
    std::vector<double> C_c_layout;
    
    // Expected outputs (Fortran column-major for comparison convenience for matrices)
    std::vector<double> A_expected_f;
    std::vector<double> E_expected_f;
    std::vector<double> B_expected_f; // if M > 0
    std::vector<double> C_expected_f; // if P > 0
    std::vector<double> LSCALE_expected;
    std::vector<double> RSCALE_expected;
    int INFO_expected;

    // Actual outputs from wrapper (matrices will be in C layout)
    // A, E, B, C are modified in-place, so we'll copy them before calling the wrapper for verification.
    std::vector<double> A_actual_c; 
    std::vector<double> E_actual_c;
    std::vector<double> B_actual_c;
    std::vector<double> C_actual_c;
    std::vector<double> LSCALE_out;
    std::vector<double> RSCALE_out;
    int INFO_out;

    // Leading dimensions for C arrays
    int LDA_c, LDE_c, LDB_c, LDC_c;

    double check_tol = 1e-4; 

    virtual void SetUpTestData() = 0;

    void SetUpBase(bool is_row_major) {
        SetUpTestData(); 

        // Size output scale vectors
        if (L_val > 0) LSCALE_out.assign(L_val, 0.0); else LSCALE_out.clear();
        if (N_val > 0) RSCALE_out.assign(N_val, 0.0); else RSCALE_out.clear();

        // Copy input layouts to actual_c versions which will be modified by the wrapper
        A_actual_c = A_c_layout;
        E_actual_c = E_c_layout;
        B_actual_c = B_c_layout;
        C_actual_c = C_c_layout;
        
        // Calculate C leading dimensions based on row_major flag
        if (is_row_major) {
            LDA_c = std::max(1, N_val); LDE_c = std::max(1, N_val);
            LDB_c = (M_val > 0) ? std::max(1, M_val) : 1;
            LDC_c = (P_val > 0) ? std::max(1, N_val) : 1;
        } else { // Column-major C
            LDA_c = std::max(1, L_val); LDE_c = std::max(1, L_val);
            LDB_c = (M_val > 0) ? std::max(1, L_val) : 1;
            LDC_c = (P_val > 0) ? std::max(1, P_val) : 1;
        }
    }

    void VerifyMatrix(const std::string& name, const std::vector<double>& actual_c, 
                      const std::vector<double>& expected_f, 
                      int rows, int cols, int ld_c, bool is_row_major_input) {
        if (rows == 0 || cols == 0) {
            EXPECT_TRUE(expected_f.empty());
            // actual_c might be empty or minimally sized, which is fine if logically empty
            if (actual_c.empty() && !expected_f.empty()){
                 FAIL() << name << " actual_c is empty but expected_f is not.";
            }
            return;
        }
        ASSERT_EQ(expected_f.size(), (size_t)rows * cols);
        if (!actual_c.empty()) {
             ASSERT_EQ(actual_c.size(), (size_t)rows * cols);
        } else {
             FAIL() << name << " actual_c is empty but expected_f is not.";
        }


        std::vector<double> actual_f(expected_f.size());
        if (is_row_major_input) {
            slicot_transpose_to_fortran_with_ld(actual_c.data(), actual_f.data(), 
                                                rows, cols, ld_c, std::max(1, rows), sizeof(double));
        } else {
            if (ld_c == rows) {
                actual_f = actual_c;
            } else {
                for (int j = 0; j < cols; ++j) {
                    for (int i = 0; i < rows; ++i) {
                        actual_f[i + j * rows] = actual_c[i + j * ld_c];
                    }
                }
            }
        }
        for (size_t i = 0; i < expected_f.size(); ++i) {
            EXPECT_NEAR(actual_f[i], expected_f[i], check_tol)
                << name << " matrix mismatch at flat index " << i 
                << ", actual_f: " << actual_f[i] << ", expected_f: " << expected_f[i];
        }
    }
    
    void VerifyVector(const std::string& name, const std::vector<double>& actual, 
                      const std::vector<double>& expected, int dim) {
        if (dim == 0) {
            EXPECT_TRUE(expected.empty());
            EXPECT_TRUE(actual.empty());
            return;
        }
        ASSERT_EQ(expected.size(), (size_t)dim);
        ASSERT_EQ(actual.size(), (size_t)dim);
        for (int i = 0; i < dim; ++i) {
            EXPECT_NEAR(actual[i], expected[i], check_tol)
                << name << " vector mismatch at index " << i;
        }
    }


    void VerifyOutputs(bool is_row_major_input) {
        ASSERT_EQ(INFO_out, INFO_expected);
        if (INFO_expected != 0) return;

        VerifyMatrix("A", A_actual_c, A_expected_f, L_val, N_val, LDA_c, is_row_major_input);
        VerifyMatrix("E", E_actual_c, E_expected_f, L_val, N_val, LDE_c, is_row_major_input);
        
        if (M_val > 0) {
            VerifyMatrix("B", B_actual_c, B_expected_f, L_val, M_val, LDB_c, is_row_major_input);
        } else { EXPECT_TRUE(B_expected_f.empty()); }
        
        if (P_val > 0) {
            VerifyMatrix("C", C_actual_c, C_expected_f, P_val, N_val, LDC_c, is_row_major_input);
        } else { EXPECT_TRUE(C_expected_f.empty()); }

        VerifyVector("LSCALE", LSCALE_out, LSCALE_expected, L_val);
        VerifyVector("RSCALE", RSCALE_out, RSCALE_expected, N_val);
    }
};

// --- Test Fixture for TG01AD Documentation Example (from Python test) ---
class Tg01adDocExampleTest : public Tg01adBaseTest {
protected:
    void SetUpTestData() override {
        L_val = 4; N_val = 4; M_val = 2; P_val = 2;
        JOB_val = 'A'; THRESH_val = 0.0;
        INFO_expected = 0;

        // Input data (Fortran column-major for A_c_layout in ColMajor test)
        A_c_layout = {-1.0, 0.0, 100.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.003, 0.02, 0.4, 0.0};
        E_c_layout = {1.0, 0.0, 300.0, 0.0, 0.2, 1.0, 90.0, 0.0, 0.0, 0.0, 6.0, 20.0, 0.0, 0.01, 0.3, 0.0};
        B_c_layout = {10.0, 0.0, 0.0, 10000.0, 0.0, 0.0, 1000.0, 10000.0};
        C_c_layout = {-0.1, 0.0, 0.0, 0.01, 0.001, -0.001, 0.0, 0.0001};

        // Expected outputs (Fortran column-major)
        A_expected_f = {-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.3, 2.0, 0.4, 0.0};
        E_expected_f = {1.0, 0.0, 3.0, 0.0, 0.2, 1.0, 0.9, 0.0, 0.0, 0.0, 0.6, 0.2, 0.0, 1.0, 0.3, 0.0};
        B_expected_f = {100.0, 0.0, 0.0, 100.0, 0.0, 0.0, 100.0, 100.0};
        C_expected_f = {-0.01, 0.0, 0.0, 0.001, 0.001, -0.001, 0.0, 0.001};
        LSCALE_expected = {10.0, 10.0, 0.1, 0.01};
        RSCALE_expected = {0.1, 0.1, 1.0, 10.0};
    }
};

TEST_F(Tg01adDocExampleTest, ColMajor) {
    SetUpBase(false); // Column-major C setup
    // A_c_layout etc. are already in Fortran column-major from SetUpTestData for this test
    
    INFO_out = slicot_tg01ad(JOB_val, L_val, N_val, M_val, P_val, THRESH_val,
                             A_actual_c.data(), LDA_c, E_actual_c.data(), LDE_c,
                             B_actual_c.data(), LDB_c, C_actual_c.data(), LDC_c,
                             LSCALE_out.data(), RSCALE_out.data(), 0);
    VerifyOutputs(false);
}

TEST_F(Tg01adDocExampleTest, RowMajor) {
    SetUpBase(true); // Row-major C setup

    // Transpose initial col-major data from SetUpTestData to row-major for A_actual_c etc.
    // A_c_layout, E_c_layout, B_c_layout, C_c_layout were set by SetUpTestData in col-major form
    // A_actual_c = A_c_layout; // This was already done in SetUpBase, now transpose it
    
    std::vector<double> temp_A = A_actual_c;
    std::vector<double> temp_E = E_actual_c;
    std::vector<double> temp_B = B_actual_c;
    std::vector<double> temp_C = C_actual_c;

    slicot_transpose_to_c_with_ld(temp_A.data(), A_actual_c.data(), L_val, N_val, L_val, N_val, sizeof(double));
    slicot_transpose_to_c_with_ld(temp_E.data(), E_actual_c.data(), L_val, N_val, L_val, N_val, sizeof(double));
    if (M_val > 0) slicot_transpose_to_c_with_ld(temp_B.data(), B_actual_c.data(), L_val, M_val, L_val, M_val, sizeof(double));
    if (P_val > 0) slicot_transpose_to_c_with_ld(temp_C.data(), C_actual_c.data(), P_val, N_val, P_val, N_val, sizeof(double));
    
    INFO_out = slicot_tg01ad(JOB_val, L_val, N_val, M_val, P_val, THRESH_val,
                             A_actual_c.data(), LDA_c, E_actual_c.data(), LDE_c,
                             B_actual_c.data(), LDB_c, C_actual_c.data(), LDC_c,
                             LSCALE_out.data(), RSCALE_out.data(), 1);
    VerifyOutputs(true);
}


// --- Test Fixture for Parameter Validation ---
class Tg01adParamValidationTest : public Tg01adBaseTest {
protected:
    void SetUpTestData() override {
        // Minimal valid setup for most validation tests
        JOB_val = 'A'; L_val = 1; N_val = 1; M_val = 0; P_val = 0; THRESH_val = 0.0;
        A_c_layout = {1.0}; E_c_layout = {1.0}; 
        B_c_layout.clear(); C_c_layout.clear(); // M=0, P=0
        // Expected outputs are not relevant for error checks, INFO_expected is key
        A_expected_f.clear(); E_expected_f.clear(); B_expected_f.clear(); C_expected_f.clear();
        LSCALE_expected.clear(); RSCALE_expected.clear();
    }
};

TEST_F(Tg01adParamValidationTest, InvalidJob) {
    INFO_expected = -1; SetUpBase(false); JOB_val = 'X';
    INFO_out = slicot_tg01ad(JOB_val,L_val,N_val,M_val,P_val,THRESH_val, A_actual_c.data(),LDA_c,E_actual_c.data(),LDE_c,B_actual_c.data(),LDB_c,C_actual_c.data(),LDC_c,LSCALE_out.data(),RSCALE_out.data(),0);
    ASSERT_EQ(INFO_out, INFO_expected);
}
TEST_F(Tg01adParamValidationTest, InvalidL) {
    INFO_expected = -2; SetUpBase(false); L_val = -1;
    INFO_out = slicot_tg01ad(JOB_val,L_val,N_val,M_val,P_val,THRESH_val, A_actual_c.data(),LDA_c,E_actual_c.data(),LDE_c,B_actual_c.data(),LDB_c,C_actual_c.data(),LDC_c,LSCALE_out.data(),RSCALE_out.data(),0);
    ASSERT_EQ(INFO_out, INFO_expected);
}
// ... Add more validation tests for N, M, P, THRESH, NULL pointers, invalid LDs ...


// --- Test Fixture for Zero Dimensions ---
class Tg01adZeroDimTest : public Tg01adBaseTest {
protected:
    void SetUpTestData() override { // Default, overridden by specific tests
        JOB_val = 'A'; THRESH_val = 0.0; INFO_expected = 0;
    }
};

TEST_F(Tg01adZeroDimTest, L0_N0_M0_P0) {
    L_val = 0; N_val = 0; M_val = 0; P_val = 0;
    SetUpBase(false);
    INFO_out = slicot_tg01ad(JOB_val,L_val,N_val,M_val,P_val,THRESH_val, nullptr,LDA_c,nullptr,LDE_c,nullptr,LDB_c,nullptr,LDC_c,nullptr,nullptr,0);
    VerifyOutputs(false); // Expects INFO=0, all output arrays empty
}

TEST_F(Tg01adZeroDimTest, L1_N1_M0_P0) {
    L_val = 1; N_val = 1; M_val = 0; P_val = 0;
    A_c_layout = {2.0}; E_c_layout = {1.0}; B_c_layout.clear(); C_c_layout.clear();
    A_expected_f = {2.0}; E_expected_f = {1.0}; // Assuming no scaling if M,P=0 and JOB='A' or 'N'
    B_expected_f.clear(); C_expected_f.clear();
    LSCALE_expected = {1.0}; RSCALE_expected = {1.0}; // Expect no scaling
    SetUpBase(false);
    INFO_out = slicot_tg01ad(JOB_val,L_val,N_val,M_val,P_val,THRESH_val, A_actual_c.data(),LDA_c,E_actual_c.data(),LDE_c,nullptr,LDB_c,nullptr,LDC_c,LSCALE_out.data(),RSCALE_out.data(),0);
    VerifyOutputs(false);
}


