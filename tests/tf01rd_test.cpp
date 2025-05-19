#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // For std::max
#include <numeric>   // For std::accumulate
#include <cctype>    // For toupper C++ style

#include "tf01rd.h"       // Wrapper header
#include "slicot_utils.h" // For transpose functions

// --- Base Test Fixture for TF01RD ---
class Tf01rdBaseTest : public ::testing::Test {
protected:
    // Parameters for slicot_tf01rd
    int NA_val, NB_val, NC_val, N_val;
    std::vector<double> A_c_layout;
    std::vector<double> B_c_layout;
    std::vector<double> C_c_layout;
    
    // Expected outputs (Fortran column-major for comparison)
    std::vector<double> H_expected_f;
    int INFO_expected;

    // Actual outputs from wrapper
    std::vector<double> H_out_c;
    int INFO_out;

    // Leading dimensions for C arrays
    int LDA_c, LDB_c, LDC_c, LDH_c;

    double check_tol = 1e-4; // Tolerance for checking double values

    virtual void SetUpTestData() = 0; // To be implemented by derived fixtures

    void SetUpBase(bool is_row_major) {
        SetUpTestData(); // Populate parameters and expected values

        // Size output H array
        // H is NC x (N*NB)
        int h_rows = NC_val;
        int h_cols = N_val * NB_val;
        if (h_rows > 0 && h_cols > 0) {
            H_out_c.assign((size_t)h_rows * h_cols, 0.0);
        } else {
            H_out_c.clear(); // If H is logically empty
        }
        
        // Calculate C leading dimensions based on row_major flag
        if (is_row_major) {
            LDA_c = std::max(1, NA_val); // cols of A
            LDB_c = std::max(1, NB_val); // cols of B
            LDC_c = std::max(1, NA_val); // cols of C
            LDH_c = std::max(1, h_cols); // cols of H
        } else { // Column-major C
            LDA_c = std::max(1, NA_val); // rows of A
            LDB_c = std::max(1, NA_val); // rows of B
            LDC_c = std::max(1, NC_val); // rows of C
            LDH_c = std::max(1, NC_val); // rows of H
        }
    }

    void VerifyOutputs(bool is_row_major_input) {
        ASSERT_EQ(INFO_out, INFO_expected);
        if (INFO_expected != 0) return;

        int expected_h_rows = NC_val;
        int expected_h_cols = N_val * NB_val;

        if (expected_h_rows == 0 || expected_h_cols == 0) {
            EXPECT_TRUE(H_expected_f.empty());
            if (expected_h_rows == 0 || expected_h_cols == 0) { 
                 EXPECT_TRUE(H_out_c.empty());
            }
            return;
        }

        ASSERT_EQ(H_expected_f.size(), (size_t)expected_h_rows * expected_h_cols);
        if(!H_out_c.empty()){ 
             ASSERT_EQ(H_out_c.size(), (size_t)expected_h_rows * expected_h_cols);
        } else {
            // This case should ideally not be hit if H_expected_f is not empty
            if (!H_expected_f.empty()) {
                 FAIL() << "H_out_c is empty but expected H is not.";
            } else {
                 // Both are empty, which is fine.
            }
        }


        std::vector<double> H_actual_f(H_expected_f.size());
         if (H_expected_f.empty()) { // If H is expected to be empty, nothing to compare
            return;
        }


        if (is_row_major_input) {
            slicot_transpose_to_fortran_with_ld(H_out_c.data(), H_actual_f.data(), 
                                                expected_h_rows, expected_h_cols, 
                                                LDH_c, std::max(1, NC_val), sizeof(double));
        } else { 
            if (LDH_c == expected_h_rows) { 
                 H_actual_f = H_out_c;
            } else { 
                 for (int j = 0; j < expected_h_cols; ++j) { // Iterate through columns
                    for (int i = 0; i < expected_h_rows; ++i) { // Iterate through rows
                        H_actual_f[i + j * expected_h_rows] = H_out_c[i + j * LDH_c];
                    }
                }
            }
        }

        for (size_t i = 0; i < H_expected_f.size(); ++i) {
            EXPECT_NEAR(H_actual_f[i], H_expected_f[i], check_tol)
                << "H matrix mismatch at flat index " << i 
                << ", actual: " << H_actual_f[i] << ", expected: " << H_expected_f[i];
        }
    }
};

// --- Test Fixture for TF01RD Documentation Example ---
class Tf01rdDocExampleTest : public Tf01rdBaseTest {
protected:
    void SetUpTestData() override {
        NA_val = 3; NB_val = 2; NC_val = 2; N_val = 5;
        INFO_expected = 0;

        // Input A, B, C from TF01RD.html "Program Data" 
        // These are defined as flat arrays in Fortran column-major order.
        // A_c_layout will be used directly for ColMajor test, transposed for RowMajor test.
        A_c_layout = {0.000, 1.000, 0.000, -0.070, 0.800, 0.000, 0.015, -0.150, 0.500};
        B_c_layout = {0.000, 2.000, 1.000, -1.000, -0.100, 1.000};
        C_c_layout = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0};

        // Expected H computed by Python script using Slycot with the above A,B,C
        // H is (NC x N*NB) = 2x10, Fortran column-major flattened.
        H_expected_f = {
            -2.0000, 0.0000, 0.1000, -1.0000, // M(1)
            -1.4500, -0.1250, 1.2300, 0.0220,  // M(2)
            -0.9600, -0.0940, 1.0370, 0.0936,  // M(3)
            -0.6365, -0.0635, 0.7735, 0.0763,  // M(4)
            -0.4270, -0.0427, 0.5612, 0.0560   // M(5)
        };
    }
};

TEST_F(Tf01rdDocExampleTest, ColMajor) {
    SetUpBase(false); // Column-major C setup
    // A_c_layout, B_c_layout, C_c_layout are already in Fortran column-major from SetUpTestData
    
    INFO_out = slicot_tf01rd(NA_val, NB_val, NC_val, N_val,
                             A_c_layout.data(), LDA_c,
                             B_c_layout.data(), LDB_c,
                             C_c_layout.data(), LDC_c,
                             H_out_c.data(), LDH_c,
                             0 /*row_major=false*/);
    VerifyOutputs(false);
}

TEST_F(Tf01rdDocExampleTest, RowMajor) {
    SetUpBase(true); // Row-major C setup

    // Original A_c_layout, B_c_layout, C_c_layout from SetUpTestData are in Fortran column-major.
    // Create temporary row-major versions for this test.
    std::vector<double> A_f_colmajor = {0.000, 1.000, 0.000, -0.070, 0.800, 0.000, 0.015, -0.150, 0.500};
    std::vector<double> B_f_colmajor = {0.000, 2.000, 1.000, -1.000, -0.100, 1.000};
    std::vector<double> C_f_colmajor = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0};

    A_c_layout.resize(A_f_colmajor.size());
    B_c_layout.resize(B_f_colmajor.size());
    C_c_layout.resize(C_f_colmajor.size());

    slicot_transpose_to_c_with_ld(A_f_colmajor.data(), A_c_layout.data(), NA_val, NA_val, NA_val, NA_val, sizeof(double));
    slicot_transpose_to_c_with_ld(B_f_colmajor.data(), B_c_layout.data(), NA_val, NB_val, NA_val, NB_val, sizeof(double));
    slicot_transpose_to_c_with_ld(C_f_colmajor.data(), C_c_layout.data(), NC_val, NA_val, NC_val, NA_val, sizeof(double));

    INFO_out = slicot_tf01rd(NA_val, NB_val, NC_val, N_val,
                             A_c_layout.data(), LDA_c,
                             B_c_layout.data(), LDB_c,
                             C_c_layout.data(), LDC_c,
                             H_out_c.data(), LDH_c,
                             1 /*row_major=true*/);
    VerifyOutputs(true);
}

// --- Test Fixture for Parameter Validation ---
class Tf01rdParamValidationTest : public Tf01rdBaseTest {
protected:
    void SetUpTestData() override {
        NA_val = 1; NB_val = 1; NC_val = 1; N_val = 1;
        A_c_layout = {1.0}; B_c_layout = {1.0}; C_c_layout = {1.0};
        H_expected_f.clear(); 
    }
};

TEST_F(Tf01rdParamValidationTest, InvalidNA) {
    INFO_expected = -1; SetUpBase(false); NA_val = -1; 
    INFO_out = slicot_tf01rd(NA_val,NB_val,NC_val,N_val, A_c_layout.data(),LDA_c, B_c_layout.data(),LDB_c, C_c_layout.data(),LDC_c, H_out_c.data(),LDH_c, 0);
    ASSERT_EQ(INFO_out, INFO_expected);
}
TEST_F(Tf01rdParamValidationTest, InvalidNB) {
    INFO_expected = -2; SetUpBase(false); NB_val = -1;
    INFO_out = slicot_tf01rd(NA_val,NB_val,NC_val,N_val, A_c_layout.data(),LDA_c, B_c_layout.data(),LDB_c, C_c_layout.data(),LDC_c, H_out_c.data(),LDH_c, 0);
    ASSERT_EQ(INFO_out, INFO_expected);
}
TEST_F(Tf01rdParamValidationTest, InvalidNC) {
    INFO_expected = -3; SetUpBase(false); NC_val = -1;
    INFO_out = slicot_tf01rd(NA_val,NB_val,NC_val,N_val, A_c_layout.data(),LDA_c, B_c_layout.data(),LDB_c, C_c_layout.data(),LDC_c, H_out_c.data(),LDH_c, 0);
    ASSERT_EQ(INFO_out, INFO_expected);
}
TEST_F(Tf01rdParamValidationTest, InvalidN) {
    INFO_expected = -4; SetUpBase(false); N_val = -1;
    INFO_out = slicot_tf01rd(NA_val,NB_val,NC_val,N_val, A_c_layout.data(),LDA_c, B_c_layout.data(),LDB_c, C_c_layout.data(),LDC_c, H_out_c.data(),LDH_c, 0);
    ASSERT_EQ(INFO_out, INFO_expected);
}


// --- Test Fixture for Zero Dimensions ---
class Tf01rdZeroDimTest : public Tf01rdBaseTest {
protected:
    void SetUpTestData() override {
        NA_val = 0; NB_val = 0; NC_val = 0; N_val = 0; // Default, overridden by tests
        A_c_layout.clear(); B_c_layout.clear(); C_c_layout.clear();
        H_expected_f.clear();
        INFO_expected = 0;
    }
};

TEST_F(Tf01rdZeroDimTest, NA_is_0) {
    NA_val = 0; NB_val = 1; NC_val = 1; N_val = 1;
    INFO_expected = 0; 
    H_expected_f.assign( (size_t)NC_val * N_val * NB_val, 0.0); 
    SetUpBase(false); 
    
    INFO_out = slicot_tf01rd(NA_val,NB_val,NC_val,N_val, nullptr,LDA_c, nullptr,LDB_c, nullptr,LDC_c, H_out_c.data(),LDH_c, 0);
    VerifyOutputs(false);
}

TEST_F(Tf01rdZeroDimTest, N_is_0) {
    NA_val = 1; NB_val = 1; NC_val = 1; N_val = 0;
    A_c_layout = {1.0}; B_c_layout = {1.0}; C_c_layout = {1.0};
    INFO_expected = 0; 
    H_expected_f.clear(); 
    SetUpBase(false);
    INFO_out = slicot_tf01rd(NA_val,NB_val,NC_val,N_val, A_c_layout.data(),LDA_c, B_c_layout.data(),LDB_c, C_c_layout.data(),LDC_c, H_out_c.data(),LDH_c, 0);
    VerifyOutputs(false);
}

TEST_F(Tf01rdZeroDimTest, NB_is_0) {
    NA_val = 1; NB_val = 0; NC_val = 1; N_val = 1;
    A_c_layout = {1.0}; B_c_layout = {}; C_c_layout = {1.0};
    INFO_expected = 0; 
    H_expected_f.clear();
    SetUpBase(false);
    INFO_out = slicot_tf01rd(NA_val,NB_val,NC_val,N_val, A_c_layout.data(),LDA_c, nullptr ,LDB_c, C_c_layout.data(),LDC_c, H_out_c.data(),LDH_c, 0);
    VerifyOutputs(false);
}

TEST_F(Tf01rdZeroDimTest, NC_is_0) {
    NA_val = 1; NB_val = 1; NC_val = 0; N_val = 1;
    A_c_layout = {1.0}; B_c_layout = {1.0}; C_c_layout = {};
    INFO_expected = 0; 
    H_expected_f.clear();
    SetUpBase(false);
    INFO_out = slicot_tf01rd(NA_val,NB_val,NC_val,N_val, A_c_layout.data(),LDA_c, B_c_layout.data(),LDB_c, nullptr ,LDC_c, H_out_c.data(),LDH_c, 0);
    VerifyOutputs(false);
}

