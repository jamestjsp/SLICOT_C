#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max, std::min
#include <iostream> 
#include <iomanip>   

#include "tb01id.h"       
#include "slicot_utils.h" 

// Helper to print matrices for debugging
void print_matrix_tb01id(const std::string& name, const std::vector<double>& mat, int r, int c, int ld_c, bool is_row_major) {
    std::cout << name << " (" << r << "x" << c << ") " << (is_row_major ? "RowMajor" : "ColMajor") << ", LD=" << ld_c << ":\n";
    if (mat.empty() || r == 0 || c == 0) {
        std::cout << "  (empty or zero dim)\n";
        return;
    }
    for (int i = 0; i < r; ++i) {
        std::cout << "  [";
        for (int j = 0; j < c; ++j) {
            double val = is_row_major ? mat[i * ld_c + j] : mat[i + j * ld_c]; // Corrected indexing for col-major
            std::cout << std::fixed << std::setprecision(8) << val << (j == c - 1 ? "" : ", ");
        }
        std::cout << "]\n";
    }
}
void print_vector_tb01id(const std::string& name, const std::vector<double>& vec) {
    std::cout << name << " (" << vec.size() << "):\n  [";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << std::fixed << std::setprecision(8) << vec[i] << (i == vec.size() - 1 ? "" : ", ");
    }
    std::cout << "]\n";
}


// --- Test Fixture Base ---
class TB01IDTest : public ::testing::Test {
protected:
    char JOB_param = 'A';
    int N_param = 0;
    int M_param = 0;
    int P_param = 0;
    double MAXRED_param_in = 0.0; // Use default in routine (10.0)

    // Input data (row-major for easy initialization)
    std::vector<double> A_data_rm, B_data_rm, C_data_rm;
    
    // Buffers for C function call (will be copies of _data_rm)
    std::vector<double> A_test_buf, B_test_buf, C_test_buf;
    double MAXRED_io_val;
    std::vector<double> SCALE_out_val;

    // Expected results (row-major)
    std::vector<double> A_expected_rm, B_expected_rm, C_expected_rm, SCALE_expected_val;
    double MAXRED_expected_out;
    int expected_info = 0;

    // Leading dimensions for C call
    int LDA, LDB, LDC;
    
    double check_tol_matrix = 1e-5; 
    double check_tol_scalar_maxred = 2.0e-4; // Tolerance for MAXRED (relative)


    void setup_slicot_example_data() {
        JOB_param = 'A';
        N_param = 5; M_param = 2; P_param = 5;
        MAXRED_param_in = 0.0; // Let routine use default

        A_data_rm = {
            0.0,  1.0000e+000,          0.0,          0.0,          0.0,
            -1.5800e+006, -1.2570e+003,          0.0,          0.0,          0.0,
            3.5410e+014,          0.0, -1.4340e+003,          0.0, -5.3300e+011,
            0.0,          0.0,          0.0,          0.0,  1.0000e+000,
            0.0,          0.0,          0.0, -1.8630e+004, -1.4820e+000
        };
        B_data_rm = {
            0.0,          0.0,
            1.1030e+002,          0.0,
            0.0,          0.0,
            0.0,          0.0,
            0.0,  8.3330e-003
        };
        C_data_rm = {
            1.0000e+000,          0.0,          0.0,          0.0,          0.0,
            0.0,          0.0,  1.0000e+000,          0.0,          0.0,
            0.0,          0.0,          0.0,  1.0000e+000,          0.0,
            6.6640e-001,          0.0, -6.2000e-013,          0.0,          0.0,
            0.0,          0.0, -1.0000e-003,  1.8960e+006,  1.5080e+002
        };

        A_expected_rm = { // From TB01ID.html example output
            0.0000000E+00,  0.1000000E+05,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,
           -0.1580000E+03, -0.1257000E+04,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,
            0.3541000E+05,  0.0000000E+00, -0.1434000E+04,  0.0000000E+00, -0.5330000E+03,
            0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.1000000E+03,
            0.0000000E+00,  0.0000000E+00,  0.0000000E+00, -0.1863000E+03, -0.1482000E+01
        };
        B_expected_rm = {
            0.0000000E+00,  0.0000000E+00,
            0.1103000E+04,  0.0000000E+00,
            0.0000000E+00,  0.0000000E+00,
            0.0000000E+00,  0.0000000E+00,
            0.0000000E+00,  0.8333000E+02
        };
        C_expected_rm = {
            0.1000000E-04,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,
            0.0000000E+00,  0.0000000E+00,  0.1000000E+06,  0.0000000E+00,  0.0000000E+00,
            0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.1000000E-05,  0.0000000E+00,
            0.6664000E-05,  0.0000000E+00, -0.6200000E-07,  0.0000000E+00,  0.0000000E+00,
            0.0000000E+00,  0.0000000E+00, -0.1000000E+03,  0.1896000E+01,  0.1508000E-01
        };
        SCALE_expected_val = {0.1000000E-04, 0.1000000E+00, 0.1000000E+06, 0.1000000E-05, 0.1000000E-03};
        MAXRED_expected_out = 0.3488E+10; 
        check_tol_matrix = 1e-4; // Example output precision for matrices
        // check_tol_scalar_maxred is 2.0e-4 from fixture default for MAXRED relative comparison

        if (N_param > 0) { 
            SCALE_out_val.resize(N_param);
        } else {
            SCALE_out_val.clear();
        }
    }
};

TEST_F(TB01IDTest, SlicotDocExampleRowMajor) {
    setup_slicot_example_data();
    LDA = N_param; LDB = M_param; LDC = N_param; 
    
    A_test_buf = A_data_rm; 
    B_test_buf = B_data_rm; 
    C_test_buf = C_data_rm;
    MAXRED_io_val = MAXRED_param_in;
    
    int info = slicot_tb01id(JOB_param, N_param, M_param, P_param,
                             &MAXRED_io_val,
                             A_test_buf.data(), LDA, 
                             B_test_buf.data(), LDB, 
                             C_test_buf.data(), LDC,
                             SCALE_out_val.data(),
                             1 /* row_major = true */);
    ASSERT_EQ(info, 0);
    if (std::abs(MAXRED_expected_out) > 1e-9) { // Use relative tolerance for MAXRED
         EXPECT_NEAR(MAXRED_io_val, MAXRED_expected_out, check_tol_scalar_maxred * std::abs(MAXRED_expected_out));
    } else { // Use absolute tolerance if expected is near zero
         EXPECT_NEAR(MAXRED_io_val, MAXRED_expected_out, check_tol_scalar_maxred);
    }

    for(size_t i=0; i < A_expected_rm.size(); ++i) EXPECT_NEAR(A_test_buf[i], A_expected_rm[i], check_tol_matrix);
    if (M_param > 0) for(size_t i=0; i < B_expected_rm.size(); ++i) EXPECT_NEAR(B_test_buf[i], B_expected_rm[i], check_tol_matrix);
    if (P_param > 0) for(size_t i=0; i < C_expected_rm.size(); ++i) EXPECT_NEAR(C_test_buf[i], C_expected_rm[i], check_tol_matrix);
    if (N_param > 0) for(size_t i=0; i < SCALE_expected_val.size(); ++i) EXPECT_NEAR(SCALE_out_val[i], SCALE_expected_val[i], check_tol_matrix);
}

TEST_F(TB01IDTest, SlicotDocExampleColMajor) {
    setup_slicot_example_data();
    LDA = std::max(1,N_param); LDB = std::max(1,N_param); LDC = std::max(1,P_param); 

    A_test_buf.resize((size_t)LDA*N_param); slicot_transpose_to_fortran_with_ld(A_data_rm.data(), A_test_buf.data(), N_param, N_param, N_param, LDA, sizeof(double));
    if (M_param > 0 && N_param > 0) { B_test_buf.resize((size_t)LDB*M_param); slicot_transpose_to_fortran_with_ld(B_data_rm.data(), B_test_buf.data(), N_param, M_param, M_param, LDB, sizeof(double));}
    else { B_test_buf.clear(); } // Clear if not used or N=0
    if (P_param > 0 && N_param > 0) { C_test_buf.resize((size_t)LDC*N_param); slicot_transpose_to_fortran_with_ld(C_data_rm.data(), C_test_buf.data(), P_param, N_param, N_param, LDC, sizeof(double));}
    else { C_test_buf.clear(); } // Clear if not used or N=0 or P=0
    MAXRED_io_val = MAXRED_param_in;
    
    int info = slicot_tb01id(JOB_param, N_param, M_param, P_param,
                             &MAXRED_io_val,
                             A_test_buf.data(), LDA, 
                             B_test_buf.empty() ? nullptr : B_test_buf.data(), LDB, 
                             C_test_buf.empty() ? nullptr : C_test_buf.data(), LDC,
                             SCALE_out_val.data(),
                             0 /* row_major = false */);
    ASSERT_EQ(info, 0);
    if (std::abs(MAXRED_expected_out) > 1e-9) { // Use relative tolerance for MAXRED
         EXPECT_NEAR(MAXRED_io_val, MAXRED_expected_out, check_tol_scalar_maxred * std::abs(MAXRED_expected_out));
    } else { // Use absolute tolerance if expected is near zero
         EXPECT_NEAR(MAXRED_io_val, MAXRED_expected_out, check_tol_scalar_maxred);
    }

    std::vector<double> A_exp_cm((size_t)N_param*N_param); slicot_transpose_to_fortran_with_ld(A_expected_rm.data(), A_exp_cm.data(), N_param,N_param,N_param,N_param,sizeof(double));
    for(size_t i=0; i < A_exp_cm.size(); ++i) EXPECT_NEAR(A_test_buf[i], A_exp_cm[i], check_tol_matrix);

    if (M_param > 0 && N_param > 0) {
        std::vector<double> B_exp_cm((size_t)N_param*M_param); slicot_transpose_to_fortran_with_ld(B_expected_rm.data(), B_exp_cm.data(), N_param,M_param,M_param,N_param,sizeof(double));
        for(size_t i=0; i < B_exp_cm.size(); ++i) EXPECT_NEAR(B_test_buf[i], B_exp_cm[i], check_tol_matrix);
    }
    if (P_param > 0 && N_param > 0) {
        std::vector<double> C_exp_cm((size_t)P_param*N_param); slicot_transpose_to_fortran_with_ld(C_expected_rm.data(), C_exp_cm.data(), P_param,N_param,N_param,P_param,sizeof(double));
        for(size_t i=0; i < C_exp_cm.size(); ++i) EXPECT_NEAR(C_test_buf[i], C_exp_cm[i], check_tol_matrix);
    }
    if (N_param > 0) for(size_t i=0; i < SCALE_expected_val.size(); ++i) EXPECT_NEAR(SCALE_out_val[i], SCALE_expected_val[i], check_tol_matrix);
}


TEST_F(TB01IDTest, ParameterValidation) {
    N_param = 1; M_param = 1; P_param = 1; 
    A_test_buf.assign(1,0.0); B_test_buf.assign(1,0.0); C_test_buf.assign(1,0.0);
    SCALE_out_val.assign(1,0.0); MAXRED_io_val = 0.0;
    LDA=1; LDB=1; LDC=1;

    int info = slicot_tb01id('X', N_param, M_param, P_param, &MAXRED_io_val, A_test_buf.data(),LDA,B_test_buf.data(),LDB,C_test_buf.data(),LDC,SCALE_out_val.data(),1);
    EXPECT_EQ(info, -1); 

    info = slicot_tb01id(JOB_param, -1, M_param, P_param, &MAXRED_io_val, nullptr,1,nullptr,1,nullptr,1,nullptr,1);
    EXPECT_EQ(info, -2); 
    
    if (N_param > 0) { 
        info = slicot_tb01id(JOB_param, N_param, M_param, P_param, &MAXRED_io_val, A_test_buf.data(),0,B_test_buf.data(),LDB,C_test_buf.data(),LDC,SCALE_out_val.data(),1);
        EXPECT_EQ(info, -7); 
    }
}

TEST_F(TB01IDTest, ZeroDimensionN) {
    N_param = 0; M_param = 1; P_param = 1; 
    A_test_buf.clear(); 
    // B is N x M. If N=0, B is 0xM (empty).
    // C is P x N. If N=0, C is Px0 (empty).
    // SCALE is N-dim, so empty.
    B_test_buf.clear(); 
    C_test_buf.clear(); 
    SCALE_out_val.clear(); 
    MAXRED_io_val = 0.0;
    LDA=1; LDB=1; LDC=1; 

    int info = slicot_tb01id('A', N_param, M_param, P_param,
                             &MAXRED_io_val,
                             nullptr, LDA, 
                             nullptr, LDB, // B is 0xM, pass nullptr
                             nullptr, LDC, // C is Px0, pass nullptr
                             nullptr,      // SCALE is 0-dim, pass nullptr
                             1 /* row_major = true */);
    EXPECT_EQ(info, 0);
}

