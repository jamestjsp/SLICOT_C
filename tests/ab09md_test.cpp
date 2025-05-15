#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max, std::min
#include <stdexcept> // For std::runtime_error
#include <iostream>  // For debugging output if needed
#include <iomanip>   // For std::fixed and std::setprecision

#include "ab09md.h"
#include "slicot_utils.h" 

// --- Column-Major Test Fixture ---
class AB09MDTestColMajor : public ::testing::Test {
protected:
    // Test parameters from AB09MD.html example
    char DICO_in = 'C';
    char JOB_in = 'N';
    char EQUIL_in = 'N'; // Matches Fortran example data and Slycot validation
    char ORDSEL_in = 'A'; // 'A' for automatic order selection
    int N_in = 7;
    int M_in = 2;
    int P_in = 3;
    int NR_io_auto_order = 0;  // Initial value for ORDSEL_in = 'A', will be updated by routine
    double ALPHA_in = -0.6;
    double TOL_in = 1.0e-1;
    
    int IWARN_out = 0;
    int info_result = -999;
    int NS_out = 0;

    double check_tol = 1e-4; // Tolerance for floating point comparisons

    std::vector<double> A_io; 
    std::vector<double> B_io; 
    std::vector<double> C_io; 
    std::vector<double> HSV_out;

    std::vector<double> A_expected_reduced;   
    std::vector<double> B_expected_reduced;   
    std::vector<double> C_expected_reduced;   
    std::vector<double> HSV_expected;     
    int expected_info = 0;
    int expected_iwarn = 0; 
    int expected_nr_auto = 5; // Expected reduced order when ORDSEL = 'A'
    int expected_ns = 5;      // Expected number of stable eigenvalues

    int LDA_io = 0;
    int LDB_io = 0;
    int LDC_io = 0;

    void SetUpBase() {
        LDA_io = std::max(1, N_in); // For A (N x N)
        LDB_io = std::max(1, N_in); // For B (N x M)
        LDC_io = std::max(1, P_in); // For C (P x N)

        // Fortran Example Program Data (A is 7x7)
        // Data file lists elements row-wise. We need to store them column-wise.
        // A(1,:): -0.04165  0.0000  4.9200  -4.9200  0.0000  0.0000  0.0000
        // A(2,:): -5.2100  -12.500  0.0000   0.0000  0.0000  0.0000  0.0000
        // A(3,:):  0.0000   3.3300 -3.3300   0.0000  0.0000  0.0000  0.0000
        // A(4,:):  0.5450   0.0000  0.0000   0.0000 -0.5450  0.0000  0.0000
        // A(5,:):  0.0000   0.0000  0.0000   4.9200 -0.04165 0.0000  4.9200
        // A(6,:):  0.0000   0.0000  0.0000   0.0000 -5.2100 -12.500  0.0000
        // A(7,:):  0.0000   0.0000  0.0000   0.0000  0.0000  3.3300 -3.3300
        A_io = {
            // Column 0 (A(:,1) from data file)
            -0.04165, -5.21000,  0.00000,  0.54500,  0.00000,  0.00000,  0.00000,
            // Column 1 (A(:,2) from data file)
             0.00000, -12.50000,  3.33000,  0.00000,  0.00000,  0.00000,  0.00000,
            // Column 2 (A(:,3) from data file)
             4.92000,  0.00000, -3.33000,  0.00000,  0.00000,  0.00000,  0.00000,
            // Column 3 (A(:,4) from data file)
            -4.92000,  0.00000,  0.00000,  0.00000,  4.92000,  0.00000,  0.00000,
            // Column 4 (A(:,5) from data file)
             0.00000,  0.00000,  0.00000, -0.54500, -0.04165, -5.21000,  0.00000,
            // Column 5 (A(:,6) from data file)
             0.00000,  0.00000,  0.00000,  0.00000,  0.00000, -12.50000,  3.33000,
            // Column 6 (A(:,7) from data file)
             0.00000,  0.00000,  0.00000,  0.00000,  4.92000,  0.00000, -3.33000
        };

        // Fortran Example Program Data for B (7x2)
        // B(1,:):  0.0000   0.0000
        // B(2,:): 12.500   0.0000
        // B(3,:):  0.0000   0.0000
        // B(4,:):  0.0000   0.0000
        // B(5,:):  0.0000   0.0000
        // B(6,:):  0.0000  12.500
        // B(7,:):  0.0000   0.0000
        B_io = {
            // Column 0 (B(:,1) from data file)
            0.00000, 12.50000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,
            // Column 1 (B(:,2) from data file)
            0.00000,  0.00000,  0.00000,  0.00000,  0.00000, 12.50000,  0.00000
        };
        
        // Fortran Example Program Data for C (3x7)
        // C(1,:):  1.0000   0.0000  0.0000   0.0000  0.0000  0.0000  0.0000
        // C(2,:):  0.0000   0.0000  0.0000   1.0000  0.0000  0.0000  0.0000
        // C(3,:):  0.0000   0.0000  0.0000   0.0000  1.0000  0.0000  0.0000
        C_io = {
            // Column 0 (C(:,1) from data file)
            1.00000,  0.00000,  0.00000,
            // Column 1 (C(:,2) from data file)
            0.00000,  0.00000,  0.00000,
            // Column 2 (C(:,3) from data file)
            0.00000,  0.00000,  0.00000,
            // Column 3 (C(:,4) from data file)
            0.00000,  1.00000,  0.00000,
            // Column 4 (C(:,5) from data file)
            0.00000,  0.00000,  1.00000,
            // Column 5 (C(:,6) from data file)
            0.00000,  0.00000,  0.00000,
            // Column 6 (C(:,7) from data file)
            0.00000,  0.00000,  0.00000
        };

        HSV_out.resize(N_in); 
        // Expected HSV from Slycot/Fortran example output
        HSV_expected = {1.9178, 0.8621, 0.7666, 0.0336, 0.0246}; 
        HSV_expected.resize(N_in, 0.0); 

        // Expected Ar (5x5) from Slycot/Fortran example output, in Column-Major Order
        // Fortran prints row-wise:
        //  -0.5181  -1.1084   0.0000   0.0000   0.0000
        //   8.8157  -0.5181   0.0000   0.0000   0.0000
        //   0.0000   0.0000   0.5124   0.0000   1.7910
        //   0.0000   0.0000   0.0000  -1.4460   0.0000
        //   0.0000   0.0000  -4.2167   0.0000  -2.9900
        A_expected_reduced = {
            // Col 0
            -0.5181,  8.8157,  0.0000,  0.0000,  0.0000,
            // Col 1
            -1.1084, -0.5181,  0.0000,  0.0000,  0.0000,
            // Col 2
             0.0000,  0.0000,  0.5124,  0.0000, -4.2167,
            // Col 3
             0.0000,  0.0000,  0.0000, -1.4460,  0.0000,
            // Col 4
             0.0000,  0.0000,  1.7910,  0.0000, -2.9900
        };

        // Expected Br (5x2) from Slycot/Fortran example output, in Column-Major Order
        // Fortran prints row-wise:
        //  -1.2837   1.2837
        //  -0.7522   0.7522
        //  -0.7447  -0.7447
        //   1.9275  -1.9275
        //  -3.6872  -3.6872
        B_expected_reduced = {
            // Col 0
            -1.2837, -0.7522, -0.7447,  1.9275, -3.6872,
            // Col 1
             1.2837,  0.7522, -0.7447, -1.9275, -3.6872
        };
        
        // Expected Cr (3x5) from Slycot/Fortran example output, in Column-Major Order
        // Fortran prints row-wise:
        //  -0.1380  -0.6445  -0.6582  -0.5771   0.2222
        //   0.6246   0.0196   0.0000   0.4131   0.0000
        //   0.1380   0.6445  -0.6582   0.5771   0.2222
        C_expected_reduced = {
            // Col 0
            -0.1380,  0.6246,  0.1380,
            // Col 1
            -0.6445,  0.0196,  0.6445,
            // Col 2
            -0.6582,  0.0000, -0.6582,
            // Col 3
            -0.5771,  0.4131,  0.5771,
            // Col 4
             0.2222,  0.0000,  0.2222
        };
    }
     void SetUp() override { 
        SetUpBase();
     }
};

// --- Row-Major Test Fixture ---
class AB09MDTestRowMajor : public AB09MDTestColMajor {
protected:
    std::vector<double> A_rm_io;
    std::vector<double> B_rm_io;
    std::vector<double> C_rm_io;
    
    std::vector<double> A_expected_reduced_rm;
    std::vector<double> B_expected_reduced_rm;
    std::vector<double> C_expected_reduced_rm;

    void SetUp() override {
        AB09MDTestColMajor::SetUpBase(); 

        int lda_c_rm = N_in; 
        int ldb_c_rm = M_in; 
        int ldc_c_rm = N_in; 

        if (N_in > 0) A_rm_io.resize((size_t)N_in * lda_c_rm); else A_rm_io.clear();
        if (N_in > 0 && M_in > 0) B_rm_io.resize((size_t)N_in * ldb_c_rm); else B_rm_io.clear();
        if (P_in > 0 && N_in > 0) C_rm_io.resize((size_t)P_in * ldc_c_rm); 
        else if (P_in > 0 && N_in == 0) C_rm_io.resize((size_t)P_in * std::max(1,ldc_c_rm)); 
        else C_rm_io.clear();

        if (N_in > 0 && !A_io.empty()) {
            slicot_transpose_to_c_with_ld(A_io.data(), A_rm_io.data(), N_in, N_in, 
                                          std::max(1,N_in) , 
                                          lda_c_rm , sizeof(double));
        }
        if (N_in > 0 && M_in > 0 && !B_io.empty()) {
            slicot_transpose_to_c_with_ld(B_io.data(), B_rm_io.data(), N_in, M_in, 
                                          std::max(1,N_in) , 
                                          ldb_c_rm , sizeof(double));
        }
        if (P_in > 0 && N_in > 0 && !C_io.empty()) {
            slicot_transpose_to_c_with_ld(C_io.data(), C_rm_io.data(), P_in, N_in, 
                                          std::max(1,P_in) , 
                                          ldc_c_rm , sizeof(double));
        }
        
        if (expected_nr_auto > 0 && !A_expected_reduced.empty()) {
            A_expected_reduced_rm.resize((size_t)expected_nr_auto * expected_nr_auto);
            slicot_transpose_to_c_with_ld(A_expected_reduced.data(), A_expected_reduced_rm.data(), 
                                          expected_nr_auto, expected_nr_auto, 
                                          std::max(1,expected_nr_auto), 
                                          std::max(1,expected_nr_auto), sizeof(double));
        }
        if (expected_nr_auto > 0 && M_in > 0 && !B_expected_reduced.empty()) {
            B_expected_reduced_rm.resize((size_t)expected_nr_auto * M_in);
            slicot_transpose_to_c_with_ld(B_expected_reduced.data(), B_expected_reduced_rm.data(), 
                                          expected_nr_auto, M_in, 
                                          std::max(1,expected_nr_auto), 
                                          std::max(1,M_in), sizeof(double));
        }
        if (P_in > 0 && expected_nr_auto > 0 && !C_expected_reduced.empty()) {
            C_expected_reduced_rm.resize((size_t)P_in * expected_nr_auto);
            slicot_transpose_to_c_with_ld(C_expected_reduced.data(), C_expected_reduced_rm.data(), 
                                          P_in, expected_nr_auto, 
                                          std::max(1,P_in), 
                                          std::max(1,expected_nr_auto), sizeof(double));
        }
    }
};


TEST_F(AB09MDTestColMajor, DocExample) {
    std::vector<double> A_run = A_io; 
    std::vector<double> B_run = B_io;
    std::vector<double> C_run = C_io;
    int nr_run = NR_io_auto_order; 
    
    info_result = slicot_ab09md(
        DICO_in, JOB_in, EQUIL_in, ORDSEL_in, 
        N_in, M_in, P_in, &nr_run, ALPHA_in,
        A_run.data(), LDA_io, 
        B_run.data(), LDB_io, 
        C_run.data(), LDC_io, 
        &NS_out, HSV_out.data(), TOL_in, 
        &IWARN_out,
        0 // col_major
    );

    ASSERT_EQ(info_result, expected_info) << "slicot_ab09md call failed with info = " << info_result;
    EXPECT_EQ(nr_run, expected_nr_auto); 
    EXPECT_EQ(NS_out, expected_ns);
    EXPECT_EQ(IWARN_out, expected_iwarn);

    for (int i = 0; i < NS_out; ++i) { 
        EXPECT_NEAR(HSV_out[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
    }

    for (int j_col = 0; j_col < nr_run; ++j_col) { 
        for (int i_row = 0; i_row < nr_run; ++i_row) { 
            EXPECT_NEAR(A_run[j_col*LDA_io + i_row], A_expected_reduced[j_col*expected_nr_auto + i_row], check_tol)
                << "A mismatch at (" << i_row << "," << j_col << ")";
        }
    }
    for (int j_col = 0; j_col < M_in; ++j_col) {
        for (int i_row = 0; i_row < nr_run; ++i_row) {
            EXPECT_NEAR(B_run[j_col*LDB_io + i_row], B_expected_reduced[j_col*expected_nr_auto + i_row], check_tol)
                << "B mismatch at (" << i_row << "," << j_col << ")";
        }
    }
    for (int j_col = 0; j_col < nr_run; ++j_col) {
        for (int i_row = 0; i_row < P_in; ++i_row) {
            EXPECT_NEAR(C_run[j_col*LDC_io + i_row], C_expected_reduced[j_col*P_in + i_row], check_tol)
                << "C mismatch at (" << i_row << "," << j_col << ")";
        }
    }
}

TEST_F(AB09MDTestRowMajor, DocExample) {
    std::vector<double> A_run_rm = A_rm_io;
    std::vector<double> B_run_rm = B_rm_io;
    std::vector<double> C_run_rm = C_rm_io;
    int nr_run = NR_io_auto_order; 
    
    int lda_c_rm_call = N_in;
    int ldb_c_rm_call = M_in;
    int ldc_c_rm_call = N_in;

    info_result = slicot_ab09md(
        DICO_in, JOB_in, EQUIL_in, ORDSEL_in, 
        N_in, M_in, P_in, &nr_run, ALPHA_in,
        A_run_rm.empty() ? nullptr : A_run_rm.data(), lda_c_rm_call, 
        B_run_rm.empty() ? nullptr : B_run_rm.data(), ldb_c_rm_call, 
        C_run_rm.empty() ? nullptr : C_run_rm.data(), ldc_c_rm_call, 
        &NS_out, HSV_out.empty() ? nullptr : HSV_out.data(), TOL_in, 
        &IWARN_out,
        1 // row_major
    );

    ASSERT_EQ(info_result, expected_info) << "slicot_ab09md call failed with info = " << info_result;
    EXPECT_EQ(nr_run, expected_nr_auto);
    EXPECT_EQ(NS_out, expected_ns);
    EXPECT_EQ(IWARN_out, expected_iwarn);

    if (NS_out > 0 && !HSV_out.empty() && !HSV_expected.empty()) { 
        for (int i = 0; i < NS_out; ++i) {
            EXPECT_NEAR(HSV_out[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
        }
    }

    if (nr_run > 0) { 
        if (!A_run_rm.empty() && !A_expected_reduced_rm.empty()) {
            for (int i_row = 0; i_row < nr_run; ++i_row) { 
                for (int j_col = 0; j_col < nr_run; ++j_col) { 
                    EXPECT_NEAR(A_run_rm[i_row*lda_c_rm_call + j_col], A_expected_reduced_rm[i_row*expected_nr_auto + j_col], check_tol)
                        << "A_rm mismatch at (" << i_row << "," << j_col << ")";
                }
            }
        }
        if (M_in > 0 && !B_run_rm.empty() && !B_expected_reduced_rm.empty()) {
            for (int i_row = 0; i_row < nr_run; ++i_row) {
                for (int j_col = 0; j_col < M_in; ++j_col) {
                    EXPECT_NEAR(B_run_rm[i_row*ldb_c_rm_call + j_col], B_expected_reduced_rm[i_row*M_in + j_col], check_tol)
                        << "B_rm mismatch at (" << i_row << "," << j_col << ")";
                }
            }
        }
        if (P_in > 0 && !C_run_rm.empty() && !C_expected_reduced_rm.empty()) {
            for (int i_row = 0; i_row < P_in; ++i_row) {
                for (int j_col = 0; j_col < nr_run; ++j_col) {
                    EXPECT_NEAR(C_run_rm[i_row*ldc_c_rm_call + j_col], C_expected_reduced_rm[i_row*expected_nr_auto + j_col], check_tol)
                        << "C_rm mismatch at (" << i_row << "," << j_col << ")";
                }
            }
        }
    }
}

TEST_F(AB09MDTestColMajor, ParameterValidation) {
    std::vector<double> dummy_A(1);
    std::vector<double> dummy_B(1);
    std::vector<double> dummy_C(1);
    std::vector<double> dummy_HSV(1);
    int dummy_nr = 0;
    int dummy_ns = 0;
    int dummy_iwarn = 0;
    int N_val = 1, M_val = 1, P_val = 1;
    int LDA_val = 1, LDB_val = 1, LDC_val = 1; 

    int result = slicot_ab09md(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, -1, M_val, P_val, &dummy_nr, ALPHA_in, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, &dummy_ns, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(result, -5); 

    result = slicot_ab09md('X', JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, ALPHA_in, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, &dummy_ns, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(result, -1); 

    dummy_nr = N_val + 1; 
    result = slicot_ab09md(DICO_in, JOB_in, EQUIL_in, 'F', N_val, M_val, P_val, &dummy_nr, ALPHA_in, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, &dummy_ns, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(result, -8); 
    dummy_nr = 0; 

    if (N_val > 0) {
        result = slicot_ab09md(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, ALPHA_in, dummy_A.data(), 0, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, &dummy_ns, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
        EXPECT_EQ(result, -11);
    }
    result = slicot_ab09md('C', JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, 0.5, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, &dummy_ns, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
    EXPECT_EQ(result, -9);
}

TEST_F(AB09MDTestColMajor, ZeroDimension) {
    int zero_n = 0;
    int zero_m = 1; 
    int zero_p = 1; 
    int zero_nr = 0; 
    int zero_ns = 0;
    int zero_iwarn = 0;
    
    int lda_z = 1;
    int ldb_z = 1; 
    int ldc_z = std::max(1, zero_p); 
    
    info_result = slicot_ab09md(
        DICO_in, JOB_in, EQUIL_in, 'F', 
        zero_n, zero_m, zero_p, &zero_nr, ALPHA_in,
        nullptr, lda_z, 
        nullptr, ldb_z, 
        nullptr, ldc_z, 
        &zero_ns, nullptr, TOL_in, 
        &zero_iwarn,
        0 // col_major
    );
    
    EXPECT_EQ(info_result, 0); 
    EXPECT_EQ(zero_nr, 0);  
    EXPECT_EQ(zero_ns, 0);
}
