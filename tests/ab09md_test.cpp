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
     char EQUIL_in = 'N';
     char ORDSEL_in = 'A';
     int N_in = 7;
     int M_in = 2;
     int P_in = 3;
     int NR_io_fixed_order = 5; 
     int NR_io_auto_order = 0;  
     double ALPHA_in = -0.6;
     double TOL_in = 1.0e-1;
     
     int IWARN_out = 0;
     int info_result = -999;
     int NS_out = 0;
 
     double check_tol = 1e-4; 
 
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
     int expected_nr_auto = 5; 
     int expected_ns = 5;      
 
     int LDA_io = 0;
     int LDB_io = 0;
     int LDC_io = 0;
 
     void SetUpBase() {
         LDA_io = std::max(1, N_in);
         LDB_io = std::max(1, N_in); 
         LDC_io = std::max(1, P_in); 
 
         A_io = { 
            -0.04165, -5.21000,  0.00000,  0.54500,  0.00000,  0.00000,  0.00000,
             0.00000, -12.50000,  3.33000,  0.00000,  0.00000,  0.00000,  0.00000,
             4.92000,  0.00000, -3.33000,  0.00000,  0.00000,  0.00000,  0.00000,
            -4.92000,  0.00000,  0.00000,  0.00000,  4.92000,  0.00000,  0.00000,
             0.00000,  0.00000,  0.00000, -0.54500, -0.04165, -5.21000,  0.00000,
             0.00000,  0.00000,  0.00000,  0.00000,  0.00000, -12.50000,  3.33000,
             0.00000,  0.00000,  0.00000,  0.00000,  4.92000,  0.00000, -3.33000
         };
         A_io.resize((size_t)LDA_io * N_in);
 
         B_io = { 
             0.00000, 12.50000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,
             0.00000,  0.00000,  0.00000,  0.00000,  0.00000, 12.50000,  0.00000
         };
         B_io.resize((size_t)LDB_io * M_in);
 
         C_io = {
             1.00000,  0.00000,  0.00000, 
             0.00000,  0.00000,  0.00000, 
             0.00000,  0.00000,  0.00000, 
             0.00000,  1.00000,  0.00000, 
             0.00000,  0.00000,  1.00000, 
             0.00000,  0.00000,  0.00000, 
             0.00000,  0.00000,  0.00000  
         };
         C_io.resize((size_t)LDC_io * N_in);
 
         HSV_out.resize(N_in);
         HSV_expected = {1.9178, 0.8621, 0.7666, 0.0336, 0.0246}; 
         HSV_expected.resize(N_in, 0.0); 
 
         A_expected_reduced = { 
             -0.5181,  8.8157,  0.0000,  0.0000,  0.0000,
             -1.1084, -0.5181,  0.0000,  0.0000,  0.0000,
              0.0000,  0.0000,  0.5124,  0.0000, -4.2167,
              0.0000,  0.0000,  0.0000, -1.4460,  0.0000,
              0.0000,  0.0000,  1.7910,  0.0000, -2.9900
         };
         A_expected_reduced.resize((size_t)expected_nr_auto * expected_nr_auto);
 
         B_expected_reduced = { 
             -1.2837, -0.7522, -0.7447,  1.9275, -3.6872,
              1.2837,  0.7522, -0.7447, -1.9275, -3.6872
         };
         B_expected_reduced.resize((size_t)expected_nr_auto * M_in);
         
         C_expected_reduced = { 
             -0.1380,  0.6246,  0.1380,
             -0.6445,  0.0196,  0.6445,
             -0.6582,  0.0000, -0.6582,
             -0.5771,  0.4131,  0.5771,
              0.2222,  0.0000,  0.2222
         };
         C_expected_reduced.resize((size_t)P_in * expected_nr_auto);
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
 
         LDA_io = std::max(1, N_in); 
         LDB_io = std::max(1, M_in); 
         LDC_io = std::max(1, N_in); 
 
         A_rm_io.resize((size_t)N_in * LDA_io);
         B_rm_io.resize((size_t)N_in * LDB_io); 
         C_rm_io.resize((size_t)P_in * LDC_io); 
 
         slicot_transpose_to_c_with_ld(A_io.data(), A_rm_io.data(), N_in, N_in, std::max(1,N_in), LDA_io, sizeof(double));
         slicot_transpose_to_c_with_ld(B_io.data(), B_rm_io.data(), N_in, M_in, std::max(1,N_in), LDB_io, sizeof(double));
         slicot_transpose_to_c_with_ld(C_io.data(), C_rm_io.data(), P_in, N_in, std::max(1,P_in), LDC_io, sizeof(double));
         
         A_expected_reduced_rm.resize((size_t)expected_nr_auto * expected_nr_auto);
         B_expected_reduced_rm.resize((size_t)expected_nr_auto * M_in);
         C_expected_reduced_rm.resize((size_t)P_in * expected_nr_auto);
 
         slicot_transpose_to_c_with_ld(A_expected_reduced.data(), A_expected_reduced_rm.data(), expected_nr_auto, expected_nr_auto, std::max(1,expected_nr_auto), std::max(1,expected_nr_auto), sizeof(double));
         slicot_transpose_to_c_with_ld(B_expected_reduced.data(), B_expected_reduced_rm.data(), expected_nr_auto, M_in, std::max(1,expected_nr_auto), std::max(1,M_in), sizeof(double));
         slicot_transpose_to_c_with_ld(C_expected_reduced.data(), C_expected_reduced_rm.data(), P_in, expected_nr_auto, std::max(1,P_in), std::max(1,expected_nr_auto), sizeof(double));
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
         0 
     );
 
     ASSERT_EQ(info_result, expected_info) << "slicot_ab09md call failed with info = " << info_result;
     EXPECT_EQ(nr_run, expected_nr_auto); 
     EXPECT_EQ(NS_out, expected_ns);
     EXPECT_EQ(IWARN_out, expected_iwarn);
 
     for (int i = 0; i < NS_out; ++i) { 
         EXPECT_NEAR(HSV_out[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
     }
 
     for (int j = 0; j < nr_run; ++j) { 
         for (int i = 0; i < nr_run; ++i) { 
             EXPECT_NEAR(A_run[j*LDA_io + i], A_expected_reduced[j*expected_nr_auto + i], check_tol)
                 << "A mismatch at (" << i << "," << j << ")";
         }
     }
     for (int j = 0; j < M_in; ++j) {
         for (int i = 0; i < nr_run; ++i) {
             EXPECT_NEAR(B_run[j*LDB_io + i], B_expected_reduced[j*expected_nr_auto + i], check_tol)
                 << "B mismatch at (" << i << "," << j << ")";
         }
     }
     for (int j = 0; j < nr_run; ++j) {
         for (int i = 0; i < P_in; ++i) {
             EXPECT_NEAR(C_run[j*LDC_io + i], C_expected_reduced[j*P_in + i], check_tol)
                 << "C mismatch at (" << i << "," << j << ")";
         }
     }
 }
 
 TEST_F(AB09MDTestRowMajor, DocExample) {
     std::vector<double> A_run_rm = A_rm_io;
     std::vector<double> B_run_rm = B_rm_io;
     std::vector<double> C_run_rm = C_rm_io;
     int nr_run = NR_io_auto_order;
     
     info_result = slicot_ab09md(
         DICO_in, JOB_in, EQUIL_in, ORDSEL_in, 
         N_in, M_in, P_in, &nr_run, ALPHA_in,
         A_run_rm.data(), LDA_io, 
         B_run_rm.data(), LDB_io, 
         C_run_rm.data(), LDC_io, 
         &NS_out, HSV_out.data(), TOL_in, 
         &IWARN_out,
         1 
     );
 
     ASSERT_EQ(info_result, expected_info) << "slicot_ab09md call failed with info = " << info_result;
     EXPECT_EQ(nr_run, expected_nr_auto);
     EXPECT_EQ(NS_out, expected_ns);
     EXPECT_EQ(IWARN_out, expected_iwarn);
 
     for (int i = 0; i < NS_out; ++i) {
         EXPECT_NEAR(HSV_out[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
     }
 
     for (int i = 0; i < nr_run; ++i) { 
         for (int j = 0; j < nr_run; ++j) { 
             EXPECT_NEAR(A_run_rm[i*LDA_io + j], A_expected_reduced_rm[i*expected_nr_auto + j], check_tol)
                 << "A_rm mismatch at (" << i << "," << j << ")";
         }
     }
     for (int i = 0; i < nr_run; ++i) {
         for (int j = 0; j < M_in; ++j) {
             EXPECT_NEAR(B_run_rm[i*LDB_io + j], B_expected_reduced_rm[i*M_in + j], check_tol)
                 << "B_rm mismatch at (" << i << "," << j << ")";
         }
     }
     for (int i = 0; i < P_in; ++i) {
         for (int j = 0; j < nr_run; ++j) {
             EXPECT_NEAR(C_run_rm[i*LDC_io + j], C_expected_reduced_rm[i*expected_nr_auto + j], check_tol)
                 << "C_rm mismatch at (" << i << "," << j << ")";
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
 
     if (N_val > 0) {
         result = slicot_ab09md(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, ALPHA_in, dummy_A.data(), 0, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, &dummy_ns, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
         EXPECT_EQ(result, -11);
     }
     result = slicot_ab09md('C', JOB_in, EQUIL_in, ORDSEL_in, N_val, M_val, P_val, &dummy_nr, 0.5, dummy_A.data(), LDA_val, dummy_B.data(), LDB_val, dummy_C.data(), LDC_val, &dummy_ns, dummy_HSV.data(), TOL_in, &dummy_iwarn, 0);
     EXPECT_EQ(result, -9);
 }
 
 TEST_F(AB09MDTestColMajor, ZeroDimension) {
     int zero_n = 0;
     int zero_m = M_in; 
     int zero_p = P_in;
     int zero_nr = 0;
     int zero_ns = 0;
     int zero_iwarn = 0;
     // For N=0: LDA, LDB can be 1. LDC must be >= MAX(1,P).
     int lda_z = 1;
     int ldb_z = 1; 
     int ldc_z = std::max(1, zero_p); // Correct LDC for P > 0, N = 0
     
     info_result = slicot_ab09md(
         DICO_in, JOB_in, EQUIL_in, 'F', 
         zero_n, zero_m, zero_p, &zero_nr, ALPHA_in,
         nullptr, lda_z, 
         nullptr, ldb_z, 
         nullptr, ldc_z, 
         &zero_ns, nullptr, TOL_in, 
         &zero_iwarn,
         0 
     );
     
     EXPECT_EQ(info_result, 0); 
     EXPECT_EQ(zero_nr, 0);  
     EXPECT_EQ(zero_ns, 0);
 }
 