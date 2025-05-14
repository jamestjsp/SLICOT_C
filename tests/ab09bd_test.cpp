#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max, std::min
#include <stdexcept> // For std::runtime_error
#include <iostream>  // For debugging output if needed
#include <iomanip>   // For std::fixed and std::setprecision

#include "ab09bd.h"
#include "slicot_utils.h" 

// --- Column-Major Test Fixture ---
class AB09BDTestColMajor : public ::testing::Test {
protected:
    // Test parameters from AB09BD.html example
    char DICO_in = 'C';    // Continuous-time system
    char JOB_in = 'N';     // Balancing-free square-root SPA method (Matches HTML example for results)
    char EQUIL_in = 'N';   // No equilibration (Matches HTML example)
    char ORDSEL_in = 'F';  // Fixed order reduction for FixedOrderReduction test
    int N_in = 7;          // Order of the original state-space
    int M_in = 2;          // Number of inputs
    int P_in = 3;          // Number of outputs
    int NR_io = 5;         // Desired reduced order (fixed for this example, matches HTML example's outcome)
    
    // Tolerances from HTML example (TOL1 used by ORDSEL='A', TOL2 for minimal realization check)
    // For ORDSEL='F', these are not strictly used for order selection by AB09BD, 
    // but passed to the routine. The example uses 1.E-1 and 1.E-14.
    double TOL1_in = 1.0e-1; 
    double TOL2_in = 1.0e-14; 
    
    int IWARN_out = 0;     // Warning indicator output
    int info_result = -999;// Result from slicot_ab09bd call

    // Verification tolerance
    double check_tol = 1e-4; // Relaxed tolerance

    // Input data vectors (column-major format)
    std::vector<double> A_io; 
    std::vector<double> B_io; 
    std::vector<double> C_io; 
    std::vector<double> D_io; 
    std::vector<double> HSV_out;

    // Expected results (from AB09BD.html documentation example for JOB='N', EQUIL='N', ORDSEL='A' giving NR=5)
    std::vector<double> A_expected_reduced;   
    std::vector<double> B_expected_reduced;   
    std::vector<double> C_expected_reduced;   
    std::vector<double> D_expected_reduced; 
    std::vector<double> HSV_expected;     
    int expected_info = 0;
    int expected_iwarn = 0; 
    int expected_nr = 5; 

    // Leading dimensions for C arrays
    int LDA_io = 0;
    int LDB_io = 0;
    int LDC_io = 0;
    int LDD_io = 0;


    void SetUp() override {
        // Set leading dimensions for column-major format (Fortran style)
        LDA_io = std::max(1, N_in);
        LDB_io = std::max(1, N_in); 
        LDC_io = std::max(1, P_in); 
        LDD_io = std::max(1, P_in);

        // Initialize matrices from AB09BD.html example (Column-Major)
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

        D_io = { 
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0
        };
        D_io.resize((size_t)LDD_io * M_in);

        HSV_out.resize(N_in);

        HSV_expected = {2.5139, 2.0846, 1.9178, 0.7666, 0.5473, 0.0253, 0.0246};
        HSV_expected.resize(N_in); 

        A_expected_reduced = { // For JOB='N', NR=5 from AB09BD.html
             1.3960, -4.1411,  0.0000,  0.0000,  1.3261,
             5.1248, -3.8605,  0.0000,  0.0000,  1.7851,
             0.0000,  0.0000,  0.5847, -4.3823,  0.0000,
             0.0000,  0.0000,  1.9230, -3.2922,  0.0000,
             4.4331, -0.6738,  0.0000,  0.0000, -0.2249
        };
        A_expected_reduced.resize((size_t)expected_nr * expected_nr);

        B_expected_reduced = { // For JOB='N', NR=5
            -0.2901, -3.4004, -0.6379, -3.9315,  1.9813,
             0.2901,  3.4004, -0.6379, -3.9315, -1.9813
        };
        B_expected_reduced.resize((size_t)expected_nr * M_in);
        
        C_expected_reduced = { // For JOB='N', NR=5
            -0.6570,  0.1094,  0.6570,
             0.2053,  0.4875, -0.2053,
            -0.6416,  0.0000, -0.6416,
             0.2526,  0.0000,  0.2526,
            -0.0364,  0.8641,  0.0364
        };
        C_expected_reduced.resize((size_t)P_in * expected_nr);

        D_expected_reduced = { // For JOB='N'
             0.0498,  0.0010, -0.0007,
            -0.0007, -0.0010,  0.0498
        };
        D_expected_reduced.resize((size_t)P_in * M_in);
    }
};

// --- Row-Major Test Fixture ---
class AB09BDTestRowMajor : public AB09BDTestColMajor {
protected:
    std::vector<double> A_rm_io;
    std::vector<double> B_rm_io;
    std::vector<double> C_rm_io;
    std::vector<double> D_rm_io;
    
    std::vector<double> A_expected_reduced_rm;
    std::vector<double> B_expected_reduced_rm;
    std::vector<double> C_expected_reduced_rm;
    std::vector<double> D_expected_reduced_rm;

    void SetUp() override {
        AB09BDTestColMajor::SetUp(); 

        LDA_io = std::max(1, N_in); 
        LDB_io = std::max(1, M_in); 
        LDC_io = std::max(1, N_in); 
        LDD_io = std::max(1, M_in);

        A_rm_io.resize((size_t)N_in * LDA_io);
        B_rm_io.resize((size_t)N_in * LDB_io); 
        C_rm_io.resize((size_t)P_in * LDC_io); 
        D_rm_io.resize((size_t)P_in * LDD_io);

        slicot_transpose_to_c_with_ld(A_io.data(), A_rm_io.data(), N_in, N_in, std::max(1,N_in), LDA_io, sizeof(double));
        slicot_transpose_to_c_with_ld(B_io.data(), B_rm_io.data(), N_in, M_in, std::max(1,N_in), LDB_io, sizeof(double));
        slicot_transpose_to_c_with_ld(C_io.data(), C_rm_io.data(), P_in, N_in, std::max(1,P_in), LDC_io, sizeof(double));
        slicot_transpose_to_c_with_ld(D_io.data(), D_rm_io.data(), P_in, M_in, std::max(1,P_in), LDD_io, sizeof(double));
        
        A_expected_reduced_rm.resize((size_t)expected_nr * expected_nr);
        B_expected_reduced_rm.resize((size_t)expected_nr * M_in);
        C_expected_reduced_rm.resize((size_t)P_in * expected_nr);
        D_expected_reduced_rm.resize((size_t)P_in * M_in);

        slicot_transpose_to_c_with_ld(A_expected_reduced.data(), A_expected_reduced_rm.data(), expected_nr, expected_nr, std::max(1,expected_nr), std::max(1,expected_nr), sizeof(double));
        slicot_transpose_to_c_with_ld(B_expected_reduced.data(), B_expected_reduced_rm.data(), expected_nr, M_in, std::max(1,expected_nr), std::max(1,M_in), sizeof(double));
        slicot_transpose_to_c_with_ld(C_expected_reduced.data(), C_expected_reduced_rm.data(), P_in, expected_nr, std::max(1,P_in), std::max(1,expected_nr), sizeof(double));
        slicot_transpose_to_c_with_ld(D_expected_reduced.data(), D_expected_reduced_rm.data(), P_in, M_in, std::max(1,P_in), std::max(1,M_in), sizeof(double));
    }
};

// --- Test Cases ---

// Test: Fixed order reduction (Column-Major) - Using Documentation Example
TEST_F(AB09BDTestColMajor, FixedOrderReduction) {
    std::vector<double> A_copy = A_io;
    std::vector<double> B_copy = B_io;
    std::vector<double> C_copy = C_io;
    std::vector<double> D_copy = D_io;
    int nr_run = NR_io; // Should be 5 for this test to match expected data
    
    // Use JOB = 'N' to match the documentation example's expected output
    info_result = slicot_ab09bd(
        DICO_in, 'N', EQUIL_in, ORDSEL_in, // JOB_in is 'N'
        N_in, M_in, P_in, &nr_run, 
        A_copy.data(), LDA_io,
        B_copy.data(), LDB_io,
        C_copy.data(), LDC_io,
        D_copy.data(), LDD_io,
        HSV_out.data(),
        TOL1_in, TOL2_in, 
        &IWARN_out,
        0 /* column-major */
    );

    ASSERT_EQ(info_result, expected_info) << "slicot_ab09bd call failed with info = " << info_result;
    EXPECT_EQ(nr_run, expected_nr); 
    EXPECT_EQ(IWARN_out, expected_iwarn);

    for (int i = 0; i < N_in; ++i) {
        EXPECT_NEAR(HSV_out[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
    }

    for (int j = 0; j < expected_nr; ++j) { 
        for (int i = 0; i < expected_nr; ++i) { 
            EXPECT_NEAR(A_copy[j*LDA_io + i], A_expected_reduced[j*expected_nr + i], check_tol)
                << "A mismatch at (" << i << "," << j << ")";
        }
    }
    for (int j = 0; j < M_in; ++j) {
        for (int i = 0; i < expected_nr; ++i) {
            EXPECT_NEAR(B_copy[j*LDB_io + i], B_expected_reduced[j*expected_nr + i], check_tol)
                << "B mismatch at (" << i << "," << j << ")";
        }
    }
    for (int j = 0; j < expected_nr; ++j) {
        for (int i = 0; i < P_in; ++i) {
            EXPECT_NEAR(C_copy[j*LDC_io + i], C_expected_reduced[j*P_in + i], check_tol)
                << "C mismatch at (" << i << "," << j << ")";
        }
    }
    for (int j = 0; j < M_in; ++j) {
        for (int i = 0; i < P_in; ++i) {
            EXPECT_NEAR(D_copy[j*LDD_io + i], D_expected_reduced[j*P_in + i], check_tol)
                << "D mismatch at (" << i << "," << j << ")";
        }
    }
}

// Test: Auto order selection (Column-Major)
TEST_F(AB09BDTestColMajor, AutoOrderSelection) {
    std::vector<double> A_copy = A_io;
    std::vector<double> B_copy = B_io;
    std::vector<double> C_copy = C_io;
    std::vector<double> D_copy = D_io;
    int nr_auto_run = 0; 
    // Use TOL1 and TOL2 from the documentation example for ORDSEL='A'
    double tol1_auto = 1.0e-1; 
    double tol2_auto = 1.0e-14; 
    
    info_result = slicot_ab09bd(
        DICO_in, 'N', EQUIL_in, 'A', // JOB='N', ORDSEL = 'A'
        N_in, M_in, P_in, &nr_auto_run, 
        A_copy.data(), LDA_io,
        B_copy.data(), LDB_io,
        C_copy.data(), LDC_io,
        D_copy.data(), LDD_io,
        HSV_out.data(),
        tol1_auto, tol2_auto, 
        &IWARN_out,
        0 /* column-major */
    );

    EXPECT_EQ(info_result, expected_info); 
    // The documentation example with these settings (JOB='N', ORDSEL='A', TOL1=0.1) results in NR=5
    EXPECT_EQ(nr_auto_run, 5); 
}

// Test: Fixed order reduction (Row-Major) - Using Documentation Example
TEST_F(AB09BDTestRowMajor, FixedOrderReduction) {
    std::vector<double> A_rm_copy = A_rm_io;
    std::vector<double> B_rm_copy = B_rm_io;
    std::vector<double> C_rm_copy = C_rm_io;
    std::vector<double> D_rm_copy = D_rm_io;
    int nr_run = NR_io; // Should be 5
    
    info_result = slicot_ab09bd(
        DICO_in, 'N', EQUIL_in, ORDSEL_in, // JOB_in is 'N'
        N_in, M_in, P_in, &nr_run,
        A_rm_copy.data(), LDA_io, 
        B_rm_copy.data(), LDB_io, 
        C_rm_copy.data(), LDC_io, 
        D_rm_copy.data(), LDD_io,
        HSV_out.data(),
        TOL1_in, TOL2_in, 
        &IWARN_out,
        1 /* row-major */
    );

    ASSERT_EQ(info_result, expected_info) << "slicot_ab09bd call failed with info = " << info_result;
    EXPECT_EQ(nr_run, expected_nr);
    EXPECT_EQ(IWARN_out, expected_iwarn);

    for (int i = 0; i < N_in; ++i) {
        EXPECT_NEAR(HSV_out[i], HSV_expected[i], check_tol) << "HSV mismatch at index " << i;
    }

    for (int i = 0; i < expected_nr; ++i) { 
        for (int j = 0; j < expected_nr; ++j) { 
            EXPECT_NEAR(A_rm_copy[i*LDA_io + j], A_expected_reduced_rm[i*expected_nr + j], check_tol)
                << "A_rm mismatch at (" << i << "," << j << ")";
        }
    }
    for (int i = 0; i < expected_nr; ++i) {
        for (int j = 0; j < M_in; ++j) {
            EXPECT_NEAR(B_rm_copy[i*LDB_io + j], B_expected_reduced_rm[i*M_in + j], check_tol)
                << "B_rm mismatch at (" << i << "," << j << ")";
        }
    }
    for (int i = 0; i < P_in; ++i) {
        for (int j = 0; j < expected_nr; ++j) {
            EXPECT_NEAR(C_rm_copy[i*LDC_io + j], C_expected_reduced_rm[i*expected_nr + j], check_tol)
                << "C_rm mismatch at (" << i << "," << j << ")";
        }
    }
     for (int i = 0; i < P_in; ++i) {
        for (int j = 0; j < M_in; ++j) {
            EXPECT_NEAR(D_rm_copy[i*LDD_io + j], D_expected_reduced_rm[i*M_in + j], check_tol)
                << "D_rm mismatch at (" << i << "," << j << ")";
        }
    }
}

// Test: Zero dimensions
TEST_F(AB09BDTestColMajor, ZeroDimensions) {
    int n_zero = 0, m_zero = 0, p_zero = 0;
    int nr_zero_io = 0;
    int lda_z = 1, ldb_z = 1, ldc_z = 1, ldd_z = 1;
    int iwarn_z = 0;
    
    info_result = slicot_ab09bd(
        DICO_in, JOB_in, EQUIL_in, 'F', 
        n_zero, m_zero, p_zero, &nr_zero_io, 
        nullptr, lda_z, 
        nullptr, ldb_z, 
        nullptr, ldc_z, 
        nullptr, ldd_z,
        nullptr,    
        TOL1_in, TOL2_in, 
        &iwarn_z,
        0 /* column-major */
    );
    
    EXPECT_EQ(info_result, 0); 
    EXPECT_EQ(nr_zero_io, 0);  
}

// Test: Parameter Validation (selected checks)
TEST_F(AB09BDTestColMajor, ParameterValidation) {
    int nr_copy = NR_io; 
    std::vector<double> a_val = A_io; 
    std::vector<double> b_val = B_io;
    std::vector<double> c_val = C_io;
    std::vector<double> d_val = D_io;
    std::vector<double> hsv_val = HSV_out;

    // Test invalid N (arg 5)
    info_result = slicot_ab09bd(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, -1, M_in, P_in, &nr_copy, a_val.data(), LDA_io, b_val.data(), LDB_io, c_val.data(), LDC_io, d_val.data(), LDD_io, hsv_val.data(), TOL1_in, TOL2_in, &IWARN_out, 0);
    EXPECT_EQ(info_result, -5);
    
    // Test invalid DICO (arg 1)
    info_result = slicot_ab09bd('X', JOB_in, EQUIL_in, ORDSEL_in, N_in, M_in, P_in, &nr_copy, a_val.data(), LDA_io, b_val.data(), LDB_io, c_val.data(), LDC_io, d_val.data(), LDD_io, hsv_val.data(), TOL1_in, TOL2_in, &IWARN_out, 0);
    EXPECT_EQ(info_result, -1);
        
    // Test invalid NR for ORDSEL='F' (arg 8)
    int nr_invalid = N_in + 1; 
    info_result = slicot_ab09bd(DICO_in, JOB_in, EQUIL_in, 'F', N_in, M_in, P_in, &nr_invalid, a_val.data(), LDA_io, b_val.data(), LDB_io, c_val.data(), LDC_io, d_val.data(), LDD_io, hsv_val.data(), TOL1_in, TOL2_in, &IWARN_out, 0);
    EXPECT_EQ(info_result, -8); 
    nr_invalid = -1;
    info_result = slicot_ab09bd(DICO_in, JOB_in, EQUIL_in, 'F', N_in, M_in, P_in, &nr_invalid, a_val.data(), LDA_io, b_val.data(), LDB_io, c_val.data(), LDC_io, d_val.data(), LDD_io, hsv_val.data(), TOL1_in, TOL2_in, &IWARN_out, 0);
    EXPECT_EQ(info_result, -8);

    // Test invalid LDA (arg 10)
    if (N_in > 0) {
        info_result = slicot_ab09bd(DICO_in, JOB_in, EQUIL_in, ORDSEL_in, N_in, M_in, P_in, &nr_copy, a_val.data(), 0, b_val.data(), LDB_io, c_val.data(), LDC_io, d_val.data(), LDD_io, hsv_val.data(), TOL1_in, TOL2_in, &IWARN_out, 0);
        EXPECT_EQ(info_result, -10);
    }
}
