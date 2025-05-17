#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max, std::min
#include <iostream> 
#include <iomanip>   

#include "sg03bd.h"       
#include "slicot_utils.h" 

// --- Test Fixture Base ---
class SG03BDTest : public ::testing::Test {
protected:
    char DICO_param = 'C';
    char FACT_param = 'N';
    char TRANS_param = 'N';
    int N_param = 0;
    int M_param = 0;

    // Input matrices (row-major by default for easy initialization)
    std::vector<double> A_data_rm, E_data_rm, Q_data_in_rm, Z_data_in_rm, B_in_data_rm;
    
    // Buffers passed to the C wrapper function - these are members of the fixture
    std::vector<double> A_test_buf, E_test_buf, Q_test_buf, Z_test_buf, B_U_test_buf;

    // Output scalars/vectors from C function
    double SCALE_out_val;
    std::vector<double> ALPHAR_out_val, ALFAI_out_val, BETA_out_val;

    // Expected results (row-major for easy comparison)
    std::vector<double> U_expected_rm;
    double SCALE_expected;
    int expected_info = 0;

    // Leading dimensions for C call
    int LDA, LDE, LDQ, LDZ, LDB;
    
    double check_tol_U = 1e-4; 

    // Sets up data based on the SLICOT SG03BD.html example
    void setup_slicot_example_data() {
        DICO_param = 'C'; FACT_param = 'N'; TRANS_param = 'N';
        N_param = 3; M_param = 1; 

        A_data_rm = {-1.0, 3.0, -4.0, 0.0, 5.0, -2.0, -4.0, 4.0, 1.0}; 
        E_data_rm = {2.0, 1.0, 3.0, 2.0, 0.0, 1.0, 4.0, 5.0, 1.0}; 
        B_in_data_rm = {2.0, -1.0, 7.0}; // M x N = 1 x 3 for TRANS='N'
        
        Q_data_in_rm.assign((size_t)N_param*N_param, 0.0); 
        Z_data_in_rm.assign((size_t)N_param*N_param, 0.0); 

        U_expected_rm = {1.6003, -0.4418, -0.1523,
                         0.0000,  0.6795, -0.2499,
                         0.0000,  0.0000,  0.2041};
        SCALE_expected = 1.0000;
        check_tol_U = 1e-4; 

        // Size output-only buffers based on current N_param
        // These are members, so they are sized here.
        if (N_param > 0) {
            ALPHAR_out_val.resize(N_param); 
            ALFAI_out_val.resize(N_param); 
            BETA_out_val.resize(N_param);
        } else { 
            ALPHAR_out_val.assign(1,0.0); ALFAI_out_val.assign(1,0.0); BETA_out_val.assign(1,0.0);
        }
        // Test buffers (A_test_buf etc.) will be sized in each test case
    }

    // Helper to prepare C call buffers for row-major tests
    void prepare_rowmajor_buffers_for_slicot_example() {
        LDA = N_param; LDE = N_param; LDQ = N_param; LDZ = N_param; 
        LDB = N_param; // C LDB (cols) for B_U_test_buf (NxN for U output) must be N_param
        
        A_test_buf = A_data_rm; 
        E_test_buf = E_data_rm; 
        Q_test_buf = Q_data_in_rm; // Will be filled by Fortran if FACT='N'
        Z_test_buf = Z_data_in_rm; // Will be filled by Fortran if FACT='N'
        
        // B_U_test_buf needs to be N_param x N_param for output U.
        // Input B_in_data_rm (M_param x N_param for TRANS='N') must be placed in its initial part.
        B_U_test_buf.assign((size_t)N_param * N_param, 0.0); // Initialize to 0
        if (M_param > 0 && N_param > 0 && TRANS_param == 'N') { 
            for(int i=0; i < M_param; ++i) { // Iterate rows of B_in (M_param)
                for(int j=0; j < N_param; ++j) { // Iterate columns of B_in (N_param)
                    // B_in_data_rm is M_param x N_param (row-major, C LDB is N_param)
                    // B_U_test_buf is N_param x N_param (row-major, C LDB is N_param)
                    if ((size_t)(i*N_param + j) < B_in_data_rm.size() && (size_t)(i*N_param+j) < B_U_test_buf.size())
                        B_U_test_buf[i*N_param + j] = B_in_data_rm[i*N_param + j];
                }
            }
        } else if (M_param > 0 && N_param > 0 && TRANS_param == 'T') { 
             // Input B_in_data_rm is N_param x M_param (row-major, C LDB is M_param)
             for(int i=0; i < N_param; ++i) { // Iterate rows of B_in (N_param)
                for(int j=0; j < M_param; ++j) { // Iterate columns of B_in (M_param)
                    // Copy B_in_data_rm[i][j] to B_U_test_buf[i][j] (conceptual)
                    // B_U_test_buf is N_param x N_param (row-major, C LDB is N_param)
                    if ((size_t)(i*M_param + j) < B_in_data_rm.size() && (size_t)(i*N_param+j) < B_U_test_buf.size()) 
                         B_U_test_buf[i*N_param + j] = B_in_data_rm[i*M_param + j]; 
                }
            }
        }
    }

    // Helper to prepare C call buffers for column-major tests
    void prepare_colmajor_buffers_for_slicot_example() {
        LDA = std::max(1,N_param); LDE = std::max(1,N_param); LDQ = std::max(1,N_param); 
        LDZ = std::max(1,N_param); 
        // Fortran LDB for B_U_test_buf. Input B is MxN (TRANS='N'), Output U is NxN.
        // So Fortran LDB must be >= max(M,N).
        LDB = std::max(1, std::max(M_param, N_param)); 

        A_test_buf.resize((size_t)LDA*N_param); slicot_transpose_to_fortran_with_ld(A_data_rm.data(), A_test_buf.data(), N_param, N_param, N_param, LDA, sizeof(double));
        E_test_buf.resize((size_t)LDE*N_param); slicot_transpose_to_fortran_with_ld(E_data_rm.data(), E_test_buf.data(), N_param, N_param, N_param, LDE, sizeof(double));
        Q_test_buf.resize((size_t)LDQ*N_param); // Will be filled by Fortran if FACT='N'
        Z_test_buf.resize((size_t)LDZ*N_param); // Will be filled by Fortran if FACT='N'
        
        // B_U_test_buf for column-major. Fortran B array is (LDB, N_param).
        B_U_test_buf.assign((size_t)LDB * N_param, 0.0); 
        if (M_param > 0 && N_param > 0) {
            std::vector<double> b_in_cm_temp; 
            if (TRANS_param == 'N') { // Input B is MxN (row-major in B_in_data_rm)
                b_in_cm_temp.resize((size_t)M_param * N_param);
                // Transpose B_in_data_rm (MxN, C_LDB=N) to b_in_cm_temp (MxN, Fortran_LDB=M)
                slicot_transpose_to_fortran_with_ld(B_in_data_rm.data(), b_in_cm_temp.data(), M_param, N_param, N_param, M_param, sizeof(double));
                // Copy b_in_cm_temp into B_U_test_buf (which is LDB x N_param Fortran)
                for(int j=0; j < N_param; ++j) { // Fortran columns (0 to N-1)
                    for(int i=0; i < M_param; ++i) { // Fortran rows (0 to M-1)
                        if ((size_t)(i + j*LDB) < B_U_test_buf.size() && (size_t)(i + j*M_param) < b_in_cm_temp.size())
                            B_U_test_buf[i + j*LDB] = b_in_cm_temp[i + j*M_param];
                    }
                }
            } else { // TRANS_param == 'T', Input B is NxM (row-major in B_in_data_rm)
                b_in_cm_temp.resize((size_t)N_param * M_param);
                // Transpose B_in_data_rm (NxM, C_LDB=M) to b_in_cm_temp (NxM, Fortran_LDB=N)
                slicot_transpose_to_fortran_with_ld(B_in_data_rm.data(), b_in_cm_temp.data(), N_param, M_param, M_param, N_param, sizeof(double));
                 for(int j=0; j < M_param; ++j) { // Fortran columns (0 to M-1)
                    for(int i=0; i < N_param; ++i) { // Fortran rows (0 to N-1)
                         if ((size_t)(i + j*LDB) < B_U_test_buf.size() && (size_t)(i + j*N_param) < b_in_cm_temp.size()) 
                            B_U_test_buf[i + j*LDB] = b_in_cm_temp[i + j*N_param]; 
                    }
                }
            }
        }
    }
};

TEST_F(SG03BDTest, SlicotDocExampleRowMajor) {
    setup_slicot_example_data();
    prepare_rowmajor_buffers_for_slicot_example();
    
    int info = slicot_sg03bd(DICO_param, FACT_param, TRANS_param, N_param, M_param,
                             A_test_buf.data(), LDA, E_test_buf.data(), LDE,
                             Q_test_buf.data(), LDQ, Z_test_buf.data(), LDZ,
                             B_U_test_buf.data(), LDB, 
                             &SCALE_out_val,
                             ALPHAR_out_val.data(), ALFAI_out_val.data(), BETA_out_val.data(),
                             1 /* row_major = true */);
    ASSERT_EQ(info, 0);
    // B_U_test_buf now contains U in row-major format
    for(size_t i=0; i < U_expected_rm.size(); ++i) {
        EXPECT_NEAR(B_U_test_buf[i], U_expected_rm[i], check_tol_U) << "U_out_rm[" << i << "]";
    }
    EXPECT_NEAR(SCALE_out_val, SCALE_expected, 1e-4);
}

TEST_F(SG03BDTest, SlicotDocExampleColMajor) {
    setup_slicot_example_data();
    prepare_colmajor_buffers_for_slicot_example();

    int info = slicot_sg03bd(DICO_param, FACT_param, TRANS_param, N_param, M_param,
                             A_test_buf.data(), LDA, E_test_buf.data(), LDE,
                             Q_test_buf.data(), LDQ, Z_test_buf.data(), LDZ,
                             B_U_test_buf.data(), LDB, 
                             &SCALE_out_val,
                             ALPHAR_out_val.data(), ALFAI_out_val.data(), BETA_out_val.data(),
                             0 /* row_major = false */);
    ASSERT_EQ(info, 0);
    
    // B_U_test_buf now contains U in column-major format.
    // Transpose expected U (row-major) to column-major for comparison.
    std::vector<double> U_expected_cm((size_t)N_param*N_param);
    slicot_transpose_to_fortran_with_ld(U_expected_rm.data(), U_expected_cm.data(), N_param,N_param,N_param,N_param,sizeof(double));

    for(int j=0; j < N_param; ++j) { // Iterate columns of U
        for(int i=0; i < N_param; ++i) { // Iterate rows of U
             // Compare elements of B_U_test_buf (Fortran LDB) with U_expected_cm (Fortran N_param LDB)
             EXPECT_NEAR(B_U_test_buf[i + j*LDB], U_expected_cm[i + j*N_param], check_tol_U) 
                 << "U_cm[" << i << "," << j << "]";
        }
    }
    EXPECT_NEAR(SCALE_out_val, SCALE_expected, 1e-4);
}


TEST_F(SG03BDTest, ParameterValidation) {
    setup_slicot_example_data(); // To get DICO_param etc.
    N_param = 1; M_param = 1; // Override N, M for this test
    
    // Re-initialize test buffers for N=1, M=1
    A_test_buf.assign(1,0.0); E_test_buf.assign(1,0.0); Q_test_buf.assign(1,0.0);
    Z_test_buf.assign(1,0.0); B_U_test_buf.assign(1,0.0);
    ALPHAR_out_val.assign(1,0.0); ALFAI_out_val.assign(1,0.0); BETA_out_val.assign(1,0.0);
    LDA=1; LDE=1; LDQ=1; LDZ=1; LDB=1;


    int info = slicot_sg03bd('X', FACT_param, TRANS_param, N_param, M_param, A_test_buf.data(),LDA,E_test_buf.data(),LDE,Q_test_buf.data(),LDQ,Z_test_buf.data(),LDZ,B_U_test_buf.data(),LDB, &SCALE_out_val,ALPHAR_out_val.data(),ALFAI_out_val.data(),BETA_out_val.data(),1);
    EXPECT_EQ(info, -1);

    info = slicot_sg03bd(DICO_param, FACT_param, TRANS_param, -1, M_param, nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1, &SCALE_out_val,nullptr,nullptr,nullptr,1);
    EXPECT_EQ(info, -4);
    
    if (N_param > 0) { 
        info = slicot_sg03bd(DICO_param, FACT_param, TRANS_param, N_param, M_param, A_test_buf.data(),0,E_test_buf.data(),LDE,Q_test_buf.data(),LDQ,Z_test_buf.data(),LDZ,B_U_test_buf.data(),LDB, &SCALE_out_val,ALPHAR_out_val.data(),ALFAI_out_val.data(),BETA_out_val.data(),1);
        EXPECT_EQ(info, -7);
    }
}

TEST_F(SG03BDTest, ZeroDimensionN) {
    setup_slicot_example_data(); // To get DICO etc.
    N_param = 0; M_param = 1; // Override N, M
    
    // Buffers for N=0
    A_test_buf.assign(1,0.0); E_test_buf.assign(1,0.0); Q_test_buf.assign(1,0.0); 
    Z_test_buf.assign(1,0.0); B_U_test_buf.assign(1,0.0); 
    ALPHAR_out_val.assign(1,0.0); ALFAI_out_val.assign(1,0.0); BETA_out_val.assign(1,0.0);
    LDA=1; LDE=1; LDQ=1; LDZ=1; LDB=1; 
    
    std::vector<double> b_input_for_n0(std::max(1, M_param * 1)); 
    if (M_param == 0 && N_param == 0) b_input_for_n0.assign(1,0.0); 

    int info = slicot_sg03bd(DICO_param, FACT_param, TRANS_param, N_param, M_param,
                             nullptr, LDA, nullptr, LDE,
                             nullptr, LDQ, nullptr, LDZ,
                             (M_param > 0 || N_param > 0 ? b_input_for_n0.data() : nullptr), LDB, 
                             &SCALE_out_val,
                             ALPHAR_out_val.data(), ALFAI_out_val.data(), BETA_out_val.data(),
                             1 /* row_major = true */);
    EXPECT_EQ(info, 0);
    EXPECT_EQ(SCALE_out_val, 1.0); 
}

