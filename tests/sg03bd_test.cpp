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

    // Input matrices (row-major)
    std::vector<double> A_data_rm, E_data_rm, Q_data_in_rm, Z_data_in_rm, B_in_data_rm;
    // Output matrices/scalars
    std::vector<double> A_out_rm, E_out_rm, Q_out_rm, Z_out_rm, U_out_rm; // U_out_rm is from b_in_u_out
    double SCALE_out_val;
    std::vector<double> ALPHAR_out_val, ALFAI_out_val, BETA_out_val;

    // Expected results
    std::vector<double> U_expected_rm;
    double SCALE_expected;
    int expected_info = 0;

    // Leading dimensions 
    int LDA, LDE, LDQ, LDZ, LDB;
    
    double check_tol_U = 1e-4; // Default tolerance for U comparison

    void setup_common_buffers(int n_val, int m_val) {
        N_param = n_val;
        M_param = m_val;
        if (n_val > 0) {
            A_out_rm.resize((size_t)n_val*n_val); E_out_rm.resize((size_t)n_val*n_val);
            Q_out_rm.resize((size_t)n_val*n_val); Z_out_rm.resize((size_t)n_val*n_val);
            U_out_rm.resize((size_t)n_val*n_val); // For output U
            ALPHAR_out_val.resize(n_val); ALFAI_out_val.resize(n_val); BETA_out_val.resize(n_val);
        } else { // N=0
            A_out_rm.assign(1,0.0); E_out_rm.assign(1,0.0); Q_out_rm.assign(1,0.0); 
            Z_out_rm.assign(1,0.0); U_out_rm.assign(1,0.0);
            ALPHAR_out_val.assign(1,0.0); ALFAI_out_val.assign(1,0.0); BETA_out_val.assign(1,0.0);
        }
    }
};

// Test based on SLICOT Documentation Example for SG03BD
TEST_F(SG03BDTest, SlicotDocExample) {
    DICO_param = 'C'; FACT_param = 'N'; TRANS_param = 'N';
    N_param = 3; M_param = 1; // M is num rows of op(B). B is MxN if TRANS='N'. Here B is 1x3.
    setup_common_buffers(N_param, M_param);

    A_data_rm = {-1.0, 3.0, -4.0, 0.0, 5.0, -2.0, -4.0, 4.0, 1.0};
    E_data_rm = {2.0, 1.0, 3.0, 2.0, 0.0, 1.0, 4.0, 5.0, 1.0};
    // B_in_data_rm is M x N = 1 x 3 for TRANS='N'
    B_in_data_rm = {2.0, -1.0, 7.0}; 
    
    // Q_data_in_rm and Z_data_in_rm are not strictly input if FACT='N', but buffers must be passed
    Q_data_in_rm.assign((size_t)N_param*N_param, 0.0); 
    Z_data_in_rm.assign((size_t)N_param*N_param, 0.0);

    U_expected_rm = {1.6003, -0.4418, -0.1523,
                     0.0000,  0.6795, -0.2499,
                     0.0000,  0.0000,  0.2041};
    SCALE_expected = 1.0000;
    check_tol_U = 1e-4; 

    // For row_major = 1 (C style)
    LDA = N_param; LDE = N_param; LDQ = N_param; LDZ = N_param; 
    // LDB for input B (MxN for TRANS='N') and output U (NxN)
    // C LDB (cols of B_in_data_rm) must be N_param.
    // C LDB for output U_out_rm (cols of U_out_rm) must be N_param.
    LDB = N_param; 
    
    A_out_rm = A_data_rm; E_out_rm = E_data_rm; 
    Q_out_rm = Q_data_in_rm; Z_out_rm = Z_data_in_rm; 
    // U_out_rm will receive the Cholesky factor U. It's an N_param x N_param matrix.
    // The input B_in_data_rm is M_param x N_param.
    // The Fortran routine expects B to be dimensioned LDB x N1, where N1 depends on TRANS.
    // On output, the first N x N part of B contains U.
    // We pass U_out_rm as the buffer for b_in_u_out, after copying B_in_data_rm to its CM equivalent.
    // For row_major=true, the wrapper will handle preparing b_cm for Fortran.
    // The b_in_u_out argument to the C wrapper should be dimensioned for the larger of input B or output U.
    // Here, U_out_rm (NxN) is used, and input B (MxN) will be copied to its CM version.
    // The C wrapper expects b_in_u_out to be a buffer that can hold the input B and receive output U.
    // If M < N, b_in_u_out needs to be at least N*N for output U.
    // If M >= N, b_in_u_out needs to be M*N for input B.
    // Let's use U_out_rm as the primary buffer and copy B_in_data_rm into it if needed for the call.
    // However, the wrapper logic for row_major will create b_cm based on b_in_nelems (input B)
    // and then copy U from b_cm back to b_in_u_out (which is U_out_rm).
    // So, we need to ensure b_in_u_out (U_out_rm) is correctly sized for output U.
    // And pass B_in_data_rm for the input part.
    // The C wrapper's b_in_u_out is used for both. So it should be sized for the larger.
    // Output U is NxN. Input B is MxN (TRANS='N').
    // So, b_in_u_out should be at least max(M*N, N*N).
    // For this test, U_out_rm is already N*N. B_in_data_rm is M*N.
    // We will pass B_in_data_rm to the wrapper, and the wrapper will handle it.
    // The output U will be written into B_in_data_rm if its LDB is sufficient.
    // This is tricky. Let's use U_out_rm for output and provide B_in_data_rm for input.
    // The wrapper needs to be smart about this.
    // The current wrapper uses b_in_u_out for both input B and output U.
    // If row_major, it creates b_cm, copies input B to b_cm, calls Fortran (which writes U to b_cm),
    // then copies U from b_cm to b_in_u_out.
    // So, b_in_u_out must be sized for U (NxN). LDB for C must be N.
    // Input B_data_rm is MxN.
    
    // For the C call, pass B_in_data_rm as the input part.
    // The output U will be written into this buffer by the wrapper if row_major=true,
    // provided LDB is appropriate for NxN output.
    // Let's make a copy for b_in_u_out to pass to the function.
    std::vector<double> b_call_buffer = B_in_data_rm; 
    // If M_param < N_param, b_call_buffer (sized M_param * N_param) might be too small for NxN output U.
    // The Fortran routine writes U into the B array.
    // LDB for Fortran must be >= N.
    // If row_major, C LDB for B_in is N. C LDB for U_out is N.
    // So, if b_in_u_out is passed as M_param x N_param (row_major),
    // and M_param < N_param, then the wrapper's copy-back of U (NxN) will be problematic.
    // The b_in_u_out for C wrapper should be dimensioned to hold the NxN output U.
    // So, we should pass U_out_rm (which is NxN) and the wrapper should handle copying B_in_data_rm into its b_cm.

    // We pass U_out_rm, which is NxN. The wrapper will create b_cm.
    // It will copy B_in_data_rm (MxN) to the relevant part of b_cm.
    // Fortran will operate on b_cm, writing U (NxN) into it.
    // Wrapper will copy U from b_cm back to U_out_rm.
    // This requires the wrapper to handle the input B dimensions correctly when populating b_cm.
    // The current wrapper copies the entire b_in_u_out to b_cm if row_major. This is wrong if b_in_u_out is input B.
    // The wrapper's b_ptr should point to a CM version of input B for the call,
    // and then the output U (which overwrites B in Fortran) should be copied back.

    // Let's assume the C wrapper correctly handles:
    // 1. Taking b_in_u_out (as a buffer for input B, M_param x N_param for TRANS='N').
    // 2. Transposing this to b_cm.
    // 3. Fortran overwrites b_cm with U (N_param x N_param).
    // 4. Transposing U from b_cm back to b_in_u_out (which must be large enough for N_param x N_param).
    // So, for the test, we pass U_out_rm (sized N_param x N_param) as b_in_u_out.
    // And we need to provide the actual B input data separately for the wrapper to use if row_major.
    // This is getting complicated. The Fortran B is truly I/O.
    // For simplicity, let's assume the C wrapper expects b_in_u_out to be the input B,
    // and that this same buffer will be overwritten by U, and LDB must be sufficient for U.
    U_out_rm = B_in_data_rm; // Start with B data
    if (TRANS_param == 'N' && M_param < N_param) { // If B_in is MxN and U_out is NxN, ensure buffer is NxN
        U_out_rm.resize((size_t)N_param * N_param); 
        // Copy B_in_data_rm (MxN) into the top MxN part of U_out_rm
        for(int i=0; i < M_param; ++i) {
            for(int j=0; j < N_param; ++j) {
                if ( (size_t)(i*N_param + j) < B_in_data_rm.size())
                     U_out_rm[i*N_param + j] = B_in_data_rm[i*N_param + j];
            }
        }
    } else if (TRANS_param == 'T' && M_param > N_param) { // If B_in is NxM and U_out is NxN
         // This case is fine as U_out_rm (NxN) will be smaller than input B (NxM)
         // The Fortran routine will write U into the top-left NxN part.
    }


    int info = slicot_sg03bd(DICO_param, FACT_param, TRANS_param, N_param, M_param,
                             A_data_rm.data(), LDA, E_data_rm.data(), LDE,
                             Q_data_in_rm.data(), LDQ, Z_data_in_rm.data(), LDZ,
                             U_out_rm.data(), LDB, // U_out_rm acts as B_in then U_out
                             &SCALE_out_val,
                             ALPHAR_out_val.data(), ALFAI_out_val.data(), BETA_out_val.data(),
                             1 /* row_major = true */);
    ASSERT_EQ(info, 0);
    for(size_t i=0; i < U_expected_rm.size(); ++i) { // Compare only the NxN part for U
        EXPECT_NEAR(U_out_rm[i], U_expected_rm[i], check_tol_U) << "U_out_rm[" << i << "]";
    }
    EXPECT_NEAR(SCALE_out_val, SCALE_expected, 1e-4);
}


TEST_F(SG03BDTest, ParameterValidation) {
    N_param = 1; M_param = 1;
    setup_common_buffers(N_param, M_param);
    LDA=1; LDE=1; LDQ=1; LDZ=1; LDB=1;

    int info = slicot_sg03bd('X', FACT_param, TRANS_param, N_param, M_param, A_out_rm.data(),LDA,E_out_rm.data(),LDE,Q_out_rm.data(),LDQ,Z_out_rm.data(),LDZ,U_out_rm.data(),LDB, &SCALE_out_val,ALPHAR_out_val.data(),ALFAI_out_val.data(),BETA_out_val.data(),1);
    EXPECT_EQ(info, -1); // Invalid DICO

    info = slicot_sg03bd(DICO_param, FACT_param, TRANS_param, -1, M_param, nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1, &SCALE_out_val,nullptr,nullptr,nullptr,1);
    EXPECT_EQ(info, -4); // Invalid N
    
    if (N_param > 0) { // LDA check only if N > 0
        info = slicot_sg03bd(DICO_param, FACT_param, TRANS_param, N_param, M_param, A_out_rm.data(),0,E_out_rm.data(),LDE,Q_out_rm.data(),LDQ,Z_out_rm.data(),LDZ,U_out_rm.data(),LDB, &SCALE_out_val,ALPHAR_out_val.data(),ALFAI_out_val.data(),BETA_out_val.data(),1);
        EXPECT_EQ(info, -7); // Invalid LDA
    }
}

TEST_F(SG03BDTest, ZeroDimensionN) {
    N_param = 0; M_param = 1; // M can be > 0 even if N=0
    setup_common_buffers(N_param, M_param);
    LDA=1; LDE=1; LDQ=1; LDZ=1; LDB=1; 

    // For N=0, A,E,Q,Z,U are 0x0. B_in is Mx0 (TRANS='N') or 0xM (TRANS='T').
    // ALPHAR etc are 0-dim.
    // Pointers for 0-dim matrices can be NULL or point to minimal buffer.
    // Wrapper should handle passing NULL to Fortran.
    
    int info = slicot_sg03bd('C', 'N', 'N', N_param, M_param,
                             nullptr, LDA, nullptr, LDE,
                             nullptr, LDQ, nullptr, LDZ,
                             nullptr, LDB, // B_in is Mx0 or 0xM (empty), U_out is 0x0
                             &SCALE_out_val,
                             ALPHAR_out_val.data(), ALFAI_out_val.data(), BETA_out_val.data(),
                             1 /* row_major = true */);
    EXPECT_EQ(info, 0);
    EXPECT_EQ(SCALE_out_val, 1.0); // Typically SCALE is 1.0 if no scaling needed or N=0
}

