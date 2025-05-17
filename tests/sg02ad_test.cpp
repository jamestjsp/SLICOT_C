#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max, std::min
#include <iostream>  // For debugging
#include <iomanip>   // For std::setprecision

#include "sg02ad.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions

// --- Test Fixture Base ---
class SG02ADTest : public ::testing::Test {
protected:
    // Default parameters for Python case, can be overridden by specific tests
    char DICO_param_py = 'D'; 
    char JOBB_param_py = 'B'; 
    char FACT_param_py = 'N'; 
    char UPLO_param_py = 'U'; 
    char JOBL_param_py = 'Z'; 
    char SCAL_param_py = 'N'; 
    char SORT_param_py = 'S'; 
    char ACC_param_py  = 'R'; 

    int N_param_py = 3;
    int M_param_py = 1;
    int P_param_py = 1; 
    double TOL_param_py = 0.0;

    double check_tol_X_py = 1e-8; 

    std::vector<double> A_data_py, E_data_py, B_data_py, Q_data_py, R_data_py, L_data_py;
    double RCONDU_out;
    std::vector<double> X_out;
    std::vector<double> ALFAR_out, ALFAI_out, BETA_out;
    std::vector<double> S_out, T_out, U_out;
    int IWARN_out;

    std::vector<double> X_expected_py_rowmajor; 
    int expected_info_py = 0; // Default expected info for Python case
    int expected_iwarn_py = 0; // Default expected iwarn for Python case


    // Leading dimensions - these will be set in specific test fixture SetUp methods
    int LDA, LDE, LDB, LDQ, LDR, LDL, LDX, LDS, LDT, LDU;
    
    void SetUpPythonCaseData() { // Renamed to avoid conflict with GTest's SetUp
        A_data_py = {
            0.63399379,  0.54906824,  0.76253406,
            0.5404729 ,  0.53745766,  0.08731853,
            0.27524045,  0.84922129,  0.46816220
        };
        E_data_py.assign(N_param_py * N_param_py, 0.0); 
        for(int i=0; i<N_param_py; ++i) E_data_py[i*N_param_py+i] = 1.0;
        B_data_py = {0.96861695, 0.05532739, 0.78934047}; 
        Q_data_py.assign(N_param_py * N_param_py, 0.0); 
        for(int i=0; i<N_param_py; ++i) Q_data_py[i*N_param_py+i] = 1.0;
        R_data_py = {1.0}; 
        L_data_py.assign(N_param_py * M_param_py, 0.0); 

        X_expected_py_rowmajor = {
            1.9127157658780856e+00, 8.9471099288711020e-01, 4.9741332606172323e-01, 
            8.9471099288711020e-01, 2.1942687126443481e+00, 5.4632712840131359e-01, 
            4.9741332606172323e-01, 5.4632712840131359e-01, 1.4427247130079717e+00
        };

        // Resize output buffers based on N_param_py
        X_out.resize(std::max(1,N_param_py * N_param_py)); // Ensure min size 1 for N=0
        if (N_param_py > 0) {
            ALFAR_out.resize(2 * N_param_py);
            ALFAI_out.resize(2 * N_param_py);
            BETA_out.resize(2 * N_param_py);
        
            int s_fortran_cols = (JOBB_param_py == 'B') ? (2 * N_param_py + M_param_py) : (2 * N_param_py);
            int t_fortran_cols = 2 * N_param_py;
            int u_fortran_cols = 2 * N_param_py;
            
            int lds_alloc_rows = (JOBB_param_py == 'B') ? std::max(1, 2*N_param_py + M_param_py) : std::max(1, 2*N_param_py);
            int ldt_alloc_rows = (JOBB_param_py == 'B') ? std::max(1, 2*N_param_py + M_param_py) : std::max(1, 2*N_param_py);
            int ldu_alloc_rows = std::max(1, 2*N_param_py);

            S_out.resize(std::max(1, lds_alloc_rows * s_fortran_cols)); 
            T_out.resize(std::max(1, ldt_alloc_rows * t_fortran_cols));
            U_out.resize(std::max(1, ldu_alloc_rows * u_fortran_cols));
        } else { // N_param_py == 0
            ALFAR_out.assign(1,0.0); ALFAI_out.assign(1,0.0); BETA_out.assign(1,0.0); // Minimal for non-null ptr
            S_out.assign(1,0.0); T_out.assign(1,0.0); U_out.assign(1,0.0);
        }
    }
};

class SG02ADTestColMajorPython : public SG02ADTest {
protected:
    std::vector<double> A_cm, E_cm, B_cm, Q_cm, R_cm, L_cm;
    std::vector<double> X_cm_out, S_cm_out, T_cm_out, U_cm_out;

    void SetUp() override {
        SetUpPythonCaseData();
        LDA = std::max(1, N_param_py); LDE = std::max(1, N_param_py); LDB = std::max(1, N_param_py);
        LDQ = (FACT_param_py == 'N' || FACT_param_py == 'D') ? std::max(1,N_param_py) : std::max(1,P_param_py);
        LDR = (JOBB_param_py == 'B' && (FACT_param_py == 'N' || FACT_param_py == 'C')) ? std::max(1,M_param_py) : ((JOBB_param_py == 'B' && (FACT_param_py == 'D' || FACT_param_py == 'B')) ? std::max(1,P_param_py) : 1);
        LDL = (JOBB_param_py == 'B' && JOBL_param_py == 'N') ? std::max(1,N_param_py) : 1;
        LDX = std::max(1, N_param_py);
        LDS = (JOBB_param_py == 'B') ? std::max(1, 2*N_param_py + M_param_py) : std::max(1, 2*N_param_py);
        LDT = (JOBB_param_py == 'B') ? std::max(1, 2*N_param_py + M_param_py) : std::max(1, 2*N_param_py); 
        LDU = std::max(1, 2*N_param_py);

        if (!A_data_py.empty()) { A_cm.resize(A_data_py.size()); slicot_transpose_to_fortran_with_ld(A_data_py.data(), A_cm.data(), N_param_py, N_param_py, N_param_py, LDA, sizeof(double));}
        if (!E_data_py.empty()) { E_cm.resize(E_data_py.size()); slicot_transpose_to_fortran_with_ld(E_data_py.data(), E_cm.data(), N_param_py, N_param_py, N_param_py, LDE, sizeof(double));}
        if (!B_data_py.empty()) { B_cm.resize(B_data_py.size()); slicot_transpose_to_fortran_with_ld(B_data_py.data(), B_cm.data(), N_param_py, M_param_py, M_param_py, LDB, sizeof(double));}
        if (!Q_data_py.empty()) { Q_cm.resize(Q_data_py.size()); slicot_transpose_to_fortran_with_ld(Q_data_py.data(), Q_cm.data(), N_param_py, N_param_py, N_param_py, LDQ, sizeof(double));}
        if (!R_data_py.empty() && JOBB_param_py == 'B') { R_cm.resize(R_data_py.size()); slicot_transpose_to_fortran_with_ld(R_data_py.data(), R_cm.data(), M_param_py, M_param_py, M_param_py, LDR, sizeof(double));}
        if (!L_data_py.empty() && JOBB_param_py == 'B' && JOBL_param_py == 'N') { L_cm.resize(L_data_py.size()); slicot_transpose_to_fortran_with_ld(L_data_py.data(), L_cm.data(), N_param_py, M_param_py, M_param_py, LDL, sizeof(double));}
        
        X_cm_out.resize(X_out.size());
        if (N_param_py > 0) {
            S_cm_out.resize(S_out.size()); T_cm_out.resize(T_out.size()); U_cm_out.resize(U_out.size());
        } else {
            S_cm_out.assign(1,0.0); T_cm_out.assign(1,0.0); U_cm_out.assign(1,0.0);
        }
    }
};

class SG02ADTestRowMajorPython : public SG02ADTest {
protected:
    void SetUp() override {
        SetUpPythonCaseData();
        LDA = N_param_py; LDE = N_param_py; LDB = M_param_py; 
        LDQ = N_param_py; LDR = M_param_py; LDL = M_param_py;
        LDX = N_param_py; 
        LDS = (JOBB_param_py == 'B') ? (2 * N_param_py + M_param_py) : (2 * N_param_py); 
        LDT = 2 * N_param_py; 
        LDU = 2 * N_param_py; 
    }
};


TEST_F(SG02ADTestColMajorPython, PythonCase1ColMajor) {
    int info_result = slicot_sg02ad(DICO_param_py, JOBB_param_py, FACT_param_py, UPLO_param_py, JOBL_param_py, SCAL_param_py, SORT_param_py, ACC_param_py,
                                   N_param_py, M_param_py, P_param_py,
                                   A_cm.data(), LDA, E_cm.data(), LDE, B_cm.data(), LDB,
                                   Q_cm.data(), LDQ, R_cm.empty() ? nullptr : R_cm.data(), LDR, 
                                   L_cm.empty() ? nullptr : L_cm.data(), LDL,
                                   &RCONDU_out, X_cm_out.data(), LDX,
                                   ALFAR_out.data(), ALFAI_out.data(), BETA_out.data(),
                                   S_cm_out.data(), LDS, T_cm_out.data(), LDT, U_cm_out.data(), LDU,
                                   TOL_param_py, &IWARN_out,
                                   0 /* row_major = false */);

    ASSERT_EQ(info_result, expected_info_py);
    ASSERT_EQ(IWARN_out, expected_iwarn_py);

    std::vector<double> X_expected_cm(X_expected_py_rowmajor.size());
    if (N_param_py > 0) { // Only transpose if N > 0
        slicot_transpose_to_fortran_with_ld(X_expected_py_rowmajor.data(), X_expected_cm.data(), N_param_py, N_param_py, N_param_py, LDX, sizeof(double));
        for(size_t i=0; i < X_cm_out.size(); ++i) {
            EXPECT_NEAR(X_cm_out[i], X_expected_cm[i], check_tol_X_py) << "X_cm_out[" << i << "]";
        }
    } else { // If N=0, X_cm_out should be empty or its contents not compared like this
        ASSERT_TRUE(X_cm_out.empty() || X_cm_out.size() == 1); // Allow for minimal allocation
    }
}

TEST_F(SG02ADTestRowMajorPython, PythonCase1RowMajor) {
    std::vector<double> a_test = A_data_py; std::vector<double> e_test = E_data_py;
    std::vector<double> b_test = B_data_py; std::vector<double> q_test = Q_data_py;
    std::vector<double> r_test = R_data_py; std::vector<double> l_test = L_data_py;

    int info_result = slicot_sg02ad(DICO_param_py, JOBB_param_py, FACT_param_py, UPLO_param_py, JOBL_param_py, SCAL_param_py, SORT_param_py, ACC_param_py,
                                   N_param_py, M_param_py, P_param_py,
                                   a_test.data(), LDA, e_test.data(), LDE, b_test.data(), LDB,
                                   q_test.data(), LDQ, r_test.empty() ? nullptr : r_test.data(), LDR, 
                                   l_test.empty() ? nullptr : l_test.data(), LDL,
                                   &RCONDU_out, X_out.data(), LDX,
                                   ALFAR_out.data(), ALFAI_out.data(), BETA_out.data(),
                                   S_out.data(), LDS, T_out.data(), LDT, U_out.data(), LDU,
                                   TOL_param_py, &IWARN_out,
                                   1 /* row_major = true */);
    ASSERT_EQ(info_result, expected_info_py);
    ASSERT_EQ(IWARN_out, expected_iwarn_py);

    if (N_param_py > 0) {
        for(size_t i=0; i < X_out.size(); ++i) {
            EXPECT_NEAR(X_out[i], X_expected_py_rowmajor[i], check_tol_X_py) << "X_out[" << i << "]";
        }
    } else {
         ASSERT_TRUE(X_out.empty() || X_out.size() == 1);
    }
}

TEST_F(SG02ADTest, SlicotExampleColMajor) {
    char DICO_ex = 'C'; char JOBB_ex = 'B'; char FACT_ex = 'B'; char UPLO_ex = 'U';
    char JOBL_ex = 'Z'; char SCAL_ex = 'N'; char SORT_ex = 'S'; char ACC_ex = 'N';
    int N_ex = 2; int M_ex = 1; int P_ex = 3; double TOL_ex = 0.0;

    std::vector<double> A_ex_rm = {0.0, 1.0, 0.0, 0.0};
    std::vector<double> E_ex_rm = {1.0, 0.0, 0.0, 1.0};
    std::vector<double> B_ex_rm = {0.0, 1.0}; 
    std::vector<double> Q_Cfact_ex_rm = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0}; 
    std::vector<double> R_Dfact_ex_rm = {0.0, 0.0, 1.0}; 
    std::vector<double> L_ex_rm; 

    std::vector<double> X_ex_expected_rm = {1.7321, 1.0000, 1.0000, 1.7321}; 

    int LDA_ex_f = std::max(1,N_ex); int LDE_ex_f = std::max(1,N_ex); int LDB_ex_f = std::max(1,N_ex);
    int LDQ_ex_f = std::max(1,P_ex); int LDR_ex_f = std::max(1,P_ex); int LDL_ex_f = 1; 
    int LDX_ex_f = std::max(1,N_ex);
    int LDS_ex_f = std::max(1, 2*N_ex + M_ex);
    int LDT_ex_f = std::max(1, 2*N_ex + M_ex);
    int LDU_ex_f = std::max(1, 2*N_ex);

    std::vector<double> A_ex_cm(N_ex*N_ex), E_ex_cm(N_ex*N_ex), B_ex_cm(N_ex*M_ex);
    std::vector<double> Q_Cfact_ex_cm(P_ex*N_ex), R_Dfact_ex_cm(P_ex*M_ex);
    
    slicot_transpose_to_fortran_with_ld(A_ex_rm.data(), A_ex_cm.data(), N_ex, N_ex, N_ex, LDA_ex_f, sizeof(double));
    slicot_transpose_to_fortran_with_ld(E_ex_rm.data(), E_ex_cm.data(), N_ex, N_ex, N_ex, LDE_ex_f, sizeof(double));
    slicot_transpose_to_fortran_with_ld(B_ex_rm.data(), B_ex_cm.data(), N_ex, M_ex, M_ex, LDB_ex_f, sizeof(double));
    slicot_transpose_to_fortran_with_ld(Q_Cfact_ex_rm.data(), Q_Cfact_ex_cm.data(), P_ex, N_ex, N_ex, LDQ_ex_f, sizeof(double));
    slicot_transpose_to_fortran_with_ld(R_Dfact_ex_rm.data(), R_Dfact_ex_cm.data(), P_ex, M_ex, M_ex, LDR_ex_f, sizeof(double));

    double rcondu_ex_out;
    std::vector<double> X_ex_cm_out(N_ex*N_ex);
    std::vector<double> alfar_ex_out(2*N_ex), alfai_ex_out(2*N_ex), beta_ex_out(2*N_ex);
    std::vector<double> S_ex_cm_out(LDS_ex_f * (2*N_ex+M_ex)); 
    std::vector<double> T_ex_cm_out(LDT_ex_f * (2*N_ex));
    std::vector<double> U_ex_cm_out(LDU_ex_f * (2*N_ex));
    int iwarn_ex_out;
    
    int info_result = slicot_sg02ad(DICO_ex, JOBB_ex, FACT_ex, UPLO_ex, JOBL_ex, SCAL_ex, SORT_ex, ACC_ex,
                                   N_ex, M_ex, P_ex,
                                   A_ex_cm.data(), LDA_ex_f, E_ex_cm.data(), LDE_ex_f, B_ex_cm.data(), LDB_ex_f,
                                   Q_Cfact_ex_cm.data(), LDQ_ex_f, R_Dfact_ex_cm.data(), LDR_ex_f, 
                                   nullptr, LDL_ex_f, 
                                   &rcondu_ex_out, X_ex_cm_out.data(), LDX_ex_f,
                                   alfar_ex_out.data(), alfai_ex_out.data(), beta_ex_out.data(),
                                   S_ex_cm_out.data(), LDS_ex_f, T_ex_cm_out.data(), LDT_ex_f, U_ex_cm_out.data(), LDU_ex_f,
                                   TOL_ex, &iwarn_ex_out,
                                   0 /* row_major = false */);

    ASSERT_EQ(info_result, 0);
    ASSERT_EQ(iwarn_ex_out, 0);

    std::vector<double> X_ex_rm_computed(N_ex*N_ex);
    slicot_transpose_to_c_with_ld(X_ex_cm_out.data(), X_ex_rm_computed.data(), N_ex, N_ex, LDX_ex_f, N_ex, sizeof(double));
    
    double slicot_example_tol = 1e-4; 
    for(size_t i=0; i < X_ex_rm_computed.size(); ++i) {
        EXPECT_NEAR(X_ex_rm_computed[i], X_ex_expected_rm[i], slicot_example_tol) << "X_ex_rm_computed[" << i << "]";
    }
}

// Test for parameter validation
TEST_F(SG02ADTest, ParameterValidation) {
    int iwarn_val; double rcondu_val;
    // Use python case parameters as a base for valid dimensions
    int n_val = N_param_py, m_val = M_param_py, p_val = P_param_py;
    double tol_val = TOL_param_py;
    char dico_val = DICO_param_py, jobb_val = JOBB_param_py, fact_val = FACT_param_py, uplo_val = UPLO_param_py;
    char jobl_val = JOBL_param_py, scal_val = SCAL_param_py, sort_val = SORT_param_py, acc_val = ACC_param_py;

    // Dummy buffers for outputs - ensure they are non-null and minimally sized
    std::vector<double> x_dummy(std::max(1,n_val*n_val)), alfar_dummy(std::max(1,2*n_val)), alfai_dummy(std::max(1,2*n_val)), beta_dummy(std::max(1,2*n_val));
    int s_cols_val = (jobb_val == 'B') ? (2 * n_val + m_val) : (2 * n_val);
    int t_cols_val = 2 * n_val;
    int u_cols_val = 2 * n_val;
    int lds_val_min = (jobb_val == 'B') ? std::max(1, 2*n_val + m_val) : std::max(1, 2*n_val);
    int ldt_val_min = (jobb_val == 'B') ? std::max(1, 2*n_val + m_val) : std::max(1, 2*n_val);
    int ldu_val_min = std::max(1, 2*n_val);
    std::vector<double> s_dummy(std::max(1, lds_val_min * s_cols_val));
    std::vector<double> t_dummy(std::max(1, ldt_val_min * t_cols_val));
    std::vector<double> u_dummy(std::max(1, ldu_val_min * u_cols_val));
    
    // Dummy input buffers (can be null if N=0, but for general validation, assume N>0)
    std::vector<double> a_dummy(std::max(1,n_val*n_val)), e_dummy(std::max(1,n_val*n_val));
    std::vector<double> b_dummy(std::max(1,n_val*((jobb_val=='B')?m_val:n_val)));
    std::vector<double> q_dummy(std::max(1,((fact_val=='N'||fact_val=='D')?n_val:p_val)*n_val));
    std::vector<double> r_dummy(std::max(1,((jobb_val=='B')?((fact_val=='N'||fact_val=='C')?m_val:p_val):1)*((jobb_val=='B')?m_val:1)));
    std::vector<double> l_dummy(std::max(1,n_val*m_val));


    int info_result;

    // Test invalid DICO
    info_result = slicot_sg02ad('X', jobb_val, fact_val, uplo_val, jobl_val, scal_val, sort_val, acc_val,
                                   n_val, m_val, p_val, a_dummy.data(),std::max(1,n_val),e_dummy.data(),std::max(1,n_val),b_dummy.data(),std::max(1,n_val),
                                   q_dummy.data(),std::max(1,n_val),r_dummy.data(),std::max(1,m_val),l_dummy.data(),std::max(1,n_val),
                                   &rcondu_val, x_dummy.data(),std::max(1,n_val),alfar_dummy.data(),alfai_dummy.data(),beta_dummy.data(),
                                   s_dummy.data(),lds_val_min,t_dummy.data(),ldt_val_min,u_dummy.data(),ldu_val_min,
                                   tol_val, &iwarn_val, 0);
    EXPECT_EQ(info_result, -1);

    // Test invalid N
     info_result = slicot_sg02ad(dico_val, jobb_val, fact_val, uplo_val, jobl_val, scal_val, sort_val, acc_val,
                                   -1, m_val, p_val, nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1,
                                   &rcondu_val, x_dummy.data(),1,alfar_dummy.data(),alfai_dummy.data(),beta_dummy.data(),
                                   s_dummy.data(),1,t_dummy.data(),1,u_dummy.data(),1,
                                   tol_val, &iwarn_val, 0);
    EXPECT_EQ(info_result, -9);

    // Test invalid M (if JOBB='B')
    if (jobb_val == 'B') {
        info_result = slicot_sg02ad(dico_val, jobb_val, fact_val, uplo_val, jobl_val, scal_val, sort_val, acc_val,
                                   n_val, -1, p_val, a_dummy.data(),std::max(1,n_val),e_dummy.data(),std::max(1,n_val),b_dummy.data(),std::max(1,n_val),
                                   q_dummy.data(),std::max(1,n_val),r_dummy.data(),1,l_dummy.data(),std::max(1,n_val),
                                   &rcondu_val, x_dummy.data(),std::max(1,n_val),alfar_dummy.data(),alfai_dummy.data(),beta_dummy.data(),
                                   s_dummy.data(),lds_val_min,t_dummy.data(),ldt_val_min,u_dummy.data(),ldu_val_min,
                                   tol_val, &iwarn_val, 0);
        EXPECT_EQ(info_result, -10);
    }

    // Test invalid LDA (if N > 0)
    if (n_val > 0) {
        info_result = slicot_sg02ad(dico_val, jobb_val, fact_val, uplo_val, jobl_val, scal_val, sort_val, acc_val,
                                   n_val, m_val, p_val, a_dummy.data(),0,e_dummy.data(),std::max(1,n_val),b_dummy.data(),std::max(1,n_val),
                                   q_dummy.data(),std::max(1,n_val),r_dummy.data(),std::max(1,m_val),l_dummy.data(),std::max(1,n_val),
                                   &rcondu_val, x_dummy.data(),std::max(1,n_val),alfar_dummy.data(),alfai_dummy.data(),beta_dummy.data(),
                                   s_dummy.data(),lds_val_min,t_dummy.data(),ldt_val_min,u_dummy.data(),ldu_val_min,
                                   tol_val, &iwarn_val, 0);
        EXPECT_EQ(info_result, -13);
    }
}

// Test for N=0
TEST_F(SG02ADTest, ZeroDimensionN) {
    char DICO_z = 'D'; char JOBB_z = 'B'; char FACT_z = 'N'; char UPLO_z = 'U';
    char JOBL_z = 'Z'; char SCAL_z = 'N'; char SORT_z = 'S'; char ACC_z = 'R';
    int N_z = 0; int M_z = 1; int P_z = 1; double TOL_z = 0.0;

    double rcondu_z_out;
    std::vector<double> x_z_out(1); // Min size for ptr
    std::vector<double> alfar_z_out(1), alfai_z_out(1), beta_z_out(1); // Min size
    std::vector<double> s_z_out(1), t_z_out(1), u_z_out(1); // Min size
    int iwarn_z_out;

    // For N=0, A, E, B, Q, L, X, S, T, U are effectively empty or not used in the same way.
    // R can be non-empty if M > 0.
    std::vector<double> R_z_data = {1.0}; // MxM = 1x1

    int LDA_z=1, LDE_z=1, LDB_z=1, LDQ_z=1, LDR_z=std::max(1,M_z), LDL_z=1;
    int LDX_z=1, LDS_z=std::max(1,M_z), LDT_z=1, LDU_z=1; // Adjusted for N=0

    int info_result = slicot_sg02ad(DICO_z, JOBB_z, FACT_z, UPLO_z, JOBL_z, SCAL_z, SORT_z, ACC_z,
                                   N_z, M_z, P_z,
                                   nullptr, LDA_z, nullptr, LDE_z, nullptr, LDB_z,
                                   nullptr, LDQ_z, R_z_data.data(), LDR_z, nullptr, LDL_z,
                                   &rcondu_z_out, x_z_out.data(), LDX_z,
                                   alfar_z_out.data(), alfai_z_out.data(), beta_z_out.data(),
                                   s_z_out.data(), LDS_z, t_z_out.data(), LDT_z, u_z_out.data(), LDU_z,
                                   TOL_z, &iwarn_z_out,
                                   0 /* row_major = false */);
    
    EXPECT_EQ(info_result, 0); // Expect successful execution
    EXPECT_EQ(iwarn_z_out, 0); // Expect no warnings
    // For N=0, RCONDU might be 1.0 or not well-defined by the routine.
    // X should be 0x0 (or not modified if x_out was passed as non-null for N=0).
    // ALFAR, etc., are 2*N, so empty.
}

