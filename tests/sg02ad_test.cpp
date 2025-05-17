#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max
#include <iostream>  // For debugging
#include <iomanip>   // For std::setprecision

#include "sg02ad.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions

// --- Test Fixture Base ---
class SG02ADTest : public ::testing::Test {
protected:
    // Parameters from Python test_sg02ad_case1
    char DICO_param_py = 'D'; // Discrete-time
    char JOBB_param_py = 'B'; // B and R are given
    char FACT_param_py = 'N'; // Q and R are not factored
    char UPLO_param_py = 'U'; // Upper triangle of Q, R stored
    char JOBL_param_py = 'Z'; // L is zero
    char SCAL_param_py = 'N'; // No scaling
    char SORT_param_py = 'S'; // Stable eigenvalues first
    char ACC_param_py  = 'R'; // Use iterative refinement for X

    int N_param_py = 3;
    int M_param_py = 1;
    int P_param_py = 1; 
    double TOL_param_py = 0.0;

    // Tolerance for comparing X matrix elements from Python case
    double check_tol_X_py = 1e-8; // Adjusted to a tighter tolerance

    // Input data vectors for Python case
    std::vector<double> A_data_py, E_data_py, B_data_py, Q_data_py, R_data_py, L_data_py;
    
    // Output data vectors
    double RCONDU_out;
    std::vector<double> X_out;
    std::vector<double> ALFAR_out, ALFAI_out, BETA_out;
    std::vector<double> S_out, T_out, U_out;
    int IWARN_out;

    // Expected X for Python case (NOW WITH HIGH PRECISION VALUES)
    std::vector<double> X_expected_py_rowmajor; 
    int expected_info = 0;
    int expected_iwarn = 0;

    // Leading dimensions
    int LDA, LDE, LDB, LDQ, LDR, LDL, LDX, LDS, LDT, LDU;
    
    void print_matrix_rm(const std::string& name, const std::vector<double>& mat, int rows, int cols) {
        std::cout << name << " (" << rows << "x" << cols << ") row-major:\n";
        if (mat.empty() && rows > 0 && cols > 0) { std::cout << "  (empty)\n"; return; }
        if (rows == 0 || cols == 0) { std::cout << "  (zero dim)\n"; return; }
        for (int i = 0; i < rows; ++i) {
            std::cout << "  [";
            for (int j = 0; j < cols; ++j) {
                std::cout << std::fixed << std::setprecision(16) << mat[i * cols + j] << (j == cols - 1 ? "" : ", ");
            }
            std::cout << "]\n";
        }
    }
    void print_matrix_cm(const std::string& name, const std::vector<double>& mat_cm, int rows, int cols, int ld_cm) {
        std::cout << name << " (" << rows << "x" << cols << ") col-major (ld=" << ld_cm << "):\n";
        if (mat_cm.empty() && rows > 0 && cols > 0) { std::cout << "  (empty)\n"; return; }
        if (rows == 0 || cols == 0) { std::cout << "  (zero dim)\n"; return; }
        for (int i = 0; i < rows; ++i) {
            std::cout << "  [";
            for (int j = 0; j < cols; ++j) {
                 std::cout << std::fixed << std::setprecision(16) << mat_cm[i + j * ld_cm] << (j == cols - 1 ? "" : ", ");
            }
            std::cout << "]\n";
        }
    }


    void SetUpPythonCase() {
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

        // Updated with high-precision values from Python script output
        X_expected_py_rowmajor = {
            1.9127157658780856e+00, 8.9471099288711020e-01, 4.9741332606172323e-01, 
            8.9471099288711020e-01, 2.1942687126443481e+00, 5.4632712840131359e-01, 
            4.9741332606172323e-01, 5.4632712840131359e-01, 1.4427247130079717e+00
        };

        X_out.resize(N_param_py * N_param_py);
        ALFAR_out.resize(2 * N_param_py);
        ALFAI_out.resize(2 * N_param_py);
        BETA_out.resize(2 * N_param_py);
        
        int s_cols = (JOBB_param_py == 'B') ? (2 * N_param_py + M_param_py) : (2 * N_param_py);
        int t_cols = 2 * N_param_py;
        int u_cols = 2 * N_param_py;
        // Use leading dimensions for allocation size calculation for S, T, U
        int lds_alloc = (JOBB_param_py == 'B') ? std::max(1, 2*N_param_py + M_param_py) : std::max(1, 2*N_param_py);
        int ldt_alloc = (JOBB_param_py == 'B') ? std::max(1, 2*N_param_py + M_param_py) : std::max(1, 2*N_param_py);
        int ldu_alloc = std::max(1, 2*N_param_py);


        S_out.resize(std::max(1, lds_alloc * s_cols)); 
        T_out.resize(std::max(1, ldt_alloc * t_cols));
        U_out.resize(std::max(1, ldu_alloc * u_cols));
    }
};

// Test fixture for Column Major data layout using Python example data
class SG02ADTestColMajorPython : public SG02ADTest {
protected:
    std::vector<double> A_cm, E_cm, B_cm, Q_cm, R_cm, L_cm;
    std::vector<double> X_cm_out, S_cm_out, T_cm_out, U_cm_out;

    void SetUp() override {
        SetUpPythonCase();
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
        }
    }
};

// Test fixture for Row Major data layout using Python example data
class SG02ADTestRowMajorPython : public SG02ADTest {
protected:
    void SetUp() override {
        SetUpPythonCase();
        LDA = N_param_py; LDE = N_param_py; LDB = M_param_py; 
        LDQ = N_param_py; LDR = M_param_py; LDL = M_param_py;
        LDX = N_param_py; 
        LDS = (JOBB_param_py == 'B') ? std::max(1, 2 * N_param_py + M_param_py) : std::max(1, 2 * N_param_py);
        LDT = (JOBB_param_py == 'B') ? std::max(1, 2 * N_param_py + M_param_py) : std::max(1, 2 * N_param_py);
        LDU = std::max(1, 2 * N_param_py);
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
                                   S_cm_out.empty() ? nullptr : S_cm_out.data(), LDS, 
                                   T_cm_out.empty() ? nullptr : T_cm_out.data(), LDT, 
                                   U_cm_out.empty() ? nullptr : U_cm_out.data(), LDU,
                                   TOL_param_py, &IWARN_out,
                                   0 /* row_major = false */);

    ASSERT_EQ(info_result, expected_info);
    ASSERT_EQ(IWARN_out, expected_iwarn);

    std::vector<double> X_expected_cm(X_expected_py_rowmajor.size());
    slicot_transpose_to_fortran_with_ld(X_expected_py_rowmajor.data(), X_expected_cm.data(), N_param_py, N_param_py, N_param_py, LDX, sizeof(double));

    for(size_t i=0; i < X_cm_out.size(); ++i) {
        EXPECT_NEAR(X_cm_out[i], X_expected_cm[i], check_tol_X_py) << "X_cm_out[" << i << "]";
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
                                   S_out.empty() ? nullptr : S_out.data(), LDS, 
                                   T_out.empty() ? nullptr : T_out.data(), LDT, 
                                   U_out.empty() ? nullptr : U_out.data(), LDU,
                                   TOL_param_py, &IWARN_out,
                                   1 /* row_major = true */);
    ASSERT_EQ(info_result, expected_info);
    ASSERT_EQ(IWARN_out, expected_iwarn);

    for(size_t i=0; i < X_out.size(); ++i) {
        EXPECT_NEAR(X_out[i], X_expected_py_rowmajor[i], check_tol_X_py) << "X_out[" << i << "]";
    }
}

// Test using data from SG02AD.html example
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
    int LDQ_ex_f = std::max(1,P_ex); 
    int LDR_ex_f = std::max(1,P_ex); 
    int LDL_ex_f = 1; 
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


TEST_F(SG02ADTestColMajorPython, ParameterValidation) {
    int iwarn_val; double rcondu_val;
    int n_val = N_param_py, m_val = M_param_py, p_val = P_param_py; // Use valid N, M, P for these
    double tol_val = TOL_param_py;
    // Buffers for outputs that are not checked but needed for the call
    std::vector<double> x_dummy(n_val*n_val), alfar_dummy(2*n_val), alfai_dummy(2*n_val), beta_dummy(2*n_val);
    std::vector<double> s_dummy(LDS*LDS), t_dummy(LDT*LDT), u_dummy(LDU*LDU); // Approx sizes
    if (n_val == 0) { // Handle N=0 for dummy outputs
        x_dummy.assign(1,0); alfar_dummy.assign(1,0); alfai_dummy.assign(1,0); beta_dummy.assign(1,0);
        s_dummy.assign(1,0); t_dummy.assign(1,0); u_dummy.assign(1,0);
    }


    // Test invalid DICO
    int info_result = slicot_sg02ad('X', JOBB_param_py, FACT_param_py, UPLO_param_py, JOBL_param_py, SCAL_param_py, SORT_param_py, ACC_param_py,
                                   n_val, m_val, p_val, nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1,
                                   &rcondu_val, x_dummy.data(),std::max(1,n_val),alfar_dummy.data(),alfai_dummy.data(),beta_dummy.data(),
                                   s_dummy.data(),LDS,t_dummy.data(),LDT,u_dummy.data(),LDU,
                                   tol_val, &iwarn_val, 0);
    EXPECT_EQ(info_result, -1);

    // Test invalid N
     info_result = slicot_sg02ad(DICO_param_py, JOBB_param_py, FACT_param_py, UPLO_param_py, JOBL_param_py, SCAL_param_py, SORT_param_py, ACC_param_py,
                                   -1, m_val, p_val, nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1,nullptr,1,
                                   &rcondu_val, x_dummy.data(),1,alfar_dummy.data(),alfai_dummy.data(),beta_dummy.data(),
                                   s_dummy.data(),1,t_dummy.data(),1,u_dummy.data(),1,
                                   tol_val, &iwarn_val, 0);
    EXPECT_EQ(info_result, -9);
}

