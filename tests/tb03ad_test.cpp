#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max, std::min
#include <iostream>
#include <iomanip>

#include "tb03ad.h"
#include "slicot_utils.h" // For transpose functions if needed for manual verification

// Helper to print matrices for debugging
void print_matrix_tb03ad(const std::string& name, const double* mat, int r, int c, int ld_c, bool is_row_major) {
    if (!mat) { std::cout << name << " is NULL\n"; return; }
    std::cout << name << " (" << r << "x" << c << ") " << (is_row_major ? "RowMajor" : "ColMajor") << ", LD=" << ld_c << ":\n";
    if (r == 0 || c == 0) { std::cout << "  (empty or zero dim)\n"; return; }
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < r; ++i) {
        std::cout << "  [";
        for (int j = 0; j < c; ++j) {
            double val = is_row_major ? mat[i * ld_c + j] : mat[i + j * ld_c];
            std::cout << val << (j == c - 1 ? "" : ", ");
        }
        std::cout << "]\n";
    }
     std::cout << std::resetiosflags(std::ios_base::fixed) << std::setprecision(6);
}

void print_3d_matrix_tb03ad(const std::string& name, const double* mat, 
                            int d1, int d2, int d3, // logical dimensions of data
                            int ld1_layout, int ld2_layout, // leading dimensions used for memory layout
                            bool is_row_major_layout) {
    if (!mat) { std::cout << name << " is NULL\n"; return; }
    std::cout << name << " (" << d1 << "x" << d2 << "x" << d3 << ") "
              << (is_row_major_layout ? "RowMajorLayout" : "ColMajorLayout")
              << ", LD1_Layout=" << ld1_layout << ", LD2_Layout=" << ld2_layout << ":\n";
    if (d1 == 0 || d2 == 0 || d3 == 0) { std::cout << "  (empty or zero dim)\n"; return; }
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < d1; ++i) { // Iterate logical dim1
        for (int j = 0; j < d2; ++j) { // Iterate logical dim2
            std::cout << "  " << name << "[" << i << "][" << j << "][:] = [";
            for (int k = 0; k < d3; ++k) { // Iterate logical dim3 (coeffs)
                double val;
                if (is_row_major_layout) { // C-style flat row-major: mat[i*ld2_layout*d3 + j*d3 + k]
                    val = mat[i * (size_t)ld2_layout * d3 + j * d3 + k];
                } else { // Fortran-style flat col-major: mat[i + j*ld1_layout + k*ld1_layout*ld2_layout]
                    val = mat[i + j * (size_t)ld1_layout + k * (size_t)ld1_layout * ld2_layout];
                }
                std::cout << val << (k == d3 - 1 ? "" : ", ");
            }
            std::cout << "]\n";
        }
    }
    std::cout << std::resetiosflags(std::ios_base::fixed) << std::setprecision(6);
}


class TB03ADTest : public ::testing::Test {
protected:
    char LERI_param;
    char EQUIL_param;
    int N_param, M_param, P_param;
    double TOL_param;

    std::vector<double> A_data_rm, B_data_rm, C_data_rm, D_data_rm;
    std::vector<double> A_test_buf, B_test_buf, C_test_buf, D_test_buf_const_copy; 
    
    int NR_out_val;
    std::vector<int> INDEX_out_val;
    std::vector<double> PCOEFF_out_val, QCOEFF_out_val, VCOEFF_out_val;

    int NR_expected;
    std::vector<int> INDEX_expected;
    std::vector<double> A_min_expected_rm, B_min_expected_rm, C_min_expected_rm;
    std::vector<double> PCOEFF_expected_rm, QCOEFF_expected_rm; 
    int info_expected;

    int LDA_c_param, LDB_c_param, LDC_c_param, LDD_c_param;
    int PCOEFF_d1_logic, PCOEFF_d2_logic;
    int QCOEFF_d1_logic, QCOEFF_d2_logic;
    int VCOEFF_d1_logic, VCOEFF_d2_logic; 
    
    int KPCOEF_expected;
    double check_tol_matrix = 1e-4;

    void setup_slicot_example_data_right() {
        LERI_param = 'R'; EQUIL_param = 'N';
        N_param = 3; M_param = 1; P_param = 2; TOL_param = 0.0;

        A_data_rm = {1.0, 2.0, 0.0, 4.0, -1.0, 0.0, 0.0, 0.0, 1.0}; 
        B_data_rm = {1.0, 0.0, 1.0}; 
        C_data_rm = {0.0, 1.0, -1.0, 0.0, 0.0, 1.0}; 
        D_data_rm = {0.0, 1.0}; 

        NR_expected = 3;
        A_min_expected_rm = {1.0000, -1.4142, 0.0000,
                            -2.8284, -1.0000, 2.8284,
                             0.0000, 1.4142, 1.0000}; 
        B_min_expected_rm = {-1.4142, 0.0000, 0.0000}; 
        C_min_expected_rm = {0.7071, 1.0000, 0.7071,
                            -0.7071, 0.0000, -0.7071}; 
        
        INDEX_expected = {3}; 
        KPCOEF_expected = (INDEX_expected.empty() || INDEX_expected[0] < 0) ? 0 : INDEX_expected[0] + 1;


        PCOEFF_expected_rm = {0.1768, -0.1768, -1.5910, 1.5910}; 
        QCOEFF_expected_rm = { 
            0.0000, -0.1768, 0.7071, 0.8839, 
            0.1768, 0.0000, -1.5910, 0.0000  
        };
        info_expected = 0;
    }

    void prepare_buffers_and_params(bool is_row_major_test) {
        A_test_buf = A_data_rm; B_test_buf = B_data_rm;
        C_test_buf = C_data_rm; D_test_buf_const_copy = D_data_rm;

        int porm_logic = (LERI_param == 'L') ? P_param : M_param;
        int porp_logic = (LERI_param == 'L') ? M_param : P_param;
        porm_logic = std::max(0, porm_logic); 
        porp_logic = std::max(0, porp_logic);
        
        INDEX_out_val.assign(porm_logic > 0 ? porm_logic : 1, 0);

        PCOEFF_d1_logic = porm_logic; PCOEFF_d2_logic = porm_logic;
        QCOEFF_d1_logic = (LERI_param == 'L') ? porm_logic : porp_logic; 
        QCOEFF_d2_logic = (LERI_param == 'L') ? porp_logic : porm_logic; 
        VCOEFF_d1_logic = porm_logic; VCOEFF_d2_logic = N_param; 
                                           
        int f_ldpco1, f_ldpco2, f_ldqco1, f_ldqco2, f_ldvco1, f_ldvco2;
        int current_max_mp = std::max(1, std::max(M_param, P_param));

        if (LERI_param == 'L') {
            f_ldpco1 = std::max(1, P_param); f_ldpco2 = std::max(1, P_param);
            f_ldqco1 = std::max(1, P_param); f_ldqco2 = std::max(1, M_param);
            f_ldvco1 = std::max(1, P_param); f_ldvco2 = std::max(1, N_param);
        } else { 
            f_ldpco1 = std::max(1, M_param); f_ldpco2 = std::max(1, M_param);
            f_ldqco1 = std::max(1, current_max_mp); f_ldqco2 = std::max(1, current_max_mp);
            f_ldvco1 = std::max(1, M_param); f_ldvco2 = std::max(1, N_param);
        }
        
        size_t pcoeff_alloc_size, qcoeff_alloc_size, vcoeff_alloc_size;
        int third_dim_coeffs = std::max(1, N_param + 1); // N+1 for coefficient arrays

        if (is_row_major_test) {
            pcoeff_alloc_size = (size_t)PCOEFF_d1_logic * PCOEFF_d2_logic * third_dim_coeffs;
            qcoeff_alloc_size = (size_t)QCOEFF_d1_logic * QCOEFF_d2_logic * third_dim_coeffs;
            vcoeff_alloc_size = (size_t)VCOEFF_d1_logic * std::max(1,VCOEFF_d2_logic) * third_dim_coeffs;
        } else { 
            pcoeff_alloc_size = (size_t)f_ldpco1 * f_ldpco2 * third_dim_coeffs;
            qcoeff_alloc_size = (size_t)f_ldqco1 * f_ldqco2 * third_dim_coeffs;
            vcoeff_alloc_size = (size_t)f_ldvco1 * f_ldvco2 * third_dim_coeffs;
        }
        if (PCOEFF_d1_logic == 0 || PCOEFF_d2_logic == 0) pcoeff_alloc_size = 0;
        if (QCOEFF_d1_logic == 0 || QCOEFF_d2_logic == 0) qcoeff_alloc_size = 0;
        if (VCOEFF_d1_logic == 0 || VCOEFF_d2_logic == 0) vcoeff_alloc_size = 0;


        PCOEFF_out_val.assign(pcoeff_alloc_size > 0 ? pcoeff_alloc_size : 1, 0.0);
        QCOEFF_out_val.assign(qcoeff_alloc_size > 0 ? qcoeff_alloc_size : 1, 0.0);
        VCOEFF_out_val.assign(vcoeff_alloc_size > 0 ? vcoeff_alloc_size : 1, 0.0);

        if (is_row_major_test) {
            LDA_c_param = N_param; LDB_c_param = M_param;
            LDC_c_param = N_param; LDD_c_param = M_param;
        } else { 
            LDA_c_param = std::max(1,N_param); LDB_c_param = std::max(1,N_param); // Fortran B is N x max_mp
            LDC_c_param = std::max(1, current_max_mp); // Fortran C is max_mp x N
            LDD_c_param = std::max(1, current_max_mp); // Fortran D is max_mp x max_mp
        }
        LDA_c_param = std::max(1, LDA_c_param); LDB_c_param = std::max(1, LDB_c_param);
        LDC_c_param = std::max(1, LDC_c_param); LDD_c_param = std::max(1, LDD_c_param);
    }
};

TEST_F(TB03ADTest, SlicotDocExampleRightRowMajor) {
    setup_slicot_example_data_right();
    prepare_buffers_and_params(true);

    int info = slicot_tb03ad(LERI_param, EQUIL_param, N_param, M_param, P_param,
                             A_test_buf.data(), LDA_c_param, B_test_buf.data(), LDB_c_param,
                             C_test_buf.data(), LDC_c_param, D_test_buf_const_copy.data(), LDD_c_param,
                             &NR_out_val, INDEX_out_val.data(),
                             PCOEFF_out_val.data(), PCOEFF_d1_logic, PCOEFF_d2_logic,
                             QCOEFF_out_val.data(), QCOEFF_d1_logic, QCOEFF_d2_logic,
                             VCOEFF_out_val.data(), VCOEFF_d1_logic, VCOEFF_d2_logic,
                             TOL_param, 1);

    ASSERT_EQ(info, info_expected);
    ASSERT_EQ(NR_out_val, NR_expected);

    for(int i=0; i<NR_expected; ++i) for(int j=0; j<NR_expected; ++j) EXPECT_NEAR(A_test_buf[i*LDA_c_param+j], A_min_expected_rm[i*NR_expected+j], check_tol_matrix);
    if (M_param > 0) for(int i=0; i<NR_expected; ++i) for(int j=0; j<M_param; ++j) EXPECT_NEAR(B_test_buf[i*LDB_c_param+j], B_min_expected_rm[i*M_param+j], check_tol_matrix);
    if (P_param > 0) for(int i=0; i<P_param; ++i) for(int j=0; j<NR_expected; ++j) EXPECT_NEAR(C_test_buf[i*LDC_c_param+j], C_min_expected_rm[i*NR_expected+j], check_tol_matrix);
    if (!INDEX_expected.empty()) for(size_t i=0; i<INDEX_expected.size(); ++i) EXPECT_EQ(INDEX_out_val[i], INDEX_expected[i]);

    if (PCOEFF_d1_logic > 0 && PCOEFF_d2_logic > 0 && KPCOEF_expected > 0)
    for(int i=0; i<PCOEFF_d1_logic; ++i) { 
        for(int j=0; j<PCOEFF_d2_logic; ++j) { 
            for(int k=0; k<KPCOEF_expected; ++k) { 
                double actual = PCOEFF_out_val[i*(size_t)PCOEFF_d2_logic*KPCOEF_expected + j*KPCOEF_expected + k];
                double expected = PCOEFF_expected_rm[i*(size_t)PCOEFF_d2_logic*KPCOEF_expected + j*KPCOEF_expected + k];
                EXPECT_NEAR(actual, expected, check_tol_matrix) << "PCOEFF_RM(" << i << "," << j << "," << k << ")";
            }
        }
    }
    if (QCOEFF_d1_logic > 0 && QCOEFF_d2_logic > 0 && KPCOEF_expected > 0)
     for(int i=0; i<QCOEFF_d1_logic; ++i) { 
        for(int j=0; j<QCOEFF_d2_logic; ++j) { 
            for(int k=0; k<KPCOEF_expected; ++k) { 
                double actual = QCOEFF_out_val[i*(size_t)QCOEFF_d2_logic*KPCOEF_expected + j*KPCOEF_expected + k];
                double expected = QCOEFF_expected_rm[i*(size_t)QCOEFF_d2_logic*KPCOEF_expected + j*KPCOEF_expected + k];
                EXPECT_NEAR(actual, expected, check_tol_matrix) << "QCOEFF_RM(" << i << "," << j << "," << k << ")";
            }
        }
    }
}

TEST_F(TB03ADTest, SlicotDocExampleRightColMajor) {
    setup_slicot_example_data_right();
    prepare_buffers_and_params(false);
    
    int current_max_mp = std::max(1, std::max(M_param, P_param));

    std::vector<double> a_cm_in((size_t)LDA_c_param * N_param);
    std::vector<double> b_cm_in((size_t)LDB_c_param * current_max_mp); // B is N x max_mp for Fortran
    std::vector<double> c_cm_in((size_t)LDC_c_param * N_param);      // C is max_mp x N for Fortran
    std::vector<double> d_cm_in((size_t)LDD_c_param * current_max_mp); // D is max_mp x max_mp for Fortran
    
    std::fill(b_cm_in.begin(), b_cm_in.end(), 0.0); // Initialize workspace parts
    std::fill(d_cm_in.begin(), d_cm_in.end(), 0.0);

    if (N_param > 0) slicot_transpose_to_fortran_with_ld(A_data_rm.data(), a_cm_in.data(), N_param, N_param, N_param, LDA_c_param, sizeof(double));
    if (N_param > 0 && M_param > 0) slicot_transpose_to_fortran_with_ld(B_data_rm.data(), b_cm_in.data(), N_param, M_param, M_param, LDB_c_param, sizeof(double));
    if (P_param > 0 && N_param > 0) slicot_transpose_to_fortran_with_ld(C_data_rm.data(), c_cm_in.data(), P_param, N_param, N_param, LDC_c_param, sizeof(double));
    if (P_param > 0 && M_param > 0) slicot_transpose_to_fortran_with_ld(D_data_rm.data(), d_cm_in.data(), P_param, M_param, M_param, LDD_c_param, sizeof(double));

    int info = slicot_tb03ad(LERI_param, EQUIL_param, N_param, M_param, P_param,
                             a_cm_in.data(), LDA_c_param, b_cm_in.data(), LDB_c_param,
                             c_cm_in.data(), LDC_c_param, d_cm_in.data(), LDD_c_param,
                             &NR_out_val, INDEX_out_val.data(),
                             PCOEFF_out_val.data(), PCOEFF_d1_logic, PCOEFF_d2_logic,
                             QCOEFF_out_val.data(), QCOEFF_d1_logic, QCOEFF_d2_logic,
                             VCOEFF_out_val.data(), VCOEFF_d1_logic, VCOEFF_d2_logic,
                             TOL_param, 0);

    ASSERT_EQ(info, info_expected);
    ASSERT_EQ(NR_out_val, NR_expected);

    std::vector<double> a_exp_cm((size_t)NR_expected*NR_expected), b_exp_cm((size_t)NR_expected*M_param), c_exp_cm((size_t)P_param*NR_expected);
    if (NR_expected>0) slicot_transpose_to_fortran_with_ld(A_min_expected_rm.data(), a_exp_cm.data(), NR_expected, NR_expected, NR_expected, NR_expected, sizeof(double));
    if (NR_expected>0 && M_param >0) slicot_transpose_to_fortran_with_ld(B_min_expected_rm.data(), b_exp_cm.data(), NR_expected, M_param, M_param, NR_expected, sizeof(double));
    if (P_param > 0 && NR_expected >0) slicot_transpose_to_fortran_with_ld(C_min_expected_rm.data(), c_exp_cm.data(), P_param, NR_expected, NR_expected, P_param, sizeof(double));

    if (NR_expected > 0) for(int j=0; j<NR_expected; ++j) for(int i=0; i<NR_expected; ++i) EXPECT_NEAR(a_cm_in[i+j*LDA_c_param], a_exp_cm[i+j*NR_expected], check_tol_matrix);
    if (M_param > 0 && NR_expected > 0) for(int j=0; j<M_param; ++j) for(int i=0; i<NR_expected; ++i) EXPECT_NEAR(b_cm_in[i+j*LDB_c_param], b_exp_cm[i+j*NR_expected], check_tol_matrix);
    if (P_param > 0 && NR_expected > 0) for(int j=0; j<NR_expected; ++j) for(int i=0; i<P_param; ++i) EXPECT_NEAR(c_cm_in[i+j*LDC_c_param], c_exp_cm[i+j*P_param], check_tol_matrix);
    if(!INDEX_expected.empty()) for(size_t i=0; i<INDEX_expected.size(); ++i) EXPECT_EQ(INDEX_out_val[i], INDEX_expected[i]);

    int f_ldpco1 = std::max(1, M_param); int f_ldpco2 = std::max(1, M_param);
    if (PCOEFF_d1_logic > 0 && PCOEFF_d2_logic > 0 && KPCOEF_expected > 0)
    for(int i1=0; i1<PCOEFF_d1_logic; ++i1) {
        for(int i2=0; i2<PCOEFF_d2_logic; ++i2) {
            for(int k=0; k<KPCOEF_expected; ++k) {
                double actual = PCOEFF_out_val[i1 + i2*(size_t)f_ldpco1 + k*(size_t)f_ldpco1*f_ldpco2];
                double expected = PCOEFF_expected_rm[i1*(size_t)PCOEFF_d2_logic*KPCOEF_expected + i2*KPCOEF_expected + k];
                EXPECT_NEAR(actual, expected, check_tol_matrix) << "PCOEFF_CM(" << i1 << "," << i2 << "," << k << ")";
            }
        }
    }
    
    int f_ldqco1 = std::max(1, current_max_mp); int f_ldqco2 = std::max(1, current_max_mp);
    if (QCOEFF_d1_logic > 0 && QCOEFF_d2_logic > 0 && KPCOEF_expected > 0)
    for(int i1=0; i1<QCOEFF_d1_logic; ++i1) { 
        for(int i2=0; i2<QCOEFF_d2_logic; ++i2) { 
            for(int k=0; k<KPCOEF_expected; ++k) {
                 double actual = QCOEFF_out_val[i1 + i2*(size_t)f_ldqco1 + k*(size_t)f_ldqco1*f_ldqco2];
                 double expected = QCOEFF_expected_rm[i1*(size_t)QCOEFF_d2_logic*KPCOEF_expected + i2*KPCOEF_expected + k];
                 EXPECT_NEAR(actual, expected, check_tol_matrix) << "QCOEFF_CM(" << i1 << "," << i2 << "," << k << ")";
            }
        }
    }
}

TEST_F(TB03ADTest, ParameterValidation) {
    N_param = 1; M_param = 1; P_param = 1; TOL_param = 0.0; LERI_param = 'L'; EQUIL_param = 'N';
    prepare_buffers_and_params(true); 
    int temp_nr; 
    INDEX_out_val.assign(std::max(1, (LERI_param == 'L' ? P_param : M_param)), 0); 
    PCOEFF_out_val.assign(std::max(1, (N_param+1)), 0.0); 
    QCOEFF_out_val.assign(std::max(1, (N_param+1)), 0.0);
    VCOEFF_out_val.assign(std::max(1, (N_param+1)), 0.0);

    EXPECT_EQ(slicot_tb03ad('X', EQUIL_param, N_param, M_param, P_param, A_test_buf.data(), LDA_c_param, B_test_buf.data(), LDB_c_param, C_test_buf.data(), LDC_c_param, D_test_buf_const_copy.data(), LDD_c_param, &temp_nr, INDEX_out_val.data(), PCOEFF_out_val.data(),1,1, QCOEFF_out_val.data(),1,1, VCOEFF_out_val.data(),1,1, TOL_param,1), -1);
    EXPECT_EQ(slicot_tb03ad(LERI_param, 'Y', N_param, M_param, P_param, A_test_buf.data(), LDA_c_param, B_test_buf.data(), LDB_c_param, C_test_buf.data(), LDC_c_param, D_test_buf_const_copy.data(), LDD_c_param, &temp_nr, INDEX_out_val.data(), PCOEFF_out_val.data(),1,1, QCOEFF_out_val.data(),1,1, VCOEFF_out_val.data(),1,1, TOL_param,1), -2);
    EXPECT_EQ(slicot_tb03ad(LERI_param, EQUIL_param, -1, M_param, P_param, nullptr, LDA_c_param, nullptr, LDB_c_param, nullptr, LDC_c_param, nullptr, LDD_c_param, &temp_nr, INDEX_out_val.data(), PCOEFF_out_val.data(),1,1, QCOEFF_out_val.data(),1,1, VCOEFF_out_val.data(),1,1, TOL_param,1), -3);
}

TEST_F(TB03ADTest, ZeroDimensionsNMP) {
    LERI_param = 'R'; EQUIL_param = 'N'; TOL_param = 0.0;
    N_param = 0; M_param = 0; P_param = 0;
    prepare_buffers_and_params(true); 
    NR_expected = 0;
    INDEX_expected.clear(); // For M=0, INDEX is empty.
    KPCOEF_expected = 1; // N=0 -> N+1 = 1. If INDEX is empty, max(INDEX)=0 effectively.
    PCOEFF_expected_rm.assign(1*1*1, 0.0); // M=0 -> porm=0. PCOEFF is 0x0x1. Test buffer might be 1.
    QCOEFF_expected_rm.assign(1*1*1, 0.0); // P=0, M=0 -> Q is 0x0x1. Test buffer might be 1.


    int info = slicot_tb03ad(LERI_param, EQUIL_param, N_param, M_param, P_param,
                             nullptr, LDA_c_param, nullptr, LDB_c_param, nullptr, LDC_c_param, nullptr, LDD_c_param,
                             &NR_out_val, INDEX_out_val.data(), 
                             PCOEFF_out_val.data(), PCOEFF_d1_logic, PCOEFF_d2_logic, 
                             QCOEFF_out_val.data(), QCOEFF_d1_logic, QCOEFF_d2_logic, 
                             VCOEFF_out_val.data(), VCOEFF_d1_logic, VCOEFF_d2_logic, 
                             TOL_param, 1);
    ASSERT_EQ(info, 0);
    ASSERT_EQ(NR_out_val, NR_expected);
    // For M=0, porm=0, INDEX_out_val is size 1 by test setup but not filled by Fortran.
    // Similar for PCOEFF, QCOEFF, VCOEFF: data dimensions are zero.
}
