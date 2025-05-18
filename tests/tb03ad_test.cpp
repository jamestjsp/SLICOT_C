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

void print_3d_matrix_tb03ad(const std::string& name, const double* mat, int d1, int d2, int d3, int ld1_c, int ld2_c, bool is_row_major) {
    if (!mat) { std::cout << name << " is NULL\n"; return; }
    std::cout << name << " (" << d1 << "x" << d2 << "x" << d3 << ") "
              << (is_row_major ? "RowMajor" : "ColMajor")
              << ", LD1_C=" << ld1_c << ", LD2_C=" << ld2_c << ":\n";
    if (d1 == 0 || d2 == 0 || d3 == 0) { std::cout << "  (empty or zero dim)\n"; return; }
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < d1; ++i) {
        for (int j = 0; j < d2; ++j) {
            std::cout << "  " << name << "[" << i << "][" << j << "][:] = [";
            for (int k = 0; k < d3; ++k) {
                // Assuming C-style flat row-major: mat[i*ld2_c*d3 + j*d3 + k]
                // This is how the wrapper is designed to output flat C arrays.
                double val = mat[i * ld2_c * d3 + j * d3 + k];
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
    std::vector<double> A_test_buf, B_test_buf, C_test_buf, D_test_buf;
    
    int NR_out_val;
    std::vector<int> INDEX_out_val;
    std::vector<double> PCOEFF_out_val, QCOEFF_out_val, VCOEFF_out_val;

    // Expected
    int NR_expected;
    std::vector<int> INDEX_expected;
    std::vector<double> A_min_expected_rm, B_min_expected_rm, C_min_expected_rm;
    std::vector<double> PCOEFF_expected_rm, QCOEFF_expected_rm; // VCOEFF not in example
    int info_expected;

    // C Leading dimensions for test buffers
    int LDA_c, LDB_c, LDC_c, LDD_c;
    int LDPCO1_c, LDPCO2_c, LDQCO1_c, LDQCO2_c, LDVCO1_c, LDVCO2_c;
    int KPCOEF_expected;


    double check_tol_matrix = 1e-4;

    void setup_slicot_example_data_right() {
        LERI_param = 'R'; EQUIL_param = 'N';
        N_param = 3; M_param = 1; P_param = 2; TOL_param = 0.0;

        A_data_rm = {1.0, 2.0, 0.0, 4.0, -1.0, 0.0, 0.0, 0.0, 1.0}; // 3x3
        B_data_rm = {1.0, 0.0, 1.0}; // 3x1
        C_data_rm = {0.0, 1.0, -1.0, 0.0, 0.0, 1.0}; // 2x3
        D_data_rm = {0.0, 1.0}; // 2x1

        NR_expected = 3;
        A_min_expected_rm = {1.0000, -1.4142, 0.0000,
                            -2.8284, -1.0000, 2.8284,
                             0.0000, 1.4142, 1.0000}; // 3x3
        B_min_expected_rm = {-1.4142, 0.0000, 0.0000}; // 3x1
        C_min_expected_rm = {0.7071, 1.0000, 0.7071,
                            -0.7071, 0.0000, -0.7071}; // 2x3
        
        INDEX_expected = {3}; // M=1 for LERI='R'
        KPCOEF_expected = INDEX_expected[0] + 1; // 3+1 = 4

        // PCOEFF is M x M x KPCOEF = 1 x 1 x 4
        PCOEFF_expected_rm = {0.1768, -0.1768, -1.5910, 1.5910};
        
        // QCOEFF is P x M x KPCOEF = 2 x 1 x 4
        QCOEFF_expected_rm = {
            0.0000, -0.1768, 0.7071, 0.8839, // Q(0,0,k)
            0.1768, 0.0000, -1.5910, 0.0000  // Q(1,0,k)
        };
        info_expected = 0;
    }

    void prepare_buffers(bool is_row_major) {
        A_test_buf = A_data_rm;
        B_test_buf = B_data_rm;
        C_test_buf = C_data_rm;
        D_test_buf = D_data_rm;

        int porm = (LERI_param == 'L') ? P_param : M_param;
        int porp = (LERI_param == 'L') ? M_param : P_param;
        
        INDEX_out_val.assign(porm > 0 ? porm : 1, 0); // Min size 1 if porm is 0

        LDPCO1_c = porm; LDPCO2_c = porm;
        PCOEFF_out_val.assign((size_t)LDPCO1_c * LDPCO2_c * (N_param + 1), 0.0);

        if (LERI_param == 'L') {
            LDQCO1_c = porm; LDQCO2_c = porp;
        } else { // LERI_param == 'R'
            LDQCO1_c = porp; LDQCO2_c = porm; // Q is P x M
        }
        QCOEFF_out_val.assign((size_t)LDQCO1_c * LDQCO2_c * (N_param + 1), 0.0);
        
        LDVCO1_c = porm; LDVCO2_c = N_param; // Use N_param for allocation, actual data uses NR
                                            // If NR_out_val becomes 0, this dim might be 0.
                                            // For safety, ensure N_param for allocation if porm > 0.
        if (porm == 0 || N_param == 0) { // If porm or N is 0, VCOEFF is effectively empty for data part
             VCOEFF_out_val.assign((size_t)LDVCO1_c * std::max(1,LDVCO2_c) * (N_param + 1), 0.0); // ensure some allocation if N_param+1 > 0
        } else {
            VCOEFF_out_val.assign((size_t)LDVCO1_c * LDVCO2_c * (N_param + 1), 0.0);
        }


        if (is_row_major) {
            LDA_c = N_param;
            LDB_c = M_param;
            LDC_c = N_param; // C is PxN, so N columns
            LDD_c = M_param; // D is PxM, so M columns
        } else { // Column Major
            LDA_c = N_param;
            LDB_c = N_param;
            // Fortran LDC for C matrix (P x N content) is MAX(1,M,P)
            LDC_c = std::max(1, std::max(M_param, P_param));
            // Fortran LDD for D matrix (P x M content) is MAX(1,M,P)
            LDD_c = std::max(1, std::max(M_param, P_param));
        }
        // Ensure minimum leading dimensions of 1
        LDA_c = std::max(1, LDA_c);
        LDB_c = std::max(1, LDB_c);
        LDC_c = std::max(1, LDC_c);
        LDD_c = std::max(1, LDD_c);
    }
};


TEST_F(TB03ADTest, SlicotDocExampleRightRowMajor) {
    setup_slicot_example_data_right();
    prepare_buffers(true /*row_major*/);

    int info = slicot_tb03ad(LERI_param, EQUIL_param, N_param, M_param, P_param,
                             A_test_buf.data(), LDA_c,
                             B_test_buf.data(), LDB_c,
                             C_test_buf.data(), LDC_c,
                             D_test_buf.data(), LDD_c,
                             &NR_out_val, INDEX_out_val.data(),
                             PCOEFF_out_val.data(), LDPCO1_c, LDPCO2_c,
                             QCOEFF_out_val.data(), LDQCO1_c, LDQCO2_c,
                             VCOEFF_out_val.data(), LDVCO1_c, LDVCO2_c, // Pass N_param for LDVCO2_c
                             TOL_param, 1 /*row_major*/);

    ASSERT_EQ(info, info_expected);
    ASSERT_EQ(NR_out_val, NR_expected);

    // Verify Amin, Bmin, Cmin (A_test_buf, B_test_buf, C_test_buf were modified)
    for(int i=0; i<NR_expected; ++i) for(int j=0; j<NR_expected; ++j) EXPECT_NEAR(A_test_buf[i*LDA_c+j], A_min_expected_rm[i*NR_expected+j], check_tol_matrix);
    for(int i=0; i<NR_expected; ++i) for(int j=0; j<M_param; ++j) EXPECT_NEAR(B_test_buf[i*LDB_c+j], B_min_expected_rm[i*M_param+j], check_tol_matrix);
    for(int i=0; i<P_param; ++i) for(int j=0; j<NR_expected; ++j) EXPECT_NEAR(C_test_buf[i*LDC_c+j], C_min_expected_rm[i*NR_expected+j], check_tol_matrix);
    
    // Verify INDEX
    for(size_t i=0; i<INDEX_expected.size(); ++i) EXPECT_EQ(INDEX_out_val[i], INDEX_expected[i]);

    // Verify PCOEFF (M x M x KPCOEF_expected for LERI='R')
    // LDPCO1_c = M_param, LDPCO2_c = M_param
    for(int i=0; i<M_param; ++i) { // dim1
        for(int j=0; j<M_param; ++j) { // dim2
            for(int k=0; k<KPCOEF_expected; ++k) { // dim3 (coeffs)
                double actual = PCOEFF_out_val[i*LDPCO2_c*KPCOEF_expected + j*KPCOEF_expected + k];
                double expected = PCOEFF_expected_rm[i*M_param*KPCOEF_expected + j*KPCOEF_expected + k];
                EXPECT_NEAR(actual, expected, check_tol_matrix) << "PCOEFF(" << i << "," << j << "," << k << ")";
            }
        }
    }
    
    // Verify QCOEFF (P x M x KPCOEF_expected for LERI='R')
    // LDQCO1_c = P_param, LDQCO2_c = M_param
     for(int i=0; i<P_param; ++i) { // dim1
        for(int j=0; j<M_param; ++j) { // dim2
            for(int k=0; k<KPCOEF_expected; ++k) { // dim3 (coeffs)
                double actual = QCOEFF_out_val[i*LDQCO2_c*KPCOEF_expected + j*KPCOEF_expected + k];
                double expected = QCOEFF_expected_rm[i*M_param*KPCOEF_expected + j*KPCOEF_expected + k];
                EXPECT_NEAR(actual, expected, check_tol_matrix) << "QCOEFF(" << i << "," << j << "," << k << ")";
            }
        }
    }
}


TEST_F(TB03ADTest, SlicotDocExampleRightColMajor) {
    setup_slicot_example_data_right();
    prepare_buffers(false /*col_major*/);
    
    // For col-major, we need to transpose input data before calling
    std::vector<double> a_cm_in((size_t)LDA_c * N_param), b_cm_in((size_t)LDB_c * std::max(1,M_param)), 
                        c_cm_in((size_t)LDC_c * N_param), d_cm_in((size_t)LDD_c * std::max(1,M_param));

    if (N_param > 0) slicot_transpose_to_fortran_with_ld(A_data_rm.data(), a_cm_in.data(), N_param, N_param, N_param, LDA_c, sizeof(double));
    if (N_param > 0 && M_param > 0) slicot_transpose_to_fortran_with_ld(B_data_rm.data(), b_cm_in.data(), N_param, M_param, M_param, LDB_c, sizeof(double));
    if (P_param > 0 && N_param > 0) slicot_transpose_to_fortran_with_ld(C_data_rm.data(), c_cm_in.data(), P_param, N_param, N_param, LDC_c, sizeof(double));
    if (P_param > 0 && M_param > 0) slicot_transpose_to_fortran_with_ld(D_data_rm.data(), d_cm_in.data(), P_param, M_param, M_param, LDD_c, sizeof(double));


    int info = slicot_tb03ad(LERI_param, EQUIL_param, N_param, M_param, P_param,
                             a_cm_in.data(), LDA_c,
                             b_cm_in.data(), LDB_c,
                             c_cm_in.data(), LDC_c,
                             d_cm_in.data(), LDD_c,
                             &NR_out_val, INDEX_out_val.data(),
                             PCOEFF_out_val.data(), LDPCO1_c, LDPCO2_c,
                             QCOEFF_out_val.data(), LDQCO1_c, LDQCO2_c,
                             VCOEFF_out_val.data(), LDVCO1_c, LDVCO2_c,
                             TOL_param, 0 /*col_major*/);

    ASSERT_EQ(info, info_expected);
    ASSERT_EQ(NR_out_val, NR_expected);

    // Verify Amin, Bmin, Cmin (a_cm_in, b_cm_in, c_cm_in were modified by wrapper)
    // Expected results are RM, so transpose them to CM for comparison
    std::vector<double> a_exp_cm((size_t)NR_expected*NR_expected), b_exp_cm((size_t)NR_expected*M_param), c_exp_cm((size_t)P_param*NR_expected);
    if (NR_expected>0) slicot_transpose_to_fortran_with_ld(A_min_expected_rm.data(), a_exp_cm.data(), NR_expected, NR_expected, NR_expected, NR_expected, sizeof(double));
    if (NR_expected>0 && M_param >0) slicot_transpose_to_fortran_with_ld(B_min_expected_rm.data(), b_exp_cm.data(), NR_expected, M_param, M_param, NR_expected, sizeof(double));
    if (P_param > 0 && NR_expected >0) slicot_transpose_to_fortran_with_ld(C_min_expected_rm.data(), c_exp_cm.data(), P_param, NR_expected, NR_expected, P_param, sizeof(double));

    for(int j=0; j<NR_expected; ++j) for(int i=0; i<NR_expected; ++i) EXPECT_NEAR(a_cm_in[i+j*LDA_c], a_exp_cm[i+j*NR_expected], check_tol_matrix);
    if (M_param > 0) for(int j=0; j<M_param; ++j) for(int i=0; i<NR_expected; ++i) EXPECT_NEAR(b_cm_in[i+j*LDB_c], b_exp_cm[i+j*NR_expected], check_tol_matrix);
    if (P_param > 0) for(int j=0; j<NR_expected; ++j) for(int i=0; i<P_param; ++i) EXPECT_NEAR(c_cm_in[i+j*LDC_c], c_exp_cm[i+j*P_param], check_tol_matrix);
    
    for(size_t i=0; i<INDEX_expected.size(); ++i) EXPECT_EQ(INDEX_out_val[i], INDEX_expected[i]);

    // PCOEFF_out_val is already column-major from the wrapper. Transpose PCOEFF_expected_rm.
    // PCOEFF_expected_rm is flat 1x1x4. PCOEFF_out_val is also flat 1x1x4.
    // For col-major C output, PCOEFF_out_val[k*LDPCO1_c*LDPCO2_c + j*LDPCO1_c + i]
    // Here LDPCO1_c=M_param=1, LDPCO2_c=M_param=1.
    // PCOEFF_out_val[k*1*1 + j*1 + i]
    std::vector<double> pcoeff_exp_cm((size_t)M_param*M_param*KPCOEF_expected);
    // transpose_3d_c2f(PCOEFF_expected_rm.data(), pcoeff_exp_cm.data(), M_param, M_param, KPCOEF_expected, M_param, M_param, M_param, M_param);
    // Since it's 1x1xKPCOEF, flat RM and CM are the same.
    for(size_t k=0; k < PCOEFF_expected_rm.size(); ++k) {
        EXPECT_NEAR(PCOEFF_out_val[k], PCOEFF_expected_rm[k], check_tol_matrix) << "PCOEFF_CM k=" << k;
    }


    // QCOEFF_out_val is P x M x KPCOEF_expected, col-major from wrapper.
    // QCOEFF_expected_rm is P x M x KPCOEF_expected, row-major.
    // LDQCO1_c = P, LDQCO2_c = M for col-major C output structure.
    std::vector<double> qcoeff_exp_cm((size_t)P_param*M_param*KPCOEF_expected);
    // transpose_3d_c2f(QCOEFF_expected_rm.data(), qcoeff_exp_cm.data(), P_param, M_param, KPCOEF_expected, P_param, M_param, P_param, M_param);
    // For QCOEFF_out_val[k*LDQCO1_c*LDQCO2_c + j*LDQCO1_c + i]
    // For QCOEFF_expected_rm[i*M_param*KPCOEF_expected + j*KPCOEF_expected + k]
    for(int i_p=0; i_p<P_param; ++i_p) { // P
        for(int j_m=0; j_m<M_param; ++j_m) { // M
            for(int k_c=0; k_c<KPCOEF_expected; ++k_c) { // K
                 double actual = QCOEFF_out_val[k_c*LDQCO1_c*LDQCO2_c + j_m*LDQCO1_c + i_p];
                 double expected = QCOEFF_expected_rm[i_p*M_param*KPCOEF_expected + j_m*KPCOEF_expected + k_c];
                 EXPECT_NEAR(actual, expected, check_tol_matrix) << "QCOEFF_CM(" << i_p << "," << j_m << "," << k_c << ")";
            }
        }
    }
}


TEST_F(TB03ADTest, ParameterValidation) {
    N_param = 1; M_param = 1; P_param = 1; TOL_param = 0.0; LERI_param = 'L'; EQUIL_param = 'N';
    prepare_buffers(true);
    int temp_nr; std::vector<int> temp_idx(1); std::vector<double> temp_p(N_param+1), temp_q(N_param+1), temp_v(N_param+1);

    EXPECT_EQ(slicot_tb03ad('X', EQUIL_param, N_param, M_param, P_param, A_test_buf.data(), LDA_c, B_test_buf.data(), LDB_c, C_test_buf.data(), LDC_c, D_test_buf.data(), LDD_c, &temp_nr, temp_idx.data(), temp_p.data(),1,1, temp_q.data(),1,1, temp_v.data(),1,1, TOL_param,1), -1);
    EXPECT_EQ(slicot_tb03ad(LERI_param, 'Y', N_param, M_param, P_param, A_test_buf.data(), LDA_c, B_test_buf.data(), LDB_c, C_test_buf.data(), LDC_c, D_test_buf.data(), LDD_c, &temp_nr, temp_idx.data(), temp_p.data(),1,1, temp_q.data(),1,1, temp_v.data(),1,1, TOL_param,1), -2);
    EXPECT_EQ(slicot_tb03ad(LERI_param, EQUIL_param, -1, M_param, P_param, A_test_buf.data(), LDA_c, B_test_buf.data(), LDB_c, C_test_buf.data(), LDC_c, D_test_buf.data(), LDD_c, &temp_nr, temp_idx.data(), temp_p.data(),1,1, temp_q.data(),1,1, temp_v.data(),1,1, TOL_param,1), -3);
    // ... more validation checks
}

TEST_F(TB03ADTest, ZeroDimensionsNMP) {
    LERI_param = 'R'; EQUIL_param = 'N'; TOL_param = 0.0;
    N_param = 0; M_param = 0; P_param = 0;
    prepare_buffers(true); // row_major
    NR_expected = 0;

    int info = slicot_tb03ad(LERI_param, EQUIL_param, N_param, M_param, P_param,
                             nullptr, LDA_c, nullptr, LDB_c, nullptr, LDC_c, nullptr, LDD_c,
                             &NR_out_val, INDEX_out_val.data(), // porm=0, index_out_val size 1
                             PCOEFF_out_val.data(), LDPCO1_c, LDPCO2_c, // porm=0, pcoeff_out_val size (N+1)
                             QCOEFF_out_val.data(), LDQCO1_c, LDQCO2_c, // porm=0, porp=0
                             VCOEFF_out_val.data(), LDVCO1_c, LDVCO2_c, // porm=0
                             TOL_param, 1);
    ASSERT_EQ(info, 0);
    ASSERT_EQ(NR_out_val, NR_expected);
    // INDEX for M=0 should be empty or not accessed. Wrapper gives size 1.
    // PCOEFF, QCOEFF, VCOEFF for M=0, P=0 should be effectively empty. Wrapper gives size N+1.
    // This case implies P(s), Q(s) are 0x0 or 1x1 constants if D was 0x0.
    // The Fortran routine should handle this gracefully.
}

