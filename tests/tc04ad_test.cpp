#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::max_element, std::max
#include <iostream>
#include <iomanip>

#include "tc04ad.h"
#include "slicot_utils.h"

// Helper to print matrices
template<typename T>
void print_test_matrix(const std::string& name, const T* data, int rows, int cols, int ld_c, bool is_row_major) {
    if (!data) {
        std::cout << name << " is NULL" << std::endl;
        return;
    }
    std::cout << name << " (" << rows << "x" << cols << "), C_LD=" << ld_c << (is_row_major ? " RM" : " CM") << ":" << std::endl;
    for (int i = 0; i < rows; ++i) {
        std::cout << "  ";
        for (int j = 0; j < cols; ++j) {
            std::cout << std::fixed << std::setw(10) << std::setprecision(4) 
                      << (is_row_major ? data[i * ld_c + j] : data[i + j * ld_c]) << " ";
        }
        std::cout << std::endl;
    }
}


class Tc04adTest : public ::testing::Test {
protected:
    // From TC04AD.html example
    char leri = 'L';
    int m = 2;
    int p = 2;
    std::vector<int> index_in;
    
    std::vector<double> pcoeff_in_flat_rm_slice; 
    std::vector<double> qcoeff_in_flat_rm_slice; 

    // Expected outputs (from example data / Python script output)
    int n_exp = 4;
    double rcond_exp = 0.25;
    std::vector<double> A_exp_cm;
    std::vector<double> B_exp_cm; 
    std::vector<double> C_exp_cm; 
    std::vector<double> D_exp_cm; 

    // Buffers for C wrapper outputs
    int n_res = 0;
    double rcond_res = 0.0;
    std::vector<double> A_res_buffer;
    std::vector<double> B_res_buffer;
    std::vector<double> C_res_buffer;
    std::vector<double> D_res_buffer;
    
    int kpcoef_val = 0;

    double tol = 1e-4; // Adjusted from 1e-4 for some values if needed

    void SetUp() override {
        index_in = {2, 2}; 
        kpcoef_val = 0;
        if (!index_in.empty()) {
            for (int deg : index_in) kpcoef_val = std::max(kpcoef_val, deg);
        }
        kpcoef_val += 1; 

        // PCOEFF data (LERI='L' -> PxP slices = 2x2 slices, kpcoef=3)
        // From Python script (matches HTML example structure)
        // P(0,0,s) = 2s^2 + 3s + 1
        // P(0,1,s) = 4s^2 - 1s - 1
        // P(1,0,s) = 5s^2 + 7s - 6
        // P(1,1,s) = 3s^2 + 2s + 2
        // Stored as: P(0,0,k0..2), P(0,1,k0..2), P(1,0,k0..2), P(1,1,k0..2) if flat
        // Or slice by slice for 3D array view in C (row-major slices)
        pcoeff_in_flat_rm_slice = {
            // k=0 (coeffs for s^2)
            2.0, 4.0,  
            5.0, 3.0,  
            // k=1 (coeffs for s^1)
            3.0, -1.0, 
            7.0, 2.0,  
            // k=2 (coeffs for s^0)
            1.0, -1.0, 
            -6.0, 2.0  
        };

        // QCOEFF data (LERI='L' -> PxM slices = 2x2 slices, kpcoef=3)
        // Q(0,0,s) = 6s^2 - 1s + 5
        // Q(0,1,s) = 1s^2 + 7s + 5
        // Q(1,0,s) = 1s^2 + 1s + 1
        // Q(1,1,s) = 4s^2 + 1s - 1
        qcoeff_in_flat_rm_slice = {
            // k=0
            6.0, 1.0, 
            1.0, 4.0,
            // k=1
            -1.0, 7.0, 
            1.0, 1.0,
            // k=2
            5.0, 5.0, 
            1.0, -1.0
        };

        n_exp = 4; 
        rcond_exp = 0.25; // From HTML example & Python output
        // Expected A, B, C, D from Python output (column-major)
        A_exp_cm = {
             0.0000,  1.0000,  0.0000,  0.0000,
             0.57142857,  1.0000, -2.0000,  0.78571429,
             0.0000,  0.0000,  0.0000,  1.0000,
            -0.42857143, -1.0000,  2.0000, -1.71428571
        };
        B_exp_cm = { 
             8.0000,  4.0000, -9.0000,  4.0000,
             3.85714286,  4.0000,  5.0000, -5.07142857
        };
        C_exp_cm = { 
             0.0000,       0.0000,
            -0.21428571,  0.35714286,
             0.0000,       0.0000,
             0.28571429, -0.14285714
        };
        D_exp_cm = { 
            -1.0000,       2.0000,
             0.92857143, -0.21428571
        };
        
        A_res_buffer.resize((size_t)n_exp * n_exp);
        B_res_buffer.resize((size_t)n_exp * m); 
        C_res_buffer.resize((size_t)p * n_exp); 
        D_res_buffer.resize((size_t)p * m);     
    }
};

TEST_F(Tc04adTest, ExampleColMajor) {
    int porm = (leri == 'L') ? p : m;
    int max_mp = std::max(1, std::max(m,p));

    std::vector<double> pcoeff_test_cm((size_t)porm * porm * kpcoef_val);
    
    size_t qcoeff_f_buf_elems;
    int q_f_ld1, q_f_ld2; // Fortran leading dimensions for QCOEFF
    if (leri == 'L') {
        q_f_ld1 = p; q_f_ld2 = m;
    } else { // LERI == 'R', QCOEFF buffer is MAXMP x MAXMP for Fortran
        q_f_ld1 = max_mp; q_f_ld2 = max_mp;
    }
    qcoeff_f_buf_elems = (size_t)q_f_ld1 * q_f_ld2 * kpcoef_val;
    std::vector<double> qcoeff_test_cm(qcoeff_f_buf_elems);


    for(int k=0; k<kpcoef_val; ++k) {
        slicot_transpose_to_fortran_with_ld(
            pcoeff_in_flat_rm_slice.data() + k * (size_t)porm * porm,
            pcoeff_test_cm.data() + k * (size_t)porm * porm, // Fortran LDs are porm, porm
            porm, porm, porm, porm, sizeof(double)
        );
        if (p > 0 && m > 0) { 
             slicot_transpose_to_fortran_with_ld(
                qcoeff_in_flat_rm_slice.data() + k * (size_t)p * m,
                qcoeff_test_cm.data() + k * (size_t)q_f_ld1 * q_f_ld2, // Use Fortran LDs for Q
                p, m, m, q_f_ld1, sizeof(double) 
            );
        }
    }
    
    // C leading dimensions for column-major call
    int ldpcoeff_c_r_cm = porm; int ldpcoeff_c_c_cm = porm;
    int ldqcoeff_c_r_cm = q_f_ld1; int ldqcoeff_c_c_cm = q_f_ld2; 


    int lda_c_test = n_exp > 0 ? n_exp : 1;
    int ldb_c_test = n_exp > 0 ? n_exp : 1; // Fortran B is N x MAXMP, so C LD for B is N
    int ldc_c_test = max_mp > 0 ? max_mp : 1; // Fortran C is MAXMP x N
    int ldd_c_test = max_mp > 0 ? max_mp : 1; // Fortran D is MAXMP x MAXMP


    int info = slicot_tc04ad(leri, m, p, index_in.data(),
                             pcoeff_test_cm.data(), ldpcoeff_c_r_cm, ldpcoeff_c_c_cm,
                             qcoeff_test_cm.data(), ldqcoeff_c_r_cm, ldqcoeff_c_c_cm,
                             &n_res, &rcond_res,
                             A_res_buffer.data(), lda_c_test,
                             B_res_buffer.data(), ldb_c_test,
                             C_res_buffer.data(), ldc_c_test,
                             D_res_buffer.data(), ldd_c_test,
                             0 /*row_major=false*/);

    ASSERT_EQ(info, 0);
    EXPECT_EQ(n_res, n_exp);
    EXPECT_NEAR(rcond_res, rcond_exp, tol);

    for(size_t i=0; i<A_exp_cm.size(); ++i) EXPECT_NEAR(A_res_buffer[i], A_exp_cm[i], tol) << "A mismatch at " << i;
    for(size_t i=0; i < (size_t)n_exp * m; ++i) EXPECT_NEAR(B_res_buffer[i], B_exp_cm[i], tol) << "B mismatch at " << i;
    for(size_t i=0; i < (size_t)p * n_exp; ++i) EXPECT_NEAR(C_res_buffer[i], C_exp_cm[i], tol) << "C mismatch at " << i;
    for(size_t i=0; i < (size_t)p * m; ++i) EXPECT_NEAR(D_res_buffer[i], D_exp_cm[i], tol) << "D mismatch at " << i;
}

TEST_F(Tc04adTest, ExampleRowMajor) {
    int porm = (leri == 'L') ? p : m;
    int max_mp = std::max(1, std::max(m,p));

    std::vector<double> pcoeff_test_rm = pcoeff_in_flat_rm_slice;
    std::vector<double> qcoeff_test_rm = qcoeff_in_flat_rm_slice; 
    
    // C leading dimensions for row-major call
    int ldpcoeff_c_r_rm = porm; int ldpcoeff_c_c_rm = porm;
    int ldqcoeff_c_r_rm = (leri == 'L') ? p : max_mp; 
    int ldqcoeff_c_c_rm = (leri == 'L') ? m : max_mp; 


    int lda_c_test = n_exp > 0 ? n_exp : 1; 
    int ldb_c_test = m > 0 ? m : 1;         
    int ldc_c_test = n_exp > 0 ? n_exp : 1; 
    int ldd_c_test = m > 0 ? m : 1;         

    int info = slicot_tc04ad(leri, m, p, index_in.data(),
                             pcoeff_test_rm.data(), ldpcoeff_c_r_rm, ldpcoeff_c_c_rm,
                             qcoeff_test_rm.data(), ldqcoeff_c_r_rm, ldqcoeff_c_c_rm,
                             &n_res, &rcond_res,
                             A_res_buffer.data(), lda_c_test,
                             B_res_buffer.data(), ldb_c_test,
                             C_res_buffer.data(), ldc_c_test,
                             D_res_buffer.data(), ldd_c_test,
                             1 /*row_major=true*/);

    ASSERT_EQ(info, 0);
    EXPECT_EQ(n_res, n_exp);
    EXPECT_NEAR(rcond_res, rcond_exp, tol);

    std::vector<double> A_exp_rm(A_exp_cm.size());
    std::vector<double> B_exp_rm(B_exp_cm.size());
    std::vector<double> C_exp_rm(C_exp_cm.size());
    std::vector<double> D_exp_rm(D_exp_cm.size());

    if (n_exp > 0) slicot_transpose_to_c(A_exp_cm.data(), A_exp_rm.data(), n_exp, n_exp, sizeof(double));
    if (n_exp > 0 && m > 0) slicot_transpose_to_c(B_exp_cm.data(), B_exp_rm.data(), n_exp, m, sizeof(double));
    if (p > 0 && n_exp > 0) slicot_transpose_to_c(C_exp_cm.data(), C_exp_rm.data(), p, n_exp, sizeof(double));
    if (p > 0 && m > 0) slicot_transpose_to_c(D_exp_cm.data(), D_exp_rm.data(), p, m, sizeof(double));

    for(size_t i=0; i<A_exp_rm.size(); ++i) EXPECT_NEAR(A_res_buffer[i], A_exp_rm[i], tol) << "A RM mismatch at " << i;
    for(size_t i=0; i<B_exp_rm.size(); ++i) EXPECT_NEAR(B_res_buffer[i], B_exp_rm[i], tol) << "B RM mismatch at " << i;
    for(size_t i=0; i<C_exp_rm.size(); ++i) EXPECT_NEAR(C_res_buffer[i], C_exp_rm[i], tol) << "C RM mismatch at " << i;
    for(size_t i=0; i<D_exp_rm.size(); ++i) EXPECT_NEAR(D_res_buffer[i], D_exp_rm[i], tol) << "D RM mismatch at " << i;
}


TEST_F(Tc04adTest, ParameterValidation) {
    std::vector<int> dummy_idx(std::max(1,std::max(m,p)), 0); 
    std::vector<double> dummy_pcoeff(1);
    std::vector<double> dummy_qcoeff(1);
    int dummy_n; double dummy_rcond;
    std::vector<double> dummy_a(1), dummy_b(1), dummy_c(1), dummy_d(1);
    int porm_val = (leri == 'L') ? p : m;
    if (porm_val == 0) porm_val = 1; 
    int max_mp_val = std::max(1, std::max(m,p));


    EXPECT_EQ(slicot_tc04ad('X', m,p,dummy_idx.data(), dummy_pcoeff.data(),porm_val,porm_val, dummy_qcoeff.data(),max_mp_val,max_mp_val, &dummy_n,&dummy_rcond, dummy_a.data(),1,dummy_b.data(),1,dummy_c.data(),1,dummy_d.data(),1,0), -1);
    EXPECT_EQ(slicot_tc04ad(leri, -1,p,dummy_idx.data(), dummy_pcoeff.data(),porm_val,porm_val, dummy_qcoeff.data(),max_mp_val,max_mp_val, &dummy_n,&dummy_rcond, dummy_a.data(),1,dummy_b.data(),1,dummy_c.data(),1,dummy_d.data(),1,0), -2);
}

TEST_F(Tc04adTest, ZeroNOutput) { 
    char leri_z = 'L';
    int m_z = 1, p_z = 1;
    std::vector<int> index_z = {0}; 
    int kpcoef_z = 1; 
    std::vector<double> pcoeff_z = {1.0}; 
    std::vector<double> qcoeff_z = {2.0}; 
    
    int n_res_z; double rcond_res_z;
    std::vector<double> d_res_z( (size_t)p_z * m_z );
    if (d_res_z.empty() && p_z > 0 && m_z > 0) d_res_z.resize(1); 

    int ldp_c_r = p_z > 0 ? p_z : 1; int ldp_c_c = p_z > 0 ? p_z : 1;
    int ldq_c_r = p_z > 0 ? p_z : 1; int ldq_c_c = m_z > 0 ? m_z : 1;


    int info = slicot_tc04ad(leri_z, m_z, p_z, index_z.data(),
                             pcoeff_z.data(), ldp_c_r, ldp_c_c, 
                             qcoeff_z.data(), ldq_c_r, ldq_c_c,
                             &n_res_z, &rcond_res_z,
                             nullptr, 1, 
                             nullptr, 1, 
                             nullptr, 1, 
                             d_res_z.empty() ? nullptr : d_res_z.data(), (p_z > 0 ? p_z : 1), 
                             0 /*row_major=false*/);
    ASSERT_EQ(info, 0);
    EXPECT_EQ(n_res_z, 0);
    if (p_z > 0 && m_z > 0 && !d_res_z.empty()) {
        EXPECT_NEAR(d_res_z[0], 2.0, tol);
    }
}
