#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // For std::max, std::min
#include <iostream>
#include <iomanip>


#include "tb04ad.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions if needed for setup/verification
#include "test_utils.h"   // For load_test_data_from_csv (not used here, data embedded)
#include "test_config.h"  // Include CMake-generated configuration for data path (not used here)

// Helper to print matrices for debugging
template<typename T>
void print_matrix(const std::string& name, const T* data, int rows, int cols, int ld, bool row_major_print) {
    std::cout << name << " (" << rows << "x" << cols << ", ld=" << ld << (row_major_print ? ", RM" : ", CM") << "):\n";
    if (!data) {
        std::cout << "  NULL\n";
        return;
    }
    for (int i = 0; i < rows; ++i) {
        std::cout << "  ";
        for (int j = 0; j < cols; ++j) {
            if (row_major_print) {
                std::cout << std::fixed << std::setw(10) << std::setprecision(4) << data[i * ld + j] << " ";
            } else {
                std::cout << std::fixed << std::setw(10) << std::setprecision(4) << data[i + j * ld] << " ";
            }
        }
        std::cout << "\n";
    }
}


// --- Test Fixture ---
class Tb04adTest : public ::testing::Test {
protected:
    // Parameters from TB04AD.html example, to be cross-verified with Python output
    char rowcol = 'R';
    int n = 3;
    int m = 2;
    int p = 2;
    double tol1 = 0.0;
    double tol2 = 0.0;

    // Input matrices (column-major for direct Fortran comparison)
    std::vector<double> A_in_cm;
    std::vector<double> B_in_cm;
    std::vector<double> C_in_cm;
    std::vector<double> D_in_cm;

    // Expected output matrices (column-major) - VALUES UPDATED FROM PYTHON SCRIPT OUTPUT
    std::vector<double> A_exp_cm;
    std::vector<double> B_exp_cm;
    std::vector<double> C_exp_cm;
    int NR_exp = 3; // From Python output
    std::vector<int> INDEX_exp; // From Python output
    std::vector<double> DCOEFF_exp_flat; 
    std::vector<double> UCOEFF_exp_flat; 

    int kmax_coeffs_exp = 0; // Will be derived from INDEX_exp

    // Output variables for the wrapper
    int NR_res = 0;
    std::vector<int> INDEX_res;
    std::vector<double> A_res_cm_test_buffer; 
    std::vector<double> B_res_cm_test_buffer; 
    std::vector<double> C_res_cm_test_buffer; 
    std::vector<double> DCOEFF_res_flat;
    std::vector<double> UCOEFF_res_flat;


    // Leading dimensions (for C arrays)
    int lda_c, ldb_c, ldc_c, ldd_c;
    int lddcoe_c, lduco1_c, lduco2_c;


    double check_tol = 1e-4; // Adjusted tolerance slightly if needed

    void SetUp() override {
        // Initialize input data (from TB04AD.html example / Python script)
        A_in_cm = {-1.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, -3.0};
        B_in_cm = {0.0, 1.0, -1.0, 1.0, 1.0, 0.0};
        C_in_cm = {0.0, 1.0,  1.0, 0.0,  1.0, 1.0}; // C(P,N) -> C(2,3)
        D_in_cm = {1.0, 0.0, 0.0, 1.0};

        // Expected outputs - UPDATED BASED ON PYTHON SCRIPT OUTPUT
        NR_exp = 3;
        // Ar_out from Python:
        // [[-2.5       -0.28867513 -0.40824829]
        //  [-0.28867513 -1.5        -0.70710678]
        //  [-0.40824829 -0.70710678 -2.        ]]
        A_exp_cm = {-2.5, -0.28867513, -0.40824829, 
                    -0.28867513, -1.5, -0.70710678, 
                    -0.40824829, -0.70710678, -2.0};
        
        // Br_out from Python:
        // [[-1.41421356 -0.70710678]
        //  [ 0.          1.22474487]
        //  [ 0.          0.        ]]
        B_exp_cm = {-1.41421356, 0.0, 0.0, 
                    -0.70710678, 1.22474487, 0.0};
        
        // Cr_out from Python:
        // [[ 0.00000000e+00  8.16496581e-01  1.15470054e+00]
        //  [ 7.07106781e-01  1.22474487e+00 -2.46633885e-16]]
        // C(P,NR) = C(2,3) column major: C(0,0),C(1,0), C(0,1),C(1,1), C(0,2),C(1,2)
        C_exp_cm = {0.0, 0.707106781, 
                    0.816496581, 1.22474487, 
                    1.15470054, -2.46633885e-16};

        INDEX_exp = {2, 2}; // From Python output
        
        kmax_coeffs_exp = 0;
        if (!INDEX_exp.empty()) {
            for(int deg : INDEX_exp) kmax_coeffs_exp = std::max(kmax_coeffs_exp, deg);
        }
        kmax_coeffs_exp +=1; // Number of coefficients is degree + 1. max(2,2)+1 = 3

        // DCOEFF_out from Python (shape: (2, 3)):
        // [[1. 5. 6.]
        //  [1. 4. 3.]]
        // Fortran DCOEFF(porm_f=2, kmax_coeffs=3)
        // Flattened column-major: D(0,0), D(1,0), D(0,1), D(1,1), D(0,2), D(1,2)
        DCOEFF_exp_flat = {
            1.0, 1.0, 
            5.0, 4.0, 
            6.0, 3.0
        };

        // UCOEFF_out from Python (shape: (2, 2, 3)):
        // [[[ 1.  5.  7.]   # U(0,0,k)
        //   [ 0.  1.  3.]]  # U(0,1,k)
        //
        //  [[ 0. -1. -1.]   # U(1,0,k)
        //   [ 1.  5.  6.]]] # U(1,1,k)
        // Fortran UCOEFF(porm_f=2, porp_f=2, kmax_coeffs=3)
        // Flattened column-major (porm varies fastest, then porp, then k_coeff_idx):
        // k_idx=0: U(0,0,0), U(1,0,0), U(0,1,0), U(1,1,0)
        // k_idx=1: U(0,0,1), U(1,0,1), U(0,1,1), U(1,1,1)
        // k_idx=2: U(0,0,2), U(1,0,2), U(0,1,2), U(1,1,2)
        UCOEFF_exp_flat = {
            // k_idx = 0 (coeffs for s^2, since max_degree is 2)
            1.0,  // U(0,0,0)
            0.0,  // U(1,0,0)
            0.0,  // U(0,1,0)
            1.0,  // U(1,1,0)
            // k_idx = 1 (coeffs for s^1)
            5.0,  // U(0,0,1)
           -1.0,  // U(1,0,1)
            1.0,  // U(0,1,1)
            5.0,  // U(1,1,1)
            // k_idx = 2 (coeffs for s^0)
            7.0,  // U(0,0,2)
           -1.0,  // U(1,0,2)
            3.0,  // U(0,1,2)
            6.0   // U(1,1,2)
        };

        int porm_val = (rowcol == 'R' ? p : m);
        int porp_val = (rowcol == 'R' ? m : p);

        A_res_cm_test_buffer.resize( (size_t)n * n); 
        B_res_cm_test_buffer.resize( (size_t)n * m); 
        C_res_cm_test_buffer.resize( (size_t)p * n); 
        if (porm_val > 0) INDEX_res.resize(porm_val); else INDEX_res.clear(); 
        if (porm_val > 0 && kmax_coeffs_exp > 0) DCOEFF_res_flat.resize( (size_t)porm_val * kmax_coeffs_exp); else DCOEFF_res_flat.clear();
        // UCOEFF must be dimensioned (porm, porp, N+1) per Fortran requirements, even though only kdcoef slices are used
        int ucoeff_third_dim = n + 1;  // Fortran requires N+1, not just kdcoef
        if (porm_val > 0 && porp_val > 0 && ucoeff_third_dim > 0) UCOEFF_res_flat.resize( (size_t)porm_val * porp_val * ucoeff_third_dim); else UCOEFF_res_flat.clear();
    }

    void copy_inputs_for_test(std::vector<double>& a_test, std::vector<double>& b_test,
                              std::vector<double>& c_test, std::vector<double>& d_test, bool is_row_major_layout) {
        if (is_row_major_layout) {
            a_test.resize(A_in_cm.size());
            b_test.resize(B_in_cm.size());
            c_test.resize(C_in_cm.size());
            d_test.resize(D_in_cm.size());
            if (n > 0 && !A_in_cm.empty()) slicot_transpose_to_c(A_in_cm.data(), a_test.data(), n, n, sizeof(double)); else a_test.clear();
            if (n > 0 && m > 0 && !B_in_cm.empty()) slicot_transpose_to_c(B_in_cm.data(), b_test.data(), n, m, sizeof(double)); else b_test.clear();
            if (p > 0 && n > 0 && !C_in_cm.empty()) slicot_transpose_to_c(C_in_cm.data(), c_test.data(), p, n, sizeof(double)); else c_test.clear();
            if (p > 0 && m > 0 && !D_in_cm.empty()) slicot_transpose_to_c(D_in_cm.data(), d_test.data(), p, m, sizeof(double)); else d_test.clear();
        } else {
            a_test = A_in_cm;
            b_test = B_in_cm;
            c_test = C_in_cm;
            d_test = D_in_cm;
        }
    }


    void calculate_c_lds(bool is_row_major_layout) { 
        if (is_row_major_layout) {
            lda_c = (n > 0) ? n : 1; 
            ldb_c = (m > 0) ? m : 1; 
            ldc_c = (n > 0) ? n : 1; 
            ldd_c = (m > 0) ? m : 1; 
        } else { 
            lda_c = (n > 0) ? n : 1;  
            ldb_c = (n > 0) ? n : 1; 
            ldc_c = (p > 0) ? p : 1; 
            ldd_c = (p > 0) ? p : 1; 
        }
        
        int porm_val = (rowcol == 'R' ? p : m);
        int porp_val = (rowcol == 'R' ? m : p);
        lddcoe_c = std::max(1,porm_val);
        lduco1_c = std::max(1,porm_val);
        lduco2_c = std::max(1,porp_val);
    }
};


TEST_F(Tb04adTest, DocExampleColMajor) {
    A_res_cm_test_buffer = A_in_cm;
    B_res_cm_test_buffer = B_in_cm;
    C_res_cm_test_buffer = C_in_cm;
    std::vector<double> d_test_const = D_in_cm;

    calculate_c_lds(false);

    int porm_val_setup = (rowcol == 'R' ? p : m);
    int porp_val_setup = (rowcol == 'R' ? m : p);
    if (porm_val_setup > 0) INDEX_res.resize(porm_val_setup); else INDEX_res.clear();
    if (porm_val_setup > 0 && kmax_coeffs_exp > 0) DCOEFF_res_flat.resize((size_t)porm_val_setup * kmax_coeffs_exp); else DCOEFF_res_flat.clear();
    // UCOEFF must be dimensioned (porm, porp, N+1) per Fortran requirements
    int ucoeff_third_dim_setup = n + 1;
    if (porm_val_setup > 0 && porp_val_setup > 0 && ucoeff_third_dim_setup > 0) UCOEFF_res_flat.resize((size_t)porm_val_setup * porp_val_setup * ucoeff_third_dim_setup); else UCOEFF_res_flat.clear();


    int info = slicot_tb04ad(&rowcol, n, m, p,
                             A_res_cm_test_buffer.data(), lda_c,
                             B_res_cm_test_buffer.data(), ldb_c,
                             C_res_cm_test_buffer.data(), ldc_c,
                             d_test_const.data(), ldd_c,
                             &NR_res, INDEX_res.empty() ? nullptr : INDEX_res.data(), 
                             DCOEFF_res_flat.empty() ? nullptr : DCOEFF_res_flat.data(), lddcoe_c,
                             UCOEFF_res_flat.empty() ? nullptr : UCOEFF_res_flat.data(), lduco1_c, lduco2_c,
                             tol1, tol2, 0 /*row_major=false*/);

    ASSERT_EQ(info, 0);
    ASSERT_EQ(NR_res, NR_exp);

    for (size_t i = 0; i < (size_t)NR_exp * NR_exp; ++i) {
        EXPECT_NEAR(A_res_cm_test_buffer[i], A_exp_cm[i], check_tol) << "A_cm mismatch at index " << i;
    }
    for (size_t i = 0; i < (size_t)NR_exp * m; ++i) {
        EXPECT_NEAR(B_res_cm_test_buffer[i], B_exp_cm[i], check_tol) << "B_cm mismatch at index " << i;
    }
    for (size_t i = 0; i < (size_t)p * NR_exp; ++i) {
        EXPECT_NEAR(C_res_cm_test_buffer[i], C_exp_cm[i], check_tol) << "C_cm mismatch at index " << i;
    }

    int porm_val = (rowcol == 'R' ? p : m);
    if (porm_val > 0) { 
        ASSERT_EQ(INDEX_res.size(), INDEX_exp.size());
        for (int i = 0; i < porm_val; ++i) {
            EXPECT_EQ(INDEX_res[i], INDEX_exp[i]) << "INDEX mismatch at index " << i;
        }
    }
    
    if (!DCOEFF_exp_flat.empty()) { 
        ASSERT_EQ(DCOEFF_res_flat.size(), DCOEFF_exp_flat.size());
        for (size_t i = 0; i < DCOEFF_exp_flat.size(); ++i) {
            EXPECT_NEAR(DCOEFF_res_flat[i], DCOEFF_exp_flat[i], check_tol) << "DCOEFF mismatch at index " << i;
        }
    } else {
        EXPECT_TRUE(DCOEFF_res_flat.empty());
    }


    if (!UCOEFF_exp_flat.empty()) {
        // UCOEFF_res_flat is sized (porm, porp, N+1) but only (porm, porp, kdcoef) elements are meaningful
        ASSERT_GE(UCOEFF_res_flat.size(), UCOEFF_exp_flat.size()) << "UCOEFF_res_flat too small";
        for (size_t i = 0; i < UCOEFF_exp_flat.size(); ++i) {
            EXPECT_NEAR(UCOEFF_res_flat[i], UCOEFF_exp_flat[i], check_tol) << "UCOEFF mismatch at index " << i;
        }
    } else {
        EXPECT_TRUE(UCOEFF_res_flat.empty());
    }
}

TEST_F(Tb04adTest, DocExampleRowMajor) {
    std::vector<double> a_test_rm, b_test_rm, c_test_rm, d_test_const_rm;
    copy_inputs_for_test(a_test_rm, b_test_rm, c_test_rm, d_test_const_rm, true);
    calculate_c_lds(true);

    int porm_val_setup = (rowcol == 'R' ? p : m);
    int porp_val_setup = (rowcol == 'R' ? m : p);
    if (porm_val_setup > 0) INDEX_res.resize(porm_val_setup); else INDEX_res.clear();
    if (porm_val_setup > 0 && kmax_coeffs_exp > 0) DCOEFF_res_flat.resize((size_t)porm_val_setup * kmax_coeffs_exp); else DCOEFF_res_flat.clear();
    // UCOEFF must be dimensioned (porm, porp, N+1) per Fortran requirements
    int ucoeff_third_dim_setup = n + 1;
    if (porm_val_setup > 0 && porp_val_setup > 0 && ucoeff_third_dim_setup > 0) UCOEFF_res_flat.resize((size_t)porm_val_setup * porp_val_setup * ucoeff_third_dim_setup); else UCOEFF_res_flat.clear();


    int info = slicot_tb04ad(&rowcol, n, m, p,
                             a_test_rm.data(), lda_c, 
                             b_test_rm.data(), ldb_c,
                             c_test_rm.data(), ldc_c,
                             d_test_const_rm.data(), ldd_c,
                             &NR_res, INDEX_res.empty() ? nullptr : INDEX_res.data(), 
                             DCOEFF_res_flat.empty() ? nullptr : DCOEFF_res_flat.data(), lddcoe_c,
                             UCOEFF_res_flat.empty() ? nullptr : UCOEFF_res_flat.data(), lduco1_c, lduco2_c,
                             tol1, tol2, 1 /*row_major=true*/);

    ASSERT_EQ(info, 0);
    ASSERT_EQ(NR_res, NR_exp);

    std::vector<double> A_exp_rm( (size_t)NR_exp * NR_exp);
    std::vector<double> B_exp_rm( (size_t)NR_exp * m);
    std::vector<double> C_exp_rm( (size_t)p * NR_exp);

    slicot_transpose_to_c(A_exp_cm.data(), A_exp_rm.data(), NR_exp, NR_exp, sizeof(double));
    slicot_transpose_to_c(B_exp_cm.data(), B_exp_rm.data(), NR_exp, m, sizeof(double));
    slicot_transpose_to_c(C_exp_cm.data(), C_exp_rm.data(), p, NR_exp, sizeof(double));

    for (size_t i = 0; i < (size_t)NR_exp * NR_exp; ++i) {
        EXPECT_NEAR(a_test_rm[i], A_exp_rm[i], check_tol) << "A_rm mismatch at index " << i;
    }
    for (size_t i = 0; i < (size_t)NR_exp * m; ++i) {
        EXPECT_NEAR(b_test_rm[i], B_exp_rm[i], check_tol) << "B_rm mismatch at index " << i;
    }
     for (size_t i = 0; i < (size_t)p * NR_exp; ++i) {
        EXPECT_NEAR(c_test_rm[i], C_exp_rm[i], check_tol) << "C_rm mismatch at index " << i;
    }

    int porm_val = (rowcol == 'R' ? p : m);
    if (porm_val > 0) {
        ASSERT_EQ(INDEX_res.size(), INDEX_exp.size());
        for (int i = 0; i < porm_val; ++i) {
            EXPECT_EQ(INDEX_res[i], INDEX_exp[i]) << "INDEX mismatch at index " << i;
        }
    }
    
    if (!DCOEFF_exp_flat.empty()) {
        ASSERT_EQ(DCOEFF_res_flat.size(), DCOEFF_exp_flat.size());
        for (size_t i = 0; i < DCOEFF_exp_flat.size(); ++i) {
            EXPECT_NEAR(DCOEFF_res_flat[i], DCOEFF_exp_flat[i], check_tol) << "DCOEFF mismatch at index " << i;
        }
    } else {
         EXPECT_TRUE(DCOEFF_res_flat.empty());
    }

    if (!UCOEFF_exp_flat.empty()) {
        // UCOEFF_res_flat is sized (porm, porp, N+1) but only (porm, porp, kdcoef) elements are meaningful
        ASSERT_GE(UCOEFF_res_flat.size(), UCOEFF_exp_flat.size()) << "UCOEFF_res_flat too small";
        for (size_t i = 0; i < UCOEFF_exp_flat.size(); ++i) {
            EXPECT_NEAR(UCOEFF_res_flat[i], UCOEFF_exp_flat[i], check_tol) << "UCOEFF mismatch at index " << i;
        }
    } else {
        EXPECT_TRUE(UCOEFF_res_flat.empty());
    }
}

TEST_F(Tb04adTest, ParameterValidation) {
    std::vector<double> dummy_A(1), dummy_B(1), dummy_C(1), dummy_D(1);
    int dummy_nr;
    std::vector<int> dummy_idx(std::max(1,std::max(this->m,this->p))); 
    
    int dummy_porm = std::max(1,std::max(this->m,this->p));
    int dummy_porp = std::max(1,std::max(this->m,this->p)); 
    int dummy_kmax = 1; 

    std::vector<double> dummy_dcoeff(dummy_porm * dummy_kmax); 
    std::vector<double> dummy_ucoeff(dummy_porm * dummy_porp * dummy_kmax);
    
    int current_n = this->n; 
    int current_m = this->m;
    int current_p = this->p;
    
    int val_lda_c = (current_n > 0) ? current_n : 1;
    int val_ldb_c = (current_n > 0) ? current_n : 1; 
    int val_ldc_c = (current_p > 0) ? current_p : 1; 
    int val_ldd_c = (current_p > 0) ? current_p : 1; 
    int val_lddcoe_c = dummy_porm;
    int val_lduco1_c = dummy_porm;
    int val_lduco2_c = dummy_porp;


    EXPECT_EQ(slicot_tb04ad("X", current_n, current_m, current_p, dummy_A.data(), val_lda_c, dummy_B.data(), val_ldb_c, dummy_C.data(), val_ldc_c, dummy_D.data(), val_ldd_c, &dummy_nr, dummy_idx.data(), dummy_dcoeff.data(), val_lddcoe_c, dummy_ucoeff.data(), val_lduco1_c, val_lduco2_c, tol1, tol2, 0), -1); 
    EXPECT_EQ(slicot_tb04ad(&rowcol, -1, current_m, current_p, dummy_A.data(), val_lda_c, dummy_B.data(), val_ldb_c, dummy_C.data(), val_ldc_c, dummy_D.data(), val_ldd_c, &dummy_nr, dummy_idx.data(), dummy_dcoeff.data(), val_lddcoe_c, dummy_ucoeff.data(), val_lduco1_c, val_lduco2_c, tol1, tol2, 0), -2); 
    
    if (current_n > 0) {
      EXPECT_EQ(slicot_tb04ad(&rowcol, current_n, current_m, current_p, nullptr, val_lda_c, dummy_B.data(), val_ldb_c, dummy_C.data(), val_ldc_c, dummy_D.data(), val_ldd_c, &dummy_nr, dummy_idx.data(), dummy_dcoeff.data(), lddcoe_c, dummy_ucoeff.data(), lduco1_c, lduco2_c, tol1, tol2, 0), -5); 
      EXPECT_EQ(slicot_tb04ad(&rowcol, current_n, current_m, current_p, dummy_A.data(), 0, dummy_B.data(), val_ldb_c, dummy_C.data(), val_ldc_c, dummy_D.data(), val_ldd_c, &dummy_nr, dummy_idx.data(), dummy_dcoeff.data(), lddcoe_c, dummy_ucoeff.data(), lduco1_c, lduco2_c, tol1, tol2, 0), -6); 
    }
}

TEST_F(Tb04adTest, ZeroDimensions) {
    int n0=0, m0=1, p0=1; 
    char rc_z = 'R'; // porm=p0=1, porp=m0=1
    std::vector<double> d_mat_z = {2.0}; 
    int nr0;
    std::vector<int> idx0(std::max(1,p0)); 
    std::vector<double> dcoeff0(std::max(1,p0) * std::max(1,(n0+1))); 
    std::vector<double> ucoeff0(std::max(1,p0) * std::max(1,m0) * std::max(1,(n0+1)));

    int lda0_c = 1; 
    int ldb0_c = 1; 
    int ldc0_c = 1; 
    int ldd0_c = std::max(1,p0); 
    int lddcoe0_c = std::max(1,p0); 
    int lduco10_c = std::max(1,p0); 
    int lduco20_c = std::max(1,m0); 


    int info = slicot_tb04ad(&rc_z, n0, m0, p0,
                             nullptr, lda0_c, nullptr, ldb0_c, nullptr, ldc0_c, d_mat_z.data(), ldd0_c,
                             &nr0, idx0.data(), dcoeff0.data(), lddcoe0_c, ucoeff0.data(), lduco10_c, lduco20_c,
                             0.0, 0.0, 0);
    ASSERT_EQ(info, 0);
    EXPECT_EQ(nr0, 0);
    if (p0 > 0) { 
       ASSERT_EQ(idx0.size(), (size_t)p0); 
       EXPECT_EQ(idx0[0], 0); 
    }
    
    if (p0 > 0 ) EXPECT_NEAR(dcoeff0[0], 1.0, 1e-9);
    if (p0 > 0 && m0 > 0 ) EXPECT_NEAR(ucoeff0[0], d_mat_z[0], 1e-9);
}

TEST_F(Tb04adTest, AllZeroDimensions) {
    int n_az=0, m_az=0, p_az=0;
    char rc_az = 'R'; // porm=0, porp=0
    int nr_az = -999; // Initialize to detect if Fortran modifies it
    
    std::vector<int> dummy_idx_az(1); 
    std::vector<double> dummy_dcoeff_az(1);
    std::vector<double> dummy_ucoeff_az(1);

    int lda_az_c = 1;
    int ldb_az_c = 1;
    int ldc_az_c = 1;
    int ldd_az_c = 1;
    int lddcoe_az_c = 1; 
    int lduco1_az_c = 1; 
    int lduco2_az_c = 1; 


    int info = slicot_tb04ad(&rc_az, n_az, m_az, p_az,
                             nullptr, lda_az_c, nullptr, ldb_az_c, nullptr, ldc_az_c, nullptr, ldd_az_c,
                             &nr_az, 
                             (p_az > 0 ? dummy_idx_az.data() : nullptr), 
                             (p_az > 0 ? dummy_dcoeff_az.data() : nullptr), lddcoe_az_c, 
                             ((p_az > 0 && m_az > 0) ? dummy_ucoeff_az.data() : nullptr), lduco1_az_c, lduco2_az_c,
                             0.0, 0.0, 0);
    ASSERT_EQ(info, 0);
    EXPECT_EQ(nr_az, 0);
}
