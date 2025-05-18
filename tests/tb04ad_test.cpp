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
    // Parameters from TB04AD.html example
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

    // Expected output matrices (column-major)
    std::vector<double> A_exp_cm;
    std::vector<double> B_exp_cm;
    std::vector<double> C_exp_cm;
    int NR_exp = 3;
    std::vector<int> INDEX_exp;
    std::vector<double> DCOEFF_exp_flat; // Flattened PORM x KMAX_D
    std::vector<double> UCOEFF_exp_flat; // Flattened PORM x PORP x KMAX_U

    int kmax_d_exp = 0; // max(INDEX_exp) + 1 for dcoeff
    int kmax_u_exp = 0; // max(INDEX_exp) + 1 for ucoeff (same, usually)

    // Output variables for the wrapper
    int NR_res = 0;
    std::vector<int> INDEX_res;
    std::vector<double> A_res_cm;
    std::vector<double> B_res_cm;
    std::vector<double> C_res_cm;
    std::vector<double> DCOEFF_res_flat;
    std::vector<double> UCOEFF_res_flat;


    // Leading dimensions (for C arrays)
    int lda_c, ldb_c, ldc_c, ldd_c;
    int lddcoe_c, lduco1_c, lduco2_c;


    double check_tol = 1e-4;

    void SetUp() override {
        // Initialize input data (from TB04AD.html example)
        // A_in = [[-1.0, 0.0, 0.0], [0.0, -2.0, 0.0], [0.0, 0.0, -3.0]]
        A_in_cm = {-1.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, -3.0};

        // B_in = [[0.0, 1.0], [1.0, 1.0], [-1.0, 0.0]] (N x M)
        B_in_cm = {0.0, 1.0, -1.0, 1.0, 1.0, 0.0};

        // C_in = [[0.0, 1.0, 1.0], [1.0, 0.0, 1.0]] (P x N)
        C_in_cm = {0.0, 1.0, 1.0, 0.0, 1.0, 1.0}; // Transposed for CM: {0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0}
                                            // Data in file: ( ( C(I,J), J = 1,N ), I = 1,P ) -> C(1,1..N), C(2,1..N)
                                            // This means C_in_cm should be {C(1,1),C(2,1), C(1,2),C(2,2), C(1,3),C(2,3)}
        C_in_cm = {0.0, 1.0,  1.0, 0.0,  1.0, 1.0};


        // D_in = [[1.0, 0.0], [0.0, 1.0]] (P x M)
        D_in_cm = {1.0, 0.0, 0.0, 1.0};

        // Expected outputs
        NR_exp = 3;
        // A_exp = [[-2.5000, -0.2887, -0.4082], [-0.2887, -1.5000, -0.7071], [-0.4082, -0.7071, -2.0000]]
        A_exp_cm = {-2.5000, -0.2887, -0.4082, -0.2887, -1.5000, -0.7071, -0.4082, -0.7071, -2.0000};
        // B_exp = [[-1.4142, -0.7071], [0.0000, 1.2247], [0.0000, 0.0000]] (NR x M)
        B_exp_cm = {-1.4142, 0.0000, 0.0000, -0.7071, 1.2247, 0.0000};
        // C_exp = [[0.0000, 0.8165, 1.1547], [0.0000, 1.6330, 0.5774]] (P x NR)
        C_exp_cm = {0.0000, 0.0000, 0.8165, 1.6330, 1.1547, 0.5774};

        INDEX_exp = {2, 3}; // For P=2 (rowcol='R')
        
        kmax_d_exp = 0;
        for(int deg : INDEX_exp) kmax_d_exp = std::max(kmax_d_exp, deg);
        kmax_d_exp +=1; // Number of coefficients
        kmax_u_exp = kmax_d_exp;


        // DCOEFF_exp for P=2, KMAX_D=4 (max(2,3)+1)
        // Row 1 (index 2 => 3 coeffs): 1.00, 5.00, 6.00 (pad with 0.00 for kmax_d_exp=4)
        // Row 2 (index 3 => 4 coeffs): 1.00, 6.00, 11.00, 6.00
        DCOEFF_exp_flat = {
            1.00, 5.00, 6.00, 0.00,
            1.00, 6.00, 11.00, 6.00
        };

        // UCOEFF_exp for P=2, M=2, KMAX_U=4
        // U(0,0): 1.00, 5.00, 7.00 (pad to 4)
        // U(0,1): 0.00, 1.00, 3.00 (pad to 4)
        // U(1,0): 0.00, 0.00, 1.00, 1.00
        // U(1,1): 1.00, 8.00, 20.00, 15.00
        // Stored as PORM slices, then PORP slices, then KMAX_U coeffs
        // PORM=P=2, PORP=M=2
        UCOEFF_exp_flat = {
            // P=0
            1.00, 5.00, 7.00, 0.00,  // M=0
            0.00, 1.00, 3.00, 0.00,  // M=1
            // P=1
            0.00, 0.00, 1.00, 1.00, // M=0
            1.00, 8.00, 20.00, 15.00 // M=1
        };

        // Resize result vectors (assuming rowcol='R' so PORM=p, PORP=m)
        int porm = (rowcol == 'R' ? p : m);
        int porp = (rowcol == 'R' ? m : p);

        A_res_cm.resize(n * n); // Max possible size, will check NR_res later
        B_res_cm.resize(n * m); // Max possible size
        C_res_cm.resize(p * n); // Max possible size
        if (porm > 0) INDEX_res.resize(porm);
        if (porm > 0 && kmax_d_exp > 0) DCOEFF_res_flat.resize(porm * kmax_d_exp);
        if (porm > 0 && porp > 0 && kmax_u_exp > 0) UCOEFF_res_flat.resize(porm * porp * kmax_u_exp);
    }

    void copy_inputs_for_test(std::vector<double>& a_test, std::vector<double>& b_test,
                              std::vector<double>& c_test, std::vector<double>& d_test, bool is_row_major) {
        if (is_row_major) {
            a_test.resize(n * n);
            b_test.resize(n * m);
            c_test.resize(p * n);
            d_test.resize(p * m);
            if (n > 0) slicot_transpose_to_c(A_in_cm.data(), a_test.data(), n, n, sizeof(double));
            if (n > 0 && m > 0) slicot_transpose_to_c(B_in_cm.data(), b_test.data(), n, m, sizeof(double));
            if (p > 0 && n > 0) slicot_transpose_to_c(C_in_cm.data(), c_test.data(), p, n, sizeof(double));
            if (p > 0 && m > 0) slicot_transpose_to_c(D_in_cm.data(), d_test.data(), p, m, sizeof(double));
        } else {
            a_test = A_in_cm;
            b_test = B_in_cm;
            c_test = C_in_cm;
            d_test = D_in_cm;
        }
    }


    void calculate_c_lds(bool is_row_major) {
        if (is_row_major) {
            lda_c = std::max(1, n); // cols of A
            ldb_c = std::max(1, m); // cols of B
            ldc_c = std::max(1, n); // cols of C
            ldd_c = std::max(1, m); // cols of D
        } else { // Column Major C
            lda_c = std::max(1, n); // rows of A
            ldb_c = std::max(1, n); // rows of B
            ldc_c = std::max(1, p); // rows of C
            ldd_c = std::max(1, p); // rows of D
        }
        // Output array leading dims (for C, these are the true first dimensions)
        int porm = (rowcol == 'R' ? p : m);
        int porp = (rowcol == 'R' ? m : p);
        lddcoe_c = std::max(1,porm);
        lduco1_c = std::max(1,porm);
        lduco2_c = std::max(1,porp);
    }
};


TEST_F(Tb04adTest, DocExampleColMajor) {
    std::vector<double> a_test, b_test, c_test, d_test_const;
    copy_inputs_for_test(a_test, b_test, c_test, d_test_const, false);
    calculate_c_lds(false);

    int info = slicot_tb04ad(&rowcol, n, m, p,
                             a_test.data(), lda_c,
                             b_test.data(), ldb_c,
                             c_test.data(), ldc_c,
                             d_test_const.data(), ldd_c,
                             &NR_res, INDEX_res.data(), DCOEFF_res_flat.data(), lddcoe_c,
                             UCOEFF_res_flat.data(), lduco1_c, lduco2_c,
                             tol1, tol2, 0 /*row_major=false*/);

    ASSERT_EQ(info, 0);
    ASSERT_EQ(NR_res, NR_exp);

    // Compare A_res (NR_res x NR_res) with A_exp_cm
    for (size_t i = 0; i < (size_t)NR_exp * NR_exp; ++i) {
        EXPECT_NEAR(a_test[i], A_exp_cm[i], check_tol) << "A_cm mismatch at index " << i;
    }
    // Compare B_res (NR_res x M) with B_exp_cm
    for (size_t i = 0; i < (size_t)NR_exp * m; ++i) {
        EXPECT_NEAR(b_test[i], B_exp_cm[i], check_tol) << "B_cm mismatch at index " << i;
    }
    // Compare C_res (P x NR_res) with C_exp_cm
    for (size_t i = 0; i < (size_t)p * NR_exp; ++i) {
        EXPECT_NEAR(c_test[i], C_exp_cm[i], check_tol) << "C_cm mismatch at index " << i;
    }

    int porm = (rowcol == 'R' ? p : m);
    for (int i = 0; i < porm; ++i) {
        EXPECT_EQ(INDEX_res[i], INDEX_exp[i]) << "INDEX mismatch at index " << i;
    }
    
    ASSERT_EQ(DCOEFF_res_flat.size(), DCOEFF_exp_flat.size());
    for (size_t i = 0; i < DCOEFF_exp_flat.size(); ++i) {
        EXPECT_NEAR(DCOEFF_res_flat[i], DCOEFF_exp_flat[i], check_tol) << "DCOEFF mismatch at index " << i;
    }

    ASSERT_EQ(UCOEFF_res_flat.size(), UCOEFF_exp_flat.size());
    for (size_t i = 0; i < UCOEFF_exp_flat.size(); ++i) {
        EXPECT_NEAR(UCOEFF_res_flat[i], UCOEFF_exp_flat[i], check_tol) << "UCOEFF mismatch at index " << i;
    }
}

TEST_F(Tb04adTest, DocExampleRowMajor) {
    std::vector<double> a_test_rm, b_test_rm, c_test_rm, d_test_const_rm;
    copy_inputs_for_test(a_test_rm, b_test_rm, c_test_rm, d_test_const_rm, true);
    calculate_c_lds(true);


    // Prepare output buffers for row major test - they will be converted back by the wrapper
    std::vector<double> A_res_rm(n * n);
    std::vector<double> B_res_rm(n * m);
    std::vector<double> C_res_rm(p * n);
    // Copy initial row-major data into these buffers as they are in/out
    A_res_rm = a_test_rm;
    B_res_rm = b_test_rm;
    C_res_rm = c_test_rm;


    int info = slicot_tb04ad(&rowcol, n, m, p,
                             A_res_rm.data(), lda_c, // Pass C LDA (cols for RM)
                             B_res_rm.data(), ldb_c,
                             C_res_rm.data(), ldc_c,
                             d_test_const_rm.data(), ldd_c,
                             &NR_res, INDEX_res.data(), DCOEFF_res_flat.data(), lddcoe_c,
                             UCOEFF_res_flat.data(), lduco1_c, lduco2_c,
                             tol1, tol2, 1 /*row_major=true*/);

    ASSERT_EQ(info, 0);
    ASSERT_EQ(NR_res, NR_exp);

    // Convert expected CM results to RM for comparison
    std::vector<double> A_exp_rm(NR_exp * NR_exp);
    std::vector<double> B_exp_rm(NR_exp * m);
    std::vector<double> C_exp_rm(p * NR_exp);

    slicot_transpose_to_c(A_exp_cm.data(), A_exp_rm.data(), NR_exp, NR_exp, sizeof(double));
    slicot_transpose_to_c(B_exp_cm.data(), B_exp_rm.data(), NR_exp, m, sizeof(double));
    slicot_transpose_to_c(C_exp_cm.data(), C_exp_rm.data(), p, NR_exp, sizeof(double));


    for (size_t i = 0; i < (size_t)NR_exp * NR_exp; ++i) {
        EXPECT_NEAR(A_res_rm[i], A_exp_rm[i], check_tol) << "A_rm mismatch at index " << i;
    }
    for (size_t i = 0; i < (size_t)NR_exp * m; ++i) {
        EXPECT_NEAR(B_res_rm[i], B_exp_rm[i], check_tol) << "B_rm mismatch at index " << i;
    }
     for (size_t i = 0; i < (size_t)p * NR_exp; ++i) {
        EXPECT_NEAR(C_res_rm[i], C_exp_rm[i], check_tol) << "C_rm mismatch at index " << i;
    }

    int porm = (rowcol == 'R' ? p : m);
    for (int i = 0; i < porm; ++i) {
        EXPECT_EQ(INDEX_res[i], INDEX_exp[i]) << "INDEX mismatch at index " << i;
    }
    
    ASSERT_EQ(DCOEFF_res_flat.size(), DCOEFF_exp_flat.size());
    for (size_t i = 0; i < DCOEFF_exp_flat.size(); ++i) {
        EXPECT_NEAR(DCOEFF_res_flat[i], DCOEFF_exp_flat[i], check_tol) << "DCOEFF mismatch at index " << i;
    }

    ASSERT_EQ(UCOEFF_res_flat.size(), UCOEFF_exp_flat.size());
    for (size_t i = 0; i < UCOEFF_exp_flat.size(); ++i) {
        EXPECT_NEAR(UCOEFF_res_flat[i], UCOEFF_exp_flat[i], check_tol) << "UCOEFF mismatch at index " << i;
    }
}

TEST_F(Tb04adTest, ParameterValidation) {
    std::vector<double> dummy_A(1), dummy_B(1), dummy_C(1), dummy_D(1);
    int dummy_nr;
    std::vector<int> dummy_idx(std::max(1,std::max(m,p)));
    std::vector<double> dummy_dcoeff(std::max(1,std::max(m,p)) * (n+1)); // Max possible sizes
    std::vector<double> dummy_ucoeff(std::max(1,std::max(m,p)) * std::max(1,std::max(m,p)) * (n+1));
    calculate_c_lds(false); // Use CM LDs for simplicity in validation calls

    EXPECT_EQ(slicot_tb04ad("X", n, m, p, dummy_A.data(), lda_c, dummy_B.data(), ldb_c, dummy_C.data(), ldc_c, dummy_D.data(), ldd_c, &dummy_nr, dummy_idx.data(), dummy_dcoeff.data(), lddcoe_c, dummy_ucoeff.data(), lduco1_c, lduco2_c, tol1, tol2, 0), -1); // Invalid rowcol
    EXPECT_EQ(slicot_tb04ad(&rowcol, -1, m, p, dummy_A.data(), lda_c, dummy_B.data(), ldb_c, dummy_C.data(), ldc_c, dummy_D.data(), ldd_c, &dummy_nr, dummy_idx.data(), dummy_dcoeff.data(), lddcoe_c, dummy_ucoeff.data(), lduco1_c, lduco2_c, tol1, tol2, 0), -2); // Invalid n
    EXPECT_EQ(slicot_tb04ad(&rowcol, n, -1, p, dummy_A.data(), lda_c, dummy_B.data(), ldb_c, dummy_C.data(), ldc_c, dummy_D.data(), ldd_c, &dummy_nr, dummy_idx.data(), dummy_dcoeff.data(), lddcoe_c, dummy_ucoeff.data(), lduco1_c, lduco2_c, tol1, tol2, 0), -3); // Invalid m
    EXPECT_EQ(slicot_tb04ad(&rowcol, n, m, -1, dummy_A.data(), lda_c, dummy_B.data(), ldb_c, dummy_C.data(), ldc_c, dummy_D.data(), ldd_c, &dummy_nr, dummy_idx.data(), dummy_dcoeff.data(), lddcoe_c, dummy_ucoeff.data(), lduco1_c, lduco2_c, tol1, tol2, 0), -4); // Invalid p

    if (n > 0) {
      EXPECT_EQ(slicot_tb04ad(&rowcol, n, m, p, nullptr, lda_c, dummy_B.data(), ldb_c, dummy_C.data(), ldc_c, dummy_D.data(), ldd_c, &dummy_nr, dummy_idx.data(), dummy_dcoeff.data(), lddcoe_c, dummy_ucoeff.data(), lduco1_c, lduco2_c, tol1, tol2, 0), -5); // Null A
      EXPECT_EQ(slicot_tb04ad(&rowcol, n, m, p, dummy_A.data(), 0, dummy_B.data(), ldb_c, dummy_C.data(), ldc_c, dummy_D.data(), ldd_c, &dummy_nr, dummy_idx.data(), dummy_dcoeff.data(), lddcoe_c, dummy_ucoeff.data(), lduco1_c, lduco2_c, tol1, tol2, 0), -6); // Invalid lda
    }
     // More specific tests for other pointers and LDs would follow
}

TEST_F(Tb04adTest, ZeroDimensions) {
    int n0=0, m0=1, p0=1; // Example: static gain system
    char rc = 'R';
    std::vector<double> d_mat = {2.0}; // 1x1 D matrix
    int nr0;
    std::vector<int> idx0(p0);
    std::vector<double> dcoeff0(p0 * (n0+1));
    std::vector<double> ucoeff0(p0 * m0 * (n0+1));

    int lda0 = std::max(1,n0);
    int ldb0 = std::max(1,n0);
    int ldc0 = std::max(1,p0);
    int ldd0 = std::max(1,p0);
    int lddcoe0 = std::max(1,p0);
    int lduco10 = std::max(1,p0);
    int lduco20 = std::max(1,m0);


    int info = slicot_tb04ad(&rc, n0, m0, p0,
                             nullptr, lda0, nullptr, ldb0, nullptr, ldc0, d_mat.data(), ldd0,
                             &nr0, idx0.data(), dcoeff0.data(), lddcoe0, ucoeff0.data(), lduco10, lduco20,
                             0.0, 0.0, 0);
    ASSERT_EQ(info, 0);
    EXPECT_EQ(nr0, 0);
    ASSERT_EQ(idx0.size(), p0);
    if (p0 > 0) EXPECT_EQ(idx0[0], 0); // Degree should be 0 for static gain
    
    // For N=0, DCOEFF should be [1.0] (monic denominator of degree 0)
    // UCOEFF should be D itself (scaled by DCOEFF, so just D)
    if (p0 > 0 && (n0+1) > 0) EXPECT_NEAR(dcoeff0[0], 1.0, 1e-9);
    if (p0 > 0 && m0 > 0 && (n0+1) > 0) EXPECT_NEAR(ucoeff0[0], d_mat[0], 1e-9);

}

// Test case for N=0, M=0, P=0
TEST_F(Tb04adTest, AllZeroDimensions) {
    int n0=0, m0=0, p0=0;
    char rc = 'R';
    int nr0;
    // Fortran expects non-NULL pointers for index, dcoeff, ucoeff even if PORM/PORP is 0,
    // but the wrapper should handle passing valid pointers to dummy 1-element arrays if needed,
    // or the Fortran might be okay with NULL if corresponding dimensions are zero.
    // The current wrapper passes what C gives, so if PORM=0, it will pass NULL for index.
    // Let's assume SLICOT handles NULL for zero-sized outputs.
    // The C wrapper itself makes lddcoe_c, lduco1_c, lduco2_c >=1 so it allocates dummy outputs.
    
    std::vector<int> dummy_idx(1);
    std::vector<double> dummy_dcoeff(1);
    std::vector<double> dummy_ucoeff(1);


    int info = slicot_tb04ad(&rc, n0, m0, p0,
                             nullptr, 1, nullptr, 1, nullptr, 1, nullptr, 1,
                             &nr0, dummy_idx.data(), dummy_dcoeff.data(), 1, dummy_ucoeff.data(), 1, 1,
                             0.0, 0.0, 0);
    ASSERT_EQ(info, 0);
    EXPECT_EQ(nr0, 0);
}