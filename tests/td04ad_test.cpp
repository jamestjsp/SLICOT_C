#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <algorithm> // For std::max, std::accumulate
#include <numeric>   // For std::accumulate
#include <cctype>    // For toupper C++ style

#include "td04ad.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions
#include "test_utils.h"   // For load_test_data_from_csv (not used here, data embedded)

// --- Base Test Fixture for TD04AD ---
class Td04adBaseTest : public ::testing::Test {
protected:
    // Parameters
    char ROWCOL_val;
    int M_val;
    int P_val;
    std::vector<int> INDEX_val;
    std::vector<double> DCOEFF_c_layout; 
    std::vector<double> UCOEFF_c_layout; 
    double TOL_val;

    // Expected outputs (always in Fortran column-major for comparison convenience)
    int NR_expected;
    std::vector<double> A_expected_f;
    std::vector<double> B_expected_f;
    std::vector<double> C_expected_f;
    std::vector<double> D_expected_f;
    int INFO_expected;

    // Actual outputs from wrapper
    int NR_out;
    std::vector<double> A_out_c; 
    std::vector<double> B_out_c;
    std::vector<double> C_out_c;
    std::vector<double> D_out_c;
    int INFO_out;

    // Dimensions
    int N_sum;
    int KDCOEF;
    int PORM;

    // Leading dimensions for C arrays passed to the wrapper
    int LDDCOE_c, LDUCO1_c, LDUCO2_c, LDA_c, LDB_c, LDC_c, LDD_c;

    double check_tol = 1e-4;

    void InitializeDimensions() {
        N_sum = 0;
        int max_idx = 0;
        char rowcol_upper_char = static_cast<char>(std::toupper(static_cast<unsigned char>(ROWCOL_val)));
        PORM = (rowcol_upper_char == 'R') ? P_val : M_val;

        if (PORM > 0 && !INDEX_val.empty()) {
            N_sum = std::accumulate(INDEX_val.begin(), INDEX_val.end(), 0);
            for (int idx_val : INDEX_val) {
                if (idx_val > max_idx) max_idx = idx_val;
            }
        }
        KDCOEF = max_idx + 1;
        
        if (PORM == 0) { 
             KDCOEF = 1; 
        } else if (INDEX_val.empty() && PORM > 0) { 
             KDCOEF = 1; 
        }
    }

    virtual void SetUpTestData() = 0; // Pure virtual, to be implemented by derived test fixtures

    void SetUpBase(bool is_row_major) {
        // Call SetUpTestData from the derived class to populate M_val, P_val, INDEX_val, etc.
        SetUpTestData(); 
        // Then initialize dimensions based on these values
        InitializeDimensions(); 

        size_t n_sum_alloc = std::max(1, N_sum);
        size_t m_alloc = std::max(1, M_val);
        size_t p_alloc = std::max(1, P_val);
        size_t kdcoef_alloc = std::max(1, KDCOEF);
        size_t porm_alloc = std::max(1, PORM);
        
        int max_mp_val = std::max(M_val, P_val);
        size_t max_mp_alloc = std::max(1, max_mp_val);


        A_out_c.assign(n_sum_alloc * n_sum_alloc, 0.0);
        if (M_val > 0 && n_sum_alloc > 0 ) B_out_c.assign(n_sum_alloc * m_alloc, 0.0); else B_out_c.clear();
        if (P_val > 0 && n_sum_alloc > 0 ) C_out_c.assign(p_alloc * n_sum_alloc, 0.0); else C_out_c.clear();
        if (P_val > 0 && M_val > 0) D_out_c.assign(p_alloc * m_alloc, 0.0); else D_out_c.clear();


        if (is_row_major) {
            LDDCOE_c = kdcoef_alloc;
            LDUCO1_c = m_alloc; 
            LDUCO2_c = kdcoef_alloc;
            LDA_c = n_sum_alloc;
            LDB_c = m_alloc;
            LDC_c = n_sum_alloc;
            LDD_c = m_alloc;
        } else { 
            LDDCOE_c = porm_alloc;
            char rowcol_upper_char = static_cast<char>(std::toupper(static_cast<unsigned char>(ROWCOL_val)));
            if (rowcol_upper_char == 'R') {
                LDUCO1_c = p_alloc;
                LDUCO2_c = m_alloc;
            } else { 
                LDUCO1_c = max_mp_alloc; 
                LDUCO2_c = max_mp_alloc;
            }
            LDA_c = n_sum_alloc;
            LDB_c = n_sum_alloc;
            LDC_c = p_alloc;
            LDD_c = p_alloc;
        }
    }

    void VerifyOutputs(bool is_row_major_input) {
        ASSERT_EQ(INFO_out, INFO_expected);
        if (INFO_expected != 0) return; 

        ASSERT_EQ(NR_out, NR_expected);

        if (NR_out == 0) { 
            if (P_val > 0 && M_val > 0) {
                 ASSERT_LE((size_t)P_val * M_val, D_out_c.size()); 
                 ASSERT_EQ(D_expected_f.size(), (size_t)P_val * M_val);
                for (int i = 0; i < P_val; ++i) {
                    for (int j = 0; j < M_val; ++j) {
                        double actual_val;
                        if (is_row_major_input) {
                            actual_val = D_out_c[i * LDD_c + j];
                        } else {
                            actual_val = D_out_c[i + j * LDD_c];
                        }
                        EXPECT_NEAR(actual_val, D_expected_f[i + j * P_val], check_tol)
                            << "D matrix mismatch at (" << i << "," << j << ")";
                    }
                }
            } else {
                EXPECT_TRUE(D_expected_f.empty()); 
            }
            EXPECT_TRUE(A_expected_f.empty()); 
            EXPECT_TRUE(B_expected_f.empty());
            EXPECT_TRUE(C_expected_f.empty());
            return; 
        }
        
        std::vector<double> A_actual_f( (size_t)NR_out * NR_out );
        std::vector<double> B_actual_f( M_val > 0 ? (size_t)NR_out * M_val : 0);
        std::vector<double> C_actual_f( P_val > 0 ? (size_t)P_val * NR_out : 0);
        std::vector<double> D_actual_f( (P_val > 0 && M_val > 0) ? (size_t)P_val * M_val : 0);

        int lda_f_comp = NR_out;
        int ldb_f_comp = NR_out; 
        int ldc_f_comp = P_val > 0 ? P_val : 1; 
        int ldd_f_comp = P_val > 0 ? P_val : 1;


        if (is_row_major_input) {
            slicot_transpose_to_fortran_with_ld(A_out_c.data(), A_actual_f.data(), NR_out, NR_out, LDA_c, lda_f_comp, sizeof(double));
            if (M_val > 0 && !B_out_c.empty()) slicot_transpose_to_fortran_with_ld(B_out_c.data(), B_actual_f.data(), NR_out, M_val, LDB_c, ldb_f_comp, sizeof(double));
            if (P_val > 0 && !C_out_c.empty()) slicot_transpose_to_fortran_with_ld(C_out_c.data(), C_actual_f.data(), P_val, NR_out, LDC_c, ldc_f_comp, sizeof(double));
            if (P_val > 0 && M_val > 0 && !D_out_c.empty()) slicot_transpose_to_fortran_with_ld(D_out_c.data(), D_actual_f.data(), P_val, M_val, LDD_c, ldd_f_comp, sizeof(double));
        } else { 
            if (!A_out_c.empty()) for(size_t j=0; j < (size_t)NR_out; ++j) for(size_t i=0; i < (size_t)NR_out; ++i) A_actual_f[i+j*(size_t)NR_out] = A_out_c[i+j*(size_t)LDA_c];
            if (M_val > 0 && !B_out_c.empty()) for(size_t j=0; j < (size_t)M_val; ++j) for(size_t i=0; i < (size_t)NR_out; ++i) B_actual_f[i+j*(size_t)NR_out] = B_out_c[i+j*(size_t)LDB_c];
            if (P_val > 0 && !C_out_c.empty()) for(size_t j=0; j < (size_t)NR_out; ++j) for(size_t i=0; i < (size_t)P_val; ++i) C_actual_f[i+j*(size_t)P_val] = C_out_c[i+j*(size_t)LDC_c];
            if (P_val > 0 && M_val > 0 && !D_out_c.empty()) for(size_t j=0; j < (size_t)M_val; ++j) for(size_t i=0; i < (size_t)P_val; ++i) D_actual_f[i+j*(size_t)P_val] = D_out_c[i+j*(size_t)LDD_c];
        }

        ASSERT_EQ(A_actual_f.size(), A_expected_f.size());
        for (size_t i = 0; i < A_expected_f.size(); ++i) EXPECT_NEAR(A_actual_f[i], A_expected_f[i], check_tol) << "A mismatch at index " << i;
        
        if (M_val > 0) {
            ASSERT_EQ(B_actual_f.size(), B_expected_f.size());
            for (size_t i = 0; i < B_expected_f.size(); ++i) EXPECT_NEAR(B_actual_f[i], B_expected_f[i], check_tol) << "B mismatch at index " << i;
        } else { EXPECT_TRUE(B_expected_f.empty()); }

        if (P_val > 0) {
            ASSERT_EQ(C_actual_f.size(), C_expected_f.size());
            for (size_t i = 0; i < C_expected_f.size(); ++i) EXPECT_NEAR(C_actual_f[i], C_expected_f[i], check_tol) << "C mismatch at index " << i;
        } else { EXPECT_TRUE(C_expected_f.empty()); }

        if (P_val > 0 && M_val > 0) {
            ASSERT_EQ(D_actual_f.size(), D_expected_f.size());
            for (size_t i = 0; i < D_expected_f.size(); ++i) EXPECT_NEAR(D_actual_f[i], D_expected_f[i], check_tol) << "D mismatch at index " << i;
        } else { EXPECT_TRUE(D_expected_f.empty());}
    }
};

class Td04adDocExampleTest : public Td04adBaseTest {
protected:
    void SetUpTestData() override {
        ROWCOL_val = 'R'; M_val = 2; P_val = 2; TOL_val = 0.0;
        INDEX_val = {3, 3}; 
        NR_expected = 3;
        A_expected_f = { 0.5000,  4.4047, -5.5541, -0.8028, -2.3380,  1.6872, 0.9387,  2.5076, -4.1620 };
        B_expected_f = { -0.2000,  0.0000,  0.0000, -1.2500, -0.6097,  2.2217 };
        C_expected_f = { 0.0000, 0.0000, -0.8679, 0.0000, 0.2119, 0.9002 };
        D_expected_f = { 1.0000, 0.0000, 0.0000, 1.0000 };
        INFO_expected = 0;
    }
};

TEST_F(Td04adDocExampleTest, ColMajor) {
    SetUpBase(false); 
    DCOEFF_c_layout = { 1.0, 1.0, 6.0, 6.0, 11.0, 11.0, 6.0, 6.0 };
    UCOEFF_c_layout = { 
        1.0, 0.0,  0.0, 1.0,   
        6.0, 0.0,  1.0, 8.0,   
        12.0,1.0,  4.0, 20.0,  
        7.0, 1.0,  3.0, 15.0   
    };

    INFO_out = slicot_td04ad(ROWCOL_val, M_val, P_val, INDEX_val.data(),
                             DCOEFF_c_layout.data(), LDDCOE_c,
                             UCOEFF_c_layout.data(), LDUCO1_c, LDUCO2_c,
                             TOL_val,
                             &NR_out, A_out_c.data(), LDA_c, B_out_c.data(), LDB_c,
                             C_out_c.data(), LDC_c, D_out_c.data(), LDD_c,
                             0);
    VerifyOutputs(false);
}

TEST_F(Td04adDocExampleTest, RowMajor) {
    SetUpBase(true); 
    DCOEFF_c_layout = { 1.0, 6.0, 11.0, 6.0,  1.0, 6.0, 11.0, 6.0 };
    UCOEFF_c_layout = { 
        1.0, 6.0, 12.0, 7.0,  0.0, 1.0, 4.0, 3.0,   
        0.0, 0.0, 1.0, 1.0,   1.0, 8.0, 20.0, 15.0  
    };

    INFO_out = slicot_td04ad(ROWCOL_val, M_val, P_val, INDEX_val.data(),
                             DCOEFF_c_layout.data(), LDDCOE_c,
                             UCOEFF_c_layout.data(), LDUCO1_c, LDUCO2_c,
                             TOL_val,
                             &NR_out, A_out_c.data(), LDA_c, B_out_c.data(), LDB_c,
                             C_out_c.data(), LDC_c, D_out_c.data(), LDD_c,
                             1);
    VerifyOutputs(true);
}

class Td04adStaticGainTest : public Td04adBaseTest {
protected:
    void SetUpTestData() override {
        ROWCOL_val = 'R'; M_val = 1; P_val = 1; TOL_val = 0.0;
        INDEX_val = {0}; 
        NR_expected = 0; A_expected_f = {}; B_expected_f = {}; C_expected_f = {};
        D_expected_f = {2.0}; INFO_expected = 0;
    }
};

TEST_F(Td04adStaticGainTest, ColMajor) {
    SetUpBase(false);
    DCOEFF_c_layout = {1.0}; UCOEFF_c_layout = {2.0};
    INFO_out = slicot_td04ad(ROWCOL_val, M_val, P_val, INDEX_val.data(),
                             DCOEFF_c_layout.data(), LDDCOE_c, UCOEFF_c_layout.data(), LDUCO1_c, LDUCO2_c,
                             TOL_val, &NR_out, A_out_c.data(), LDA_c, B_out_c.data(), LDB_c,
                             C_out_c.data(), LDC_c, D_out_c.data(), LDD_c, 0);
    VerifyOutputs(false);
}

TEST_F(Td04adStaticGainTest, RowMajor) {
    SetUpBase(true);
    DCOEFF_c_layout = {1.0}; UCOEFF_c_layout = {2.0};
    INFO_out = slicot_td04ad(ROWCOL_val, M_val, P_val, INDEX_val.data(),
                             DCOEFF_c_layout.data(), LDDCOE_c, UCOEFF_c_layout.data(), LDUCO1_c, LDUCO2_c,
                             TOL_val, &NR_out, A_out_c.data(), LDA_c, B_out_c.data(), LDB_c,
                             C_out_c.data(), LDC_c, D_out_c.data(), LDD_c, 1);
    VerifyOutputs(true);
}

class Td04adParamValidationTest : public Td04adBaseTest { 
protected:
    void SetUpTestData() override { // Called by SetUpBase
        // Common defaults for validation tests, specific tests will override parameters
        TOL_val = 0.0;
        // INFO_expected is set by each test case directly
    }
};

TEST_F(Td04adParamValidationTest, InvalidRowCol) {
    ROWCOL_val = 'X'; M_val = 1; P_val = 1; INDEX_val = {0}; 
    DCOEFF_c_layout = {1.0}; UCOEFF_c_layout = {1.0};
    INFO_expected = -1; // Expect error for invalid ROWCOL
    SetUpBase(false); 
    INFO_out = slicot_td04ad(ROWCOL_val,M_val,P_val,INDEX_val.data(),DCOEFF_c_layout.data(),LDDCOE_c,UCOEFF_c_layout.data(),LDUCO1_c,LDUCO2_c,TOL_val,&NR_out,A_out_c.data(),LDA_c,B_out_c.data(),LDB_c,C_out_c.data(),LDC_c,D_out_c.data(),LDD_c,0); 
    ASSERT_EQ(INFO_out, INFO_expected);
}
TEST_F(Td04adParamValidationTest, InvalidM) {
    ROWCOL_val = 'R'; M_val = -1; P_val = 1; INDEX_val = {0}; 
    DCOEFF_c_layout = {1.0}; UCOEFF_c_layout = {1.0};
    INFO_expected = -2; // Expect error for invalid M
    SetUpBase(false);
    INFO_out = slicot_td04ad(ROWCOL_val,M_val,P_val,INDEX_val.data(),DCOEFF_c_layout.data(),LDDCOE_c,UCOEFF_c_layout.data(),LDUCO1_c,LDUCO2_c,TOL_val,&NR_out,A_out_c.data(),LDA_c,B_out_c.data(),LDB_c,C_out_c.data(),LDC_c,D_out_c.data(),LDD_c,0); 
    ASSERT_EQ(INFO_out, INFO_expected);
}
TEST_F(Td04adParamValidationTest, InvalidP) {
    ROWCOL_val = 'R'; M_val = 1; P_val = -1; INDEX_val = {0}; 
    DCOEFF_c_layout = {1.0}; UCOEFF_c_layout = {1.0};
    INFO_expected = -3; // Expect error for invalid P
    SetUpBase(false);
    INFO_out = slicot_td04ad(ROWCOL_val,M_val,P_val,INDEX_val.data(),DCOEFF_c_layout.data(),LDDCOE_c,UCOEFF_c_layout.data(),LDUCO1_c,LDUCO2_c,TOL_val,&NR_out,A_out_c.data(),LDA_c,B_out_c.data(),LDB_c,C_out_c.data(),LDC_c,D_out_c.data(),LDD_c,0); 
    ASSERT_EQ(INFO_out, INFO_expected);
}
TEST_F(Td04adParamValidationTest, NullIndexWhenNeeded) {
    ROWCOL_val = 'R'; P_val = 1; M_val = 0; 
    INDEX_val.clear(); 
    DCOEFF_c_layout = {1.0}; UCOEFF_c_layout.clear(); 
    INFO_expected = -4; // Expect error for NULL index when PORM > 0
    SetUpBase(false); 
    INFO_out = slicot_td04ad(ROWCOL_val,M_val,P_val,nullptr,DCOEFF_c_layout.data(),LDDCOE_c,UCOEFF_c_layout.empty()? nullptr : UCOEFF_c_layout.data(),LDUCO1_c,LDUCO2_c,TOL_val,&NR_out,A_out_c.data(),LDA_c,B_out_c.data(),LDB_c,C_out_c.data(),LDC_c,D_out_c.data(),LDD_c,0); 
    ASSERT_EQ(INFO_out, INFO_expected); 
}
TEST_F(Td04adParamValidationTest, InvalidIndexContent) {
    ROWCOL_val = 'R'; P_val = 1; M_val = 1; 
    INDEX_val = {-1}; 
    DCOEFF_c_layout = {1.0}; 
    UCOEFF_c_layout = {1.0}; 
    INFO_expected = -4; // Expect error for negative value in INDEX
    SetUpBase(false); 
    
    INFO_out = slicot_td04ad(ROWCOL_val,M_val,P_val,INDEX_val.data(),DCOEFF_c_layout.data(),LDDCOE_c,UCOEFF_c_layout.data(),LDUCO1_c,LDUCO2_c,TOL_val,&NR_out,A_out_c.data(),LDA_c,B_out_c.data(),LDB_c,C_out_c.data(),LDC_c,D_out_c.data(),LDD_c,0); 
    ASSERT_EQ(INFO_out, INFO_expected); 
}

class Td04adZeroDimTest : public Td04adBaseTest {
protected:
    void SetUpTestData() override {
        ROWCOL_val = 'R'; M_val = 0; P_val = 0; TOL_val = 0.0;
        INDEX_val.clear(); NR_expected = 0; 
        A_expected_f.clear(); B_expected_f.clear(); C_expected_f.clear(); D_expected_f.clear();
        INFO_expected = 0;
    }
};

TEST_F(Td04adZeroDimTest, ColMajor) {
    SetUpBase(false);
    DCOEFF_c_layout.clear(); UCOEFF_c_layout.clear(); 
    INFO_out = slicot_td04ad(ROWCOL_val, M_val, P_val, 
                             nullptr, nullptr, LDDCOE_c, nullptr, LDUCO1_c, LDUCO2_c,
                             TOL_val, &NR_out, A_out_c.data(), LDA_c, B_out_c.data(), LDB_c,
                             C_out_c.data(), LDC_c, D_out_c.data(), LDD_c, 0);
    VerifyOutputs(false);
}

TEST_F(Td04adZeroDimTest, RowMajor) {
    SetUpBase(true);
    DCOEFF_c_layout.clear(); UCOEFF_c_layout.clear();
     INFO_out = slicot_td04ad(ROWCOL_val, M_val, P_val, 
                             nullptr, nullptr, LDDCOE_c, nullptr, LDUCO1_c, LDUCO2_c,
                             TOL_val, &NR_out, A_out_c.data(), LDA_c, B_out_c.data(), LDB_c,
                             C_out_c.data(), LDC_c, D_out_c.data(), LDD_c, 1);
    VerifyOutputs(true);
}

