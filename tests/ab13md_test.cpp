#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm> // For std::count, std::max, std::accumulate
#include <numeric>   // For std::accumulate
#include <iostream>  // For debugging
#include <iomanip>   // For std::setprecision


#include "ab13md.h"
#include "slicot_utils.h" // Assumed to provide slicot_complex_double, transpose functions, etc.

// --- Column-Major Test Fixture ---
class AB13MDTestColMajor : public ::testing::Test {
protected:
    // Test parameters from AB13MD.html example
    int N_default = 6;
    int M_default = 5;
    double check_tol = 1e-7; // Tightened tolerance

    // Variables for individual tests, initialized in SetUp or test body
    int N;
    int M;
    std::vector<int> nblock_in;
    std::vector<int> itype_in;

    std::vector<slicot_complex_double> Z_in;
    std::vector<double> X_io;
    double bound_out;
    std::vector<double> D_out;
    std::vector<double> G_out;

    double expected_bound;
    int expected_info_success = 0;

    int info_result = -999; // Default to indicate not run
    int LDZ_in = 0;

    // Setup for the main example data
    void SetupExampleData() {
        N = N_default;
        M = M_default;
        nblock_in = {1, 1, 2, 1, 1};
        itype_in = {1, 1, 2, 2, 2};

        if (N > 0) {
            D_out.resize(N);
            G_out.resize(N);
        } else {
            D_out.clear();
            G_out.clear();
        }
        
        int num_real_blocks = 0;
        if (M > 0 && !itype_in.empty()) {
            for(int type : itype_in) {
                if (type == 1) {
                    num_real_blocks++;
                }
            }
        }
        
        int x_dim = 0;
        if (M > 0) { 
           x_dim = M + num_real_blocks -1;
        }
        if (x_dim > 0) {
            X_io.resize(x_dim);
            std::fill(X_io.begin(), X_io.end(), 0.0); 
        } else {
             X_io.clear();
        }

        // Data from AB13MD.html example (Z is 6x6, column-major)
        using std::complex;
        std::vector<complex<double>> Z_complex_col_major_from_html = {
            // Column 1
            complex<double>(-1.0, 6.0), complex<double>(2.0, -3.0), complex<double>(3.0, 8.0),
            complex<double>(3.0, 8.0), complex<double>(-5.0, -9.0), complex<double>(-6.0, 2.0),
            // Column 2
            complex<double>(4.0, 2.0), complex<double>(-2.0, 5.0), complex<double>(-6.0, -7.0),
            complex<double>(-4.0, 11.0), complex<double>(8.0, -7.0), complex<double>(12.0, -1.0),
            // Column 3
            complex<double>(5.0, -4.0), complex<double>(-4.0, -8.0), complex<double>(1.0, -3.0),
            complex<double>(-6.0, 14.0), complex<double>(2.0, -5.0), complex<double>(4.0, 16.0),
            // Column 4 (repeated from col 1 in HTML example)
            complex<double>(-1.0, 6.0), complex<double>(2.0, -3.0), complex<double>(3.0, 8.0),
            complex<double>(3.0, 8.0), complex<double>(-5.0, -9.0), complex<double>(-6.0, 2.0),
            // Column 5 (repeated from col 2 in HTML example)
            complex<double>(4.0, 2.0), complex<double>(-2.0, 5.0), complex<double>(-6.0, -7.0),
            complex<double>(-4.0, 11.0), complex<double>(8.0, -7.0), complex<double>(12.0, -1.0),
            // Column 6 (repeated from col 3 in HTML example)
            complex<double>(5.0, -4.0), complex<double>(-4.0, -8.0), complex<double>(1.0, -3.0),
            complex<double>(-6.0, 14.0), complex<double>(2.0, -5.0), complex<double>(4.0, 16.0)
        };

        if (N > 0) {
            Z_in.resize((size_t)N * N);
            size_t elements_to_copy = std::min(Z_in.size(), Z_complex_col_major_from_html.size());
            for (size_t i = 0; i < elements_to_copy; ++i) {
                 Z_in[i] = reinterpret_cast<slicot_complex_double&>(Z_complex_col_major_from_html[i]);
            }
        } else {
            Z_in.clear();
        }
        
        expected_bound = 41.74753408; 
        LDZ_in = std::max(1, N); // For column major C, LD is number of rows
    }

    void SetUp() override {
        SetupExampleData();
    }
};

// --- Row-Major Test Fixture ---
class AB13MDTestRowMajor : public AB13MDTestColMajor {
protected:
    std::vector<slicot_complex_double> Z_rm_in;

    void SetUp() override {
        AB13MDTestColMajor::SetupExampleData(); 
        
        if (N > 0 && !Z_in.empty()) {
            Z_rm_in.resize((size_t)N * N);
            slicot_transpose_to_c_with_ld(Z_in.data(), Z_rm_in.data(), 
                                          N, N, 
                                          LDZ_in, // Source LD (rows for col-major Z_in)
                                          N,      // Dest. LD (cols for row-major Z_rm_in)
                                          sizeof(slicot_complex_double));
        } else {
            Z_rm_in.clear();
        }
        LDZ_in = (N > 0) ? N : 1; // For row major C, LD is number of columns
    }
};


TEST_F(AB13MDTestColMajor, BasicFunctionality) {
    info_result = slicot_ab13md(
        'N', N, Z_in.empty() ? nullptr : Z_in.data(), LDZ_in,
        M, nblock_in.data(), itype_in.data(),
        X_io.empty() ? nullptr : X_io.data(),
        &bound_out, D_out.empty() ? nullptr : D_out.data(), G_out.empty() ? nullptr : G_out.data(),
        0 // col_major
    );
    ASSERT_EQ(info_result, expected_info_success) << "Call with FACT='N' failed.";
    EXPECT_NEAR(bound_out, expected_bound, check_tol);

    double bound2;
    std::vector<double> D2(N > 0 ? N : 0);
    std::vector<double> G2(N > 0 ? N : 0);
    // X_io is input/output, so its content from the first call is used here.

    info_result = slicot_ab13md(
        'F', N, Z_in.empty() ? nullptr : Z_in.data(), LDZ_in,
        M, nblock_in.data(), itype_in.data(),
        X_io.empty() ? nullptr : X_io.data(), 
        &bound2, D2.empty() ? nullptr : D2.data(), G2.empty() ? nullptr : G2.data(),
        0 // col_major
    );
    ASSERT_EQ(info_result, expected_info_success) << "Call with FACT='F' failed.";
    // The bound with FACT='F' might be slightly different or refined.
    // For this test, we expect it to be very close to the FACT='N' result if the state X converges quickly.
    // Or, if the problem is simple enough, it might be identical.
    // Comparing to `expected_bound` (from FACT='N' documentation) is a strong check.
    EXPECT_NEAR(bound2, expected_bound, check_tol); 
}

TEST_F(AB13MDTestRowMajor, BasicFunctionality) {
    info_result = slicot_ab13md(
        'N', N, Z_rm_in.empty() ? nullptr : Z_rm_in.data(), LDZ_in, 
        M, nblock_in.data(), itype_in.data(),
        X_io.empty() ? nullptr : X_io.data(),
        &bound_out, D_out.empty() ? nullptr : D_out.data(), G_out.empty() ? nullptr : G_out.data(),
        1 // row_major
    );
    ASSERT_EQ(info_result, expected_info_success);
    EXPECT_NEAR(bound_out, expected_bound, check_tol);
}

TEST_F(AB13MDTestColMajor, DifferentBlockStructures) {
    // Test with NBLOCK[0] == N, so X is not used.
    int n_alt = 3; 
    int m_alt = 1;
    std::vector<int> nblock_alt = {n_alt}; 
    std::vector<int> itype_alt = {2};  
    
    int ldz_alt = std::max(1, n_alt);
    std::vector<slicot_complex_double> z_alt_data;
    if (n_alt > 0) {
        z_alt_data.assign((size_t)n_alt * n_alt, {1.0, 1.0}); // Dummy data
    }
    
    std::vector<double> d_alt(n_alt > 0 ? n_alt : 0); // Size 0 if n_alt is 0
    std::vector<double> g_alt(n_alt > 0 ? n_alt : 0);
    double bound_alt_val;

    // X should not be used as NBLOCK[0] == N. Pass nullptr.
    info_result = slicot_ab13md(
        'N', n_alt, z_alt_data.empty() ? nullptr : z_alt_data.data(), ldz_alt,
        m_alt, nblock_alt.data(), itype_alt.data(),
        nullptr, 
        &bound_alt_val, d_alt.empty() ? nullptr : d_alt.data(), g_alt.empty() ? nullptr : g_alt.data(),
        0 // col_major
    );
    EXPECT_EQ(info_result, 0); 
    if (n_alt > 0) {
      EXPECT_GT(bound_alt_val, -1.0); // Bound should be non-negative, check against -1 for some margin with 0
    } else {
      EXPECT_DOUBLE_EQ(bound_alt_val, 0.0);
    }
}

TEST_F(AB13MDTestColMajor, ParameterValidation) {
    int n_val = 2, m_val = 1;
    int ldz_val_ok = std::max(1, n_val);
    std::vector<slicot_complex_double> z_val( (size_t)ldz_val_ok * n_val );
    std::vector<int> nblock_val_ok = {n_val}; // Sum is N, positive
    std::vector<int> itype_val_ok = {2};      // Complex
    
    int mr_val = 0; for(int type : itype_val_ok) if(type==1) mr_val++;
    int x_dim_val = (m_val > 0) ? (m_val + mr_val -1) : 0;
    std::vector<double> x_val(x_dim_val > 0 ? x_dim_val : 0);

    double bound_val_dummy;
    std::vector<double> d_val_dummy(n_val > 0 ? n_val : 0); 
    std::vector<double> g_val_dummy(n_val > 0 ? n_val : 0);

    // FACT invalid
    info_result = slicot_ab13md('X', n_val, z_val.data(), ldz_val_ok, m_val, nblock_val_ok.data(), itype_val_ok.data(), x_val.empty() ? nullptr : x_val.data(), &bound_val_dummy, d_val_dummy.empty()?nullptr:d_val_dummy.data(), g_val_dummy.empty()?nullptr:g_val_dummy.data(), 0);
    EXPECT_EQ(info_result, -1);
    
    // N invalid
    info_result = slicot_ab13md('N', -1, z_val.data(), ldz_val_ok, m_val, nblock_val_ok.data(), itype_val_ok.data(), x_val.empty() ? nullptr : x_val.data(), &bound_val_dummy, d_val_dummy.empty()?nullptr:d_val_dummy.data(), g_val_dummy.empty()?nullptr:g_val_dummy.data(), 0);
    EXPECT_EQ(info_result, -2);

    // Z invalid (NULL when N > 0)
    if (n_val > 0) {
        info_result = slicot_ab13md('N', n_val, nullptr, ldz_val_ok, m_val, nblock_val_ok.data(), itype_val_ok.data(), x_val.empty() ? nullptr : x_val.data(), &bound_val_dummy, d_val_dummy.empty()?nullptr:d_val_dummy.data(), g_val_dummy.empty()?nullptr:g_val_dummy.data(), 0);
        EXPECT_EQ(info_result, -3);
    }
    
    // LDZ invalid (col-major)
    if (n_val > 0) {
        info_result = slicot_ab13md('N', n_val, z_val.data(), 0, m_val, nblock_val_ok.data(), itype_val_ok.data(), x_val.empty() ? nullptr : x_val.data(), &bound_val_dummy, d_val_dummy.empty()?nullptr:d_val_dummy.data(), g_val_dummy.empty()?nullptr:g_val_dummy.data(), 0);
        EXPECT_EQ(info_result, -4); 
    }
     // LDZ invalid (row-major)
    if (n_val > 0) { // ldz (cols) < n_val
        info_result = slicot_ab13md('N', n_val, z_val.data(), (n_val - 1 > 0 ? n_val -1 : 0) , m_val, nblock_val_ok.data(), itype_val_ok.data(), x_val.empty() ? nullptr : x_val.data(), &bound_val_dummy, d_val_dummy.empty()?nullptr:d_val_dummy.data(), g_val_dummy.empty()?nullptr:g_val_dummy.data(), 1);
        EXPECT_EQ(info_result, -4);
    }

    // M invalid (<1)
    info_result = slicot_ab13md('N', n_val, z_val.data(), ldz_val_ok, 0, nullptr, nullptr, nullptr, &bound_val_dummy, d_val_dummy.empty()?nullptr:d_val_dummy.data(), g_val_dummy.empty()?nullptr:g_val_dummy.data(), 0);
    EXPECT_EQ(info_result, -5); 
    
    // NBLOCK invalid (NULL when M > 0)
    if (m_val > 0) {
         info_result = slicot_ab13md('N', n_val, z_val.data(), ldz_val_ok, m_val, nullptr, itype_val_ok.data(), x_val.empty() ? nullptr : x_val.data(), &bound_val_dummy, d_val_dummy.empty()?nullptr:d_val_dummy.data(), g_val_dummy.empty()?nullptr:g_val_dummy.data(), 0);
         EXPECT_EQ(info_result, -6);
    }

    // ITYPE invalid (NULL when M > 0)
     if (m_val > 0) {
         info_result = slicot_ab13md('N', n_val, z_val.data(), ldz_val_ok, m_val, nblock_val_ok.data(), nullptr, x_val.empty() ? nullptr : x_val.data(), &bound_val_dummy, d_val_dummy.empty()?nullptr:d_val_dummy.data(), g_val_dummy.empty()?nullptr:g_val_dummy.data(), 0);
         EXPECT_EQ(info_result, -7);
    }
    
    // X invalid (NULL when FACT='F' and its dim > 0)
    if (x_dim_val > 0) { // Only if X is supposed to have elements
        info_result = slicot_ab13md('F', n_val, z_val.data(), ldz_val_ok, m_val, nblock_val_ok.data(), itype_val_ok.data(), nullptr, &bound_val_dummy, d_val_dummy.empty()?nullptr:d_val_dummy.data(), g_val_dummy.empty()?nullptr:g_val_dummy.data(), 0);
        EXPECT_EQ(info_result, -8);
    }
    
    // Fortran error: Block sizes must be positive (INFO = 1)
    std::vector<int> nblock_has_zero = {0}; // NBLOCK(1)=0
    // For this test, N must be sum of NBLOCK. If NBLOCK={0}, then N should be 0.
    info_result = slicot_ab13md('N', 0, nullptr, 1, m_val, nblock_has_zero.data(), itype_val_ok.data(), nullptr, &bound_val_dummy, nullptr, nullptr, 0);
    EXPECT_EQ(info_result, 1);

    // Fortran error: Sum of block sizes must be N (INFO = 2)
    std::vector<int> nblock_sum_neq_n = {n_val + 1}; // Sum = N+1, but N=N
    if (n_val > 0) { // Ensure N is positive for this specific test of sum
        info_result = slicot_ab13md('N', n_val, z_val.data(), ldz_val_ok, m_val, nblock_sum_neq_n.data(), itype_val_ok.data(), x_val.empty()? nullptr : x_val.data(), &bound_val_dummy, d_val_dummy.empty()?nullptr:d_val_dummy.data(), g_val_dummy.empty()?nullptr:g_val_dummy.data(), 0);
        EXPECT_EQ(info_result, 2);
    }

    // Fortran error: Real block size must be 1 (INFO = 3)
    std::vector<int> nblock_real_bad = {2}; // Real block of size 2
    std::vector<int> itype_real = {1};      // Real block type
    int n_for_real_bad = 2; // N must be sum of NBLOCK
    int m_for_real_bad = 1;
    std::vector<slicot_complex_double> z_for_real_bad( (size_t)n_for_real_bad * n_for_real_bad );
    std::vector<double> d_for_real_bad(n_for_real_bad > 0 ? n_for_real_bad : 0);
    std::vector<double> g_for_real_bad(n_for_real_bad > 0 ? n_for_real_bad : 0);
    int mr_for_real_bad = 1;
    int x_dim_for_real_bad = m_for_real_bad + mr_for_real_bad - 1; // 1+1-1 = 1
    std::vector<double> x_for_real_bad(x_dim_for_real_bad > 0 ? x_dim_for_real_bad : 0);

    info_result = slicot_ab13md('N', n_for_real_bad, z_for_real_bad.data(), std::max(1,n_for_real_bad), m_for_real_bad, nblock_real_bad.data(), itype_real.data(), x_for_real_bad.empty()?nullptr:x_for_real_bad.data(), &bound_val_dummy, d_for_real_bad.empty()?nullptr:d_for_real_bad.data(), g_for_real_bad.empty()?nullptr:g_for_real_bad.data(), 0);
    EXPECT_EQ(info_result, 3);
    
    // Fortran error: Block type must be 1 or 2 (INFO = 4)
    std::vector<int> itype_invalid_val = {3}; // Invalid type
    if (n_val > 0 && m_val > 0) { // Ensure N, M are valid for this
        info_result = slicot_ab13md('N', n_val, z_val.data(), ldz_val_ok, m_val, nblock_val_ok.data(), itype_invalid_val.data(), x_val.empty() ? nullptr : x_val.data(), &bound_val_dummy, d_val_dummy.empty()?nullptr:d_val_dummy.data(), g_val_dummy.empty()?nullptr:g_val_dummy.data(), 0);
        EXPECT_EQ(info_result, 4);
    }
}

TEST_F(AB13MDTestColMajor, ZeroDimensions) {
    int n_zero = 0;
    double bound_zero_out = 0.0; // Expected bound for N=0 is 0.0
    int m_zd = 1; // M must be >= 1
    std::vector<int> itype_zd = {2}; // Dummy type
    
    // Case 1: N=0, M=0 (Wrapper error: M must be >= 1)
    info_result = slicot_ab13md('N', n_zero, nullptr, 1, 0, nullptr, nullptr, nullptr, &bound_zero_out, nullptr, nullptr, 0 );
    EXPECT_EQ(info_result, -5);

    // Case 2: N=0, M=1, NBLOCK={0} (Fortran INFO = 1: block size must be positive)
    std::vector<int> nblock_zd_zero_val = {0}; 
    std::vector<double> x_zd_case2; // x_dim = 1 + 0 - 1 = 0
    info_result = slicot_ab13md('N', n_zero, nullptr, 1, m_zd, nblock_zd_zero_val.data(), itype_zd.data(), x_zd_case2.empty() ? nullptr : x_zd_case2.data(), &bound_zero_out, nullptr, nullptr, 0 );
    EXPECT_EQ(info_result, 1); 
    EXPECT_DOUBLE_EQ(bound_zero_out, 0.0); 

    // Case 3: N=0, M=1, NBLOCK={1} (Fortran INFO = 2: Sum NBLOCK != N. Sum=1, N=0)
    std::vector<int> nblock_zd_sum_neq_n = {1}; 
    std::vector<double> x_zd_case3; // x_dim = 1 + 0 - 1 = 0
    info_result = slicot_ab13md('N', n_zero, nullptr, 1, m_zd, nblock_zd_sum_neq_n.data(), itype_zd.data(), x_zd_case3.empty() ? nullptr : x_zd_case3.data(), &bound_zero_out, nullptr, nullptr, 0);
    EXPECT_EQ(info_result, 2); 
    EXPECT_DOUBLE_EQ(bound_zero_out, 0.0);
}
