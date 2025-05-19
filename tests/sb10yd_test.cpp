#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <numeric>   // For std::iota
#include <algorithm> // For std::max, std::min
#include <iostream>

// Define M_PI if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "sb10yd.h"       // Include the wrapper header
#include "slicot_utils.h" // For transpose functions if needed for A

// Helper function to generate simple frequency response data for testing
void generate_dummy_freq_data(int lendat, 
                              std::vector<double>& rfrdat, 
                              std::vector<double>& ifrdat, 
                              std::vector<double>& omega_data, 
                              int discfl) {
    if (lendat <= 0) { // Guard against zero or negative size
        rfrdat.clear();
        ifrdat.clear();
        omega_data.clear();
        return;
    }
    rfrdat.resize(lendat);
    ifrdat.resize(lendat);
    omega_data.resize(lendat);

    for (int i = 0; i < lendat; ++i) {
        if (discfl == 0) { // Continuous
            omega_data[i] = static_cast<double>(i + 1) * 0.1; 
        } else { // Discrete
            // Ensure lendat-1 is not zero if lendat is 1 (though lendat >= 2 is a precondition for SB10YD)
            if (lendat > 1) {
                omega_data[i] = (static_cast<double>(i) / (lendat - 1)) * M_PI; 
            } else { // lendat is 1 (though invalid for SB10YD, handle for safety)
                omega_data[i] = 0.0;
            }
        }
        // Generic dummy response, avoid division by zero if omega_data[i] is 0
        if (std::abs(omega_data[i]) < 1e-9 && omega_data[i]*omega_data[i] < 1e-9) { // if omega is effectively zero
             rfrdat[i] = 1.0;
             ifrdat[i] = 0.0;
        } else {
            rfrdat[i] = 1.0 / (1.0 + omega_data[i]*omega_data[i]); 
            ifrdat[i] = -omega_data[i] / (1.0 + omega_data[i]*omega_data[i]); 
        }
    }
    if (lendat > 0 && discfl == 1) omega_data[0] = 0.0; 
}


// --- Test Fixture Base ---
class SB10YDTest : public ::testing::Test {
protected:
    int LENDAT_param = 50; 
    double TOL_param = 0.0;   

    std::vector<double> RFRDAT_data, IFRDAT_data, OMEGA_data;
    std::vector<double> A_out, B_out, C_out, D_out_scalar; 
    int N_io_param; 

    int expected_info = 0;
    double check_tol_abcd = 1e-1; 

    int LDA;

    void SetUpTest(int n_in, int discfl, int flag) {
        N_io_param = n_in;
        // LENDAT_param is a fixture member, typically initialized.
        // Ensure it's what's expected for the test if SetUpTest is the primary setup path.
        // For ParameterValidation, LENDAT_param is used directly from fixture.
        // Here, we ensure it's 50 if this helper is used for other tests.
        if (LENDAT_param <= 1) LENDAT_param = 50;


        RFRDAT_data.resize(LENDAT_param);
        IFRDAT_data.resize(LENDAT_param);
        OMEGA_data.resize(LENDAT_param);

        if (n_in == 0 && discfl == 0) { 
            for(int i=0; i<LENDAT_param; ++i) {
                OMEGA_data[i] = static_cast<double>(i+1) * 0.1; 
                RFRDAT_data[i] = 1.0; 
                IFRDAT_data[i] = 0.0;
            }
        } else if (n_in == 2 && discfl == 0) { 
            for(int i=0; i<LENDAT_param; ++i) {
                double w = static_cast<double>(i+1) * (10.0 / LENDAT_param); 
                OMEGA_data[i] = w;
                double den_real = -w*w + 0.5;
                double den_imag = 0.1*w;
                double common_den = den_real*den_real + den_imag*den_imag;
                if (std::abs(common_den) < 1e-9) common_den = 1e-9; 
                RFRDAT_data[i] = (0.1*w*w) / common_den;
                IFRDAT_data[i] = (w*(-w*w+0.5)) / common_den;
            }
        }
        else { 
            generate_dummy_freq_data(LENDAT_param, RFRDAT_data, IFRDAT_data, OMEGA_data, discfl);
        }

        if (N_io_param > 0) {
            A_out.resize(N_io_param * N_io_param);
            B_out.resize(N_io_param);
            C_out.resize(N_io_param);
        } else { 
            A_out.clear(); 
            B_out.clear(); 
            C_out.clear(); 
        }
        D_out_scalar.resize(1); 
        LDA = std::max(1, N_io_param); 
    }
};


TEST_F(SB10YDTest, ContN0Flag0) {
    int n_initial = 0;
    int discfl = 0; 
    int flag = 0;   
    SetUpTest(n_initial, discfl, flag);
    
    int n_val_returned = n_initial; 

    int info_result = slicot_sb10yd(discfl, flag, LENDAT_param,
                                   RFRDAT_data.data(), IFRDAT_data.data(), OMEGA_data.data(),
                                   &n_val_returned,
                                   A_out.empty() ? nullptr : A_out.data(), LDA,
                                   B_out.empty() ? nullptr : B_out.data(), 
                                   C_out.empty() ? nullptr : C_out.data(), 
                                   D_out_scalar.data(),
                                   TOL_param, 0 /*row_major=false*/);

    ASSERT_EQ(info_result, expected_info);
    if (flag == 0) {
        ASSERT_EQ(n_val_returned, n_initial);
    }
}

TEST_F(SB10YDTest, ContN2Flag0) {
    int n_initial = 2;
    int discfl = 0; 
    int flag = 0;   
    SetUpTest(n_initial, discfl, flag);
    
    int n_val_returned = n_initial;

    int info_result = slicot_sb10yd(discfl, flag, LENDAT_param,
                                   RFRDAT_data.data(), IFRDAT_data.data(), OMEGA_data.data(),
                                   &n_val_returned,
                                   A_out.data(), LDA, B_out.data(), C_out.data(), D_out_scalar.data(),
                                   TOL_param, 0 /*row_major=false*/);

    ASSERT_EQ(info_result, expected_info);
    if (flag == 0) {
        ASSERT_EQ(n_val_returned, n_initial);
    }
    if (info_result == 0 && n_val_returned > 0) {
        for(size_t i=0; i < (size_t)n_val_returned * n_val_returned; ++i) EXPECT_FALSE(std::isnan(A_out[i]));
        for(int i=0; i < n_val_returned; ++i) EXPECT_FALSE(std::isnan(B_out[i]));
        for(int i=0; i < n_val_returned; ++i) EXPECT_FALSE(std::isnan(C_out[i]));
        EXPECT_FALSE(std::isnan(D_out_scalar[0]));
    }
}

TEST_F(SB10YDTest, DiscN0Flag0) {
    int n_initial = 0;
    int discfl = 1; 
    int flag = 0;   
    SetUpTest(n_initial, discfl, flag);
    
    int n_val_returned = n_initial;

    int info_result = slicot_sb10yd(discfl, flag, LENDAT_param,
                                   RFRDAT_data.data(), IFRDAT_data.data(), OMEGA_data.data(),
                                   &n_val_returned,
                                   A_out.empty() ? nullptr : A_out.data(), LDA,
                                   B_out.empty() ? nullptr : B_out.data(), 
                                   C_out.empty() ? nullptr : C_out.data(), 
                                   D_out_scalar.data(),
                                   TOL_param, 0 /*row_major=false*/);

    ASSERT_EQ(info_result, expected_info);
    if (flag == 0) {
        ASSERT_EQ(n_val_returned, n_initial);
    }
}


TEST_F(SB10YDTest, ContN2Flag1Stability) {
    int n_initial = 2;
    int discfl = 0; 
    int flag = 1;   
    SetUpTest(n_initial, discfl, flag);
    
    int n_val_returned = n_initial;

    int info_result = slicot_sb10yd(discfl, flag, LENDAT_param,
                                   RFRDAT_data.data(), IFRDAT_data.data(), OMEGA_data.data(),
                                   &n_val_returned,
                                   A_out.data(), LDA, B_out.data(), C_out.data(), D_out_scalar.data(),
                                   TOL_param, 0 /*row_major=false*/);

    ASSERT_EQ(info_result, expected_info);
    EXPECT_GE(n_val_returned, 0); 
    EXPECT_LE(n_val_returned, n_initial); 
}


TEST_F(SB10YDTest, ParameterValidation) {
    int n_val = 0;
    int info_result;
    // LENDAT_param is 50 by default from fixture. TOL_param is 0.0.

    // Test invalid DISCFL
    info_result = slicot_sb10yd(2, 0, LENDAT_param, nullptr,nullptr,nullptr, &n_val, nullptr,1,nullptr,nullptr,nullptr,TOL_param,0);
    EXPECT_EQ(info_result, -1);

    // Test invalid LENDAT
    info_result = slicot_sb10yd(0, 0, 1, nullptr,nullptr,nullptr, &n_val, nullptr,1,nullptr,nullptr,nullptr,TOL_param,0);
    EXPECT_EQ(info_result, -3);
    
    // Test NULL n_io
    // Ensure RFRDAT_data etc. are valid for this call, even if n_io is NULL
    RFRDAT_data.resize(LENDAT_param); 
    IFRDAT_data.resize(LENDAT_param);
    OMEGA_data.resize(LENDAT_param);
    generate_dummy_freq_data(LENDAT_param, RFRDAT_data, IFRDAT_data, OMEGA_data, 0); 
    info_result = slicot_sb10yd(0, 0, LENDAT_param, RFRDAT_data.data(),IFRDAT_data.data(),OMEGA_data.data(), nullptr, nullptr,1,nullptr,nullptr,nullptr,TOL_param,0);
    EXPECT_EQ(info_result, -7);

    // Test N > LENDAT-1
    n_val = LENDAT_param; // n_val = 50 (LENDAT_param is 50 from fixture)
    
    // Manual setup for this specific case to ensure RFRDAT_data and other inputs are valid.
    // This bypasses SetUpTest for this specific sub-case to isolate the issue.
    RFRDAT_data.resize(LENDAT_param); 
    IFRDAT_data.resize(LENDAT_param);
    OMEGA_data.resize(LENDAT_param);
    generate_dummy_freq_data(LENDAT_param, RFRDAT_data, IFRDAT_data, OMEGA_data, 0); 

    // Size A_out, B_out, C_out based on n_val for this call
    if (n_val > 0) {
        A_out.resize((size_t)n_val * n_val); // Use size_t for safety
        B_out.resize(n_val);
        C_out.resize(n_val);
    } else {
        A_out.clear(); B_out.clear(); C_out.clear();
    }
    D_out_scalar.resize(1); // D is always scalar
    LDA = std::max(1, n_val); // LDA for A
    
    info_result = slicot_sb10yd(0, 0, LENDAT_param, 
                                   RFRDAT_data.data(), IFRDAT_data.data(), OMEGA_data.data(),
                                   &n_val, 
                                   A_out.empty() ? nullptr : A_out.data(), LDA,
                                   B_out.empty() ? nullptr : B_out.data(), 
                                   C_out.empty() ? nullptr : C_out.data(), 
                                   D_out_scalar.data(),
                                   TOL_param, 0 /*row_major=false*/);
    EXPECT_EQ(info_result, -7); 
}

// Row-major specific test
TEST_F(SB10YDTest, ContN2Flag0RowMajor) {
    int n_initial = 2;
    int discfl = 0; 
    int flag = 0;   
    SetUpTest(n_initial, discfl, flag); 
    
    LDA = n_initial; 
    if (n_initial == 0) LDA = 1;

    // Create copies for row-major call because SetUpTest populates A_out etc.
    // which might be sized for column-major if other tests ran first (though GTest isolates fixtures)
    // Best to use fresh vectors or ensure correct sizing for row-major A.
    std::vector<double> a_rm_out;
    if (n_initial > 0) a_rm_out.resize((size_t)LDA * n_initial); // LDA is cols, N_initial is rows for RM A
    else a_rm_out.clear();
    
    std::vector<double> b_rm_out = B_out; // B, C, D are 1D or scalar-like, direct copy is fine
    std::vector<double> c_rm_out = C_out;
    std::vector<double> d_rm_out = D_out_scalar;

    int n_val_returned = n_initial;

    int info_result = slicot_sb10yd(discfl, flag, LENDAT_param,
                                   RFRDAT_data.data(), IFRDAT_data.data(), OMEGA_data.data(),
                                   &n_val_returned,
                                   a_rm_out.empty() ? nullptr : a_rm_out.data(), LDA, 
                                   b_rm_out.empty() ? nullptr : b_rm_out.data(), 
                                   c_rm_out.empty() ? nullptr : c_rm_out.data(), 
                                   d_rm_out.data(),
                                   TOL_param, 1 /*row_major=true*/);

    ASSERT_EQ(info_result, expected_info);
    if (flag == 0) {
        ASSERT_EQ(n_val_returned, n_initial);
    }
    if (info_result == 0 && n_val_returned > 0 && !a_rm_out.empty()) {
        EXPECT_FALSE(std::isnan(a_rm_out[0])); 
    }
}

