#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // For std::max, std::min
#include <iostream>  // For debugging prints
#include <iomanip>   // For std::fixed, std::setprecision

#include "tc01od.h"
#include "slicot_utils.h" 

// Helper to compare 3D arrays (coefficient by coefficient)
// Assumes C arrays are flat, and slices are laid out contiguously.
// Accesses elements assuming row-major C layout for each slice.
void EXPECT_3D_ARRAY_NEAR_RM_SLICE(const double* expected_flat_c_rm_slice_order,
                                   const double* actual_flat_c_rm_slice_order,
                                   int slice_dim1, int slice_dim2, int num_slices, 
                                   int c_ld_rows_alloc, int c_ld_cols_alloc, // Allocated strides for actual_flat_c
                                   double tol, const std::string& arr_name) {
    if (slice_dim1 == 0 || slice_dim2 == 0 || num_slices == 0) { 
        if (!((expected_flat_c_rm_slice_order == nullptr || expected_flat_c_rm_slice_order[0] == 0) && // crude check for empty
              (actual_flat_c_rm_slice_order == nullptr || actual_flat_c_rm_slice_order[0] == 0) )) {
             // Allow if both are effectively empty or NULL.
        }
        return; 
    }
    ASSERT_NE(expected_flat_c_rm_slice_order, nullptr) << arr_name << " expected data is NULL.";
    ASSERT_NE(actual_flat_c_rm_slice_order, nullptr) << arr_name << " actual data is NULL.";

    for (int k = 0; k < num_slices; ++k) {     // Slice index
        for (int i = 0; i < slice_dim1; ++i) { // Row index of current logical slice
            for (int j = 0; j < slice_dim2; ++j) { // Col index of current logical slice
                // Accessing actual C flat array (row-major slice order)
                size_t actual_idx = (size_t)k * c_ld_rows_alloc * c_ld_cols_alloc + 
                                    (size_t)i * c_ld_cols_alloc + j;
                // Accessing expected C flat array (assumed tightly packed row-major slices)
                size_t expected_idx = (size_t)k * slice_dim1 * slice_dim2 +
                                      (size_t)i * slice_dim2 + j;
                
                ASSERT_LT(actual_idx, (size_t)num_slices * c_ld_rows_alloc * c_ld_cols_alloc) << "Actual index out of bounds";
                ASSERT_LT(expected_idx, (size_t)num_slices * slice_dim1 * slice_dim2) << "Expected index out of bounds";


                EXPECT_NEAR(actual_flat_c_rm_slice_order[actual_idx], expected_flat_c_rm_slice_order[expected_idx], tol)
                    << arr_name << " mismatch at slice " << k << ", row " << i << ", col " << j;
            }
        }
    }
}


class Tc01odTest : public ::testing::Test {
protected:
    char leri = 'L';
    int m = 2;
    int p = 2;
    int indlim = 3;

    std::vector<double> pcoeff_in_flat_rm_slice; 
    std::vector<double> qcoeff_in_flat_rm_slice_PxM; 

    std::vector<double> pcoeff_expected_flat_rm_slice;
    std::vector<double> qcoeff_expected_flat_rm_slice_MxP; 

    std::vector<double> pcoeff_test_buffer;
    std::vector<double> qcoeff_test_buffer;

    double tol = 1e-5;

    // C leading dimensions for row-major layout
    int ldpcoeff_c_rm_rows, ldpcoeff_c_rm_cols;
    int ldqcoeff_c_rm_rows_alloc, ldqcoeff_c_rm_cols_alloc; 
    
    // C leading dimensions for column-major layout (used for Fortran LDs)
    int ldpcoeff_f_ld1, ldpcoeff_f_ld2;
    int ldqcoeff_f_ld1, ldqcoeff_f_ld2;


    void SetUp() override {
        int porm_val = (leri == 'L') ? p : m;

        pcoeff_in_flat_rm_slice = {
            2.0, 4.0, 5.0, 3.0,  // k=0 (s^2)
            3.0, -1.0, 7.0, 2.0, // k=1 (s^1)
            1.0, -1.0, -6.0, 2.0 // k=2 (s^0)
        };
        pcoeff_expected_flat_rm_slice = {
            2.0, 5.0, 4.0, 3.0,   // k=0
            3.0, 7.0, -1.0, 2.0,  // k=1
            1.0, -6.0, -1.0, 2.0  // k=2
        };

        qcoeff_in_flat_rm_slice_PxM = { // PxM = 2x2 slices
            6.0, 1.0, 1.0, 4.0,    // k=0
            -1.0, 7.0, 1.0, 1.0,   // k=1
            5.0, 5.0, 1.0, -1.0    // k=2
        };
        qcoeff_expected_flat_rm_slice_MxP = { // MxP = 2x2 slices
            6.0, 1.0, 1.0, 4.0,    // k=0
            -1.0, 1.0, 7.0, 1.0,   // k=1
            5.0, 1.0, 5.0, -1.0    // k=2
        };

        // For row-major C, dimensions of slices for allocation
        ldpcoeff_c_rm_rows = porm_val;
        ldpcoeff_c_rm_cols = porm_val;
        ldqcoeff_c_rm_rows_alloc = std::max(p,m); 
        ldqcoeff_c_rm_cols_alloc = std::max(p,m); 
        
        // Fortran LDs (these are the number of rows for the first dimension, etc.)
        ldpcoeff_f_ld1 = MAX(1, porm_val);
        ldpcoeff_f_ld2 = MAX(1, porm_val);
        ldqcoeff_f_ld1 = MAX(1, MAX(m,p));
        ldqcoeff_f_ld2 = MAX(1, MAX(m,p));
    }
};

TEST_F(Tc01odTest, ExampleColMajor) {
    int porm_val = (leri == 'L') ? p : m;
    
    // Prepare column-major input buffers for the test
    pcoeff_test_buffer.resize((size_t)ldpcoeff_f_ld1 * ldpcoeff_f_ld2 * indlim);
    qcoeff_test_buffer.resize((size_t)ldqcoeff_f_ld1 * ldqcoeff_f_ld2 * indlim);
    
    // Convert C row-major slice-flat input to Fortran column-major flat for pcoeff_test_buffer
    for(int k=0; k<indlim; ++k) {
        double* c_rm_slice_start = pcoeff_in_flat_rm_slice.data() + k * (size_t)porm_val * porm_val;
        double* f_cm_slice_start = pcoeff_test_buffer.data() + k * (size_t)ldpcoeff_f_ld1 * ldpcoeff_f_ld2;
        slicot_transpose_to_fortran_with_ld(c_rm_slice_start, f_cm_slice_start, 
                                            porm_val, porm_val, 
                                            porm_val, ldpcoeff_f_ld1, sizeof(double));
    }

    // Convert C row-major slice-flat input (PxM) to Fortran column-major flat for qcoeff_test_buffer (MAXMP x MAXMP buffer)
    memset(qcoeff_test_buffer.data(), 0, qcoeff_test_buffer.size() * sizeof(double)); // Zero out buffer
    if (p > 0 && m > 0) {
        for(int k=0; k<indlim; ++k) {
            double* c_rm_slice_start = qcoeff_in_flat_rm_slice_PxM.data() + k * (size_t)p * m;
            double* f_cm_slice_start = qcoeff_test_buffer.data() + k * (size_t)ldqcoeff_f_ld1 * ldqcoeff_f_ld2;
            slicot_transpose_to_fortran_with_ld(c_rm_slice_start, f_cm_slice_start, 
                                                p, m, 
                                                m, ldqcoeff_f_ld1, sizeof(double));
        }
    }

    int info = slicot_tc01od(leri, m, p, indlim,
                             pcoeff_test_buffer.data(), ldpcoeff_f_ld1, ldpcoeff_f_ld2, // Pass Fortran LDs for CM C
                             qcoeff_test_buffer.data(), ldqcoeff_f_ld1, ldqcoeff_f_ld2, // Pass Fortran LDs for CM C
                             0 /*row_major=false*/);
    ASSERT_EQ(info, 0);

    // Prepare expected results in Fortran column-major flat order
    std::vector<double> pcoeff_expected_flat_cm(pcoeff_test_buffer.size());
    std::vector<double> qcoeff_expected_flat_cm(qcoeff_test_buffer.size());

    for(int k=0; k<indlim; ++k) {
        double* c_rm_slice_start = pcoeff_expected_flat_rm_slice.data() + k * (size_t)porm_val * porm_val;
        double* f_cm_slice_start = pcoeff_expected_flat_cm.data() + k * (size_t)ldpcoeff_f_ld1 * ldpcoeff_f_ld2;
        slicot_transpose_to_fortran_with_ld(c_rm_slice_start, f_cm_slice_start, 
                                            porm_val, porm_val, 
                                            porm_val, ldpcoeff_f_ld1, sizeof(double));
    }
    
    memset(qcoeff_expected_flat_cm.data(), 0, qcoeff_expected_flat_cm.size() * sizeof(double));
    if (m > 0 && p > 0) { // Output is MxP
        for(int k=0; k<indlim; ++k) {
            double* c_rm_slice_start = qcoeff_expected_flat_rm_slice_MxP.data() + k * (size_t)m * p;
            double* f_cm_slice_start = qcoeff_expected_flat_cm.data() + k * (size_t)ldqcoeff_f_ld1 * ldqcoeff_f_ld2;
            slicot_transpose_to_fortran_with_ld(c_rm_slice_start, f_cm_slice_start, 
                                                m, p, 
                                                p, ldqcoeff_f_ld1, sizeof(double));
        }
    }

    for(size_t i=0; i < pcoeff_test_buffer.size(); ++i) {
        EXPECT_NEAR(pcoeff_test_buffer[i], pcoeff_expected_flat_cm[i], tol) << "PCOEFF CM mismatch at flat index " << i;
    }
    
    // Compare the relevant MxP part of QCOEFF output
    for(int k=0; k<indlim; ++k) {
        for(int j_col_out=0; j_col_out<p; ++j_col_out) { // Output has P columns
            for(int i_row_out=0; i_row_out<m; ++i_row_out) { // Output has M rows
                size_t idx = (size_t)i_row_out + (size_t)j_col_out*ldqcoeff_f_ld1 + (size_t)k*ldqcoeff_f_ld1*ldqcoeff_f_ld2;
                EXPECT_NEAR(qcoeff_test_buffer[idx], qcoeff_expected_flat_cm[idx], tol) 
                    << "QCOEFF CM mismatch at k=" << k << ", out_row=" << i_row_out << ", out_col=" << j_col_out;
            }
        }
    }
}


TEST_F(Tc01odTest, ExampleRowMajor) {
    int porm_val = (leri == 'L') ? p : m;
    pcoeff_test_buffer = pcoeff_in_flat_rm_slice; 
    
    qcoeff_test_buffer.assign((size_t)ldqcoeff_c_rm_rows_alloc * ldqcoeff_c_rm_cols_alloc * indlim, 0.0);
    if (p > 0 && m > 0) { // Copy PxM input into the allocated buffer
        for(int k=0; k<indlim; ++k) {
            for(int i_row_in=0; i_row_in<p; ++i_row_in) { 
                for(int j_col_in=0; j_col_in<m; ++j_col_in) { 
                    qcoeff_test_buffer[k*((size_t)ldqcoeff_c_rm_rows_alloc * ldqcoeff_c_rm_cols_alloc) + i_row_in*ldqcoeff_c_rm_cols_alloc + j_col_in] =
                        qcoeff_in_flat_rm_slice_PxM[k*((size_t)p*m) + i_row_in*m + j_col_in];
                }
            }
        }
    }


    int info = slicot_tc01od(leri, m, p, indlim,
                             pcoeff_test_buffer.data(), ldpcoeff_c_rm_rows, ldpcoeff_c_rm_cols,
                             qcoeff_test_buffer.data(), ldqcoeff_c_rm_rows_alloc, ldqcoeff_c_rm_cols_alloc,
                             1 /*row_major=true*/);
    ASSERT_EQ(info, 0);

    EXPECT_3D_ARRAY_NEAR_RM_SLICE(pcoeff_expected_flat_rm_slice.data(), pcoeff_test_buffer.data(),
                                  porm_val, porm_val, indlim, 
                                  ldpcoeff_c_rm_rows, ldpcoeff_c_rm_cols, tol, "PCOEFF RM");
    
    EXPECT_3D_ARRAY_NEAR_RM_SLICE(qcoeff_expected_flat_rm_slice_MxP.data(), qcoeff_test_buffer.data(),
                                  m, p, indlim, // Logical output dimensions MxP
                                  ldqcoeff_c_rm_rows_alloc, ldqcoeff_c_rm_cols_alloc, tol, "QCOEFF RM");
}


TEST_F(Tc01odTest, ParameterValidation) {
    std::vector<double> dummy_p(1);
    std::vector<double> dummy_q(1);
    int porm_val = (leri == 'L') ? p : m;
    int max_mp = std::max(m,p);
    if (max_mp == 0) max_mp = 1; // Ensure non-zero for dummy allocation
    if (porm_val == 0) porm_val = 1;


    EXPECT_EQ(slicot_tc01od('X', m, p, indlim, dummy_p.data(), porm_val, porm_val, dummy_q.data(), max_mp, max_mp, 0), -1); 
    EXPECT_EQ(slicot_tc01od(leri, -1, p, indlim, dummy_p.data(), porm_val, porm_val, dummy_q.data(), max_mp, max_mp, 0), -2); 
    EXPECT_EQ(slicot_tc01od(leri, m, -1, indlim, dummy_p.data(), porm_val, porm_val, dummy_q.data(), max_mp, max_mp, 0), -3); 
    EXPECT_EQ(slicot_tc01od(leri, m, p, 0, dummy_p.data(), porm_val, porm_val, dummy_q.data(), max_mp, max_mp, 0), -4);    
    
    int current_porm = (leri == 'L') ? this->p : this->m;
    int current_max_mp = std::max(this->m, this->p);

    if (current_porm > 0) {
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, nullptr, current_porm, current_porm, dummy_q.data(), current_max_mp, current_max_mp, 0), -5);    
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, dummy_p.data(), 0, current_porm, dummy_q.data(), current_max_mp, current_max_mp, 0), -6); 
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, dummy_p.data(), current_porm, 0, dummy_q.data(), current_max_mp, current_max_mp, 0), -6); 
    }
    if (this->m > 0 || this->p > 0) {
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, dummy_p.data(), current_porm, current_porm, nullptr, current_max_mp, current_max_mp, 0), -8);    
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, dummy_p.data(), current_porm, current_porm, dummy_q.data(), 0, current_max_mp, 0), -9); 
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, dummy_p.data(), current_porm, current_porm, dummy_q.data(), current_max_mp, 0, 0), -9); 
    }
}
