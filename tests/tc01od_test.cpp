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
        // If logically empty, we expect actual_flat_c_rm_slice_order to not have been modified
        // or to be NULL if it was passed as NULL.
        // This helper might not be the best for zero-dim arrays if complex NULL handling is needed.
        // For TC01OD, if input dims are 0, output dims are 0, so no elements to compare.
        // However, if actual is not NULL and expected is (or vice versa for non-empty logical), it's an issue.
        if (!((expected_flat_c_rm_slice_order == nullptr && actual_flat_c_rm_slice_order == nullptr) ||
              (expected_flat_c_rm_slice_order != nullptr && actual_flat_c_rm_slice_order != nullptr && num_slices*slice_dim1*slice_dim2 == 0 )
            )) {
             // This condition is tricky. If dims are 0, ideally both pointers are effectively NULL or point to zero-sized data.
             // For this test, if any dim is 0, we just return.
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
                
                ASSERT_LT(actual_idx, (size_t)num_slices * c_ld_rows_alloc * c_ld_cols_alloc) << "Actual index out of bounds for " << arr_name;
                ASSERT_LT(expected_idx, (size_t)num_slices * slice_dim1 * slice_dim2) << "Expected index out of bounds for " << arr_name;


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

        ldpcoeff_c_rm_rows = porm_val;
        ldpcoeff_c_rm_cols = porm_val;
        ldqcoeff_c_rm_rows_alloc = std::max(p,m); 
        ldqcoeff_c_rm_cols_alloc = std::max(p,m); 
        
        ldpcoeff_f_ld1 = MAX(1, porm_val);
        ldpcoeff_f_ld2 = MAX(1, porm_val);
        ldqcoeff_f_ld1 = MAX(1, MAX(m,p));
        ldqcoeff_f_ld2 = MAX(1, MAX(m,p));
    }
};

TEST_F(Tc01odTest, ExampleColMajor) {
    int porm_val = (leri == 'L') ? p : m;
    
    pcoeff_test_buffer.resize((size_t)ldpcoeff_f_ld1 * ldpcoeff_f_ld2 * indlim);
    qcoeff_test_buffer.resize((size_t)ldqcoeff_f_ld1 * ldqcoeff_f_ld2 * indlim);
    
    for(int k=0; k<indlim; ++k) {
        double* c_rm_slice_start = pcoeff_in_flat_rm_slice.data() + k * (size_t)porm_val * porm_val;
        double* f_cm_slice_start = pcoeff_test_buffer.data() + k * (size_t)ldpcoeff_f_ld1 * ldpcoeff_f_ld2;
        slicot_transpose_to_fortran_with_ld(c_rm_slice_start, f_cm_slice_start, 
                                            porm_val, porm_val, 
                                            porm_val, ldpcoeff_f_ld1, sizeof(double));
    }

    memset(qcoeff_test_buffer.data(), 0, qcoeff_test_buffer.size() * sizeof(double)); 
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
                             pcoeff_test_buffer.data(), ldpcoeff_f_ld1, ldpcoeff_f_ld2, 
                             qcoeff_test_buffer.data(), ldqcoeff_f_ld1, ldqcoeff_f_ld2, 
                             0 /*row_major=false*/);
    ASSERT_EQ(info, 0);

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
    
    for(int k=0; k<indlim; ++k) {
        for(int j_col_out=0; j_col_out<p; ++j_col_out) { 
            for(int i_row_out=0; i_row_out<m; ++i_row_out) { 
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
    if (p > 0 && m > 0) { 
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
                                  m, p, indlim, 
                                  ldqcoeff_c_rm_rows_alloc, ldqcoeff_c_rm_cols_alloc, tol, "QCOEFF RM");
}


TEST_F(Tc01odTest, ParameterValidation) {
    std::vector<double> dummy_p(1);
    std::vector<double> dummy_q(1);
    int porm_val = (leri == 'L') ? p : m;
    int max_mp = std::max(m,p);
    if (max_mp == 0) max_mp = 1; 
    if (porm_val == 0) porm_val = 1;

    EXPECT_EQ(slicot_tc01od('X', m, p, indlim, dummy_p.data(), porm_val, porm_val, dummy_q.data(), max_mp, max_mp, 0), -1); 
    EXPECT_EQ(slicot_tc01od(leri, -1, p, indlim, dummy_p.data(), porm_val, porm_val, dummy_q.data(), max_mp, max_mp, 0), -2); 
    EXPECT_EQ(slicot_tc01od(leri, m, -1, indlim, dummy_p.data(), porm_val, porm_val, dummy_q.data(), max_mp, max_mp, 0), -3); 
    EXPECT_EQ(slicot_tc01od(leri, m, p, 0, dummy_p.data(), porm_val, porm_val, dummy_q.data(), max_mp, max_mp, 0), -4);    
    
    int current_porm = (leri == 'L') ? this->p : this->m;
    int current_max_mp = std::max(this->m, this->p);
    if (current_porm == 0) current_porm = 1; // for dummy LDs
    if (current_max_mp == 0) current_max_mp = 1;


    if (this->p > 0 || (this->leri == 'R' && this->m >0) ) { // If PCOEFF is not logically empty
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, nullptr, current_porm, current_porm, dummy_q.data(), current_max_mp, current_max_mp, 0), -5);    
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, dummy_p.data(), 0, current_porm, dummy_q.data(), current_max_mp, current_max_mp, 0), -6); 
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, dummy_p.data(), current_porm, 0, dummy_q.data(), current_max_mp, current_max_mp, 0), -6); 
    }
    if (this->m > 0 || this->p > 0) { // If QCOEFF is not logically empty
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, dummy_p.data(), current_porm, current_porm, nullptr, current_max_mp, current_max_mp, 0), -8);    
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, dummy_p.data(), current_porm, current_porm, dummy_q.data(), 0, current_max_mp, 0), -9); 
        EXPECT_EQ(slicot_tc01od(leri, this->m, this->p, indlim, dummy_p.data(), current_porm, current_porm, dummy_q.data(), current_max_mp, 0, 0), -9); 
    }
}

TEST_F(Tc01odTest, ZeroDimensionM) {
    char leri_z = 'L'; 
    int m_z = 0;
    int p_z = 2;
    int indlim_z = 1;
    int porm_z = p_z; 

    std::vector<double> pcoeff_z((size_t)porm_z * porm_z * indlim_z); 
    for(size_t i=0; i<pcoeff_z.size(); ++i) pcoeff_z[i] = (double)i+1.0; // Fill with some data
    
    // QCOEFF is PxM -> Px0. Logically empty.
    // Fortran expects QCOEFF buffer based on MAX(M,P) for LDs.
    // LDQCO1 >= MAX(0,2)=2, LDQCO2 >= MAX(0,2)=2.
    // The C wrapper will allocate qcoeff_cm if (p_c > 0 || m_c > 0). Here p_c=2.
    // The test should pass a valid buffer for qcoeff_c if it's not logically empty for the wrapper's validation.
    // The wrapper's validation: `if (qcoeff_c == NULL && (p_c > 0 || m_c > 0) )`
    // Here (p_c > 0 || m_c > 0) is true (2 > 0 || 0 > 0). So qcoeff_c cannot be NULL.
    int ldq_alloc_rows = std::max(1, std::max(p_z, m_z));
    int ldq_alloc_cols = std::max(1, std::max(p_z, m_z));
    std::vector<double> qcoeff_z_buffer((size_t)ldq_alloc_rows * ldq_alloc_cols * indlim_z, 0.0);


    int ldp_rows_rm = porm_z;
    int ldp_cols_rm = porm_z;
    
    int info = slicot_tc01od(leri_z, m_z, p_z, indlim_z,
                             pcoeff_z.data(), ldp_rows_rm, ldp_cols_rm,
                             qcoeff_z_buffer.data(), ldq_alloc_rows, ldq_alloc_cols, 
                             1 /*row_major=true*/);
    ASSERT_EQ(info, 0);
    // PCOEFF should be transposed. QCOEFF output is MxP -> 0xP, effectively no data to check beyond success.
}

TEST_F(Tc01odTest, ZeroDimensionP) {
    char leri_z = 'L'; 
    int m_z = 2;
    int p_z = 0;
    int indlim_z = 1;
    int porm_z = p_z; // PCOEFF is 0x0

    // PCOEFF is logically 0x0. Wrapper expects NULL if porm_f is 0.
    // QCOEFF is PxM -> 0xM. Logically empty.
    // Fortran expects QCOEFF buffer based on MAX(M,P) for LDs.
    // LDQCO1 >= MAX(2,0)=2, LDQCO2 >= MAX(2,0)=2.
    // Wrapper validation: `if (qcoeff_c == NULL && (p_c > 0 || m_c > 0) )`
    // Here (0 > 0 || 2 > 0) is true. So qcoeff_c cannot be NULL.
    int ldq_alloc_rows = std::max(1, std::max(p_z, m_z));
    int ldq_alloc_cols = std::max(1, std::max(p_z, m_z));
    std::vector<double> qcoeff_z_buffer((size_t)ldq_alloc_rows * ldq_alloc_cols * indlim_z, 0.0);

    int ldp_rows_rm = std::max(1,porm_z); 
    int ldp_cols_rm = std::max(1,porm_z);
    
    int info = slicot_tc01od(leri_z, m_z, p_z, indlim_z,
                             nullptr, ldp_rows_rm, ldp_cols_rm, 
                             qcoeff_z_buffer.data(), ldq_alloc_rows, ldq_alloc_cols, 
                             1 /*row_major=true*/);
    ASSERT_EQ(info, 0);
    // PCOEFF is 0x0. QCOEFF output is MxP -> Mx0, effectively no data to check.
}

TEST_F(Tc01odTest, AllZeroDimensions) {
    char leri_z = 'L';
    int m_z = 0;
    int p_z = 0;
    int indlim_z = 1;
    int porm_z = p_z; // 0

    // PCOEFF is 0x0. QCOEFF is 0x0.
    // Wrapper validation: `if (pcoeff_c == NULL && porm_f > 0)` -> false
    // Wrapper validation: `if (qcoeff_c == NULL && (p_c > 0 || m_c > 0) )` -> false
    // So, NULL can be passed.

    int ldp_rows_rm = 1; 
    int ldp_cols_rm = 1;
    int ldq_rows_alloc_rm = 1;
    int ldq_cols_alloc_rm = 1;

    int info = slicot_tc01od(leri_z, m_z, p_z, indlim_z,
                             nullptr, ldp_rows_rm, ldp_cols_rm,
                             nullptr, ldq_rows_alloc_rm, ldq_cols_alloc_rm,
                             1 /*row_major=true*/);
    ASSERT_EQ(info, 0);
}
