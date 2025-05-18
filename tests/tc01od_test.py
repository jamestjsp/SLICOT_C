import numpy as np

def print_poly_matrix_coeffs_py(name, coeffs_3d_list_of_slices, num_rows, num_cols, indlim_val):
    """
    Prints polynomial matrix coefficients from a list of 2D numpy slices.
    """
    print(f"{name} (slices for s^({indlim_val-1}) down to s^0):")
    if not coeffs_3d_list_of_slices:
        print("  None or empty.")
        return
    for k_idx in range(indlim_val):
        print(f"  Coefficients for s^({indlim_val-1-k_idx}) (Slice k_idx={k_idx}):")
        if k_idx < len(coeffs_3d_list_of_slices):
            slice_k = coeffs_3d_list_of_slices[k_idx]
            if slice_k.shape == (num_rows, num_cols):
                 for i_row in range(num_rows):
                    row_str = "    ["
                    for j_col in range(num_cols):
                        row_str += f"{slice_k[i_row, j_col]:8.2f}" # Changed to .2f to match example output
                        if j_col < num_cols - 1:
                            row_str += ", "
                    row_str += "]"
                    print(row_str)
            else:
                 print(f"    Unexpected slice shape: {slice_k.shape}, expected ({num_rows},{num_cols})")
        else:
            print(f"    Slice k_idx={k_idx} not available (indlim_val might be > actual num_coeffs).")
    print("-" * 30)

def tc01od_slycot_equivalent(leri_char, m_val, p_val, indlim_val, pcoeff_3d_fortran, qcoeff_3d_fortran_PxM):
    """
    Simulates TC01OD by transposing coefficient matrices.
    pcoeff_3d_fortran: (porm, porm, indlim) Fortran-ordered numpy array for P.
    qcoeff_3d_fortran_PxM: (P, M, indlim) Fortran-ordered numpy array for Q.
    Returns P_dual, Q_dual as Fortran-ordered numpy arrays.
    """
    pcoeff_dual = np.zeros_like(pcoeff_3d_fortran)
    for k in range(indlim_val):
        pcoeff_dual[:,:,k] = pcoeff_3d_fortran[:,:,k].T
    
    # QCOEFF input is (P, M, indlim)
    # Q_dual is (M, P, indlim)
    qcoeff_dual = np.zeros((m_val, p_val, indlim_val)) # Initialize with correct output dimensions
    if p_val > 0 and m_val > 0 : # only transpose if both dimensions are non-zero
        for k in range(indlim_val):
            qcoeff_dual[:,:,k] = qcoeff_3d_fortran_PxM[:,:,k].T
    elif p_val == 0 and m_val > 0: # Input Q is 0xM, Output Q' is Mx0
        pass # qcoeff_dual is already zeros of shape (M,0,indlim) effectively
    elif m_val == 0 and p_val > 0: # Input Q is Px0, Output Q' is 0xP
        pass # qcoeff_dual is already zeros of shape (0,P,indlim) effectively
    # if both m_val and p_val are 0, qcoeff_dual is (0,0,indlim)

    return pcoeff_dual, qcoeff_dual


def main():
    leri = 'L'
    m = 2
    p = 2
    indlim = 3 

    porm = p if leri == 'L' else m

    pcoeff_f_order = np.zeros((porm, porm, indlim))
    pcoeff_f_order[:,:,0] = np.array([[2.0, 4.0], [5.0, 3.0]]) 
    pcoeff_f_order[:,:,1] = np.array([[3.0, -1.0], [7.0, 2.0]]) 
    pcoeff_f_order[:,:,2] = np.array([[1.0, -1.0], [-6.0, 2.0]]) 
    
    qcoeff_f_order_PxM = np.zeros((p, m, indlim))
    qcoeff_f_order_PxM[:,:,0] = np.array([[6.0, 1.0], [1.0, 4.0]])
    qcoeff_f_order_PxM[:,:,1] = np.array([[-1.0, 7.0], [1.0, 1.0]])
    qcoeff_f_order_PxM[:,:,2] = np.array([[5.0, 5.0], [1.0, -1.0]])

    print("--- Input Parameters (Fortran Order Slices) ---")
    print(f"LERI: '{leri}'")
    print(f"M: {m}")
    print(f"P: {p}")
    print(f"INDLIM: {indlim}")
    
    pcoeff_in_slices_for_print = [pcoeff_f_order[:,:,k] for k in range(indlim)]
    qcoeff_in_slices_for_print = [qcoeff_f_order_PxM[:,:,k] for k in range(indlim)]
    print_poly_matrix_coeffs_py("PCOEFF_in", pcoeff_in_slices_for_print, porm, porm, indlim)
    print_poly_matrix_coeffs_py("QCOEFF_in (PxM)", qcoeff_in_slices_for_print, p, m, indlim)
    print("=" * 50)

    pcoeff_dual_f_order, qcoeff_dual_f_order_MxP = tc01od_slycot_equivalent(
        leri, m, p, indlim, pcoeff_f_order, qcoeff_f_order_PxM
    )

    print("\n--- Output Results (Fortran Order Slices) ---")
    pcoeff_out_slices_for_print = [pcoeff_dual_f_order[:,:,k] for k in range(indlim)]
    qcoeff_out_slices_for_print = [qcoeff_dual_f_order_MxP[:,:,k] for k in range(indlim)]
    print_poly_matrix_coeffs_py("PCOEFF_out (P')", pcoeff_out_slices_for_print, porm, porm, indlim)
    print_poly_matrix_coeffs_py("QCOEFF_out (Q' - MxP)", qcoeff_out_slices_for_print, m, p, indlim)

    
    if pcoeff_dual_f_order is not None and pcoeff_dual_f_order.size > 0 :
        pcoeff_flat_c_rm_slice = []
        for k_idx in range(indlim):
            pcoeff_flat_c_rm_slice.extend(pcoeff_dual_f_order[:,:,k_idx].flatten(order='C').tolist())
            
        print("\nPCOEFF_out flattened for C++ (Row-Major Slices Style):")
        print("{", end="")
        for i, val in enumerate(pcoeff_flat_c_rm_slice):
            print(f"{val:.2f}", end="") 
            if i < len(pcoeff_flat_c_rm_slice) - 1:
                print(", ", end="")
        print("};")
        print("-" * 30)

    if qcoeff_dual_f_order_MxP is not None and qcoeff_dual_f_order_MxP.size > 0:
        qcoeff_flat_c_rm_slice = []
        for k_idx in range(indlim):
            qcoeff_flat_c_rm_slice.extend(qcoeff_dual_f_order_MxP[:,:,k_idx].flatten(order='C').tolist())

        print("\nQCOEFF_out (MxP) flattened for C++ (Row-Major Slices Style):")
        print("{", end="")
        for i, val in enumerate(qcoeff_flat_c_rm_slice):
            print(f"{val:.2f}", end="") 
            if i < len(qcoeff_flat_c_rm_slice) - 1:
                print(", ", end="")
        print("};")
        print("-" * 30)

if __name__ == "__main__":
    main()
