import numpy as np
# Assuming _tc01od_l and _tc01od_r are low-level wrappers if direct transform.tc01od isn't available
# For this example, we'll simulate the transpose operation directly as TC01OD is simple.

def print_poly_matrix_coeffs(name, coeffs_3d_list_of_slices, p, m, indlim):
    """
    Prints polynomial matrix coefficients.
    coeffs_3d_list_of_slices: A list of 2D numpy arrays, where each 2D array is a coefficient slice P(k) or Q(k).
                              Order is P(0)/Q(0), P(1)/Q(1), ..., P(indlim-1)/Q(indlim-1)
                              (corresponding to s^(indlim-1), s^(indlim-2), ..., s^0)
    """
    print(f"{name} (slices for s^({indlim-1}) down to s^0):")
    if not coeffs_3d_list_of_slices:
        print("  None or empty.")
        return
    for k in range(indlim):
        print(f"  Coefficients for s^({indlim-1-k}) (Slice k={k}):")
        # Check if slice exists, handle cases where indlim might be larger than actual degree
        if k < len(coeffs_3d_list_of_slices):
            slice_k = coeffs_3d_list_of_slices[k]
            if slice_k.shape == (p,m) : # Check if shape matches expected PxM or MxP etc.
                 for i_row in range(p):
                    row_str = "    ["
                    for j_col in range(m):
                        row_str += f"{slice_k[i_row, j_col]:8.4f}"
                        if j_col < m - 1:
                            row_str += ", "
                    row_str += "]"
                    print(row_str)
            else:
                 print(f"    Unexpected slice shape: {slice_k.shape}, expected ({p},{m})")
        else:
            print(f"    Slice k={k} not available (indlim might be > actual num_coeffs).")
    print("-" * 30)

def tc01od_manual_transpose(m_in, p_in, indlim_in, pcoeff_in_slices, qcoeff_in_slices, leri_in):
    """
    Manually performs the transpose operation of TC01OD.
    pcoeff_in_slices: list of porm x porm numpy arrays (slices of P)
    qcoeff_in_slices: list of P x M numpy arrays (slices of Q)
    Returns:
    pcoeff_out_slices: list of porm x porm numpy arrays (P')
    qcoeff_out_slices: list of M x P numpy arrays (Q')
    """
    pcoeff_out_slices = []
    for k_slice in range(indlim_in):
        if k_slice < len(pcoeff_in_slices):
            p_k = pcoeff_in_slices[k_slice]
            pcoeff_out_slices.append(p_k.T.copy()) # P'(k) = P(k)^T
        else: # Should not happen if indlim matches actual data
            pcoeff_out_slices.append(np.zeros_like(pcoeff_in_slices[0]).T)


    qcoeff_out_slices = []
    for k_slice in range(indlim_in):
        if k_slice < len(qcoeff_in_slices):
            q_k = qcoeff_in_slices[k_slice]
            qcoeff_out_slices.append(q_k.T.copy()) # Q'(k) = Q(k)^T
        else:
             qcoeff_out_slices.append(np.zeros((m_in, p_in))) # MxP zero matrix

    return pcoeff_out_slices, qcoeff_out_slices


def main():
    # Parameters from TC01OD.html example
    leri = 'L'
    m = 2
    p = 2
    indlim = 3 # Number of coefficient matrices (max_degree + 1)

    # PCOEFF input data (for LERI='L', PCOEFF is PxP)
    # Slices for K=1 (s^2), K=2 (s^1), K=3 (s^0)
    # P(s) = P0*s^2 + P1*s^1 + P2*s^0 (if indlim=3, degrees are indlim-k)
    # Fortran K=1 -> s^(indlim-1) = s^2
    # Fortran K=2 -> s^(indlim-2) = s^1
    # Fortran K=3 -> s^(indlim-3) = s^0
    pcoeff_slices_in = [
        np.array([[2.0, 4.0], [5.0, 3.0]]),    # P_coeffs for s^2 (Fortran K=1)
        np.array([[3.0, -1.0], [7.0, 2.0]]),   # P_coeffs for s^1 (Fortran K=2)
        np.array([[1.0, -1.0], [-6.0, 2.0]])   # P_coeffs for s^0 (Fortran K=3)
    ]

    # QCOEFF input data (PxM = 2x2)
    qcoeff_slices_in = [
        np.array([[6.0, 1.0], [1.0, 4.0]]),    # Q_coeffs for s^2
        np.array([[-1.0, 7.0], [1.0, 1.0]]),   # Q_coeffs for s^1
        np.array([[5.0, 5.0], [1.0, -1.0]])    # Q_coeffs for s^0
    ]

    print("--- Input Parameters ---")
    print(f"LERI: '{leri}'")
    print(f"M: {m}")
    print(f"P: {p}")
    print(f"INDLIM: {indlim}")
    
    porm_dim = p if leri == 'L' else m
    print_poly_matrix_coeffs("PCOEFF_in", pcoeff_slices_in, porm_dim, porm_dim, indlim)
    print_poly_matrix_coeffs("QCOEFF_in (PxM)", qcoeff_slices_in, p, m, indlim)
    print("=" * 50)

    # Perform the dual operation (manual transpose for TC01OD)
    pcoeff_slices_out, qcoeff_slices_out = tc01od_manual_transpose(
        m, p, indlim, pcoeff_slices_in, qcoeff_slices_in, leri
    )

    print("\n--- Output Results (Manually Transposed) ---")
    print_poly_matrix_coeffs("PCOEFF_out (P')", pcoeff_slices_out, porm_dim, porm_dim, indlim)
    # Q' is MxP
    print_poly_matrix_coeffs("QCOEFF_out (Q' - MxP)", qcoeff_slices_out, m, p, indlim)


    # --- Helper for C++ flat array generation (Row-Major Slice Order) ---
    # This is how the C++ test expects the flat arrays (concatenation of row-major slices)
    
    # PCOEFF_out (porm x porm slices)
    if pcoeff_slices_out:
        pcoeff_flat_c_rm_slice = np.concatenate([s.flatten(order='C') for s in pcoeff_slices_out]).tolist()
        print("\nPCOEFF_out flattened for C++ (Row-Major Slices Style):")
        print("{", end="")
        for i, val in enumerate(pcoeff_flat_c_rm_slice):
            print(f"{val:.4f}", end="")
            if i < len(pcoeff_flat_c_rm_slice) - 1:
                print(", ", end="")
        print("};")
        print("-" * 30)

    # QCOEFF_out (MxP slices)
    if qcoeff_slices_out:
        qcoeff_flat_c_rm_slice = np.concatenate([s.flatten(order='C') for s in qcoeff_slices_out]).tolist()
        print("\nQCOEFF_out (MxP) flattened for C++ (Row-Major Slices Style):")
        print("{", end="")
        for i, val in enumerate(qcoeff_flat_c_rm_slice):
            print(f"{val:.4f}", end="")
            if i < len(qcoeff_flat_c_rm_slice) - 1:
                print(", ", end="")
        print("};")
        print("-" * 30)


if __name__ == "__main__":
    main()
