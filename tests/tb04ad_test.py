import numpy as np
from slycot import transform
from slycot.exceptions import SlycotParameterError # For error handling

def print_array_details(name, arr):
    """Helper function to print array details."""
    if arr is None:
        print(f"{name}: None")
        return
    print(f"{name} (shape: {arr.shape}, dtype: {arr.dtype}):")
    print(arr)
    print("-" * 30)

def main():
    # Parameters from TB04AD.html example
    # ROWCOL = 'R' # This parameter is implicit in the slycot.transform.tb04ad wrapper,
                   # likely defaulting to 'R' as per the example.
    n = 3        # Order of the state-space representation
    m = 2        # Number of system inputs
    p = 2        # Number of system outputs
    tol1 = 0.0   # Tolerance for determining transfer function coefficients
    tol2 = 0.0   # Tolerance for separating controllable/observable subsystem

    # Input matrices (as NumPy arrays)
    # A_in: State dynamics matrix (N x N)
    A_in = np.array([
        [-1.0,  0.0,  0.0],
        [ 0.0, -2.0,  0.0],
        [ 0.0,  0.0, -3.0]
    ])

    # B_in: Input matrix (N x M)
    B_in = np.array([
        [ 0.0,  1.0],
        [ 1.0,  1.0],
        [-1.0,  0.0]
    ])

    # C_in: Output matrix (P x N)
    C_in = np.array([
        [0.0, 1.0, 1.0],
        [1.0, 0.0, 1.0]
    ])

    # D_in: Direct transmission matrix (P x M)
    D_in = np.array([
        [1.0, 0.0],
        [0.0, 1.0]
    ])

    print("--- Input Parameters (Assuming ROWCOL='R' for slycot.transform.tb04ad) ---")
    print(f"N: {n}")
    print(f"M: {m}")
    print(f"P: {p}")
    print_array_details("A_in", A_in)
    print_array_details("B_in", B_in)
    print_array_details("C_in", C_in)
    print_array_details("D_in", D_in)
    print(f"TOL1: {tol1}")
    print(f"TOL2: {tol2}")
    print("=" * 50)

    try:
        # Call tb04ad from slycot.transform
        # The slycot wrapper tb04ad returns:
        # Ar, Br, Cr, nr_out, index_out, dcoeff_out, ucoeff_out
        # Note: The slycot wrapper for tb04ad handles the internal workspace.
        # The input matrices A, B, C are modified by the Fortran routine.
        # We pass copies to preserve the original inputs for printing.
        A_f = A_in.copy()
        B_f = B_in.copy()
        C_f = C_in.copy()
        # D_in is not modified by Fortran if ROWCOL='R'.
        
        # Corrected call: removed the 'rowcol' keyword argument.
        Ar_out, Br_out, Cr_out, nr_out, index_out, dcoeff_out, ucoeff_out = \
            transform.tb04ad(n, m, p, A_f, B_f, C_f, D_in.copy(), tol1=tol1, tol2=tol2)

        print("\n--- Output Results ---")
        print(f"NR_out: {nr_out}")
        print_array_details("Ar_out (transformed A)", Ar_out)
        print_array_details("Br_out (transformed B)", Br_out)
        print_array_details("Cr_out (transformed C)", Cr_out)
        print_array_details("INDEX_out", index_out)
        
        # DCOEFF_out is (porm, kmax_coeffs)
        # For ROWCOL='R', porm = p
        print_array_details("DCOEFF_out", dcoeff_out)
        
        # UCOEFF_out is (porm, porp, kmax_coeffs)
        # For ROWCOL='R', porm = p, porp = m
        print_array_details("UCOEFF_out", ucoeff_out)

        # --- Helper for C++ flat array generation (Column-Major) ---
        # This demonstrates how to flatten for C++ assuming Fortran's column-major order.
        
        if dcoeff_out is not None and dcoeff_out.size > 0:
            porm_f_val = dcoeff_out.shape[0]
            kmax_coeffs = dcoeff_out.shape[1]
            dcoeff_flat_cm = []
            # Fortran DCOEFF(PORM, KMAX_COEFFS)
            # C++ flat: D(0,0), D(1,0)...D(PORM-1,0), D(0,1), D(1,1)...
            for k_coeff_idx in range(kmax_coeffs): # Iterate over columns (coefficient index)
                for i_porm_idx in range(porm_f_val):    # Iterate over rows (porm index)
                    dcoeff_flat_cm.append(dcoeff_out[i_porm_idx, k_coeff_idx])
            print("\nDCOEFF_out flattened for C++ (Column-Major Style):")
            print("{", end="")
            for i, val in enumerate(dcoeff_flat_cm):
                print(f"{val:.4f}", end="")
                if i < len(dcoeff_flat_cm) - 1:
                    print(", ", end="")
            print("};")
            print("-" * 30)

        if ucoeff_out is not None and ucoeff_out.size > 0:
            porm_f_val = ucoeff_out.shape[0]
            porp_f_val = ucoeff_out.shape[1]
            kmax_coeffs = ucoeff_out.shape[2]
            ucoeff_flat_cm = []
            # Fortran UCOEFF(PORM, PORP, KMAX_COEFFS)
            # C++ flat: U(0,0,0), U(1,0,0)...U(PORM-1,0,0), U(0,1,0), U(1,1,0)... then next k_coeff_idx
            for k_coeff_idx in range(kmax_coeffs):      # Outermost loop for coefficient index
                for j_porp_idx in range(porp_f_val):    # Middle loop for PORP dimension
                    for i_porm_idx in range(porm_f_val):# Innermost loop for PORM dimension
                        ucoeff_flat_cm.append(ucoeff_out[i_porm_idx, j_porp_idx, k_coeff_idx])
            
            print("\nUCOEFF_out flattened for C++ (Column-Major Style - PORM varies fastest, then PORP, then K_COEFF_IDX):")
            print("{", end="")
            for i, val in enumerate(ucoeff_flat_cm):
                print(f"{val:.4f}", end="")
                if i < len(ucoeff_flat_cm) - 1:
                    print(", ", end="")
            print("};")
            print("-" * 30)

    except SlycotParameterError as e:
        print(f"\nSlycotParameterError: {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
