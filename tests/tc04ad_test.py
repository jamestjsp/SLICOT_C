import numpy as np
from slycot import transform
from slycot.exceptions import SlycotParameterError, SlycotArithmeticError

def print_array_details_py(name, arr):
    """Helper function to print array details."""
    if arr is None:
        print(f"{name}: None")
        return
    print(f"{name} (shape: {arr.shape}, dtype: {arr.dtype}):")
    # Use a more controlled print for potentially large arrays
    with np.printoptions(precision=4, suppress=True, threshold=100):
        print(arr)
    print("-" * 30)

def print_poly_matrix_coeffs_py(name, coeffs_3d_fortran_order, dim1, dim2, kpcoef_val):
    """
    Prints polynomial matrix coefficients from a Fortran-ordered 3D numpy array.
    coeffs_3d_fortran_order: (dim1, dim2, kpcoef_val)
    """
    print(f"{name} (Fortran order: {dim1}x{dim2}x{kpcoef_val}), slices for s^(kpcoef-1-k) down to s^0:")
    if coeffs_3d_fortran_order is None or coeffs_3d_fortran_order.size == 0:
        print("  None or empty.")
        return
    for k_idx in range(kpcoef_val): # Iterate through coefficient slices
        print(f"  Coeffs for s^({kpcoef_val-1-k_idx}) (Fortran K={k_idx+1}):")
        # Check if slice exists for this k_idx
        if k_idx < coeffs_3d_fortran_order.shape[2]:
            slice_k = coeffs_3d_fortran_order[:,:,k_idx]
            for i_row in range(dim1):
                # Check if row exists for this i_row (relevant if dim1 > slice_k.shape[0])
                if i_row < slice_k.shape[0]:
                    row_str = "    ["
                    for j_col in range(dim2):
                         # Check if col exists for this j_col
                        if j_col < slice_k.shape[1]:
                            row_str += f"{slice_k[i_row, j_col]:8.4f}"
                            if j_col < dim2 - 1:
                                row_str += ", "
                        else:
                            row_str += " OOB_COL " # Out of bounds column
                    row_str += "]"
                    print(row_str)
                else:
                    print(f"    Row {i_row} OOB_ROW") # Out of bounds row
        else:
            print(f"    Slice k_idx={k_idx} OOB_SLICE_K")
    print("-" * 30)


def main():
    leri = 'L'
    m_py = 2
    p_py = 2
    index_py = np.array([2, 2], dtype=np.int32)
    
    kpcoef_py = 0
    if index_py.size > 0 :
        kpcoef_py = np.max(index_py) + 1
    else: 
        kpcoef_py = 1

    pcoeff_py = np.zeros((p_py, p_py, kpcoef_py))
    if p_py > 0: 
        pcoeff_py[0,0,0]=2.0; pcoeff_py[0,0,1]=3.0; pcoeff_py[0,0,2]=1.0; 
        pcoeff_py[0,1,0]=4.0; pcoeff_py[0,1,1]=-1.0;pcoeff_py[0,1,2]=-1.0; 
        pcoeff_py[1,0,0]=5.0; pcoeff_py[1,0,1]=7.0; pcoeff_py[1,0,2]=-6.0; 
        pcoeff_py[1,1,0]=3.0; pcoeff_py[1,1,1]=2.0; pcoeff_py[1,1,2]=2.0;  

    qcoeff_py = np.zeros((p_py, m_py, kpcoef_py))
    if p_py > 0 and m_py > 0: 
        qcoeff_py[0,0,0]=6.0; qcoeff_py[0,0,1]=-1.0;qcoeff_py[0,0,2]=5.0; 
        qcoeff_py[0,1,0]=1.0; qcoeff_py[0,1,1]=7.0; qcoeff_py[0,1,2]=5.0; 
        qcoeff_py[1,0,0]=1.0; qcoeff_py[1,0,1]=1.0; qcoeff_py[1,0,2]=1.0; 
        qcoeff_py[1,1,0]=4.0; qcoeff_py[1,1,1]=1.0; qcoeff_py[1,1,2]=-1.0;

    print("--- Input Parameters (Python/Slycot) ---")
    print(f"LERI: '{leri}'")
    print(f"M: {m_py}")
    print(f"P: {p_py}")
    print_array_details_py("INDEX", index_py)
    print(f"kpcoef calculated: {kpcoef_py}")
    print_poly_matrix_coeffs_py("PCOEFF_in", pcoeff_py, p_py, p_py, kpcoef_py)
    print_poly_matrix_coeffs_py("QCOEFF_in", qcoeff_py, p_py, m_py, kpcoef_py)
    print("=" * 50)

    try:
        if (leri == 'L' and p_py == 0) or (leri == 'R' and m_py == 0):
             print("\nSkipping SLYCOT call due to zero primary dimension for INDEX.")
             print("N_out: 0")
             print("RCOND_out: N/A (not computed)")
             print("A_out, B_out, C_out, D_out are effectively empty.")
             return


        n_out_slycot, rcond_out_slycot, A_out_slycot, B_out_slycot, C_out_slycot, D_out_slycot = \
            transform.tc04ad(m_py, p_py, index_py, 
                             pcoeff_py.copy(), 
                             qcoeff_py.copy(), 
                             leri=leri) 

        print("\n--- Output Results (from Slycot) ---")
        print(f"N_out: {n_out_slycot}")
        print(f"RCOND_out: {rcond_out_slycot:.4f}")
        print_array_details_py("A_out", A_out_slycot)
        print_array_details_py("B_out", B_out_slycot) 
        print_array_details_py("C_out", C_out_slycot) 
        print_array_details_py("D_out", D_out_slycot) 

    except (SlycotParameterError, SlycotArithmeticError) as e:
        print(f"\nSLYCOT Error: {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
