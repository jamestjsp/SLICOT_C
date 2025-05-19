import numpy as np
import slycot.transform

def print_matrix_details(name, matrix_data, shape, order='F'):
    """Prints matrix details and its representation."""
    print(f"--- {name} ---")
    print(f"Shape: {shape}")
    print(f"Input flat data (Fortran column-major order): {matrix_data.tolist()}")
    matrix = matrix_data.reshape(shape, order=order)
    print(f"Matrix as used by Slycot (numpy default print, effectively row-major view of {shape} matrix):")
    print(matrix)
    print("-" * 20)

# Parameters from TF01RD.html "Program Data"
N_param = 5
NA_param = 3
NB_param = 2
NC_param = 2

print("=== Input Parameters ===")
print(f"N : {N_param}")
print(f"NA: {NA_param}")
print(f"NB: {NB_param}")
print(f"NC: {NC_param}\n")

# Matrix A: (NA x NA) = 3x3
# Fortran reads: ( ( A(I,J), I = 1,NA ), J = 1,NA )
# Data lines from HTML:
#   0.000 -0.070  0.015
#   1.000  0.800 -0.150
#   0.000  0.000  0.500
# This implies flat column-major order:
A_flat_f = np.array([0.000, 1.000, 0.000, -0.070, 0.800, 0.000, 0.015, -0.150, 0.500])
A_slycot = A_flat_f.reshape((NA_param, NA_param), order='F')
print_matrix_details("Matrix A", A_flat_f, (NA_param, NA_param))

# Matrix B: (NA x NB) = 3x2
# Fortran reads: ( ( B(I,J), I = 1,NA ), J = 1,NB )
# Data lines from HTML:
#   0.000  2.000  1.000  (Column 1 of B)
#  -1.000 -0.100  1.000  (Column 2 of B)
# This implies flat column-major order:
B_flat_f = np.array([0.000, 2.000, 1.000, -1.000, -0.100, 1.000])
B_slycot = B_flat_f.reshape((NA_param, NB_param), order='F')
print_matrix_details("Matrix B", B_flat_f, (NA_param, NB_param))

# Matrix C: (NC x NA) = 2x3
# Fortran reads: ( ( C(I,J), I = 1,NC ), J = 1,NA )
# Data lines from HTML (this part is tricky, assuming numbers are read sequentially to fill columns):
#   0.000  1.000  0.000
#   0.000  1.000  0.000
# C(1,1)=0.000, C(2,1)=0.000 (from first two numbers of "0.000  1.000  0.000")
# C(1,2)=1.000, C(2,2)=1.000 (from last of first line, first of second line if numbers are split, or next two if all on one line)
# The example data file for TF01RD.f has C as:
# 0.000  1.000
# -1.000  0.000
# 0.000  0.000
# This corresponds to C = [[0, -1, 0], [1, 0, 0]] (standard display)
# Flat column major: 0.0, 1.0,  -1.0, 0.0,  0.0, 0.0
C_flat_f = np.array([0.0, 1.0, -1.0, 0.0, 0.0, 0.0])
C_slycot = C_flat_f.reshape((NC_param, NA_param), order='F')
print_matrix_details("Matrix C", C_flat_f, (NC_param, NA_param))

# Call Slycot's tf01rd
# ldwork >= max(1, 2*na*nc) = max(1, 2*3*2) = 12
try:
    H_slycot_output = slycot.transform.tf01rd(NA_param, NB_param, NC_param, N_param, 
                                              A_slycot, B_slycot, C_slycot, ldwork=12)
    
    print("\n=== Output H Matrix (from Slycot) ===")
    print(f"Shape: ({H_slycot_output.shape[0]} x {H_slycot_output.shape[1]})")
    print("H (Slycot returns it as NC x (N*NB), effectively row-major of concatenated M_k):")
    with np.printoptions(precision=4, suppress=True, linewidth=120):
        print(H_slycot_output)
    
    print("\nDecomposed Markov Parameters M(k) from H_slycot_output:")
    for k_idx in range(N_param):
        mk = H_slycot_output[:, k_idx*NB_param:(k_idx+1)*NB_param]
        print(f"M({k_idx+1}):\n{mk}")

    # Convert H_slycot_output to the flattened Fortran column-major order
    # for the C++ test's H_expected_f vector.
    # H_expected_f stores [M(1)(col1), M(1)(col2), ..., M(N)(col1), M(N)(col2)]
    # where each M(k)(col_j) is a column vector of length NC.
    H_fortran_flat_list = []
    for k_loop in range(N_param):      # For each Markov parameter M(k)
        for j_loop in range(NB_param): # For each column in M(k)
            for i_loop in range(NC_param): # For each row in M(k)
                H_fortran_flat_list.append(H_slycot_output[i_loop, k_loop*NB_param + j_loop])
    
    print("\n-------------------------------------------------------------")
    print("H_expected_f for C++ test (Fortran column-major flattened):")
    print("-------------------------------------------------------------")
    print("{", end="")
    for i, val in enumerate(H_fortran_flat_list):
        print(f"{val:.4f}", end="")
        if i < len(H_fortran_flat_list) - 1:
            print(", ", end="")
        if (i + 1) % (NC_param * NB_param) == 0 and i < len(H_fortran_flat_list) - 1 : # Newline after each M(k) block
             print(f", // M({(i+1)//(NC_param*NB_param)})")
             print(" ", end="")
        elif (i + 1) % NC_param == 0 and NB_param > 1 and (i+1) % (NC_param * NB_param) != 0  and i < len(H_fortran_flat_list) - 1 : # Space between columns of M(k)
             print(", ", end="")


    print("};")
    print("-------------------------------------------------------------")

except Exception as e:
    print(f"\nError calling slycot.tf01rd: {e}")

