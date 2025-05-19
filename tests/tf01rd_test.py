import numpy as np
import slycot.transform

def print_matrix_f(name, matrix, rows, cols):
    """Helper to print matrices in Fortran column-major style (transposed)."""
    print(f"{name} ({rows}x{cols}, as stored Fortran column-major, printed row-wise):")
    if matrix is None or matrix.size == 0:
        print("  None or Empty")
    else:
        # Print row by row, but it's column data
        temp_matrix = matrix.reshape((cols, rows)).T # Reshape to (rows, cols) for printing
        with np.printoptions(precision=4, suppress=True, linewidth=120):
            print(temp_matrix)
    print("-" * 20)

def print_slycot_h(name, h_matrix, nc, nb, N):
    print(f"{name} ({nc}x{N*nb}) from Slycot (effectively row-major of concatenated M_k):")
    if h_matrix is None or h_matrix.size == 0:
        print("  None or Empty")
    else:
        with np.printoptions(precision=4, suppress=True, linewidth=120):
            print(h_matrix) # Slycot returns H as p x (N*m)
    print("Decomposed Markov Parameters M(k):")
    for k in range(N):
        mk = h_matrix[:, k*nb:(k+1)*nb]
        print(f"M({k+1}):\n{mk}")
    print("-" * 20)


# Data from TF01RD.html "Program Data"
na = 3
nb = 2
nc = 2
N_in = 5

# A (Fortran column-major order, then reshaped to 3x3 for Slycot)
# SLICOT example reads: ( ( A(I,J), I = 1,NA ), J = 1,NA )
#   0.000 -0.070  0.015
#   1.000  0.800 -0.150
#   0.000  0.000  0.500
a_f_flat = np.array([0.0, 1.0, 0.0, -0.070, 0.800, 0.000, 0.015, -0.150, 0.500])
a_slycot = a_f_flat.reshape((na, na), order='F')
# print_matrix_f("A_slycot (input to slycot.tf01rd)", a_slycot, na, na)


# B (Fortran column-major order, then reshaped to NAxNB for Slycot)
# SLICOT example reads: ( ( B(I,J), I = 1,NA ), J = 1,NB )
#   0.000  2.000  1.000  (col 1 of B)
#  -1.000 -0.100  1.000  (col 2 of B)
b_f_flat = np.array([0.0, 2.0, 1.0, -1.0, -0.100, 1.000])
b_slycot = b_f_flat.reshape((na, nb), order='F')
# print_matrix_f("B_slycot (input to slycot.tf01rd)", b_slycot, na, nb)


# C (Fortran column-major order, then reshaped to NCxNA for Slycot)
# SLICOT example reads: ( ( C(I,J), I = 1,NC ), J = 1,NA )
#   0.000  1.000  (col 1 of C)
#  -1.000  0.000  (col 2 of C)
#   0.000  0.000  (col 3 of C)
c_f_flat = np.array([0.0, 1.0, -1.0, 0.0, 0.0, 0.0])
c_slycot = c_f_flat.reshape((nc, na), order='F')
# print_matrix_f("C_slycot (input to slycot.tf01rd)", c_slycot, nc, na)

# Call Slycot's tf01rd
# ldwork >= max(1, 2*na*nc) = max(1, 2*3*2) = max(1,12) = 12
try:
    H_slycot = slycot.transform.tf01rd(na, nb, nc, N_in, a_slycot, b_slycot, c_slycot, ldwork=12)
    print_slycot_h("H_slycot (output)", H_slycot, nc, nb, N_in)

    # Convert H_slycot (which is NC x N*NB, row-major like) to Fortran column-major flat string for C++ test
    # H_fortran_flat will be [M(1)(:,0), M(1)(:,1), M(2)(:,0), M(2)(:,1), ...]
    H_fortran_flat_list = []
    for k in range(N_in): # Iterate through M_k
        for j in range(nb): # Iterate through columns of M_k
            for i in range(nc): # Iterate through rows of M_k
                H_fortran_flat_list.append(H_slycot[i, k*nb + j])
    
    print("H_expected_f for C++ test (Fortran column-major flattened):")
    print("{", end="")
    for i, val in enumerate(H_fortran_flat_list):
        print(f"{val:.4f}", end="")
        if i < len(H_fortran_flat_list) - 1:
            print(", ", end="")
    print("};")

except Exception as e:
    print(f"Error calling slycot.tf01rd: {e}")

