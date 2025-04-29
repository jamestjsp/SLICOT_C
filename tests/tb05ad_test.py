import numpy as np
import scipy.linalg # Using SciPy for solve is generally recommended

# --- Define Matrices (using NumPy arrays) ---
A = np.array([
    [1.0, 2.0, 0.0],
    [4.0, -1.0, 0.0],
    [0.0, 0.0, 1.0]
])

B = np.array([
    [1.0],
    [0.0],
    [1.0]
])

C = np.array([
    [1.0, 0.0, -1.0],
    [0.0, 0.0, 1.0]
])

D = np.array([
    [0.0],
    [0.0]
])

# --- 1. Calculate Eigenvalues of A ---
eigenvalues_A = np.linalg.eigvals(A)

print("Eigenvalues of A from np.linalg.eigvals:")
# Sort them for easier comparison (optional)
print(np.sort(eigenvalues_A))
print("-" * 30)
print("Fortran example Eigenvalues:")
print("   3.00 + 0.00j")
print("  -3.00 + 0.00j")
print("   1.00 + 0.00j")
print("-" * 30)

# --- 2. Calculate H_inv_B ---
# Define the complex frequency
s = complex(0.0, 0.5)
N = A.shape[0]

# Construct H = sI - A
H = s * np.identity(N) - A

# Solve HX = B for X = H_inv_B
try:
    H_inv_B = scipy.linalg.solve(H, B)
    print("H_inv_B = (sI - A)^-1 * B from scipy.linalg.solve (Python):")
    print(f"[[{H_inv_B[0,0]: .6f}{H_inv_B[0,0].imag:+.6f}j]")
    print(f" [{H_inv_B[1,0]: .6f}{H_inv_B[1,0].imag:+.6f}j]")
    print(f" [{H_inv_B[2,0]: .6f}{H_inv_B[2,0].imag:+.6f}j]]")


    # Optional Check: Recalculate G = C*X + D
    G_check = C @ H_inv_B + D
    print("\nCheck: G = C @ H_inv_B + D:")
    print(G_check)


except (np.linalg.LinAlgError, scipy.linalg.LinAlgError):
     print(f"Matrix H=sI-A is singular at s={s}. Cannot solve for H_inv_B.")
except ValueError as e:
     print(f"Error during solve for H_inv_B: {e}")

print("-" * 30)
print("Fortran example H(inverse)*B:")
print(" (-0.11, -0.05j)")
print(" (-0.43,  0.00j)")
print(" (-0.80, -0.40j)")
print("-" * 30)