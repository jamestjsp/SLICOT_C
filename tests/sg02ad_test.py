import numpy as np
from slycot import synthesis # Assuming slycot is installed

def get_expected_X_for_sg02ad_case1():
    """
    Runs the sg02ad case1 and prints the resulting X matrix
    with full precision for C++ test comparison.
    """
    n = 3
    m = 1
    p = 1 # P is not strictly used by Fortran when FACT='N', but slycot wrapper might use it.
          # The slycot wrapper for sg02ad passes p to the Fortran routine.

    A_np = np.array([[ 0.63399379,  0.54906824,  0.76253406],
                     [ 0.5404729 ,  0.53745766,  0.08731853],
                     [ 0.27524045,  0.84922129,  0.4681622 ]])
    E_np = np.eye(3)
    B_np = np.array([[ 0.96861695],
                     [ 0.05532739],
                     [ 0.78934047]]) # Shape (3,1)
    Q_np = np.eye(3)
    R_np = np.ones((1,1), dtype=float)
    L_np = np.zeros((n,m), dtype=float) # Shape (3,1)

    # Mode parameters from the Python test
    dico = 'D'
    jobb = 'B'
    fact = 'N'
    uplo = 'U'
    jobl = 'Z' # L is zero, so L_np is not strictly used by Fortran but passed by slycot
    scal = 'N'
    sort = 'S'
    acc  = 'R'
    tol  = 0.0

    print(f"Running sg02ad with parameters:")
    print(f"DICO='{dico}', JOBB='{jobb}', FACT='{fact}', UPLO='{uplo}', JOBL='{jobl}'")
    print(f"SCAL='{scal}', SORT='{sort}', ACC='{acc}'")
    print(f"N={n}, M={m}, P={p}, TOL={tol}\n")

    rcondu, X, alphar, alphai, beta, S_out, T_out, U_out, iwarn = \
        synthesis.sg02ad(dico, jobb, fact, uplo, jobl, scal, sort, acc,
                         n, m, p, # p is passed by slycot wrapper
                         A_np, E_np, B_np, Q_np, R_np, L_np, tol=tol)

    print(f"INFO from sg02ad: Not directly returned by slycot's sg02ad, but iwarn is: {iwarn}")
    print(f"RCONDU: {rcondu}\n")

    print("Expected X matrix (row-major for C++ vector initialization):")
    print("{\n    ", end="")
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            print(f"{X[i, j]:.16e}", end="") # Print with high precision
            if not (i == X.shape[0] - 1 and j == X.shape[1] - 1):
                print(", ", end="")
            if j == X.shape[1] - 1 and i != X.shape[0] -1 :
                print("\n    ", end="") # Newline for next row
    print("\n};")

    # For verification, the Python test checks the Riccati residual:
    # LATXB = L_np + A_np.T @ X @ B_np
    # residual = A_np.T @ X @ A_np - E_np.T @ X @ E_np - LATXB @ np.linalg.solve(R_np + B_np.T @ X @ B_np, LATXB.T) + Q_np
    # print("\nRiccati Residual (should be close to zero matrix):")
    # print(residual)

if __name__ == '__main__':
    get_expected_X_for_sg02ad_case1()

