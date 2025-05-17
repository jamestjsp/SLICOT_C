import numpy as np
import slycot.synthesis
import slycot.exceptions

def run_sb10fd_test():
    """
    Runs the SB10FD test using slycot with the example data and prints inputs and outputs.
    """
    # Define input parameters
    n = 6
    m = 5
    np_val = 5 # 'np' is a reserved keyword in numpy, so using np_val
    ncon = 2
    nmeas = 2
    gamma = 15.0
    tol = 1e-8 # As per example data (0.00000001)

    print("--- SB10FD Test using SLYCOT ---")
    print("\n--- Input Parameters ---")
    print(f"N (Order of system): {n}")
    print(f"M (Column size of B): {m}")
    print(f"NP (Row size of C): {np_val}")
    print(f"NCON (Number of control inputs): {ncon}")
    print(f"NMEAS (Number of measurements): {nmeas}")
    print(f"GAMMA: {gamma}")
    print(f"TOL: {tol}")

    # Define input matrices (row-major, as typical for NumPy)
    # These are taken from the Python test data you provided initially
    A = np.array(((-1.0,  0.0,  4.0,  5.0, -3.0, -2.0),
                  (-2.0,  4.0, -7.0, -2.0,  0.0,  3.0),
                  (-6.0,  9.0, -5.0,  0.0,  2.0, -1.0),
                  (-8.0,  4.0,  7.0, -1.0, -3.0,  0.0),
                  ( 2.0,  5.0,  8.0, -9.0,  1.0, -4.0),
                  ( 3.0, -5.0,  8.0,  0.0,  2.0, -6.0)))

    B = np.array(((-3.0, -4.0, -2.0,  1.0,  0.0),
                  ( 2.0,  0.0,  1.0, -5.0,  2.0),
                  (-5.0, -7.0,  0.0,  7.0, -2.0),
                  ( 4.0, -6.0,  1.0,  1.0, -2.0),
                  (-3.0,  9.0, -8.0,  0.0,  5.0),
                  ( 1.0, -2.0,  3.0, -6.0, -2.0)))

    C = np.array(((1.0, -1.0,  2.0, -4.0,  0.0, -3.0),
                  (-3.0,  0.0,  5.0, -1.0,  1.0,  1.0),
                  (-7.0,  5.0,  0.0, -8.0,  2.0, -2.0),
                  ( 9.0, -3.0,  4.0,  0.0,  3.0,  7.0),
                  ( 0.0,  1.0, -2.0,  1.0, -6.0, -2.0)))

    D = np.array((( 1.0, -2.0, -3.0,  0.0,  0.0),
                  ( 0.0,  4.0,  0.0,  1.0,  0.0),
                  ( 5.0, -3.0, -4.0,  0.0,  1.0),
                  ( 0.0,  1.0,  0.0,  1.0, -3.0),
                  ( 0.0,  0.0,  1.0,  7.0,  1.0)))

    print("\n--- Input Matrix A ---")
    print(A)
    print("\n--- Input Matrix B ---")
    print(B)
    print("\n--- Input Matrix C ---")
    print(C)
    print("\n--- Input Matrix D ---")
    print(D)

    # Make copies of A, B, C to check if slycot modifies them
    # (slycot's f2py wrapper with intent(in) should prevent modification of originals)
    A_orig = A.copy()
    B_orig = B.copy()
    C_orig = C.copy()

    # Call SB10FD
    # The slycot wrapper handles ldwork automatically if set to None or 0.
    # The Python test calls it with ldwork=None (implicitly by not passing it)
    # or ldwork=0 for minimum, or a specific large value.
    # We'll use the default (ldwork=None).
    try:
        # Corrected: Expect 5 return values from slycot.synthesis.sb10fd
        Ak, Bk, Ck, Dk, rcond = slycot.synthesis.sb10fd(
            n, m, np_val, ncon, nmeas, gamma, A, B, C, D, tol=tol
        )

        print("\n--- Output Results ---")
        # For a successful call in slycot, info is implicitly 0.
        # If there was an error, an exception would have been raised.
        print(f"INFO code from sb10fd: 0 (Successful SLYCOT call)")

        print("\n--- Output Controller Matrix Ak ---")
        print(Ak)
        print("\n--- Output Controller Matrix Bk ---")
        print(Bk)
        print("\n--- Output Controller Matrix Ck ---")
        print(Ck)
        print("\n--- Output Controller Matrix Dk ---")
        print(Dk)
        print("\n--- Output RCOND ---")
        print(rcond)

        # Check if original A, B, C were modified by the slycot call
        print("\n--- Check for Modification of Input Matrices (A, B, C) by slycot ---")
        if np.array_equal(A, A_orig):
            print("Input matrix A was NOT modified by slycot call.")
        else:
            print("Input matrix A WAS MODIFIED by slycot call.")
            # print("Original A:\n", A_orig)
            # print("A after call:\n", A)


        if np.array_equal(B, B_orig):
            print("Input matrix B was NOT modified by slycot call.")
        else:
            print("Input matrix B WAS MODIFIED by slycot call.")
            # print("Original B:\n", B_orig)
            # print("B after call:\n", B)

        if np.array_equal(C, C_orig):
            print("Input matrix C was NOT modified by slycot call.")
        else:
            print("Input matrix C WAS MODIFIED by slycot call.")
            # print("Original C:\n", C_orig)
            # print("C after call:\n", C)


    except slycot.exceptions.SlycotArithmeticError as e:
        print(f"\nSLYCOT Arithmetic Error (INFO code usually embedded in message): {e}")
    except slycot.exceptions.SlycotParameterError as e:
        print(f"\nSLYCOT Parameter Error (INFO code usually embedded in message): {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    run_sb10fd_test()
