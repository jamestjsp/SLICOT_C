import numpy as np
import slycot

def run_slycot_sb02mt():
    # Inputs based on sb02mt_test.cpp SB02MTTestColMajor::SetUp()
    N_in = 3
    M_in = 2

    # JOBG_in = 'G' (this is handled by slycot wrapper not needing it as direct input)
    JOBL_in = 'Z'
    FACT_in = 'N' # Slycot sb02mt only supports 'N' or 'C' effectively. 'U' is not directly exposed.
    UPLO_in = 'U'

    # Matrices (ensure they are Fortran contiguous for slycot if passed directly,
    # or let slycot handle conversion from C-order numpy arrays)
    # Python/NumPy default is row-major (C order)
    A_py = np.array([[0.5, 0.0, 0.1],
                     [0.2, 0.5, 0.1],
                     [0.0, 0.0, 0.8]], order='C') # Transposed for row-major definition
    A_py = np.asfortranarray(A_py.T) # Convert to column-major as typically used in Fortran examples

    B_py = np.array([[1.0, 0.0],
                     [0.0, 1.0],
                     [0.5, 0.0]], order='C') # N_in x M_in
    B_py = np.asfortranarray(B_py)


    Q_py = np.array([[1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0],
                     [0.0, 0.0, 1.0]], order='C')
    Q_py = np.asfortranarray(Q_py.T)

    R_py = np.array([[2.0, 0.0],
                     [0.0, 2.0]], order='C')
    R_py = np.asfortranarray(R_py.T) # R is M_in x M_in

    # L is zero as JOBL_in = 'Z'. Slycot handles optional args.
    L_py = None
    if JOBL_in == 'N': # For completeness, though not this test case
        L_py = np.zeros((N_in, M_in), order='F')


    print(f"N_in: {N_in}, M_in: {M_in}")
    print(f"JOBL_in: {JOBL_in}, FACT_in: {FACT_in}, UPLO_in: {UPLO_in}")
    print("\nInput A (Fortran order):\n", A_py)
    print("\nInput B (Fortran order):\n", B_py)
    print("\nInput Q (Fortran order):\n", Q_py)
    print("\nInput R (Fortran order):\n", R_py)
    if L_py is not None:
        print("\nInput L (Fortran order):\n", L_py)

    # Call slycot.sb02mt
    # The slycot sb02mt wrapper determines JOBG based on whether G is requested as output implicitly.
    # It seems the direct _wrapper.sb02mt_x calls are used internally.
    # We want the behavior matching JOBG='G'. The Python function returns G.
    try:
        A_b, B_b, Q_b, R_b, L_b, ipiv, oufact, G = \
            slycot.sb02mt(n=N_in, m=M_in, B=B_py, R=R_py, A=A_py, Q=Q_py, L=L_py,
                          fact=FACT_in, jobl=JOBL_in, uplo=UPLO_in)

        print("\n--- Outputs from slycot.sb02mt ---")
        if A_b is not None:
            print("\nA_b (Fortran order):\n", A_b)
        else:
            print("\nA_b: None (as JOBL='Z')")

        print("\nB_b (Fortran order):\n", B_b)

        if Q_b is not None:
            print("\nQ_b (Fortran order):\n", Q_b)
        else:
            print("\nQ_b: None (as JOBL='Z')")

        print("\nR_b (Fortran order - Cholesky factor if oufact=1):\n", R_b)

        if L_b is not None:
            print("\nL_b (Fortran order):\n", L_b)
        else:
            print("\nL_b: None (as JOBL='Z')")

        if ipiv is not None:
            print("\nIPIV:\n", ipiv)
        else:
            print("\nIPIV: None")

        print("\nOUFACT:\n", oufact)
        print("\nG (Fortran order):\n", G)

        print("\n--- For C++ Test G_expected (Column Major Flattened) ---")
        # Ensure G is Fortran contiguous before flattening for column-major
        if G is not None:
            g_f_cont = np.asfortranarray(G)
            print("{", ", ".join(map(str, g_f_cont.flatten(order='F'))), "};")

        print("\n--- For C++ Test A_expected (if modified, Column Major Flattened) ---")
        if A_b is not None: # A_b is A_io modified
            a_b_f_cont = np.asfortranarray(A_b)
            print("{", ", ".join(map(str, a_b_f_cont.flatten(order='F'))), "};")
        elif A_py is not None and JOBL_in == 'Z': # A is not modified by L terms
             a_py_f_cont = np.asfortranarray(A_py)
             printComment = "// A should be unchanged as JOBL='Z'\n"
             print(printComment + "{", ", ".join(map(str, a_py_f_cont.flatten(order='F'))), "};")


        print("\n--- For C++ Test Q_expected (if modified, Column Major Flattened) ---")
        if Q_b is not None: # Q_b is Q_io modified
            q_b_f_cont = np.asfortranarray(Q_b)
            print("{", ", ".join(map(str, q_b_f_cont.flatten(order='F'))), "};")
        elif Q_py is not None and JOBL_in == 'Z': # Q is not modified by L terms
             q_py_f_cont = np.asfortranarray(Q_py)
             printComment = "// Q should be unchanged as JOBL='Z'\n"
             print(printComment + "{", ", ".join(map(str, q_py_f_cont.flatten(order='F'))), "};")

        print(f"\nEXPECTED_OUFACT = {oufact};")


    except Exception as e:
        print("\nError calling slycot.sb02mt:", e)

if __name__ == '__main__':
    run_slycot_sb02mt()