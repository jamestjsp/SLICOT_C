import numpy as np
import slycot.transform

def print_matrix(name, matrix):
    """Helper to print matrices."""
    print(f"{name}:")
    if matrix is None or matrix.size == 0:
        print("  None or Empty")
    else:
        # Limit precision for cleaner output, matching Fortran example somewhat
        with np.printoptions(precision=4, suppress=True, linewidth=120):
            print(matrix)
    print("-" * 20)

def run_td04ad_slycot(rowcol, m, p, index, dcoeff, ucoeff, tol=0.0, test_name="Test Case"):
    """
    Runs slycot.transform.td04ad and prints inputs and outputs.
    Assumes dcoeff and ucoeff are already in the correct shape for slycot.
    """
    print(f"\n======= {test_name} (Using Slycot) =======")
    print("INPUTS:")
    print(f"ROWCOL: '{rowcol}'")
    print(f"M: {m}")
    print(f"P: {p}")
    print_matrix("INDEX", np.array(index))
    print_matrix("DCOEFF", dcoeff)
    print_matrix("UCOEFF", ucoeff)
    print(f"TOL: {tol}")
    print("----------------------------------------")

    try:
        # Slycot's td04ad might return fewer elements if info is non-zero
        # The Fortran routine TD04AD has INFO as the last output of the F77 call.
        # Slycot's td04ad wrapper returns nr, A, B, C, D. It raises an exception on error.
        # We'll try to mimic how one might call it and catch errors.
        
        # Note: Slycot's td04ad does not directly return the 'info' code from the Fortran routine
        # in the same way. It raises exceptions for errors.
        # We are primarily interested in the state-space matrices for successful calls.

        nr, A, B, C, D = slycot.transform.td04ad(rowcol, m, p, np.array(index, dtype=np.int32), dcoeff, ucoeff, tol=tol)
        slycot_info = 0 # Assume 0 if no exception
    except Exception as e:
        print(f"Slycot raised an exception: {e}")
        nr, A, B, C, D = -1, None, None, None, None # Or some other error indication
        slycot_info = -999 # Placeholder for error

    print("\nOUTPUTS (from Slycot):")
    print(f"NR: {nr}")
    print_matrix("A", A)
    print_matrix("B", B)
    print_matrix("C", C)
    print_matrix("D", D)
    print(f"Slycot Info (derived): {slycot_info}") # Not the Fortran INFO
    print("========================================")

    return nr, A, B, C, D, slycot_info

if __name__ == "__main__":
    # 1. Data from Td04adDocExampleTest
    # ROWCOL='R', M=2, P=2, INDEX={3,3}, TOL=0.0
    # DCOEFF (PxK = 2x4) for slycot (row-major numpy array)
    # Original Fortran example data:
    # D(1,:)=[1,6,11,6]
    # D(2,:)=[1,6,11,6]
    dcoeff_doc = np.array([
        [1.0, 6.0, 11.0, 6.0],
        [1.0, 6.0, 11.0, 6.0]
    ])

    # UCOEFF (PxMxK = 2x2x4) for slycot (row-major numpy array)
    # U(i,j,k)
    # U(0,0,:) = [1,6,12,7]
    # U(0,1,:) = [0,1,4,3]
    # U(1,0,:) = [0,0,1,1]
    # U(1,1,:) = [1,8,20,15]
    ucoeff_doc = np.array([
        [ # P=0
            [1.0, 6.0, 12.0, 7.0], # M=0
            [0.0, 1.0,  4.0, 3.0]  # M=1
        ],
        [ # P=1
            [0.0, 0.0,  1.0, 1.0], # M=0
            [1.0, 8.0, 20.0, 15.0] # M=1
        ]
    ])
    run_td04ad_slycot('R', 2, 2, [3, 3], dcoeff_doc, ucoeff_doc, test_name="Doc Example")

    # 2. Data from Td04adStaticGainTest
    # ROWCOL='R', M=1, P=1, INDEX={0}, TOL=0.0
    # DCOEFF (PxK = 1x1)
    dcoeff_static = np.array([[1.0]])
    # UCOEFF (PxMxK = 1x1x1)
    ucoeff_static = np.array([[[2.0]]])
    run_td04ad_slycot('R', 1, 1, [0], dcoeff_static, ucoeff_static, test_name="Static Gain Example")

    # 3. Data from Td04adZeroDimTest
    # ROWCOL='R', M=0, P=0, INDEX={}, TOL=0.0
    # For M=0 or P=0, INDEX, DCOEFF, UCOEFF handling can be tricky.
    # If M=0, PORM for ROWCOL='C' is 0. If P=0, PORM for ROWCOL='R' is 0.
    # Slycot might expect empty arrays or handle it internally.
    
    # Case 3a: P=0, M=any (ROWCOL='R') -> PORM = 0
    # INDEX should be empty for PORM=0. DCOEFF for P=0 rows. UCOEFF for P=0, M=any.
    index_zero_p = []
    dcoeff_zero_p = np.empty((0,1)) # P=0, K=1 (from index max 0 + 1)
    ucoeff_zero_p = np.empty((0,1,1))# P=0, M=1, K=1
    run_td04ad_slycot('R', 1, 0, index_zero_p, dcoeff_zero_p, ucoeff_zero_p, test_name="Zero Dim Example (P=0, M=1, ROWCOL='R')")

    # Case 3b: M=0, P=any (ROWCOL='R')
    # PORM = P. DCOEFF for P rows. UCOEFF for P, M=0.
    index_zero_m_r = [0] # P=1, so index has 1 element
    dcoeff_zero_m_r = np.array([[1.0]]) # P=1, K=1
    ucoeff_zero_m_r = np.empty((1,0,1)) # P=1, M=0, K=1
    run_td04ad_slycot('R', 0, 1, index_zero_m_r, dcoeff_zero_m_r, ucoeff_zero_m_r, test_name="Zero Dim Example (P=1, M=0, ROWCOL='R')")
    
    # Case 3c: M=0, P=0 (ROWCOL='R') -> PORM = 0
    index_zero_mp = []
    dcoeff_zero_mp = np.empty((0,1))
    ucoeff_zero_mp = np.empty((0,0,1))
    run_td04ad_slycot('R', 0, 0, index_zero_mp, dcoeff_zero_mp, ucoeff_zero_mp, test_name="Zero Dim Example (P=0, M=0, ROWCOL='R')")

    # Case 3d: M=0, P=any (ROWCOL='C') -> PORM = 0
    index_zero_m_c = []
    dcoeff_zero_m_c = np.empty((0,1)) # M=0, K=1
    # Slycot expects ucoeff for ROWCOL='C' to be (max(1,M,P), max(1,M,P), K)
    # If M=0, P=1, then max(1,0,1)=1. So ucoeff_f_shape = (1,1,1)
    # The actual PxM part is 1x0, which is empty.
    ucoeff_zero_m_c_slycot_shape = np.empty((1,1,1)) # P=1, M=0, K=1 (this is for slycot's internal needs)
    run_td04ad_slycot('C', 0, 1, index_zero_m_c, dcoeff_zero_m_c, ucoeff_zero_m_c_slycot_shape, test_name="Zero Dim Example (P=1, M=0, ROWCOL='C')")

