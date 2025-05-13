import numpy as np

# Attempt to import slycot, handle import errors
try:
    import slycot
    slycot_found = True
except ImportError:
    slycot_found = False
    print("WARNING: Slycot library not found. Install with 'pip install slycot'")

# System dimensions from the AB09MD example
N = 7  # Number of states
M = 2  # Number of inputs
P = 3  # Number of outputs
ALPHA = -0.6  # Stability boundary
TOL = 0.1  # Tolerance

# Input matrices from the HTML documentation example
A = np.array([
    [-0.04165,  0.0000,  4.9200, -4.9200,  0.0000,  0.0000,  0.0000],
    [-5.2100, -12.500,   0.0000,  0.0000,  0.0000,  0.0000,  0.0000],
    [ 0.0000,   3.3300, -3.3300,  0.0000,  0.0000,  0.0000,  0.0000],
    [ 0.5450,   0.0000,  0.0000,  0.0000, -0.5450,  0.0000,  0.0000],
    [ 0.0000,   0.0000,  0.0000,  4.9200, -0.04165, 0.0000,  4.9200],
    [ 0.0000,   0.0000,  0.0000,  0.0000, -5.2100, -12.500,  0.0000],
    [ 0.0000,   0.0000,  0.0000,  0.0000,  0.0000,  3.3300, -3.3300]
])

B = np.array([
    [ 0.0000,   0.0000],
    [12.500,    0.0000],
    [ 0.0000,   0.0000],
    [ 0.0000,   0.0000],
    [ 0.0000,   0.0000],
    [ 0.0000,  12.500],
    [ 0.0000,   0.0000]
])

C = np.array([
    [1.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000],
    [0.0000,  0.0000,  0.0000,  1.0000,  0.0000,  0.0000,  0.0000],
    [0.0000,  0.0000,  0.0000,  0.0000,  1.0000,  0.0000,  0.0000]
])

# Expected results from SLICOT HTML documentation
expected_NR = 5
expected_HSV = np.array([1.9178, 0.8621, 0.7666, 0.0336, 0.0246])

expected_AR = np.array([
    [-0.5181, -1.1084,  0.0000,  0.0000,  0.0000],
    [ 8.8157, -0.5181,  0.0000,  0.0000,  0.0000],
    [ 0.0000,  0.0000,  0.5124,  0.0000,  1.7910],
    [ 0.0000,  0.0000,  0.0000, -1.4460,  0.0000],
    [ 0.0000,  0.0000, -4.2167,  0.0000, -2.9900]
])

expected_BR = np.array([
    [-1.2837,  1.2837],
    [-0.7522,  0.7522],
    [-0.7447, -0.7447],
    [ 1.9275, -1.9275],
    [-3.6872, -3.6872]
])

expected_CR = np.array([
    [-0.1380, -0.6445, -0.6582, -0.5771,  0.2222],
    [ 0.6246,  0.0196,  0.0000,  0.4131,  0.0000],
    [ 0.1380,  0.6445, -0.6582,  0.5771,  0.2222]
])

print("AB09MD Validation using Slycot")
print("==============================")
print("\nInput Parameters:")
print(f"N={N}, M={M}, P={P}, ALPHA={ALPHA}, TOL={TOL}")
print("DICO='C', JOB='N', EQUIL='N', ORDSEL='A'")

# Run Slycot and compare results
if slycot_found:
    try:
        # Call Slycot ab09md function with the same parameters as our C++ test
        nr, ar, br, cr, ns, hsv = slycot.ab09md(
            dico='C', job='N', equil='N', 
            n=N, m=M, p=P,
            A=A, B=B, C=C, 
            alpha=ALPHA, tol=TOL
        )
        
        print("\n--- Results from Slycot ---")
        print(f"Reduced order (nr): {nr}")
        print(f"Number of stable eigenvalues (ns): {ns}")
        print(f"Hankel singular values: {hsv[:ns]}")
        print("\nReduced Ar matrix:")
        print(np.round(ar, 4))
        print("\nReduced Br matrix:")
        print(np.round(br, 4))
        print("\nReduced Cr matrix:")
        print(np.round(cr, 4))
        
        print("\n--- Comparing with expected results from SLICOT docs ---")
        print(f"Order matches: {nr == expected_NR}")
        
        # Use a reasonable tolerance for numerical comparisons
        hsv_match = np.allclose(hsv[:ns], expected_HSV, atol=1e-2)
        ar_match = np.allclose(ar, expected_AR, atol=1e-2)
        br_match = np.allclose(br, expected_BR, atol=1e-2)
        cr_match = np.allclose(cr, expected_CR, atol=1e-2)
        
        print(f"HSV close to expected: {hsv_match}")
        print(f"Ar close to expected: {ar_match}")
        print(f"Br close to expected: {br_match}")
        print(f"Cr close to expected: {cr_match}")
        
        # Compare our C++ test results with Slycot results
        print("\n--- Comparing C++ test values with Slycot results ---")
        # Values observed from C++ test failure
        cpp_hsv = np.array([1.669422400872794, 0.48155211383813412, 0.36689329366395884, 0.06805670480642699, 0.00076658086133285193])
        
        print(f"HSV (C++ test vs Slycot):")
        for i in range(min(len(cpp_hsv), len(hsv[:ns]))):
            print(f"  HSV[{i}]: C++={cpp_hsv[i]:.6f}, Slycot={hsv[i]:.6f}, Diff={cpp_hsv[i]-hsv[i]:.6f}")
        
        # If results don't match, show detailed differences
        if not hsv_match or not ar_match or not br_match or not cr_match:
            print("\n--- Detailed differences between Slycot and SLICOT docs ---")
            
            if not hsv_match:
                print("\nHSV differences:")
                for i in range(min(len(hsv[:ns]), len(expected_HSV))):
                    print(f"HSV[{i}]: Slycot={hsv[i]:.6f}, Expected={expected_HSV[i]:.6f}, Diff={hsv[i]-expected_HSV[i]:.6f}")
        
        print("\n--- Recommendation for C++ Test ---")
        if np.allclose(cpp_hsv, hsv[:ns], atol=1e-2):
            print("C++ test results match Slycot. Use the Slycot results as expected values with tolerance=1e-2.")
        else:
            print("C++ test results differ from Slycot. Check C++ implementation for issues.")
        
    except Exception as e:
        print(f"\nError during Slycot ab09md call: {e}")
else:
    print("\nSkipping Slycot test as library is not available.")
