import numpy as np
import matplotlib.pyplot as plt

# Attempt to import control and harold, handle import errors
try:
    import control
    control_found = True
except ImportError:
    control_found = False
    print("WARNING: python-control library not found. Install with 'pip install control'")

try:
    from harold import State, simulate_linear_system
    harold_found = True
except ImportError:
    harold_found = False
    print("WARNING: harold library not found. Install with 'pip install harold'")


# System Dimensions
N = 3  # Number of states
M = 2  # Number of inputs
P = 2  # Number of outputs
NY = 10 # Number of time steps

# State-Space Matrices (from the C++ test fixture, matching TF01MD example)
# A (N x N)
A = np.array([
    [0.0,  -0.07,  0.015],
    [1.0,   0.8,  -0.15],
    [0.0,   0.0,   0.5]
])

# B (N x M)
B = np.array([
    [0.0,  2.0],
    [-1.0, -0.1],
    [1.0,  1.0]
])

# C (P x N)
C = np.array([
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 0.5]
])

# D (P x M)
D = np.array([
    [0.0, 0.5],
    [1.0, 0.0]
])

# Initial State Vector X(1) used by TF01MD
# Note: Python libraries often use X(0) as the initial state.
# We will use this X_initial as X(0) for the Python simulations.
X_initial = np.array([1.0, 1.0, 1.0])

# Input Sequence U
# Data is taken from the Fortran column-major vector in the C++ test
# Fortran layout: [u1(1), u2(1), u1(2), u2(2), ..., u1(10), u2(10)]
U_fortran_flat = np.array([
    -0.6922, 0.2614,  # k=1
    -1.4934, -0.9160, # k=2
     0.3081, -0.6030, # k=3
    -2.7726, 1.2556,  # k=4
     2.0039, 0.2951,  # k=5
    -1.5734, -0.9942, # k=6
     1.5639, 1.8957,  # k=7
    -0.9942, 0.8988,  # k=8
     1.8957, -0.0701, # k=9
     0.8988, 0.4118   # k=10
])
# Reshape to Python convention (M rows, NY columns)
# We need [[u1(1), u1(2), ...], [u2(1), u2(2), ...]]
U = U_fortran_flat.reshape((NY, M)).T  # Reshape to NYxM then transpose to MxNY

# Expected Output Y (from SLICOT Documentation Example)
# Fortran layout: [y1(1), y2(1), y1(2), y2(2), ...]
Y_expected_fortran_flat = np.array([
     0.3078, -0.0928, # k=1
    -1.5125, 1.2611,  # k=2
    -1.2577, 3.4002,  # k=3
    -0.2947, -0.7060, # k=4
    -0.5632, 5.4532,  # k=5
    -1.0846, 1.1846,  # k=6
    -1.2427, 2.2286,  # k=7
     1.8097, -1.9534, # k=8
     0.6685, -4.4965, # k=9
    -0.0896, 1.1654   # k=10
])
# Reshape to Python convention (P rows, NY columns)
Y_expected = Y_expected_fortran_flat.reshape((NY, P)).T # Reshape to NYxP then transpose to PxNY

# --- Print System Info ---
print("System Matrices:")
print("A:\n", A)
print("\nB:\n", B)
print("\nC:\n", C)
print("\nD:\n", D)
print("\nInitial State (X0 for Python libs):\n", X_initial)
print("\nInput Sequence U (M x NY):\n", U)
print("\nExpected Output Y (from SLICOT docs) (P x NY):\n", Y_expected)

# --- Simulation ---
# Generate time vector (0 to NY-1)
T = np.arange(NY) # Time steps k = 0, 1, ..., NY-1

# Initialize results
yout_control = None
yout_harold = None

# --- Python Control Simulation ---
if control_found:
    try:
        sys_control = control.ss(A, B, C, D, dt=True) # dt=True indicates discrete time
        timepts_control, yout_control, xout_control = control.forced_response(
            sys_control, T=T, U=U, X0=X_initial, return_x=True
        )
        print("\nCalculated Output yout (P x NY) using python-control:\n", yout_control)
    except Exception as e:
        print(f"\nError during python-control simulation: {e}")
else:
     print("\nSkipping python-control simulation.")

# --- Harold Simulation ---
if harold_found:
    try:
        sys_harold = State(A, B, C, D, dt=1)  # dt=1 for discrete time
        # Transpose U to time x inputs format that harold expects
        U_harold_fmt = U.T
        yout_harold_raw, tout_harold = simulate_linear_system(sys_harold, U_harold_fmt, t=T, x0=X_initial)
        # Transpose harold output back to (outputs, time) format
        yout_harold = yout_harold_raw.T
        print("\nCalculated Output using harold (P x NY):\n", yout_harold)
    except Exception as e:
        print(f"\nError during harold simulation: {e}")
else:
    print("\nSkipping harold simulation.")


# --- Comparison ---
print("\n--- Comparisons (atol=1e-4) ---")
if yout_control is not None:
    print("Control outputs close to expected (SLICOT docs):", np.allclose(yout_control, Y_expected, atol=1e-4))
if yout_harold is not None:
    print("Harold outputs close to expected (SLICOT docs): ", np.allclose(yout_harold, Y_expected, atol=1e-4))
if yout_control is not None and yout_harold is not None:
    print("Harold and Control outputs close:               ", np.allclose(yout_harold, yout_control, atol=1e-8)) # Use tighter tolerance for self-comparison


# --- Plotting (Optional) ---
if (yout_control is not None or yout_harold is not None):
    plt.figure(figsize=(10, 8))
    # Plot Y1
    plt.subplot(P, 1, 1) # Adjust subplot based on P
    plt.plot(T, Y_expected[0,:], 'ko-', label='Expected Y1 (SLICOT Docs)')
    if yout_control is not None:
        plt.plot(T, yout_control[0,:], 'r.--', label='Control Y1')
    if yout_harold is not None:
        plt.plot(T, yout_harold[0,:], 'b^:', label='Harold Y1')
    plt.ylabel('Output Y1')
    plt.title('System Response Comparison')
    plt.legend()
    plt.grid(True)

    # Plot Y2
    if P > 1:
        plt.subplot(P, 1, 2) # Adjust subplot based on P
        plt.plot(T, Y_expected[1,:], 'ko-', label='Expected Y2 (SLICOT Docs)')
        if yout_control is not None:
            plt.plot(T, yout_control[1,:], 'r.--', label='Control Y2')
        if yout_harold is not None:
            plt.plot(T, yout_harold[1,:], 'b^:', label='Harold Y2')
        plt.ylabel('Output Y2')
        plt.xlabel('Time step k')
        plt.legend()
        plt.grid(True)

    # Add more subplots if P > 2

    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust layout
    try:
        plt.show()
    except Exception as e:
        print(f"\nPlotting failed (display not available?): {e}")
else:
    print("\nNo simulation results available for plotting.")

