import casadi as ca
import numpy as np

# Define the black line parameters
A_middle = ca.DM([-0.0645709, -0.997913])
b_middle = 0.0524145
A_norm = ca.norm_2(A_middle)

# Define blue line points (example points, replace with your data)
blue_points = np.array([
    [0.0714334, 0.00638967, 1.66001],
    [0.171036, 0.015299, 1.66001],
    [0.270638, 0.0242084, 1.66001],
    [0.37024, 0.0331177, 1.66001],
    [0.469843, 0.0420271, 1.66001],
    [0.569445, 0.0509365, 1.66001],
    [0.669047, 0.0598458, 1.66001],
    [0.76865, 0.0687552, 1.66001],
    [0.868252, 0.0776645, 1.66001],
    [0.967854, 0.0865739, 1.66001],
    [1.06746, 0.0954833, 1.66001],
    [1.16706, 0.104393, 1.66001],
    [1.26666, 0.113302, 1.66001],
    [1.36626, 0.122211, 1.66001],
    [1.46587, 0.131121, 1.66001],
    [1.56547, 0.14003, 1.66001],
    [1.66507, 0.148939, 1.66001],
    [1.76467, 0.157849, 1.66001],
    [1.86428, 0.166758, 1.66001],
    [1.96388, 0.175668, 1.66001],
    [2.06348, 0.184577, 1.66001],
    [2.16308, 0.193486, 1.66001],
    [2.26268, 0.202396, 1.66001],
    [2.36229, 0.211305, 1.66001],
    [2.46189, 0.220214, 1.66001],
    [2.56149, 0.229124, 1.66001],
    [2.66109, 0.238033, 1.66001],
    [2.7607, 0.246942, 1.66001],
    [2.8603, 0.255852, 1.66001],
    [2.9599, 0.264761, 1.66001],
    [3.0595, 0.27367, 1.66001],
    [3.15911, 0.28258, 1.66001],
    [3.25871, 0.291489, 1.66001],
    [3.35831, 0.300399, 1.66001],
    [3.45791, 0.309308, 1.66001],
    [3.55751, 0.318217, 1.66001],
    [3.65712, 0.327127, 1.66001],
    [3.75672, 0.336036, 1.66001],
    [3.85632, 0.344945, 1.66001],
    [3.95592, 0.353855, 1.66001],
    [4.05553, 0.362764, 1.66001],
    [4.15513, 0.371673, 1.66001]
])

x = ca.MX.sym('x', blue_points.shape[0], 2)  # Decision variables (optimized positions for blue points)

# Objective function: minimize squared distances to the black line
objective = 0
for i in range(blue_points.shape[0]):
    point = x[i, :]
    distance = (A_middle.T @ point + b_middle) / A_norm
    objective += ca.sqr(distance)

# Define optimization problem
opts = {'ipopt.print_level': 0, 'print_time': 0}
nlp = {'x': x.reshape((-1,)), 'f': objective}
solver = ca.nlpsol('solver', 'ipopt', nlp, opts)

# Solve the optimization
blue_initial = blue_points.flatten()  # Flatten initial guess
solution = solver(x0=blue_initial)
optimized_points = np.array(solution['x']).reshape((-1, 2))

# Print results
print("Optimized points on the black line:", optimized_points)
