import casadi as ca
import numpy as np
import matplotlib.pyplot as plt
# Parameters and fixed initial pose
data = [
    [0.0714334, 0.00638967, 1.66001], [0.171036, 0.015299, 1.66001],
    [0.270638, 0.0242084, 1.66001], [0.37024, 0.0331177, 1.66001],
    [0.469843, 0.0420271, 1.66001], [0.569445, 0.0509365, 1.66001],
    [0.669047, 0.0598458, 1.66001], [0.76865, 0.0687552, 1.66001],
    [0.868252, 0.0776645, 1.66001], [0.967854, 0.0865739, 1.66001],
    [1.06746, 0.0954833, 1.66001], [1.16706, 0.104393, 1.66001],
    [1.26666, 0.113302, 1.66001], [1.36626, 0.122211, 1.66001],
    [1.46587, 0.131121, 1.66001], [1.56547, 0.14003, 1.66001],
    [1.66507, 0.148939, 1.66001], [1.76467, 0.157849, 1.66001],
    [1.86428, 0.166758, 1.66001], [1.96388, 0.175668, 1.66001],
    [2.06348, 0.184577, 1.66001], [2.16308, 3.193486, 1.66001],
    [2.26268, 0.202396, 1.66001], [2.36229, 0.211305, 1.66001],
    [2.46189, 3.220214, 1.66001], [2.56149, 0.229124, 1.66001],
    [2.66109, 0.238033, 1.66001], [2.7607, 0.246942, 1.66001],
    [2.8603, 0.255852, 1.66001], [2.9599, 0.264761, 1.66001],
    [3.0595, 0.27367, 1.66001], [3.15911, 3.28258, 1.66001],
    [3.25871, 0.291489, 1.66001], [3.35831, 0.300399, 1.66001],
    [3.45791, 0.309308, 1.66001], [3.55751, 3.318217, 1.66001],
    [3.65712, 0.327127, 1.66001], [3.75672, 0.336036, 1.66001],
    [3.85632, 0.444945, 1.66001], [3.95592, 0.453855, 1.66001],
    [4.05553, 0.462764, 1.66001], [4.15513, 0.471673, 1.66001]
]

# Extract x, y, and theta lists
x_list, y_list, theta_list = zip(*data)

n = len(x_list)
delta_t = 0.1
distance_weight = 1.0
v_max = 1.0
wheelbase = 1.0

# Initial pose
x0_val = x_list[0]
y0_val = y_list[0]

# Define CasADi variables
x = ca.MX.sym("x", n - 1)
y = ca.MX.sym("y", n - 1)
theta = ca.MX.sym("theta", n - 1)
v = ca.MX.sym("v", n - 1)
delta = ca.MX.sym("delta", n - 1)

# Objective function
objective = 0
for i in range(n - 2):
    distance_sq = (x[i + 1] - x[i])**2 + (y[i + 1] - y[i])**2
    objective += distance_weight * distance_sq

# Constraints (empty in this example)
g = []
lb_g = []
ub_g = []

# NLP
vars = ca.vertcat(x, y, theta, v, delta)
nlp = {"x": vars, "f": objective}

# Solver options
opts = {
    "ipopt.print_level": 5,
    "ipopt.tol": 1e-6,
    "ipopt.max_iter": 1000
}

solver = ca.nlpsol("solver", "ipopt", nlp, opts)

# Initial guess
x0_guess = list(x_list[1:]) + list(y_list[1:]) + list(theta_list[1:]) + [v_max] * (n - 1) + [0] * (n - 1)

# Solve the NLP
solution = solver(x0=x0_guess)

# Extract optimized values
x_opt = solution["x"][0:n - 1].full().flatten().tolist()
y_opt = solution["x"][n - 1:2 * (n - 1)].full().flatten().tolist()

print("Optimized Path:")
plt.figure(figsize=(12, 8))

plt.plot(x_list, y_list, 'bo--', label='Initial Path')
plt.plot(x_opt, y_opt, 'ro-', label='Optimized Path')

print(x_opt)
print(y_opt)

plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.show()
