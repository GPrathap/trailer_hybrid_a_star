import casadi as ca
import numpy as np
import matplotlib.pyplot as plt
import json 
# Parameters and fixed initial pose
data = np.array([
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
    [2.16308, 3.193486, 1.66001],
    [2.26268, 0.202396, 1.66001],
    [2.36229, 0.211305, 1.66001],
    [2.46189, 3.220214, 1.66001],
    [2.56149, 0.229124, 1.66001],
    [2.66109, 0.238033, 1.66001],
    [2.7607, 0.246942, 1.66001],
    [2.8603, 0.255852, 1.66001],
    [2.9599, 0.264761, 1.66001],
    [3.0595, 0.27367, 1.66001],
    [3.15911, 3.28258, 1.66001],
    [3.25871, 0.291489, 1.66001],
    [3.35831, 0.300399, 1.66001],
    [3.45791, 0.309308, 1.66001],
    [3.55751, 3.318217, 1.66001],
    [3.65712, 0.327127, 1.66001],
    [3.75672, 0.336036, 1.66001],
    [3.85632, 0.444945, 1.66001],
    [3.95592, 0.453855, 1.66001],
    [4.05553, 0.462764, 1.66001],
    [4.15513, 0.471673, 1.66001]
])

A_p = np.array([[-0.0645709, -0.997913], [0.0645709, 0.997913]])
b_p = np.array([0.737278, 0.632449])

x_list = data[:, 0]
y_list = data[:, 1]
theta_list = data[:, 2]  # Sample orientation values

curvature_weight = 1.0
curvature_rate_weight = 1.0
distance_weight = 1.0
v_max = 1.0
delta_t = 0.1
wheelbase = 1.0  # Length of the wheelbase for car-like dynamics
n = len(x_list)

# Initial fixed pose
x0_val = x_list[0]
y0_val = y_list[0]
theta0_val = theta_list[0]

# Define optimization variables for remaining poses, velocities, and steering angles
x = ca.MX.sym('x', n - 1)
y = ca.MX.sym('y', n - 1)
theta = ca.MX.sym('theta', n - 1)
v = ca.MX.sym('v', n - 1)        # Velocity variable for each timestep
delta = ca.MX.sym('delta', n - 1)  # Steering angle for each timestep

# Objective function
objective = 0
# for i in range(n - 3):  # Adjust indices since initial pose is excluded
#     dds_x = x[i+2] - 2 * x[i+1] + x[i]
#     dds_y = y[i+2] - 2 * y[i+1] + y[i]
#     objective += curvature_weight * (dds_x**2 + dds_y**2)
    
# for i in range(n - 4):
#     ddds_x = x[i+3] - 3 * x[i+2] + 3 * x[i+1] - x[i]
#     ddds_y = y[i+3] - 3 * y[i+2] + 3 * y[i+1] - y[i]
#     objective += curvature_rate_weight * (ddds_x**2 + ddds_y**2)

# Add penalty for distance between consecutive points
# for i in range(n - 2):
#     distance_sq = (x[i+1] - x[i])**2 + (y[i+1] - y[i])**2
#     print(distance_sq)
#     objective += distance_weight * distance_sq  # Penalize large distances between consecutive points

# # Constraints
g = []
lb_g = []
ub_g = []


# Kinematic constraints for car-like dynamics
for i in range(n - 2):
    # x_t+1 = x_t + delta_t * v * cos(theta)
    g.append(x[i+1] - (x[i] + delta_t * v[i] * ca.cos(theta[i])))
    # print(x[i+1] - (x[i] + delta_t * v[i] * ca.cos(theta[i])))
    lb_g.append(0)
    ub_g.append(0)
    
    # # y_t+1 = y_t + delta_t * v * sin(theta)
    g.append(y[i+1] - (y[i] + delta_t * v[i] * ca.sin(theta[i])))
    lb_g.append(0)
    ub_g.append(0)
    
    # # theta_t+1 = theta_t + delta_t * v * tan(delta) / L
    # g.append(theta[i+1] - (theta[i] + delta_t * v[i] * ca.tan(delta[i]) / wheelbase))
    g.append(theta[i+1] - (theta[i] + delta_t * v[i] * ca.tan(delta[i]) / wheelbase))
    print("   ", theta[i+1] - (theta[i] + delta_t * v[i] * ca.tan(delta[i]) / wheelbase))
    lb_g.append(0)
    ub_g.append(0)
    


# # Velocity and steering angle constraints
# max_steering_angle = np.deg2rad(25)
max_steering_angle = 0.1
for i in range(n - 1):
    # Velocity constraints
    g.append(v[i])
    lb_g.append(0.0)  # Minimum velocity (no reversing)
    ub_g.append(v_max)  # Maximum velocity
    # Steering angle constraints (example: max steering angle of 30 degrees)
    g.append(delta[i])
    lb_g.append(-max_steering_angle)
    ub_g.append(max_steering_angle)
    

# for i in range(0, n-1):
#     for j in range(0, 2):
#         g.append(x[i]*A_p[j, 0] + y[i]*A_p[j, 1])
#         lb_g.append(-ca.inf)
#         ub_g.append(b_p[j])
        
# Define NLP problem without the fixed initial pose in the variables
nlp = {'x': ca.vertcat(x, y, theta, v, delta), 'f': objective, 'g': ca.vertcat(*g)}

# Set IPOPT solver options
opts = {
    'ipopt.print_level': 5,
    'print_time': False,
    'ipopt.tol': 1e-6,
    'ipopt.max_iter': 1000
}
solver = ca.nlpsol('solver', 'ipopt', nlp, opts)

# Initial guess (excluding the fixed initial pose)
x0_guess = np.hstack((x_list[1:], y_list[1:], theta_list[1:], np.ones(n - 1) * v_max, np.zeros(n - 1)))

print(x0_guess)

# Solve the problem
# solution = solver(x0=x0_guess, lbx=-ca.inf, ubx=ca.inf)
solution = solver(x0=x0_guess, lbx=-ca.inf, ubx=ca.inf, lbg=lb_g, ubg=ub_g)
# solution = solver(x0=x0_guess)


problem_data = {
    'variables': nlp['x'].size1(),
    'objective': str(objective),
    'constraints': [str(nlp['g'])],
    'bounds': {
        'lb_vars': -ca.inf,
        'ub_vars': ca.inf,
        'lb_constraints': lb_g,
        'ub_constraints': ub_g,
    }
}

# Write to JSON file
with open("nlp_debug_python.json", "w") as f:
    json.dump(problem_data, f, indent=4)

print("NLP problem saved to nlp_debug.json")

# Extract optimized values, adding back fixed initial pose
x_opt = solution['x'][:n - 1]
y_opt = solution['x'][n - 1:2 * (n - 1)]
theta_opt = np.vstack((np.array([[theta0_val]]), solution['x'][2 * (n - 1):]))
v_opt = solution['x'][3 * (n - 1):4 * (n - 1)]
delta_opt = solution['x'][4 * (n - 1):]

# Visualization
plt.figure(figsize=(12, 8))
plt.plot(x_list, y_list, 'bo--', label='Initial Path')
plt.plot(x_opt, y_opt, 'ro-', label='Optimized Path')

# Boundary line plots
x_vals = np.linspace(-5, 5, 400)
y1_vals = (b_p[0] - A_p[0, 0] * x_vals) / A_p[0, 1]
y2_vals = (b_p[1] - A_p[1, 0] * x_vals) / A_p[1, 1]
plt.plot(x_vals, y1_vals, label=r'$-0.368x - 0.929y \leq 2.24$')
plt.plot(x_vals, y2_vals, label=r'$-0.368x - 0.929y \leq 0.87$')

plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.show()