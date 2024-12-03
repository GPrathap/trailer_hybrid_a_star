import numpy as np
import matplotlib.pyplot as plt

def plot_line_and_distance(A, b, point):
    """
    Plot a line, a point, and the perpendicular distance from the point to the line.
    
    Parameters:
    A (array): Coefficients of the line (e.g., [a1, a2]).
    b (float): Intercept of the line.
    point (array): Point coordinates (e.g., [x, y]).
    """
    # Create a grid for plotting the line
    x_vals = np.linspace(-10, 10, 500)
    y_vals = (b - A[0] * x_vals) / A[1]
    
    # Compute the closest point on the line
    A = np.array(A)
    point = np.array(point)
    t = (b - np.dot(A, point)) / np.dot(A, A)
    closest_point = point + t * A

    # Plot the line
    plt.plot(x_vals, y_vals, label=f'Line: {A[0]}x + {A[1]}y = {b}', color='blue')

    # Plot the original point
    plt.scatter(point[0], point[1], color='red', label='Point')
    
    # Plot the closest point on the line
    plt.scatter(closest_point[0], closest_point[1], color='green', label='Closest point on line')
    
    # Draw the perpendicular distance
    plt.plot([point[0], closest_point[0]], [point[1], closest_point[1]], 'k--', label='Perpendicular distance')

    # Set up the plot
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(True)
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Line and Perpendicular Distance')
    plt.show()


def plot_two_lines(A1, b1, A2, b2):
    """
    Plot two lines given their coefficients and intercepts.
    
    Parameters:
    A1, A2 (array): Coefficients of the lines (e.g., [a1, a2]).
    b1, b2 (float): Intercepts of the lines.
    """
    # Create a grid of x values for plotting
    x_vals = np.linspace(-10, 10, 500)
    
    # Calculate y values for both lines
    y1_vals = (b1 - A1[0] * x_vals) / A1[1]
    y2_vals = (b2 - A2[0] * x_vals) / A2[1]
    
    # Plot the lines
    plt.plot(x_vals, y1_vals, label=f'Line 1: {A1[0]}x + {A1[1]}y = {b1}', color='blue')
    plt.plot(x_vals, y2_vals, label=f'Line 2: {A2[0]}x + {A2[1]}y = {b2}', color='orange')

    # Set up the plot
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    


# A1 = np.array([2, 3])
# b1 = 5

# A2 = np.array([2.1, 3.1])
# b2 = 5.2
A1 = np.array([-0.0645709, -0.997913])
b1 = 0.737278

A2 = np.array([0.0645709,  0.997913])
b2 = 0.632449

if np.allclose(A1, -A2):
    # Compute the middle line
    A_middle = A1  # Either A1 or A2 works, as they are normalized and opposites
    b_middle = (b1 - b2) / 2
else:
    A_middle = (A1 + A2) / 2.0
    b_middle = (b1 + b2) / 2.0

print("Middle line coefficients (A):", A_middle)
print("Middle line intercept (b):", b_middle)

# # Example usage
# A = [2.05, 3.05]
# b = 5.1
point = [1, 2]
plot_two_lines(A1, b1, A2, b2)
plot_line_and_distance(A_middle, b_middle, point)
