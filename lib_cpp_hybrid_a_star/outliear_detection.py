import pandas as pd
import numpy as np 

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from geopy.distance import geodesic
import numpy as np
import matplotlib.pyplot as plt


# Function to calculate 3D distance
def compute_3d_distance(lat1, lon1, alt1, lat2, lon2, alt2):
    # Calculate 2D distance (latitude and longitude)
    point1 = (lat1, lon1)
    point2 = (lat2, lon2)
    horizontal_distance = geodesic(point1, point2).kilometers

    # Include altitude for 3D distance
    vertical_distance = abs(alt2 - alt1) / 1000.0  # Convert altitude to kilometers
    distance_3d = np.sqrt(horizontal_distance**2 + vertical_distance**2)
    return distance_3d

def compute_2d_distance(lat1, lon1, lat2, lon2):
    # Calculate 2D geodesic distance
    point1 = (lat1, lon1)
    point2 = (lat2, lon2)
    distance_2d = geodesic(point1, point2).kilometers
    return distance_2d


# Load the CSV file
csv_filename = "/home/op/projects/wp-tt-elite/navsatfix_data_20241230_165640.csv"  # Replace with your actual CSV file name
data = pd.read_csv(csv_filename)

# data = data[0:2600]
from geopy.distance import geodesic  # To calculate geographic distance

# Initial ground height and parameters
GROUND_HEIGHT = 8.2
THRESHOLD = 0.3  # Adjust threshold as needed
alpha = 0.001  # Learning rate for updating ground height

# Initialize dynamic ground height
dynamic_ground_height = GROUND_HEIGHT
outlier_periods = []
in_outlier = False

# Iterate through data to find outlier periods
for i, altitude in enumerate(data['Altitude']):
    is_outlier = np.abs(altitude - dynamic_ground_height) > THRESHOLD

    if is_outlier:
        if not in_outlier:
            start_index = i
            in_outlier = True
    else:
        if in_outlier:
            end_index = i - 1
            outlier_periods.append((start_index, end_index))
            in_outlier = False

        # Update dynamic ground height
        dynamic_ground_height = (1 - alpha) * dynamic_ground_height + alpha * altitude

# Handle case where the data ends in an outlier period
if in_outlier:
    outlier_periods.append((start_index, len(data) - 1))

# Extract indices of outliers
outlier_indices = []
for start, end in outlier_periods:
    outlier_indices.extend(range(start, end + 1))



# Parameters
Z_SCORE_THRESHOLD = 8  # Z-score threshold
IQR_MULTIPLIER = 4.0   # IQR multiplier

# Initialize variables
historical_data = []  # Holds the non-outlier historical data
outliers = []         # Records indices of outliers
altitudes = data['Altitude']

# Real-time outlier detection
for index, altitude in enumerate(altitudes):
    if len(historical_data) > 1:
        # Calculate mean and standard deviation
        mean_altitude = np.mean(historical_data)
        std_altitude = np.std(historical_data)
        
        # Calculate IQR
        Q1 = np.percentile(historical_data, 25)
        Q3 = np.percentile(historical_data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_MULTIPLIER * IQR
        upper_bound = Q3 + IQR_MULTIPLIER * IQR

        # Detect outliers
        z_score = (altitude - mean_altitude) / std_altitude if std_altitude > 0 else 0
        if (np.abs(z_score) > Z_SCORE_THRESHOLD) or (altitude < lower_bound) or (altitude > upper_bound):
            outliers.append(index)
            continue  # Skip updating historical stats for outliers
    
    # Update historical data
    historical_data.append(altitude)


# final_outliers = []  # Indices of confirmed outliers
# DISTANCE_THRESHOLD = 3.0  # Maximum distance (in kilometers) to validate outliers

# for idx in outliers:
#     current_point = (data.loc[idx, 'Latitude'], data.loc[idx, 'Longitude'])
#     is_isolated = True

#     # Compare to neighbors within a window
#     for j in range(max(0, idx - 5), min(len(data), idx + 5)):  # Adjust window size as needed
#         if j == idx:
#             continue
#         neighbor_point = (data.loc[j, 'Latitude'], data.loc[j, 'Longitude'])
#         distance = geodesic(current_point, neighbor_point).kilometers
#         if distance < DISTANCE_THRESHOLD:
#             is_isolated = False
#             break

#     if is_isolated:
#         final_outliers.append(idx)


# Calculate velocities
velocities = []
velocities_2d = []

timestamps = data['Time (s)'].to_numpy()

for i in range(1, len(data)):
    lat1, lon1, alt1 = data.loc[i - 1, ['Latitude', 'Longitude', 'Altitude']]
    lat2, lon2, alt2 = data.loc[i, ['Latitude', 'Longitude', 'Altitude']]
    distance_3d = compute_3d_distance(lat1, lon1, alt1, lat2, lon2, alt2)
    
    # Calculate time difference in hours
    time_diff = (timestamps[i] - timestamps[i - 1]) / 3600.0  # Convert seconds to hours
    if time_diff > 0:
        velocity = distance_3d / time_diff  # Velocity in km/h
    else:
        velocity = 0  # Handle cases where time difference is 0

    distance_2d = compute_2d_distance(lat1, lon1, lat2, lon2)
    if time_diff > 0:
        velocity_2d = distance_2d / time_diff  # Velocity in km/h
    else:
        velocity_2d = 0  # Handle cases where time difference is 0

    velocities_2d.append(velocity_2d)
    velocities.append(velocity)


VELOCITY_THRESHOLD = 50  # Example threshold in km/h
outlier_indices_3d = [i for i, v in enumerate(velocities) if v > VELOCITY_THRESHOLD]

# Plot velocity profile
plt.figure(figsize=(15, 8))
plt.plot(range(1, len(data)), velocities, label='3D Velocity (km/h)', color='green')
plt.scatter(outlier_indices_3d, [velocities[i] for i in outlier_indices_3d], color='red', label='High Velocity Outliers', s=20)
plt.plot(data.index, data['Altitude'], label='Altitude', color='blue', linewidth=2)
plt.axhline(GROUND_HEIGHT, color='orange', linestyle='--', label='Initial Ground Height', linewidth=2)
plt.xlabel('Index')
plt.ylabel('Altitude')
plt.title('Velocity Outliers')
plt.legend()
plt.grid()

# # Plot altitude with marked outliers
plt.figure(figsize=(10, 6))
plt.plot(data.index, data['Altitude'], label='Altitude', color='blue')
plt.scatter(data.index[outliers], data.loc[outliers, 'Altitude'], color='red', label='Outliers', s=20)
plt.axhline(GROUND_HEIGHT, color='orange', linestyle='--', label='Initial Ground Height', linewidth=2)
plt.xlabel('Index')
plt.ylabel('Altitude')
plt.title('Z-score threshold and IQR multiplier')
plt.legend()
plt.grid()

plt.figure(figsize=(10, 6))
plt.scatter(data['Longitude'], data['Latitude'], c='blue', label='Normal Points', s=10)
plt.scatter(data.loc[outliers, 'Longitude'], data.loc[outliers, 'Latitude'], c='red', label='Outliers Z-score threshold and IQR multiplier', s=20,)
plt.legend()
plt.grid()

plt.figure(figsize=(10, 6))
plt.scatter(data['Longitude'], data['Latitude'], c='blue', label='Normal Points', s=10)
plt.scatter(data.loc[outlier_indices_3d, 'Longitude'], data.loc[outlier_indices_3d, 'Latitude'], c='red', label='Outliers Velocity', s=20,)
plt.legend()
plt.grid()

plt.show()


