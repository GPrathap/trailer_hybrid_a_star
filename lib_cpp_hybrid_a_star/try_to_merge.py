import pandas as pd
import numpy as np 

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from geopy.distance import geodesic
import numpy as np
import matplotlib.pyplot as plt
# data = data[0:2600]
from geopy.distance import geodesic  # To calculate geographic distance
import math

class OnlineStats:
    def __init__(self):
        self.count = 0       # Number of data points
        self.mean = 0.0      # Current mean
        self.m2 = 0.0        # Sum of squared differences from the mean

    def update(self, value):
        """
        Update the statistics with a new data point.
        """
        self.count += 1
        delta = value - self.mean
        self.mean += delta / self.count
        delta2 = value - self.mean
        self.m2 += delta * delta2

    def get_mean(self):
        """
        Return the current mean.
        """
        return self.mean

    def get_stddev(self):
        """
        Return the current standard deviation.
        """
        if self.count > 1:
            return math.sqrt(self.m2 / self.count)
        return 0.0


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

# # Load the CSV file
# csv_filename = "/home/op/projects/trailer_hybrid_a_star/lib_cpp_hybrid_a_star/navsatfix_data_20250102_102601.csv"  # Replace with your actual CSV file name
# IQR_MULTIPLIER = Z_SCORE_THRESHOLD = GROUND_HEIGHT = 9.0 # Initial ground height and parameters

# csv_filename = "/home/op/projects/trailer_hybrid_a_star/lib_cpp_hybrid_a_star/navsatfix_data_20241230_165640.csv"  # Replace with your actual CSV file name
# IQR_MULTIPLIER = Z_SCORE_THRESHOLD = GROUND_HEIGHT = 8.3 # Initial ground height and parameters

# csv_filename = "/home/op/projects/trailer_hybrid_a_star/lib_cpp_hybrid_a_star/navsatfix_data_20250102_121117.csv"  # Replace with your actual CSV file name
# IQR_MULTIPLIER = Z_SCORE_THRESHOLD = GROUND_HEIGHT = 8.0 # Initial ground height and parameters


# csv_filename = "/home/op/projects/trailer_hybrid_a_star/lib_cpp_hybrid_a_star/test6/combined_data.csv"  # Update with your actual file path
csv_filename = "/home/op/projects/trailer_hybrid_a_star/lib_cpp_hybrid_a_star/test7/combined_data.csv"
IQR_MULTIPLIER = Z_SCORE_THRESHOLD = GROUND_HEIGHT = 8.0 # Initial ground height and parameters



# data = pd.read_csv(csv_file)


data = pd.read_csv(csv_filename)
altitudes = data['altitude']
latitudes = data['latitude']
longitudes = data['longitude']
heading_stdev = data['heading_stdev']
# data = data[300:-1]



# Parameters
MAX_ALT_THRESHOLD = 3.5
MAX_HEADING_STD_THRESHOLD = 10.0

THRESHOLD = 0.3  # Adjust threshold as needed
alpha = 0.001  # Learning rate for updating ground height

# Initialize variables
historical_data = []  # Holds the non-outlier historical data
historical_data_diff = []  # Holds the non-outlier historical data
outliers = []         # Records indices of outliers
outlier_indices_3d = []
outlier_indices_2d = []


# Calculate velocities
velocities = []
velocities_2d = []
previous_lat_log = []
# Real-time outlier detection
detect_outlier =  False
outlier_periods = []
in_outlier = False
start_index = 0


z_score_outlier = False
ground_height_outlier = False 
velocity_outlier = False 
alt_outlier_extrem = False 
velocity_2d_outlier = False 
iqr_outlier = False 

outlier_jumps =  []
final_outlier =  []

stats = OnlineStats()
mean_altitude_all = []
last_alt = 0.0

all_history_data = []
data_between_peaks = []
peak_heading_last = 0.0
window_size = 3  # Sliding window for peak detection

# Helper function to check if a point is a peak
def is_peak(window):
    if len(window) < 3:
        return False  
    return window[0] < window[1] and window[1] > window[2] and (window[1] >=0.5) # Simple peak condition

capturing = False
capturing_peak = True
first_peak_index = 0 
second_peak_index = 0
is_outlier = False 
cumulative_std = []
for index, (alt, lat, log, hd_stddev) in enumerate(zip(altitudes, latitudes, longitudes, heading_stdev)):
    if len(historical_data) > 1:
        lat1, lon1, alt1 = previous_lat_log[0], previous_lat_log[1], previous_lat_log[2]
        lat2, lon2, alt2 = lat, log, alt
        distance_3d = compute_3d_distance(lat1, lon1, alt1, lat2, lon2, alt2)
        cumulative_std.append(distance_3d)
        all_history_data.append(hd_stddev)

        
        # Calculate mean and standard deviation
        mean_altitude = np.mean(historical_data)
        std_altitude = np.std(historical_data)
        
        # stats.update(last_alt)
        # mean_altitude_for_all = stats.get_mean()
        # mean_altitude_all.append(mean_altitude_for_all)
        # std_altitude = stats.get_stddev()
        
        # Calculate IQR
        Q1 = np.percentile(historical_data, 25)
        Q3 = np.percentile(historical_data, 60)
        # print(f"historical_data: {historical_data}")
        

        IQR = Q3 - Q1
        lower_bound = Q1 - IQR_MULTIPLIER * IQR
        upper_bound = Q3 + IQR_MULTIPLIER * IQR

       

        # Detect outliers
        z_score = (alt - mean_altitude) / std_altitude if std_altitude > 0 else 0
        if (alt < lower_bound) or (alt > upper_bound):
            outliers.append(index)
            detect_outlier = True  # Skip updating historical stats for outliers
            iqr_outlier = True
        else:
            detect_outlier = False
            iqr_outlier = False

        # print(f"Q1: {Q1} Q3: {Q3} meanAltitude: {mean_altitude} stdAltitude: {std_altitude} iqr_outlier: {iqr_outlier}")

        if (np.abs(z_score) > Z_SCORE_THRESHOLD):
            outliers.append(index)
            detect_outlier = True  # Skip updating historical stats for outliers
            z_score_outlier = True 
        else:
            detect_outlier = False
            z_score_outlier = False 

       

        # check_out = np.array(all_history_data)
        # index_outlier = np.diff(check_out, 2)
        
        alt_outlier_extrem = (abs(GROUND_HEIGHT - alt) > MAX_ALT_THRESHOLD)
        
        is_outlier = False

        

        if((peak_heading_last >= MAX_HEADING_STD_THRESHOLD) and (hd_stddev < MAX_HEADING_STD_THRESHOLD) and abs(second_peak_index -index)>100):
            capturing = True
            first_peak_index = index 
            # is_outlier = True
            # capturing_peak = True 
            print(f"First peak detected at index {first_peak_index} index {index}")
        elif((peak_heading_last < MAX_HEADING_STD_THRESHOLD) and (hd_stddev >= MAX_HEADING_STD_THRESHOLD) and capturing== True):
            second_peak_index = index  # Middle of the sliding window
            print(f"Second peak detected at index {second_peak_index}  index {index}")
            data_between_peaks.append([first_peak_index, second_peak_index])
            capturing_peak = True
            capturing = False 


        if(MAX_HEADING_STD_THRESHOLD < hd_stddev):
            is_outlier = True
        elif(capturing):
             is_outlier = True
            
        peak_heading_last = hd_stddev
        
        # if(len(index_outlier)>0):
        #     mean_altitude_all.append(index_outlier[-1])
        #     window = mean_altitude_all[-4:-1]
        #     if is_peak(window):  # Check if the current point is a peak
        #         if not capturing:
        #             # First peak detected
        #             capturing = True
        #             first_peak_index = index  # Middle of the sliding window
        #             # is_outlier = True
        #             # data_between_peaks.append(first_peak_index)
        #             print(f"First peak detected at index {first_peak_index}")
        #             # all_history_data = []
        #         else:
        #             # Second peak detected
        #             second_peak_index = index  # Middle of the sliding window
        #             print(f"Second peak detected at index {second_peak_index}")
        #             data_between_peaks.append([first_peak_index, second_peak_index])
        #             capturing = False
        #             is_outlier = False
        #             # all_history_data = []
        #             # Stop capturing data after finding the second peak

        # if(alt_outlier_extrem):
        #     capturing = False
        #     is_outlier = False
        # if(capturing):
        #     is_outlier = True 

        is_outlier_final = (z_score_outlier + iqr_outlier ) >= 1 
        # is_outlier_final = is_outlier
        # is_outlier_final = (z_score_outlier + iqr_outlier + velocity_outlier) >= 2  
        # is_outlier_final = (z_score_outlier + iqr_outlier) >= 2  
        # is_outlier_final = is_outlier 
        
        
        if(alt_outlier_extrem):
            final_outlier.append(index)   
        elif(is_outlier_final):
            final_outlier.append(index)   
    
    previous_lat_log = [lat, log, alt]

    if (abs(GROUND_HEIGHT - alt) < MAX_ALT_THRESHOLD):
        historical_data.append(alt)
    
    # if(len(historical_data_diff)>0.0):
    #     incre_alt = alt - last_alt + historical_data_diff[-1]
    # else:
    #     incre_alt = last_alt

    # historical_data_diff.append(incre_alt)
    # last_alt =  alt 
    
    if len(historical_data) > 500:  # Maintain the last 200 entries
            historical_data.pop(0)
            
    
    


outlier_indices = []
for start, end in data_between_peaks:
    print(start, end)
    outlier_indices.extend(range(start, end + 1))




print("--------------------------------------------", len(cumulative_std))
# # # # Plot locations of outliers
plt.figure(figsize=(10, 6))
plt.scatter(data['longitude'], data['latitude'], c='blue', label='Normal Points', s=10)
plt.scatter(data.loc[outlier_indices, 'longitude'], data.loc[outlier_indices, 'latitude'], c='red', label='Outliers', s=20)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Location of Outliers')
plt.legend()
plt.grid()

# Plot velocity profile
plt.figure(figsize=(15, 8))
# plt.plot(range(0, len(velocities)), velocities, label='3D Velocity (km/h)', color='green')
# plt.scatter(outlier_indices_3d, [velocities[i] for i in outlier_indices_3d], color='red', label='High Velocity Outliers 3d', s=20)
print("----------------------d1----------------------")
# plt.plot(range(0, len(cumulative_std)), cumulative_std, label='3d distance', color='blue', linewidth=2)
print("----------------------d2----------------------")
# plt.plot(data.index, data['altitude'], label='Altitude', color='blue', linewidth=2)
plt.plot(data.index, data['heading_stdev'], label="heading_stdev")
# mean_altitude_all = np.array(mean_altitude_all)
# print(mean_altitude_all.shape, data.index.shape)
# plt.plot(historical_data_diff)

# plt.axhline(GROUND_HEIGHT, color='orange', linestyle='--', label='Initial Ground Height', linewidth=2)
# plt.axhline(MAX_HEADING_STD_THRESHOLD, color='red', linestyle='--', label='Heading Dev', linewidth=2)

# diff_alt = np.diff(data['Altitude'])
# diff_alt_1 = np.diff(data['Altitude'], 2)
# diff_alt = np.append(diff_alt, 0)

# diff_alt_1 = np.append(diff_alt_1, 0)
# diff_alt_1 = np.append(diff_alt_1, 0)

# plt.plot(indexes, diff_alt, label='rtrgh', color='blue', linewidth=2)
# plt.plot(data.index, diff_alt_1, label='Altitude', color='red', linewidth=1)

plt.xlabel('Index')
plt.ylabel('Altitude')
plt.title('Velocity Outliers')
plt.legend()
plt.grid()

# # # Plot altitude with marked outliers
plt.figure(figsize=(10, 6))
plt.plot(data.index, data['altitude'], label='altitude', color='blue')
# plt.scatter(data.index[outliers], data.loc[outliers, 'Altitude'], color='red', label='Outliers', s=20)
plt.scatter(data.index[outlier_indices], data.loc[outlier_indices, 'altitude'], color='red', label='Outliers', s=20)
plt.axhline(GROUND_HEIGHT, color='orange', linestyle='--', label='Initial Ground Height', linewidth=2)
plt.xlabel('Index')
plt.ylabel('Altitude')
plt.title('Z-score threshold and IQR multiplier')
plt.legend()
plt.grid()

# plt.figure(figsize=(10, 6))
# plt.scatter(data['Longitude'], data['Latitude'], c='blue', label='Normal Points', s=10)
# plt.scatter(data.loc[outliers, 'Longitude'], data.loc[outliers, 'Latitude'], c='red', label='Outliers Z-score threshold and IQR multiplier', s=20,)
# plt.legend()
# plt.grid()

# plt.figure(figsize=(10, 6))
# plt.scatter(data['Longitude'], data['Latitude'], c='blue', label='Normal Points', s=10)
# plt.scatter(data.loc[outlier_indices_3d, 'Longitude'], data.loc[outlier_indices_3d, 'Latitude'], c='red', label='Outliers Velocity 3d', s=20,)
# plt.scatter(data.loc[outlier_indices_2d, 'Longitude'], data.loc[outlier_indices_2d, 'Latitude'], c='black', label='Outliers Velocity 2d', s=5,)

# plt.legend()
# plt.grid()

plt.figure(figsize=(10, 6))
plt.scatter(data['longitude'], data['latitude'], c='blue', label='Normal Points', s=10)
plt.scatter(data.loc[final_outlier, 'longitude'], data.loc[final_outlier, 'latitude'], c='red', label='Outliers Final', s=20,)
plt.legend()
plt.grid()

plt.show()


