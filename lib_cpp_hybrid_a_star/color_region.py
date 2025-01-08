import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


csv_filename = "/home/op/projects/trailer_hybrid_a_star/lib_cpp_hybrid_a_star/navsatfix_data_20250102_102601.csv"  # Replace with your actual CSV file name
IQR_MULTIPLIER = Z_SCORE_THRESHOLD = GROUND_HEIGHT = 9.0 # Initial ground height and parameters

# csv_filename = "/home/op/projects/trailer_hybrid_a_star/lib_cpp_hybrid_a_star/navsatfix_data_20241230_165640.csv"  # Replace with your actual CSV file name
# IQR_MULTIPLIER = Z_SCORE_THRESHOLD = GROUND_HEIGHT = 8.3 # Initial ground height and parameters

# csv_filename = "/home/op/projects/trailer_hybrid_a_star/lib_cpp_hybrid_a_star/navsatfix_data_20250102_121117.csv"  # Replace with your actual CSV file name
# IQR_MULTIPLIER = Z_SCORE_THRESHOLD = GROUND_HEIGHT = 8.0 # Initial ground height and parameters


data = pd.read_csv(csv_filename)

# Sample data (replace with your actual data)
# data = np.array([i for i in range(1000)])  # Example: array from 0 to 999

# Apply np.diff iteratively for 20th difference on the last element only
data = data['Altitude'][0:-1]
data_mean = data.mean()
# data = data - data_mean

# Calculate the np.diff for every 20 points
indexi = 100
diff_array = []
for i in range(len(data) - indexi):
    diff_array.append((data[i+indexi] - data[i]))


plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(data, label="Original Data", marker='o')
plt.title("Original Data")
plt.xlabel("Index")
plt.ylabel("Value")
plt.legend()

# Plot the difference array
plt.subplot(2, 1, 2)
plt.plot(diff_array, label="Difference (np.diff of 20 points)", color='orange', marker='o')
plt.title("Difference Array")
plt.xlabel("Index")
plt.ylabel("Difference")
plt.legend()

# Show the plots
plt.tight_layout()
plt.show()
# # Calculate the third difference of Altitude
# diff_alt_3 = np.diff(data['Altitude'], n=3)



# def custom_diff(data, n=20):
#     def calculate_diff_20(data, n):
#         for i in range(len(data)):
#             if i < n:
#                 yield 0.0  # No diff for the first 20 elements
#             else:
#                 yield data[i] - data[i - n]
    
#     diff_3_generator = calculate_diff_20(data, n)
#     diff_3_updated = []
#     for i, value in enumerate(diff_3_generator):
#         diff_3_updated.append(value)
#     return np.array(diff_3_updated)
    

# diff_3_updated = custom_diff(data['Altitude'], n=20)
# # diff_3_generator = np.array(diff_3_generator)
# # Define a threshold for outliers (example: 0.5 meters)
# threshold = 0.5

# # Identify outlier indices based on the third difference
# outlier_indices = np.where(np.abs(diff_alt_3) > threshold)[0]
# outlier_indices_custom = np.where(np.abs(diff_3_updated) > threshold)[0]

# outlier_indices = outlier_indices_custom

# # Group consecutive outlier indices into regions
# outlier_regions = []
# start_idx = None

# for i in range(len(outlier_indices) - 1):
#     if start_idx is None:
#         start_idx = outlier_indices[i]
#     if outlier_indices[i + 1] - outlier_indices[i] > 1:
#         outlier_regions.append((start_idx, outlier_indices[i]))
#         start_idx = None
# if start_idx is not None:
#     outlier_regions.append((start_idx, outlier_indices[-1]))

# # Prepare Altitude + diff_alt_3 for visualization
# shifted_diff_alt_3 = np.pad(diff_alt_3, (3, 0), 'constant', constant_values=np.nan) + data['Altitude']

# # Plot Altitude and 3rd Diff
# fig, ax = plt.subplots(figsize=(14, 8))

# # Original Altitude
# ax.plot(data['Time (s)'], data['Altitude'], label='Altitude', color='blue', zorder=1)

# # Third difference of Altitude added to Altitude (shifted up for visibility)
# ax.plot(data['Time (s)'], shifted_diff_alt_3, label='3rd Diff + Altitude', color='green', linewidth=2, zorder=2)

# # Highlight Outlier Regions using fill_between
# for start, end in outlier_regions:
#     start_time = data.iloc[start]['Time (s)']
#     end_time = data.iloc[min(end + 1, len(data) - 1)]['Time (s)']
#     ax.fill_between(
#         [start_time, end_time],
#         data['Altitude'].min(),
#         data['Altitude'].max(),
#         color='red',
#         alpha=1.0,
#         label='Outlier Region' if 'Outlier Region' not in ax.get_legend_handles_labels()[1] else "",
#         zorder=0
#     )

# # Highlight Sign Changes
# sign_changes = np.where(np.diff(np.sign(diff_alt_3)) != 0)[0]
# sign_change_times = data.iloc[sign_changes]['Time (s)']
# sign_change_altitudes = data.iloc[sign_changes]['Altitude']
# ax.scatter(sign_change_times, sign_change_altitudes, color='orange', label='Sign Changes', zorder=3)

# # Add Labels, Title, and Legend
# ax.set_xlabel('Time (s)')
# ax.set_ylabel('Altitude')
# ax.set_title('Altitude with Highlighted Outlier Regions and Sign Changes')
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels)
# ax.grid(True)

# plt.show()
