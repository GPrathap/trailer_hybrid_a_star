import pandas as pd
import numpy as np 
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt

# Load the CSV file
csv_filename = "/home/op/projects/trailer_hybrid_a_star/lib_cpp_hybrid_a_star/navsatfix_data_20250102_121117.csv"  # Replace with your actual CSV file name
# csv_filename = "/home/op/projects/trailer_hybrid_a_star/lib_cpp_hybrid_a_star/navsatfix_data_20250102_102601.csv"  # Replace with your actual CSV file name
# csv_filename = "/home/op/projects/trailer_hybrid_a_star/lib_cpp_hybrid_a_star/navsatfix_data_20241230_165640.csv"  # Replace with your actual CSV file name
data = pd.read_csv(csv_filename)

# Ensure the columns exist in the data
required_columns = ["Altitude", "Latitude", "Longitude", "Time (s)"]
if not all(column in data.columns for column in required_columns):
    raise ValueError("The CSV file must contain the following columns: Altitude, Latitude, Longitude, Time (s)")

# Extract data
altitude = data["Altitude"]
latitude = data["Latitude"]
longitude = data["Longitude"]
time = data["Time (s)"]

altitude_diff = altitude.diff().fillna(0) 
altitude_diff_cummulative_sum = np.cumsum(altitude_diff)


# Create plots
plt.figure(figsize=(12, 8))

# Plot Altitude
plt.subplot(3, 1, 1)
plt.plot(time, altitude, marker='o', label='Altitude')
plt.plot(time, altitude_diff, label="Altitude Difference", color="green")
plt.plot(time, altitude_diff_cummulative_sum, label="Altitude Difference", color="blue")
plt.title("Altitude over Time")
plt.xlabel("Time (s)")
plt.ylabel("Altitude")
plt.grid(True)
plt.legend()

# Plot Latitude
plt.subplot(3, 1, 2)
plt.plot(time, latitude, marker='o', color='orange', label='Latitude')
plt.title("Latitude over Time")
plt.xlabel("Time (s)")
plt.ylabel("Latitude")
plt.grid(True)
plt.legend()

# Plot Longitude
plt.subplot(3, 1, 3)
plt.plot(time, longitude, marker='o', color='green', label='Longitude')
plt.title("Longitude over Time")
plt.xlabel("Time (s)")
plt.ylabel("Longitude")
plt.grid(True)
plt.legend()

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
