#ifndef OUTLIER_DETECTOR_H_
#define OUTLIER_DETECTOR_H_

#include <deque>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

namespace outliear_detector {

  // Encapsulates parameters for outlier detection.
  struct OutlierDetectionParams{

    double iqr_multiplier;                 // Multiplier for IQR-based outlier detection
    double z_score_threshold;              // Threshold for Z-score-based outlier detection
    double ground_height;                  // Initial ground height
    double altitude_threshold;             // Threshold for altitude-based outlier detection
    double alpha;                       // Learning rate for dynamic ground height
    double velocity_threshold;           // Maximum velocity threshold (km/h)
    size_t max_history_size;               // Maximum size of the historical data buffer
    double max_altitude_threshold;         // Hard threshold for extreme altitude outliers
    double max_altitude_threshold_for_z_score;  // Soft threshold for Z-score adjustments

    // Constructor with default values
    OutlierDetectionParams(
        double ground_height = 8.0,
        double altitude_threshold = 0.3,
        double alpha = 0.001,
        double velocity_threshold = 100.0,
        size_t max_history_size = 500,
        double max_altitude_threshold = 1.5,
        double max_altitude_threshold_for_z_score = 1.2)
        : iqr_multiplier(ground_height),
          z_score_threshold(ground_height),
          ground_height(ground_height),
          altitude_threshold(altitude_threshold),
          alpha(alpha),
          velocity_threshold(velocity_threshold),
          max_history_size(max_history_size),
          max_altitude_threshold(max_altitude_threshold),
          max_altitude_threshold_for_z_score(max_altitude_threshold_for_z_score) {};

  };

   
  // Loads data from a CSV file into a 2D vector.
  std::vector<std::vector<double>> LoadCSV(const std::string& filename);

  // Class to detect outliers based on multiple criteria.
  class OutlierDetector {
    public:
      // Constructor to initialize the detector with user-defined parameters.
      OutlierDetector() = default;
      ~OutlierDetector() = default;

      void initParam(const OutlierDetectionParams& params);

      // Processes a single row of data and updates the state.
      void ProcessRow(int index, const std::vector<double>& row);

    // private:
      // Computes the 3D distance between two geographical points.
      double Compute3DDistance(double lat1, double lon1, double alt1,
                              double lat2, double lon2, double alt2);

      OutlierDetectionParams params_;  // Parameters for outlier detection.

      // Data buffers and state variables.
      std::deque<double> historical_data_;  // Buffer to store recent altitude data.
      std::vector<int> outliers_;           // Indices of identified outliers.
      std::vector<int> final_outliers_;     // Finalized list of outliers.
      std::vector<int> outlier_indices_3d_; // Indices of velocity-based outliers.
      std::vector<int> outlier_indices_; // Indices of velocity-based outliers.
      std::vector<double> velocities_;      // List of computed velocities.
      std::vector<std::vector<int>> outlier_periods_;  // Periods of consecutive outliers.

      double dynamic_ground_height_;  // Adjusted ground height based on recent data.
      double previous_time_index_ = 0.0;  // Timestamp of the last processed row.
      std::vector<double> previous_lat_lon_alt_;  // Previous latitude, longitude, and altitude.

      bool in_outlier_ = false;  // Tracks whether we are currently in an outlier period.
      int start_index_ = 0;      // Start index of the current outlier period.
      int end_index_ = 0;        // End index of the current outlier period.
  };


}
#endif  // OUTLIER_DETECTOR_H_
