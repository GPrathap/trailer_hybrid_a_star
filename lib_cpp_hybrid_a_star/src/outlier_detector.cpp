#include "lib_cpp_hybrid_a_star/outlier_detector.hpp"

namespace outliear_detector {

    std::vector<std::vector<double>> LoadCSV(const std::string& filename) {
        std::vector<std::vector<double>> data;
        std::ifstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return data;
        }

        std::string line;
        int row_number = 0;
        while (std::getline(file, line)) {
            row_number++;
            std::stringstream ss(line);
            std::string value;
            std::vector<double> row;

            // Parse each value in the CSV row.
            while (std::getline(ss, value, ',')) {
            try {
                row.push_back(std::stod(value));
            } catch (const std::invalid_argument&) {
                std::cerr << "Error: Invalid data at row " << row_number << " - " << value << std::endl;
                row.clear();
                break;  // Skip this row if invalid data is encountered.
            }
            }

            // Ensure the row has the expected number of columns.
            if (!row.empty() && row.size() == 4) {
            data.push_back(row);
            }
        }

        file.close();
        return data;
    }

    void OutlierDetector::initParam(const OutlierDetectionParams& params){
        params_ = params;
        dynamic_ground_height_ = params.ground_height;
    };

    double OutlierDetector::Compute3DDistance(double lat1, double lon1, double alt1,
                                                double lat2, double lon2, double alt2) {
        const double kEarthRadiusKm = 6371.0;  // Earth's radius in kilometers.
        double d_lat = (lat2 - lat1) * M_PI / 180.0;
        double d_lon = (lon2 - lon1) * M_PI / 180.0;

        // Haversine formula for horizontal distance.
        double a = sin(d_lat / 2) * sin(d_lat / 2) +
                    cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) *
                    sin(d_lon / 2) * sin(d_lon / 2);
        double c = 2 * atan2(sqrt(a), sqrt(1 - a));
        double horizontal_distance = kEarthRadiusKm * c;

        // Calculate the vertical distance in kilometers.
        double vertical_distance = std::fabs(alt2 - alt1) / 1000.0;

        // Return the total 3D distance.
        return std::sqrt(horizontal_distance * horizontal_distance +
                        vertical_distance * vertical_distance);
    }

    void OutlierDetector::ProcessRow(int index, const std::vector<double>& row) {
        double time_index = row[3];
        double lat = row[1];
        double lon = row[2];
        double alt = row[0];

        bool z_score_outlier = false;
        bool iqr_outlier = false;
        bool ground_height_outlier = false;
        bool velocity_outlier = false;

        if (!previous_lat_lon_alt_.empty()) {
            // Calculate 3D distance and velocity.
            double distance_3d = Compute3DDistance(
                previous_lat_lon_alt_[0], previous_lat_lon_alt_[1], previous_lat_lon_alt_[2],
                lat, lon, alt);
            double time_diff = (time_index - previous_time_index_) / 3600.0;  // Convert seconds to hours.
            double velocity = (time_diff > 0) ? distance_3d / time_diff : 0.0;
            velocities_.push_back(velocity);

            // Check for velocity-based outliers.
            if (velocity > params_.velocity_threshold) {
            outlier_indices_3d_.push_back(index);
            velocity_outlier = true;
            }

            // Statistical outlier detection using IQR and Z-score.
            if (historical_data_.size() > 1) {
            Eigen::VectorXd data(historical_data_.size());
            for (size_t i = 0; i < historical_data_.size(); ++i) {
                data(i) = historical_data_[i];
            }

            double mean_altitude = data.mean();
            double std_altitude = std::sqrt((data.array() - mean_altitude).square().mean());

            std::nth_element(data.data(), data.data() + data.size() / 4, data.data() + data.size());
            double q1 = data(data.size() / 4);
            std::nth_element(data.data(), data.data() + 3 * data.size() / 4, data.data() + data.size());
            double q3 = data(3 * data.size() / 4);

            double iqr = q3 - q1;
            double lower_bound = q1 - params_.iqr_multiplier * iqr;
            double upper_bound = q3 + params_.iqr_multiplier * iqr;

            if (alt < lower_bound || alt > upper_bound) {
                outliers_.push_back(index);
                iqr_outlier = true;
            }

            if (std_altitude > 0) {
                double z_score = (alt - mean_altitude) / std_altitude;
                if (std::fabs(z_score) > params_.z_score_threshold) {
                outliers_.push_back(index);
                z_score_outlier = true;
                }
            }
            }

            // Altitude-based dynamic ground height adjustment.
            if (std::fabs(alt - dynamic_ground_height_) > params_.altitude_threshold) {
                if (!in_outlier_) {
                    start_index_ = index;
                    in_outlier_ = true;
                }
                ground_height_outlier = true;
            } else {
                if (in_outlier_) {
                    end_index_ = index - 1;
                    outlier_periods_.push_back({start_index_, end_index_});
                    in_outlier_ = false;
                }
                dynamic_ground_height_ =
                    (1 - params_.alpha) * dynamic_ground_height_ + params_.alpha * alt;
            }

            // Finalize outlier decisions.
            bool extreme_altitude_outlier = (std::fabs(params_.ground_height - alt) > params_.max_altitude_threshold);
            bool can_be_outlier = (z_score_outlier + iqr_outlier +
                                ground_height_outlier + velocity_outlier) >= 2;
            if (extreme_altitude_outlier || can_be_outlier) {
                final_outliers_.push_back(index);
            }
        }

        // Update state with the current row.
        previous_time_index_ = time_index;
        previous_lat_lon_alt_ = {lat, lon, alt};

        // Update historical data buffer.
        if (std::fabs(params_.ground_height - alt) < params_.max_altitude_threshold_for_z_score) {
            historical_data_.push_back(alt);
        }

        if (historical_data_.size() >= params_.max_history_size) {
            historical_data_.pop_front();
        }
    }

}

