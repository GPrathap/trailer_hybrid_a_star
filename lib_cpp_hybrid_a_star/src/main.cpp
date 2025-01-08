#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>
#include <cmath>
#include "lib_cpp_hybrid_a_star/rs_paths.hpp"
#include "lib_cpp_hybrid_a_star/grid_a_star.hpp"
#include "lib_cpp_hybrid_a_star/trailerlib.hpp"
#include "lib_cpp_hybrid_a_star/trailer_hybrid_a_star.hpp"
#include "matplotlibcpp.h"

#include "lib_cpp_hybrid_a_star/outlier_detector.hpp"

namespace plt = matplotlibcpp;

namespace plt = matplotlibcpp;
using namespace std;
using namespace Eigen;

// int main(int argc, char *argv[]) {

//     rs_paths::RSPaths rs_path;
//     Eigen::Vector3d s(14.0, 10.0, math_utility::deg2rad(0.0));
//     Eigen::Vector3d g(0.0, 0.0, math_utility::deg2rad(90.0));
//     double max_curvature = 0.1;
    
    
    
//     rs_paths::Path path = rs_path.calc_shortest_path(s, g, max_curvature);
//     std::vector<double> rc;
//     std::vector<double> rds;
//     std::cout<< "========1main 0" << path.poses.rows() << std::endl;
//     rs_path.calc_curvature(path.poses, rc, rds);
//     std::vector<double> path_short_x, path_short_y;
//     for(int i=0; i<path.poses.rows(); i++){
//         path_short_x.push_back(path.poses.row(i)[0]);
//         path_short_y.push_back(path.poses.row(i)[1]);
//     }
  

//     std::vector<std::vector<double>> bpath_x, bpath_y;
//     std::vector<rs_paths::Path> paths;
//     rs_path.calc_paths(s, g, max_curvature, paths);
//     for(auto path : paths){
//         std::vector<double> path_info_x, path_info_y;
//         for(int i=0; i<path.poses.rows(); i++){
//             path_info_x.push_back(path.poses.row(i)[0]);
//             path_info_y.push_back(path.poses.row(i)[1]);
//         }
//         bpath_x.push_back(path_info_x);
//         bpath_y.push_back(path_info_y);
//     }
//     // std::cout<< "========1main 1" << path.poses.rows() << std::endl;

//     // // First subplot
//     plt::figure();
//     for(int i=0; i< bpath_x.size(); i++){
//         plt::plot(bpath_x[i], bpath_y[i]);
//     }

//     plt::plot({s.x()}, {s.y()}, "bo"); // Start point
//     plt::plot({g.x()}, {g.y()}, "go");    // End point

//     plt::legend();
//     plt::grid(true);
//     plt::axis("equal");

//     // // // Second subplot for curvature
//     plt::figure();
//     plt::plot(path_short_x, path_short_y, "-r");
//     // plt::plot(rc, ".r");
//     plt::grid(true);
//     // plt::title("Curvature");

//     // // // Show all plots
//     plt::show();




//     return 0;
// }


struct PointCloud {
    Eigen::MatrixXd points;

    inline size_t kdtree_get_point_count() const { return points.rows(); }

    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return points(idx, dim);
    }

    template <class BBOX> bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>;

// int main(int arc, char *argv[]){

//     // Eigen::Vector4d s(-10.0, 6.0, math_utility::deg2rad(0.0), math_utility::deg2rad(0.0));
//     Eigen::Vector4d s(20.0, 6.0, math_utility::deg2rad(0.0), math_utility::deg2rad(0.0));
//     // Eigen::Vector4d s(14.0, 10.0, math_utility::deg2rad(0.0), math_utility::deg2rad(0.0));
//     Eigen::Vector4d g(0.0, 0.0, math_utility::deg2rad(90.0), math_utility::deg2rad(90.0));
    
//     Eigen::MatrixXd obss(139, 2);
//     int obs_index = 0;
//     for(int i=-25; i<= 25; i++){
//         obss.row(obs_index) << i*1.0, 15.0;
//         obs_index++;
//     }
//     for(int i=-25; i< -4; i++){
//         obss.row(obs_index) << i*1.0, 4.0;
//         obs_index++;
//     }
//     for(int i=-15; i< 4; i++){
//         obss.row(obs_index) << -4.0, i*1.0;
//         obs_index++;
//     }
//     for(int i=-15; i< 4; i++){
//         obss.row(obs_index) << 4.0, i*1.0;
//         obs_index++;
//     }
//     for(int i=4; i< 25; i++){
//         obss.row(obs_index) << i*1.0, 4.0;
//         obs_index++;
//     }
//     for(int i=-4; i< 4; i++){
//         obss.row(obs_index) << i*1.0, -15.0;
//         obs_index++;
//     }

//     // double radius = 3.1;
//     // PointCloud cloud;
//     // cloud.points =  obss;

//     // // // Create KDTree index
//     // KDTree kd_tree(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
//     // kd_tree.buildIndex();
//     // double cx = 5.810300232532859;
//     // double cy = 6.347042978647656;
//     //  std::cout<< " cx "<< cx<< " cy "<< cy << std::endl;
//     // std::vector<nanoflann::ResultItem<size_t, double>> indices_dists;
//     // // std::vector<double> dists;
//     // nanoflann::RadiusResultSet<double, size_t> result_set(radius*radius, indices_dists);
//     // // // resultSet.init(ret_indexes.data(), out_dists.data());
//     // Eigen::Vector2d query_point(cx, cy);
//     // kd_tree.findNeighbors(result_set, query_point.data(), nanoflann::SearchParameters());
//     // // std::cout << "points within radius " << radius << " from point (" << query_point[0] << ", " << query_point[1] << " " << ret_indexes.size() << "):\n";
//     // std::cout << "Points within radius " << radius << " of (" 
//     //           << query_point[0] << ", " << query_point[1] << "):\n";

//     // for (const auto& item : indices_dists) {
//     //     size_t idx = item.first;
//     //     double dist = item.second;
//     //     std::cout << "Index: " << idx << ", Distance: " << dist 
//     //               << ", Point: (" << cloud.points(idx, 0) << ", " 
//     //               << cloud.points(idx, 1) << ")\n";
//     // }


//     // Eigen::VectorXd col_min = obss.colwise().minCoeff(); 
//     // Eigen::VectorXd col_max = obss.colwise().maxCoeff();

//     // std::cout<< col_min << std::endl;
//     // std::cout<< col_max << std::endl;


//     planning::HybridPath path;
//     planning::TrailerHybridAStar trailer_hybrid_astar;
//     planning::TrailerLib trailer_lib;
//     planning::PlannerParams planner_params;
//     bool find_path = trailer_hybrid_astar.calc_hybrid_astar_path(s, g, obss, path);

//     plt::figure();
//     std::vector<double> ox, oy;
//     for(int i=0; i< obss.rows(); i++){
//         ox.push_back(obss.row(i)[0]);
//         oy.push_back(obss.row(i)[1]);
//     }

//     std::vector<double> path_x, path_y;
//     for(int i=0; i< path.poses.rows(); i++){
//         path_x.push_back(path.poses.row(i)[0]);
//         path_y.push_back(path.poses.row(i)[1]);
//     }
//     // plt::plot(ox, oy, ".r");
//     // plt::plot(path_x, path_y, ".y");
//     // plt::plot({s.x()}, {s.y()}, "bo"); // Start point
//     // plt::plot({g.x()}, {g.y()}, "go");    // End point
//     // // trailer_lib.plot_trailer(s.x(), s.y(), s[2], s[3], 0.0);
//     // // trailer_lib.plot_trailer(g.x(), g.y(), g[2], g[3], 0.0);
//     // plt::show();

//     Eigen::VectorXd x = path.poses.col(0);
//     Eigen::VectorXd y = path.poses.col(1);
//     Eigen::VectorXd yaw = path.poses.col(2);
//     Eigen::VectorXd yaw1 = path.poses.col(3);
//     Eigen::VectorXd direction = path.poses.col(4);
//     // std::cout<< direction.transpose() << std::endl;
//     double steer = 0.0;
//     for (size_t ii = 0; ii < x.size(); ++ii) {
//         plt::clf();
//         plt::plot(ox, oy, ".r");
//         plt::plot(path_x, path_y, ".y");

//         if (ii < x.size() - 1) {
//             double k = (yaw[ii + 1] - yaw[ii]) / planner_params.MOTION_RESOLUTION;
//             // std::cout<< " k: " << direction[ii] << "," << ii << std::endl;
//             if (direction[ii]< 0.0) {
//                 // std::cout<< " k: " << direction[ii] << "," << ii << std::endl;
//                 k *= -1.0;
//             }
//             steer = std::atan2(planner_params.WB * k, 1.0);  // Equivalent of Julia's `Base.atan(WB*k, 1.0)`
//         } else {
//             // std::cout<< " index: " << ii << std::endl;
//             steer = 0.0;
//         }
//         // std::cout<< " steer: " << steer << std::endl;
//         // std::cout<< x[ii]<< "," << y[ii]<< ","<< yaw[ii]<< ","<<yaw1[ii]<< "," <<steer<< std::endl;

//         trailer_lib.plot_trailer(x[ii], y[ii], yaw[ii], yaw1[ii], steer);

//         plt::grid(true);
//         plt::axis("equal");
//         plt::pause(0.0001);  // Small pause for animation
//     }
// }



// int main(int argc, char *argv[]) {
//     double x = 0.0;
//     double y = 0.0;
//     double yaw0 = math_utility::deg2rad(10.0);
//     double yaw1 = math_utility::deg2rad(-10.0);
//     planning::TrailerLib trailer_lib;
//     trailer_lib.plot_trailer(x, y, yaw0, yaw1, 0.0);
// }

// int main(int argc, char *argv[]) {

//     rs_paths::RSPaths rs_path;
//     Eigen::Vector3d s(3.0, 10.0, math_utility::deg2rad(40.0));
//     Eigen::Vector3d g(0.0, 1.0, math_utility::deg2rad(0.0));
//     double max_curvature = 0.1;
    
    
    
//     rs_paths::Path path = rs_path.calc_shortest_path(s, g, max_curvature);
//     std::vector<double> rc;
//     std::vector<double> rds;
//     std::cout<< "========1main 0" << path.poses.rows() << std::endl;
//     rs_path.calc_curvature(path.poses, rc, rds);
//     std::vector<double> path_short_x, path_short_y;
//     for(int i=0; i<path.poses.rows(); i++){
//         path_short_x.push_back(path.poses.row(i)[0]);
//         path_short_y.push_back(path.poses.row(i)[1]);
//     }
  

//     std::vector<std::vector<double>> bpath_x, bpath_y;
//     std::vector<rs_paths::Path> paths;
//     rs_path.calc_paths(s, g, max_curvature, paths);
//     for(auto path : paths){
//         std::vector<double> path_info_x, path_info_y;
//         for(int i=0; i<path.poses.rows(); i++){
//             path_info_x.push_back(path.poses.row(i)[0]);
//             path_info_y.push_back(path.poses.row(i)[1]);
//         }
//         bpath_x.push_back(path_info_x);
//         bpath_y.push_back(path_info_y);
//     }
//     // std::cout<< "========1main 1" << path.poses.rows() << std::endl;

//     // // First subplot
//     plt::figure();
//     for(int i=0; i< bpath_x.size(); i++){
//         plt::plot(bpath_x[i], bpath_y[i]);
//     }

//     plt::plot({s.x()}, {s.y()}, "bo"); // Start point
//     plt::plot({g.x()}, {g.y()}, "go");    // End point

//     plt::legend();
//     plt::grid(true);
//     plt::axis("equal");

//     // // // Second subplot for curvature
//     plt::figure();
//     plt::plot(path_short_x, path_short_y, "-r");
//     // plt::plot(rc, ".r");
//     plt::grid(true);
//     // plt::title("Curvature");

//     // // // Show all plots
//     plt::show();




//     return 0;
// }




// Parameters
constexpr double IQR_MULTIPLIER =  8.0;
constexpr double Z_SCORE_THRESHOLD =  8.0;
constexpr double GROUND_HEIGHT =  8.0; // Initial ground height
constexpr double THRESHOLD = 0.3;     // Altitude threshold
constexpr double ALPHA = 0.001;       // Learning rate for dynamic ground height
constexpr double VELOCITY_THRESHOLD = 100.0; // km/h
constexpr size_t MAX_HISTORY_SIZE = 500;
constexpr double MAX_ALT_THRESHOLD = 1.5; //if this exceed this value, it is consider as outlier without any checking 
constexpr double MAX_ALT_THRESHOLD_FOR_Z_SCORE = 1.2; //if this exceed this value, it is consider as soft-outlier that only for reduce distarbanuce for Z-score estimaiton



// Function to load data from a CSV file
vector<vector<double>> loadCSV(const string& filename) {
    vector<vector<double>> data;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return data;
    }
    string line;
    int rowNumber = 0;
    while (getline(file, line)) {
        rowNumber++;
        stringstream ss(line);
        string value;
        vector<double> row;

        // Parse each line as a row of doubles
        while (getline(ss, value, ',')) {
            try {
                row.push_back(stod(value));
            } catch (const std::invalid_argument& e) {
                cerr << "Error: Invalid data at row " << rowNumber << " - " << value << endl;
                row.clear();
                break; // Skip this row
            }
        }

        if (!row.empty() && row.size() == 4) { // Ensure valid data with 4 columns
            data.push_back(row);
        }
    }
    file.close();
    return data;
}

class OutlierDetector {
public:
    deque<double> historicalData;
    vector<int> outliers;
    vector<int> finalOutliers;
    vector<int> outlierIndices3D;
    vector<int> outlier_indices;
    vector<double> velocities;
    vector<vector<int>> outlier_periods;

    double dynamicGroundHeight = GROUND_HEIGHT;
    double previousTimeIndex = 0.0;
    vector<double> previousLatLonAlt;

    bool in_outlier{false};
    int start_index{0};
    int end_index{0};

    // Helper to compute 3D distance
    double compute3DDistance(double lat1, double lon1, double alt1, double lat2, double lon2, double alt2) {
        double R = 6371.0; // Earth radius in kilometers
        double dLat = (lat2 - lat1) * M_PI / 180.0;
        double dLon = (lon2 - lon1) * M_PI / 180.0;

        double a = sin(dLat / 2) * sin(dLat / 2) +
                cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) *
                sin(dLon / 2) * sin(dLon / 2);
        double c = 2 * atan2(sqrt(a), sqrt(1 - a));
        double horizontalDistance = R * c;

        double verticalDistance = fabs(alt2 - alt1) / 1000.0; // Convert meters to kilometers
        return sqrt(horizontalDistance * horizontalDistance + verticalDistance * verticalDistance);
    }

    void processRow(int index, const vector<double>& row) {
        double timeIndex = row[3];
        double lat = row[1];
        double lon = row[2];
        double alt = row[0];

        bool zScoreOutlier = false, iqrOutlier = false, groundHeightOutlier = false, velocityOutlier = false;

        if (!previousLatLonAlt.empty()) {
            // Compute 3D distance and velocity
            double distance3D = compute3DDistance(previousLatLonAlt[0], previousLatLonAlt[1], previousLatLonAlt[2],
                                                  lat, lon, alt);
            double timeDiff = (timeIndex - previousTimeIndex) / 3600.0; // Convert to hours
            double velocity = (timeDiff > 0) ? distance3D / timeDiff : 0.0;
            velocities.push_back(velocity);

            if (velocity > VELOCITY_THRESHOLD) {
                outlierIndices3D.push_back(index);
                velocityOutlier = true;
            }

            // Statistical computations
            double meanAltitude = 0, stdAltitude = 0, Q1 = 0, Q3 = 0;
            if (historicalData.size()>1) {
                VectorXd data(historicalData.size());
                for (size_t i = 0; i < historicalData.size(); ++i)
                    data(i) = historicalData[i];

                meanAltitude = data.mean();
                stdAltitude = sqrt((data.array() - meanAltitude).square().mean());

                std::nth_element(data.data(), data.data() + data.size() / 4, data.data() + data.size());
                Q1 = data(data.size() / 4);
                std::nth_element(data.data(), data.data() + 3 * data.size() / 4, data.data() + data.size());
                Q3 = data(3 * data.size() / 4);

                // std::cout<< " data: " << data.transpose() << std::endl; 

                double IQR = Q3 - Q1;
                double lowerBound = Q1 - IQR_MULTIPLIER * IQR;
                double upperBound = Q3 + IQR_MULTIPLIER * IQR;

                if (alt < lowerBound || alt > upperBound) {
                    outliers.push_back(index);
                    iqrOutlier = true;
                }
                // std::cout<< "q1 "<< Q1 << " q2: " << Q3 << " meanAltitude " << meanAltitude << " stdAltitude: " << stdAltitude << " iqrOutlier " << iqrOutlier << std::endl;
                if (stdAltitude > 0) {
                    double zScore = (alt - meanAltitude) / stdAltitude;
                    if (fabs(zScore) > Z_SCORE_THRESHOLD) {
                        outliers.push_back(index);
                        zScoreOutlier = true;
                    }
                }
            }

            bool is_outlier = (fabs(alt - dynamicGroundHeight) > THRESHOLD);
            if(is_outlier){
                if(in_outlier==false){
                    start_index = index;
                    in_outlier = true;
                }
                groundHeightOutlier = true;
            }else{
                if(in_outlier){
                    end_index = index-1;
                    vector<int> section = {start_index, end_index};
                    outlier_periods.push_back(section);
                    in_outlier = false;
                }
                groundHeightOutlier = false;
                dynamicGroundHeight = (1 - ALPHA) * dynamicGroundHeight + ALPHA * alt;
            }

            bool alt_outlier_extrem = (abs(GROUND_HEIGHT - alt) > MAX_ALT_THRESHOLD);
            bool can_be_consider_as_outlier = (zScoreOutlier + iqrOutlier + groundHeightOutlier + velocityOutlier) >= 2;
            if(alt_outlier_extrem){
                finalOutliers.push_back(index);
            }else if(can_be_consider_as_outlier){
                finalOutliers.push_back(index);
            }
        }

        // Update state
        previousTimeIndex = timeIndex;
        previousLatLonAlt = {lat, lon, alt};

        bool soft_outlier = abs(GROUND_HEIGHT - alt) < MAX_ALT_THRESHOLD_FOR_Z_SCORE;
        if (soft_outlier) {
            historicalData.push_back(alt);
        }

        if (historicalData.size() >= MAX_HISTORY_SIZE) {
            historicalData.pop_front();  // Remove oldest entry
        }

    }
};


void plotOutliers(
    const std::vector<double>& longitudes,
    const std::vector<double>& latitudes,
    const std::vector<int>& outlierIndices,
    const std::vector<int>& outlierIndices3D,
    const std::vector<int>& finalOutliers,
    const std::vector<double>& altitudes,
    const std::vector<double>& velocities,
    double groundHeight,
    const std::vector<int>& altitudeOutliers) {

    // Convert indices to outlier coordinates
    std::vector<double> outlierLongitudes, outlierLatitudes, velocityOutlierLongitudes, velocityOutlierLatitudes, finalOutlierLongitudes, finalOutlierLatitudes;
    for (int idx : altitudeOutliers) {
        outlierLongitudes.push_back(longitudes[idx]);
        outlierLatitudes.push_back(latitudes[idx]);
    }
    for (int idx : outlierIndices3D) {
        velocityOutlierLongitudes.push_back(longitudes[idx]);
        velocityOutlierLatitudes.push_back(latitudes[idx]);
    }
    for (int idx : finalOutliers) {
        finalOutlierLongitudes.push_back(longitudes[idx]);
        finalOutlierLatitudes.push_back(latitudes[idx]);
    }

    // Plot locations of outliers
    plt::figure(1);
    plt::scatter(longitudes, latitudes, 10, {{"c", "blue"}, {"label", "Normal Points"}});
    plt::scatter(outlierLongitudes, outlierLatitudes, 20, {{"c", "red"}, {"label", "Outliers"}});
    plt::xlabel("Longitude");
    plt::ylabel("Latitude");
    plt::title("Location of Outliers");
    plt::legend();
    plt::save("v1.png");

    // // Plot velocity profile
    plt::figure(2);
    plt::plot(velocities, {{"color", "green"}, {"label", "3D Velocity (km/h)"}});
    vector<double> updated_vel;
    for (int idx : outlierIndices3D) {
        updated_vel.push_back(velocities[idx]);
    }
    plt::scatter(outlierIndices3D, updated_vel, 20, {{"c", "red"}, {"label", "High Velocity Outliers"}});
    // plt::axhline(groundHeight, 0, velocities.size(), {{"color", "orange"}, {"linestyle", "--"}, {"label", "Initial Ground Height"}});
    plt::plot(altitudes, {{"color", "red"}, {"label", "Altitude"}});
    plt::xlabel("Index");
    plt::ylabel("Velocity / Altitude");
    plt::title("Velocity Outliers");
    plt::legend();
    plt::save("v2.png");


    // // Plot altitude with marked outliers
    plt::figure(3);
    plt::plot(altitudes, {{"color", "blue"}, {"label", "Altitude"}});
    vector<double> altitudeOutliers_vel;
    for (int idx : altitudeOutliers) {
        altitudeOutliers_vel.push_back(altitudes[idx]);
    }
    plt::scatter(altitudeOutliers, altitudeOutliers_vel, 20, {{"c", "red"}, {"label", "Altitude Outliers"}});
    // plt::axhline(groundHeight, 0, altitudes.size(), {{"color", "orange"}, {"linestyle", "--"}, {"label", "Initial Ground Height"}});
    plt::xlabel("Index");
    plt::ylabel("Altitude");
    plt::title("Z-score threshold and IQR multiplier");
    plt::legend();
    plt::save("v3.png");
    

    // Plot locations of Z-score/IQR outliers
    plt::figure(4);
    plt::scatter(longitudes, latitudes, 10, {{"c", "blue"}, {"label", "Normal Points"}});
    plt::scatter(outlierLongitudes, outlierLatitudes, 20, {{"c", "red"}, {"label", "Outliers Z-score threshold and IQR multiplier"}});
    plt::legend();
    plt::save("v4.png");

    // Plot velocity outliers on map
    plt::figure(5);
    plt::scatter(longitudes, latitudes, 10, {{"c", "blue"}, {"label", "Normal Points"}});
    plt::scatter(velocityOutlierLongitudes, velocityOutlierLatitudes, 20, {{"c", "red"}, {"label", "Outliers Velocity"}});
    plt::legend();
    plt::save("v5.png");

    // // Plot final outliers
    plt::figure(6);
    plt::scatter(longitudes, latitudes, 10, {{"c", "blue"}, {"label", "Normal Points"}});
    plt::scatter(finalOutlierLongitudes, finalOutlierLatitudes, 20, {{"c", "red"}, {"label", "Outliers Final"}});
    plt::legend();
    plt::save("v6.png");

    plt::show();
}

// void plotResults(const vector<double>& altitudes, const vector<int>& outliers, const vector<double>& velocities) {
//     vector<double> indices(altitudes.size());
//     for (size_t i = 0; i < indices.size(); ++i)
//         indices[i] = i;

//     // // Plot Altitudes
//     plt::figure();
//     plt::plot(indices, altitudes, "r");
//     vector<double> outlierAlts;
//     vector<double> outlierIndices;
//     for (int idx : outliers) {
//         outlierAlts.push_back(altitudes[idx]);
//         outlierIndices.push_back(idx);
//     }
//     std::vector<double> colors = {0.1, 0.3, 0.5, 0.7, 0.9};
//     plt::scatter(outlierIndices, outlierAlts, 20.0);
//     plt::xlabel("Index");
//     plt::ylabel("Altitude");
//     plt::title("Altitude Data with Outliers");
//     plt::legend();
//     plt::grid(true);

    

//     // if (indices.size() != velocities.size()) {
//     //     std::cerr << "Error: x and y vectors must have the same size for plotting." << std::endl;
//     //     std::cerr << "Size of x: " << indices.size() << ", Size of y: " << velocities.size() << std::endl;
//     //     return;
//     // }

//     try {
//         // Plot Velocities
//         plt::figure();
//         plt::plot(indices, velocities, "g-");
//         plt::xlabel("Index");
//         plt::ylabel("Velocity (km/h)");
//         plt::title("Velocity Profile");
//         plt::grid(true);
//     } catch (const std::exception& e) {
//         std::cerr << "Exception during plotting: " << e.what() << std::endl;
//     }
//     plt::show();
// }

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <CSV file>" << endl;
        return 1;
    }

    string filename = argv[1];
    vector<vector<double>> data = loadCSV(filename);

    if (data.empty()) {
        cerr << "Error: No data loaded from file " << filename << endl;
        return 1;
    }

    // Define parameters for the outlier detector.
    outliear_detector::OutlierDetectionParams params;
    params.iqr_multiplier = 9.0;
    params.z_score_threshold = 9.0;
    params.ground_height = 9.0;
    params.altitude_threshold = 0.3;
    params.alpha = 0.001;
    params.velocity_threshold = 100.0;  // km/h
    params.max_history_size = 500;
    params.max_altitude_threshold = 1.5;
    params.max_altitude_threshold_for_z_score = 1.2;

    // Instantiate the outlier detector with the parameters.
    outliear_detector::OutlierDetector detector;
    detector.initParam(params);

    vector<double> longitudes;
    vector<double> latitudes;
    vector<double> altitudes;
    for (size_t i = 0; i < data.size(); ++i) {
    // for (size_t i = 0; i < 20; ++i) {
        altitudes.push_back(data[i][0]);
        latitudes.push_back(data[i][1]);
        longitudes.push_back(data[i][2]);
        detector.ProcessRow(i, data[i]);
    }

    // Loop through the start and end indices of outlier periods
    for (const auto& period : detector.outlier_periods_) {
        size_t start = period[0];
        size_t end = period[1];
        for (size_t i = start; i <= end; ++i) {
            detector.outlier_indices_.push_back(i);
        }
    }

    detector.velocities_.push_back(0.0);

    plotOutliers(longitudes, latitudes, detector.outlier_indices_, detector.outlier_indices_3d_, detector.final_outliers_
            , altitudes, detector.velocities_, params.ground_height, detector.outliers_);

    // plotResults(altitudes, detector.finalOutliers, detector.velocities);

    return 0;
}


