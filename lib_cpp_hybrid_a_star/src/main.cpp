#include <iostream>
#include <Eigen/Dense>

#include "lib_cpp_hybrid_a_star/rs_paths.hpp"
#include "lib_cpp_hybrid_a_star/grid_a_star.hpp"
#include "lib_cpp_hybrid_a_star/trailerlib.hpp"
#include "lib_cpp_hybrid_a_star/trailer_hybrid_a_star.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

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

int main(int arc, char *argv[]){

    // Eigen::Vector4d s(-10.0, 6.0, math_utility::deg2rad(0.0), math_utility::deg2rad(0.0));
    Eigen::Vector4d s(20.0, 6.0, math_utility::deg2rad(0.0), math_utility::deg2rad(0.0));
    // Eigen::Vector4d s(14.0, 10.0, math_utility::deg2rad(0.0), math_utility::deg2rad(0.0));
    Eigen::Vector4d g(0.0, 0.0, math_utility::deg2rad(90.0), math_utility::deg2rad(90.0));
    
    Eigen::MatrixXd obss(139, 2);
    int obs_index = 0;
    for(int i=-25; i<= 25; i++){
        obss.row(obs_index) << i*1.0, 15.0;
        obs_index++;
    }
    for(int i=-25; i< -4; i++){
        obss.row(obs_index) << i*1.0, 4.0;
        obs_index++;
    }
    for(int i=-15; i< 4; i++){
        obss.row(obs_index) << -4.0, i*1.0;
        obs_index++;
    }
    for(int i=-15; i< 4; i++){
        obss.row(obs_index) << 4.0, i*1.0;
        obs_index++;
    }
    for(int i=4; i< 25; i++){
        obss.row(obs_index) << i*1.0, 4.0;
        obs_index++;
    }
    for(int i=-4; i< 4; i++){
        obss.row(obs_index) << i*1.0, -15.0;
        obs_index++;
    }

    // double radius = 3.1;
    // PointCloud cloud;
    // cloud.points =  obss;

    // // // Create KDTree index
    // KDTree kd_tree(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    // kd_tree.buildIndex();
    // double cx = 5.810300232532859;
    // double cy = 6.347042978647656;
    //  std::cout<< " cx "<< cx<< " cy "<< cy << std::endl;
    // std::vector<nanoflann::ResultItem<size_t, double>> indices_dists;
    // // std::vector<double> dists;
    // nanoflann::RadiusResultSet<double, size_t> result_set(radius*radius, indices_dists);
    // // // resultSet.init(ret_indexes.data(), out_dists.data());
    // Eigen::Vector2d query_point(cx, cy);
    // kd_tree.findNeighbors(result_set, query_point.data(), nanoflann::SearchParameters());
    // // std::cout << "points within radius " << radius << " from point (" << query_point[0] << ", " << query_point[1] << " " << ret_indexes.size() << "):\n";
    // std::cout << "Points within radius " << radius << " of (" 
    //           << query_point[0] << ", " << query_point[1] << "):\n";

    // for (const auto& item : indices_dists) {
    //     size_t idx = item.first;
    //     double dist = item.second;
    //     std::cout << "Index: " << idx << ", Distance: " << dist 
    //               << ", Point: (" << cloud.points(idx, 0) << ", " 
    //               << cloud.points(idx, 1) << ")\n";
    // }


    // Eigen::VectorXd col_min = obss.colwise().minCoeff(); 
    // Eigen::VectorXd col_max = obss.colwise().maxCoeff();

    // std::cout<< col_min << std::endl;
    // std::cout<< col_max << std::endl;


    planning::HybridPath path;
    planning::TrailerHybridAStar trailer_hybrid_astar;
    planning::TrailerLib trailer_lib;
    planning::PlannerParams planner_params;
    bool find_path = trailer_hybrid_astar.calc_hybrid_astar_path(s, g, obss, path);

    plt::figure();
    std::vector<double> ox, oy;
    for(int i=0; i< obss.rows(); i++){
        ox.push_back(obss.row(i)[0]);
        oy.push_back(obss.row(i)[1]);
    }

    std::vector<double> path_x, path_y;
    for(int i=0; i< path.poses.rows(); i++){
        path_x.push_back(path.poses.row(i)[0]);
        path_y.push_back(path.poses.row(i)[1]);
    }
    // plt::plot(ox, oy, ".r");
    // plt::plot(path_x, path_y, ".y");
    // plt::plot({s.x()}, {s.y()}, "bo"); // Start point
    // plt::plot({g.x()}, {g.y()}, "go");    // End point
    // // trailer_lib.plot_trailer(s.x(), s.y(), s[2], s[3], 0.0);
    // // trailer_lib.plot_trailer(g.x(), g.y(), g[2], g[3], 0.0);
    // plt::show();

    Eigen::VectorXd x = path.poses.col(0);
    Eigen::VectorXd y = path.poses.col(1);
    Eigen::VectorXd yaw = path.poses.col(2);
    Eigen::VectorXd yaw1 = path.poses.col(3);
    Eigen::VectorXd direction = path.poses.col(4);
    // std::cout<< direction.transpose() << std::endl;
    double steer = 0.0;
    for (size_t ii = 0; ii < x.size(); ++ii) {
        plt::clf();
        plt::plot(ox, oy, ".r");
        plt::plot(path_x, path_y, ".y");

        if (ii < x.size() - 1) {
            double k = (yaw[ii + 1] - yaw[ii]) / planner_params.MOTION_RESOLUTION;
            // std::cout<< " k: " << direction[ii] << "," << ii << std::endl;
            if (direction[ii]< 0.0) {
                // std::cout<< " k: " << direction[ii] << "," << ii << std::endl;
                k *= -1.0;
            }
            steer = std::atan2(planner_params.WB * k, 1.0);  // Equivalent of Julia's `Base.atan(WB*k, 1.0)`
        } else {
            // std::cout<< " index: " << ii << std::endl;
            steer = 0.0;
        }
        // std::cout<< " steer: " << steer << std::endl;
        // std::cout<< x[ii]<< "," << y[ii]<< ","<< yaw[ii]<< ","<<yaw1[ii]<< "," <<steer<< std::endl;

        trailer_lib.plot_trailer(x[ii], y[ii], yaw[ii], yaw1[ii], steer);

        plt::grid(true);
        plt::axis("equal");
        plt::pause(0.0001);  // Small pause for animation
    }
}



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