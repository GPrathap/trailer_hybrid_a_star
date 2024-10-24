#include <iostream>
#include <Eigen/Dense>

#include "lib_cpp_hybrid_a_star/rs_paths.hpp"
#include "lib_cpp_hybrid_a_star/grid_a_star.hpp"
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;


int main(int argc, char *argv[]){
    Eigen::Vector2d s(10.0, 10.0);
    Eigen::Vector2d g(50.0, 50.0);

    Eigen::MatrixXd obss(320, 2);
    int obs_index = 0;
    for(int i=0; i< 60; i++){
        obss.row(obs_index) << i*1.0, 0.0;
        obs_index++;
    }
    for(int i=0; i< 60; i++){
        obss.row(obs_index) << 60.0, i*1.0;
        obs_index++;
    }
    for(int i=0; i< 60; i++){
        obss.row(obs_index) << i*1.0, 60.0;
        obs_index++;
    }
    for(int i=0; i< 60; i++){
        obss.row(obs_index) << 0.0, i*1.0;
        obs_index++;
    }
    for(int i=0; i< 40; i++){
        obss.row(obs_index) << 20.0, i*1.0;
        obs_index++;
    }
    for(int i=0; i< 40; i++){
        obss.row(obs_index) << 40.0, 60.0-i*1.0;
        obs_index++;
    }

    double VEHICLE_RADIUS = 5.0;
    double GRID_RESOLUTION = 1.0;

    grid_search::GridAStar grid_a_star;
    grid_a_star.calc_dist_policy(s, g, obss, GRID_RESOLUTION, VEHICLE_RADIUS);
}



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