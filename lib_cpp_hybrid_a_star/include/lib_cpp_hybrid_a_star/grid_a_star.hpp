#ifndef GRID_SEARCH_A_STAR
#define GRID_SEARCH_A_STAR

#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <initializer_list>
#include <Eigen/Dense>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <thread>
#include <iostream>
#include <cstdlib> 
#include <execution>
#include "math_utility.hpp"
#include <nanoflann.hpp>
#include <queue>

namespace grid_search
{
    using namespace math_utility;
    template <typename T>
    void printVector(const std::vector<T>& vec) {
        std::cout << "[ ";
        for (const auto& elem : vec) {
            std::cout << elem << " ";
        }
        std::cout << "]" << std::endl;
    }

    typedef struct Node
    {
        Eigen::Vector2i pose;
        double cost;
        int pind;
        
        Node(Eigen::Vector2i position, double cost_, int id){
            pose = position;
            cost = cost_;
            pind = id;
        }

    };

    struct CompareNode {
        bool operator()(const Node& a, const Node& b) {
            return a.cost > b.cost;
        }
    };
    

    struct PointCloud {
        Eigen::MatrixXd points;

        inline size_t kdtree_get_point_count() const { return points.rows(); }

        inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
            return points(idx, dim); // accessing the coordinate
        }

        template <class BBOX> bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
    };

    using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>;


    class GridAStar{
        public:
            GridAStar();

            void calc_obstacle_map(Eigen::MatrixXd& obses, double& reso, double& vr);
            void calc_dist_policy(Eigen::Vector2d s, Eigen::Vector2d g
                            , Eigen::MatrixXd obses, double reso, double vr);
            int calc_index(Node node);
            double h(int x, int y);
            double calc_cost(Node n, Node ngoal);
            Eigen::MatrixXd get_motion_model();
            bool verify_node(Node node);

        private:
            Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> obs_map;
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> pmap;
            int minx_, miny_, maxx_, maxy_;
            int xwidth_, ywidth_;
    };   
}


#endif 