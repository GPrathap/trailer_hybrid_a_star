#include "lib_cpp_hybrid_a_star/grid_a_star.hpp"


namespace grid_search
{

    GridAStar::GridAStar(){

    }


    void GridAStar::calc_obstacle_map(Eigen::MatrixXd& obses, double& reso, double& vr){

        Eigen::VectorXd col_min = obses.colwise().minCoeff(); 
        Eigen::VectorXd col_max = obses.colwise().maxCoeff();

        minx_ = std::round(col_min[0]);
        miny_ = std::round(col_min[1]);
        maxx_ = std::round(col_max[0]);
        maxy_ = std::round(col_max[1]);
        std::cout<< "min: " << col_min.transpose() << std::endl;
        std::cout<< "max: " << col_max.transpose() << std::endl;
        xwidth_ = maxx_ - minx_;
        ywidth_ = maxy_ - miny_;
        obs_map.resize(xwidth_, ywidth_);
        obs_map.setZero();

        PointCloud cloud;
        cloud.points = obses;

        KDTree kd_tree(2 , cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
        kd_tree.buildIndex();
        const size_t num_results = 1;
        std::vector<size_t> ret_indexes(num_results); 
        std::vector<double> out_dist_sqr(num_results);

        nanoflann::KNNResultSet<double> resultSet(num_results);
        resultSet.init(&ret_indexes[0], &out_dist_sqr[0]);

        for(int ix=0; ix<xwidth_; ix++ ){
            int x = ix + minx_;
            for(int iy=0; iy<ywidth_; iy++ ){
                int y = iy + miny_;
                const double query_pt[2] = {x*1.0, y*1.0};
                kd_tree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(10));
                // std::cout << "K nearest neighbors for point (" << query_pt[0] << ", " << query_pt[1] << "):" << std::endl;
                double dis = std::sqrt(out_dist_sqr[0]);
                if (dis <= vr/reso){
                    obs_map(ix,iy) = 1;
                } 
            }
        }
    }

    int GridAStar::calc_index(Node node){
        return (int)((node.pose.y() - miny_)*xwidth_ + (node.pose.x() - minx_));
    }

    double GridAStar::calc_cost(Node n, Node ngoal){
        return (n.cost + h(n.pose.x() - ngoal.pose.x(), n.pose.y() - ngoal.pose.y()));
    }

    double GridAStar::h(int x, int y){
        return sqrt( pow(x,2) + pow(y,2));
    }

    Eigen::MatrixXd GridAStar::get_motion_model(){
        Eigen::MatrixXd motion(24, 3);
        motion<< 1, 0, 1, 0, 1, 1, -1, 0, 1, 0, -1, 1, -1, -1, sqrt(2), -1, 1, sqrt(2), 1, -1, sqrt(2), 1, 1, sqrt(2);
        return motion;
    }


    void GridAStar::calc_dist_policy(Eigen::Vector2d s, Eigen::Vector2d g
                                                , Eigen::MatrixXd obses, double reso, double vr){
        
        Eigen::Vector2i goal_pose (g.x()/reso, g.y()/reso);
        Node ngoal(goal_pose, 0.0, -1);  
        for(int i=0; i<obses.rows(); i++){
            obses.row(i) = obses.row(i)/reso;
        }   
        calc_obstacle_map(obses, reso, vr);

        std::unordered_map<int, Node*> open;
        std::unordered_map<int, Node*> closed;
        open[calc_index(ngoal)] = &ngoal;

        Eigen::MatrixXd motion = get_motion_model();
        int nmotion = motion.rows();
        // std::priority_queue<int> pq;
        std::priority_queue<Node, std::vector<Node>, CompareNode> pq;
        pq.push(ngoal);
        while( true ){
            if(open.empty()){
                std::cerr << " Error: No open set " << std::endl;
                break;
            }
            Node next_node = pq.top();
            int c_id = next_node.pind;
            Node* current = open[c_id];
        }
    }

    bool GridAStar::verify_node(Node node){
        if ((node.pose.x() - minx_) >= xwidth_){
            return false;
        } else if (( node.pose.x() - minx_) <= 0 ){
            return false;
        } 
        if ((node.pose.y() - miny_) >= ywidth_){
            return false;
        } else if (( node.pose.y() - miny_) <= 0 ){
            return false;
        }       
        if (obs_map(node.pose.x()-minx_, node.pose.y()-miny_) > 0){
            return false;
        }
        return true;                    
    }

//     void GridAStar::calc_policy_map(closed){
//         pmap.resize(xwidth_, ywidth_);
//         pmap.setConstant(std::numeric_limits<double>::infinity());
//          for n in values(closed)
//         pmap[n.x-minx, n.y-miny] = n.cost
//     end
//     # println(pmap)

//     return pmap
// end
//     }

}