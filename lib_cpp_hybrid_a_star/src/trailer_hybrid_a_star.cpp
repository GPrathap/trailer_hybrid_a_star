#include "lib_cpp_hybrid_a_star/trailer_hybrid_a_star.hpp"

namespace planning
{
    TrailerHybridAStar::TrailerHybridAStar(){
    }

    void TrailerHybridAStar::calc_config(Eigen::MatrixXd& obses, double xyreso, double yawreso){
        
        // Eigen::VectorXd col_min = obses.colwise().minCoeff(); 
        // Eigen::VectorXd col_max = obses.colwise().maxCoeff();

        // col_min = col_min.array() - PlannerParams::EXTEND_AREA;
        // col_max = col_max.array() + PlannerParams::EXTEND_AREA;

        // int minx = std::round(col_min[0]/xyreso);
        // int miny = std::round(col_min[1]/xyreso);
        // int maxx = std::round(col_max[0]/xyreso);
        // int maxy = std::round(col_max[1]/xyreso);
        int minx = -15;
        int miny = -10;
        int maxx = 15;
        int maxy = 10;
        
        int xw = maxx - minx;
        int yw = maxy - miny;

        int minyaw = std::round(-M_PI/yawreso) - 1;
        int maxyaw = std::round(M_PI/yawreso);
        int yaww = std::round((maxyaw - minyaw));

        int minyawt = minyaw;
        int maxyawt = maxyaw;
        int yawtw = yaww;

        config_ = Config(minx, miny, minyaw, minyawt, maxx, maxy, maxyaw, maxyawt, xw, yw, yaww, yawtw, xyreso, yawreso);

    }

    void TrailerHybridAStar::calc_holonomic_with_obstacle_heuristic(HybridNode& gnode
                                    , Eigen::MatrixXd& obses, double xyreso){
        Eigen::Vector2d goal_pose = gnode.poses.row(gnode.poses.rows()-1).head(2);
        
        grid_a_star_.calc_dist_policy(goal_pose, obses, xyreso, 1.0, h_dp);
    }

    int TrailerHybridAStar::calc_index(HybridNode& node){
        
        int ind = (node.yawind - config_.minyaw)*config_.xw*config_.yw
                    +(node.yind - config_.miny)*config_.xw + (node.xind - config_.minx);
        // 4D grid
        int rows = node.poses.rows();
        if( rows>0 ){
            double yaw1 = node.poses.row(node.poses.rows()-1)[3];
            int yaw1ind = std::round(yaw1/config_.yawreso);
            ind += (yaw1ind - config_.minyawt) *config_.xw*config_.yw*config_.yaww;
        } else{
            std::cout<< " Error(calc_index): " << std::endl;
        }
        if(ind <= 0){
            std::cout<< " Error(calc_index): " << std::endl;
        }
        return ind;
    }

    void TrailerHybridAStar::calc_motion_inputs(std::vector<double>& u, std::vector<double>& d){
        std::vector<double> up;
        double increment = PlannerParams::MAX_STEER/PlannerParams::N_STEER;
        for (double i=increment; i <= PlannerParams::MAX_STEER+0.002; i+=increment) {
            up.push_back(i);
        }
        u.clear();
        d.clear();

        u.push_back(0.0);
        u.insert(u.end(), up.begin(), up.end());
        std::transform(up.begin(), up.end(), std::back_inserter(u), [](double i) { return -i; });

        int len_u = u.size();
        d.resize(2*len_u);
        double* d_ptr = d.data();
        double* mid_ptr = d_ptr + len_u;
        for (int i = 0; i < len_u; ++i) {
            *(d_ptr++) = 1.0;
            *(mid_ptr++) = -1.0;
        }
        u.reserve(2 * u.size());
        u.insert(u.end(), u.begin(), u.end());
    }

    double TrailerHybridAStar::calc_rs_path_cost(rs_paths::Path& rspath, Eigen::VectorXd& yaw1){
        double cost = 0.0;
        for(int i=0; i<rspath.lengths.size(); i++){
            double l = rspath.lengths[i];
            cost += (rspath.lengths[i] >= 0) ? l : std::abs(l) * PlannerParams::BACK_COST;
        }
        // std::cout<< " 1cost " << cost << std::endl;
        // swich back penalty
        for(int i=0; i<rspath.lengths.size()-1; i++){
            if(rspath.lengths[i] * rspath.lengths[i+1] < 0.0){
                cost += PlannerParams::SB_COST;
            }
        }
        // std::cout<< " 2cost " << cost << std::endl;
        // steer penalyty
        for(int i=0; i<rspath.ctypes.size(); i++){
            if(rspath.ctypes[i] != "S"){
                cost += PlannerParams::STEER_COST*std::abs(PlannerParams::MAX_STEER);
            }
        }
        // std::cout<< " 3cost " << cost << std::endl;
        // steer change penalty
        // calc steer profile
        int nctypes = rspath.ctypes.size();
        std::vector<double> ulist(nctypes, 0.0);
        for(int i=0; i<nctypes; i++){
            if(rspath.ctypes[i] == "R"){
                ulist[i] = - PlannerParams::MAX_STEER;
            }else if(rspath.ctypes[i] == "L"){
                ulist[i] = PlannerParams::MAX_STEER;
            }
        }
        // std::cout<< " 4cost " << cost << std::endl;
        for(int i=0; i<rspath.ctypes.size()-1; i++){
            cost += PlannerParams::STEER_CHANGE_COST*std::abs(ulist[i+1] - ulist[i]);
        }
        // std::cout<< " 5cost " << cost << std::endl;
        // std::cout<< " 4cost yaw1 " << rspath.poses.col(2).transpose() << std::endl;
        // std::cout<< " 4cost yaw1 " << yaw1.transpose() << std::endl;

        Eigen::VectorXd yaw_diff = rspath.poses.col(2) - yaw1;

        // wrap the differences to the range [-π, π]
        // std::cout<< " 5cost " << yaw_diff.transpose() << std::endl;
        // std::cout<< " rspath.poses.col(2) " << rspath.poses.col(2) << std::endl;
        // std::cout<< " rspath.poses.col(2) " << yaw1.col(2) << std::endl;
        Eigen::VectorXd wrapped_yaw_diff = yaw_diff.unaryExpr([](double x) {
            return pi_to_pi(x);
        });

        // calculate the absolute values
        Eigen::VectorXd abs_yaw_diff = wrapped_yaw_diff.array().abs();

        // sum the absolute values
        double sum_abs_yaw_diff = abs_yaw_diff.sum();
        cost += PlannerParams::JACKKNIF_COST * sum_abs_yaw_diff;
        // std::cout<< " 6cost " << cost << " " << sum_abs_yaw_diff << std::endl;
        return cost;
    }

    bool TrailerHybridAStar::update_node_with_analystic_expantion(HybridNode& current
                , HybridNode& ngoal,  Eigen::MatrixXd& obses, grid_search::KDTree& kdtree
                , double gyaw1, HybridNode& updated_node){
            
            rs_paths::Path apath;
            
            bool find_path = analystic_expantion(current, ngoal, obses, kdtree, apath);
            if(find_path){
                // std::cout<< "----analystic_expantion---- " << find_path << std::endl;
                Eigen::VectorXd steps = apath.poses.col(3).array()*PlannerParams::MOTION_RESOLUTION;
                Eigen::VectorXd current_pose = current.poses.row(current.poses.rows()-1);
                Eigen::VectorXd yaw1;
                // std::cout<< "steps: " << steps.transpose() << std::endl;
                // std::cout<< "yaw: " << apath.poses.col(2).transpose() << std::endl;
                bool can_estimate = trailerlib_.calc_trailer_yaw_from_xyyaw(apath.poses, current_pose[3], steps, yaw1);
                // if(!can_estimate){
                //     return false;
                // }
                // std::cout<< "id: "<< current.pind << " yaw1 " << yaw1[yaw1.size()-1] << " gyaw1 " << gyaw1  << " pi_to_pi " << math_utility::pi_to_pi(yaw1[yaw1.size()-1] - gyaw1) << " yaw1 "<< current_pose[3] << std::endl;
                if (std::abs(math_utility::pi_to_pi(yaw1[yaw1.size()-1] - gyaw1)) >= PlannerParams::GOAL_TYAW_TH){
                    return false;
                }
                double fcost = current.cost + calc_rs_path_cost(apath, yaw1);
                int fpind = calc_index(current);

                // for(int i=1; i<apath.poses.rows(); i++){
                //     double d = apath.poses.row(i)[3];
                //     if(d > 0){

                //     }
                // }
                // std::cout<< "----analystic_expantion---- fcost " << fcost << std::endl;
                
                // to get rows from index 1 to the end (2nd to last rows)
                Eigen::MatrixXd updated_poses(apath.poses.rows() - 1, apath.poses.cols()+1);

                std::cout<< "----updated_poses " << updated_poses.rows() << "," << updated_poses.cols() << " " << yaw1.size() << std::endl;

                updated_poses.block(0, 0, apath.poses.rows() - 1, apath.poses.cols()) = apath.poses.block(1, 0, apath.poses.rows() - 1, apath.poses.cols());
                updated_poses.col(4) = updated_poses.col(3);
                updated_poses.col(3) = yaw1.head(yaw1.size()-1);
                
                // std::cout<< updated_poses.col(3).transpose() << std::endl;
                
                double  fsteer = 0.0;
                HybridNode est_node(current.xind, current.yind, current.yawind
                            , current.direction, updated_poses, fsteer, fcost, fpind);
                updated_node = est_node;
                return true;
            }
        return false;
    }


    void TrailerHybridAStar::get_final_path(std::unordered_map<int, HybridNode>& closed
                                                    , HybridNode& ngoal, HybridNode& nstart, HybridPath& path){
        
        std::cout<< "directions goal "<< ngoal.poses.rows() << " " <<  ngoal.poses.cols() << std::endl;
        Eigen::MatrixXd g_poses = ngoal.poses.colwise().reverse();
        // std::cout<< " g_poses directions "<< g_poses << std::endl;
        std::cout<< "directions "<< g_poses.rows() << " " <<  g_poses.cols() << std::endl;

        int nid = ngoal.pind;
        int finalcost = ngoal.cost;
        std::vector<Eigen::MatrixXd> ref_path;
        ref_path.push_back(g_poses);
        int total_rows = g_poses.rows();
        int total_cols = 5;
        // std::cout<< "Final path info "<< nid << " finalcost " << finalcost << std::endl;
        
        
        while(true){
            // std::cout<< "------------d1------------" << std::endl;
            if (closed.find(nid) == closed.end()) {
                std::cout<< nid << " cant find the requested node" << std::endl;
                break;
            }
            // std::cout<< "------------d2------------" << std::endl;
            HybridNode n = closed[nid];
            // std::cout<< " inter directions "<< n.poses << std::endl;
            // std::cout<< "------------d3------------" << n.poses << std::endl;
            Eigen::MatrixXd n_poses = n.poses.colwise().reverse();
            // std::cout<< "------------d4------------" << std::endl;
            ref_path.push_back(n_poses);
            total_rows += n_poses.rows();
            nid = n.pind;
            if (is_same_grid(n, nstart)){
                break;
            }
        }
        // std::cout<< "------------d5------------ " << total_rows << std::endl;
        Eigen::MatrixXd final_path(total_rows, total_cols);
        int current_row = 0;
        for (const auto& mat : ref_path) {
            // std::cout<< "------------d5----1-------- " << mat.rows() << " " << mat.cols() <<  " " << total_cols << std::endl;
            final_path.block(current_row, 0, mat.rows(), mat.cols()) = mat;
            // std::cout<< "------------d5-----2------" << current_row << std::endl;
            current_row += mat.rows();
        }
        // std::cout<< "------------d6------------" << final_path << " " << final_path.cols() << std::endl;
        // adjuct first direction
        // direction[1] = direction[2];
        Eigen::MatrixXd final_path_up  = final_path.colwise().reverse();
        // std::cout<< "------------d6------------" << final_path << " " << final_path.cols() << std::endl;
        if(final_path_up.rows() > 2 && final_path_up.cols() > 3 ){
            final_path_up.row(0)[3] = final_path_up.row(1)[3];
            HybridPath path_final(final_path_up, finalcost);
            path = path_final;
        }else{
            std::cout << "Path can not be found..." << std::endl;
        }

        // std::cout<< " final_path length " << final_path_up.rows() << std::endl;
        // std::cout<< " final_path length " << final_path_up.col(2).transpose() << std::endl;
        
    }

    double TrailerHybridAStar::calc_cost(HybridNode& n, HybridNode& ngoal){
        // std::cout<< h_dp << std::endl;
        // int index_x = ((n.xind - config_.minx) < 0) ? 0 : n.xind - config_.minx;
        // int index_y = ((n.yind - config_.miny) < 0) ? 0 : n.yind - config_.miny;
        int index_x = (n.xind - config_.minx) - 1;
        int index_y = (n.yind - config_.miny) - 1;
        // std::cout<< h_dp.rows() << " " << h_dp.cols() << " " << " " << n.xind << " "<< n.yind <<"    index: " <<  index_x << "," << index_y  << std::endl;
        double total_cost = n.cost + PlannerParams::H_COST*h_dp(index_x, index_y);
        // std::cout<< " n.cost " << n.cost << " hp " << h_dp(index_x, index_y) << " total: "<< total_cost << std::endl;
        return total_cost;
    }


    bool TrailerHybridAStar::calc_hybrid_astar_path(Eigen::Vector4d s, Eigen::Vector4d g
                        , Eigen::MatrixXd& obses, HybridPath& path){
            double syaw = math_utility::pi_to_pi(s[2]);
            double gyaw = math_utility::pi_to_pi(g[2]);

            double xyreso = PlannerParams::XY_GRID_RESOLUTION;
            double yawreso = PlannerParams::YAW_GRID_RESOLUTION;
            grid_search::PointCloud cloud;
            cloud.points = obses;

            // std::cout<< " ============= "<< obses.rows()*obses.cols() << std::endl;
            // std::cout<< " ============= "<< obses << std::endl;

            grid_search::KDTree kdtree(2 , cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
            kdtree.buildIndex();
            calc_config(obses, xyreso, yawreso);
            
            // std::cout<< " ==========s=== "<< s.transpose() << " " << xyreso << std::endl;
            // std::cout<< " ==========g=== "<< g.transpose() << " " << xyreso  << " " << yawreso << std::endl;
            int xind_s = std::round(s[0]/xyreso);
            int yind_s = std::round(s[1]/xyreso);
            int yawind_s = std::round(s[2]/yawreso);

            int xind_g = std::round(g[0]/xyreso);
            int yind_g = std::round(g[1]/xyreso);
            int yawind_g = std::round(g[2]/yawreso);
            Eigen::MatrixXd poses_s(1,5);
            poses_s << s[0], s[1], s[2], s[3], 1.0;

            Eigen::MatrixXd poses_g(1,5);
            poses_g << g[0], g[1], g[2], g[3], 1.0;

            HybridNode nstart(xind_s, yind_s, yawind_s, true, poses_s, 0.0, 0.0, -1);
            HybridNode ngoal(xind_g, yind_g, yawind_g, true, poses_g, 0.0, 0.0, -1);
            calc_holonomic_with_obstacle_heuristic(ngoal, obses, xyreso);

            // std::cout<< "nstart: " << nstart << std::endl;
            // std::cout<< "ngoal: " << nstart << std::endl;
            
            std::unordered_map<int, HybridNode> open;
            std::unordered_map<int, HybridNode> closed;

            open[calc_index(nstart)] = nstart;

            std::priority_queue<CostNode, std::vector<CostNode>, CompareCostNode> pq;
            CostNode cost_start(calc_cost(nstart, ngoal), calc_index(nstart));
            pq.push(cost_start);

            std::vector<double> u;
            std::vector<double> d;
            calc_motion_inputs(u, d);

            int nmotion = u.size();
            HybridNode fnode;
            // std::cout<< "   u  d " << nmotion << " " << d.size() << std::endl;
            // print_vec(u);
            // print_vec(d);
            // std::cout<< "   ============================== " << std::endl;
            int counter = 0;

            // HybridNode n_tmp(xind_s, yind_s, yawind_s, true, poses_s, 0.0, 0.0, -1);
            // double inityaw1 = 0.0;
            // verify_index(n_tmp, obses, inityaw1, kdtree);
            // return false;
               



            while (true){
                if(open.empty()){
                    std::cout<< "Error: Cannot find path, No open set" << std::endl;
                    break;
                }
                if(pq.empty()){
                    std::cerr << " Can not find a path " << std::endl;
                    break;
                }
                CostNode next_node = pq.top();
                pq.pop();
                int c_id = next_node.id;
                HybridNode current = open[c_id];
                // move current node from open to closed
                // if( (current->pose.x() == ngoal.pose.x())&& (current->pose.y() == ngoal.pose.y())){
                //     std::cout << " Goal!! " << std::endl;
                //     closed[c_id] = current;
                //     break;
                // }

                open.erase(c_id);
                closed[c_id] = current;
                HybridNode fpath;
                // std::cout << " c_id " << c_id << std::endl; 
                // std::cout<< "------------------------node info------------------------" << std::endl;
                // std::cout << " current " << current << std::endl; 
                // std::cout << " ngoal " << ngoal << std::endl; 
                // std::cout << " gyaw1 " << g[3] << std::endl; 
                
                bool isupdated = update_node_with_analystic_expantion(current, ngoal, obses, kdtree, g[3], fpath);
                // std::cout << " isupdated " << isupdated << std::endl;
                // std::cout << " fpath " << fpath << std::endl;
                // break;
                // counter++;
                // if(counter == 4){
                //     break;
                // }
                if (isupdated){
                    // std::cout << " find the path " << fpath << std::endl;
                    fnode = fpath;
                    break;
                }
                // std::cout<< " -------current poses----------- " << current.poses << std::endl;
                double inityaw1 = current.poses.row(0)[3];
                // std::cout<< " --------------------------------- " << std::endl;
                for(int i=0; i<nmotion; i++){
                    HybridNode node = calc_next_node(current, c_id, u[i], d[i]);
                    if (!verify_index(node, obses, inityaw1, kdtree)){
                        // std::cout<< " verify node " << std::endl;
                        continue;
                    }
                    // if(node.pind == 203278){
                    //     std::cout<< "fiood"<< std::endl;
                    //     return false;
                    // }
                    int node_ind = calc_index(node);
                    // std::cout<< "------------------------node info------------------------" << std::endl;
                    // std::cout<< " node_ind " << node_ind << std::endl;
                    // std::cout<< " node " << node << std::endl;
                    // std::cout<< "------------------------end----------------------" << std::endl;
                    //  If it is already in the closed set, skip it
                    if (closed.find(node_ind) != closed.end()) {
                        continue;
                    }
                    if (open.find(node_ind) == open.end()) {
                        open[node_ind] = node;
                        double cost = calc_cost(node, ngoal);
                        // std::cout<< " node_ind " << node_ind << " cost " << cost << std::endl;
                        // break;
                        CostNode cost_node(cost, node_ind);
                        pq.push(cost_node);
                    }else{
                        if(open[node_ind].cost > node.cost){
                            // If so, update the node to have a new parent
                            // std::cout<< " open " << open[node_ind].cost << " cost " << node.cost << " id: " << node.pind << std::endl;
                            // std::cout<< " open " << open[node_ind].poses.col(4).transpose() << std::endl;
                            open[node_ind] = node;
                        }
                    }

                    // counter++;
                    // if(counter == 2){
                    //     return false;
                    // }
                }
                // std::cout<< " -------------------end-------------- " << std::endl;
                // break;
                
            }

            
            std::cout<< "final expand node:" << open.size() + closed.size() << std::endl;
            get_final_path(closed, fnode, nstart, path);
            return true;
    }


    bool TrailerHybridAStar::analystic_expantion(HybridNode& node, HybridNode& ngoal
                        ,  Eigen::MatrixXd& obses, grid_search::KDTree& kdtree
                        , rs_paths::Path& selected_path){

        Eigen::VectorXd current_pose = node.poses.row(node.poses.rows()-1);
        Eigen::VectorXd target_pose = ngoal.poses.row(ngoal.poses.rows()-1);
        double sx = current_pose[0];
        double sy = current_pose[1];
        double syaw = current_pose[2];
        double max_curvature = tan(PlannerParams::MAX_STEER)/PlannerParams::WB;

        std::vector<rs_paths::Path> paths;
        rs_path_.calc_paths(current_pose.head(3), target_pose.head(3)
                            , max_curvature, paths, PlannerParams::MOTION_RESOLUTION);

        // std::cout << " 1analystic_expantion paths.size() "<< paths.size() << std::endl;
        if(paths.size() == 0){
            return false;
        }
        // std::cout << " 1analystic_expantion paths "<< paths.size() << std::endl;

        std::priority_queue<rs_paths::Path, std::vector<rs_paths::Path>, rs_paths::CompareNode> pathqueue;
        

        for(auto path : paths){
            Eigen::VectorXd steps = path.poses.col(3).array()*PlannerParams::MOTION_RESOLUTION;
            Eigen::VectorXd yaw1;
            // std::cout <<  "======d 1" << std::endl;
            bool can_estimate = trailerlib_.calc_trailer_yaw_from_xyyaw(path.poses, current_pose[3], steps, yaw1);
            // if(!can_estimate){
            //     return false;
            // }
            // std::cout <<  "======d 2" << std::endl;
            path.cost = calc_rs_path_cost(path, yaw1);
            // std::cout << " 1analystic_expantion path.cost "<< path.cost << std::endl;
            pathqueue.push(path);
        }

        while (!pathqueue.empty()) {
            rs_paths::Path path = pathqueue.top();
            pathqueue.pop();

            Eigen::VectorXd steps = path.poses.col(3).array()*PlannerParams::MOTION_RESOLUTION;
            Eigen::VectorXd yaws1;
            bool can_estimate = trailerlib_.calc_trailer_yaw_from_xyyaw(node.poses, current_pose[3], steps, yaws1);
            // if(!can_estimate){
            //     return false;
            // }
            std::vector<int> indices;
            for(int i=0; i<path.poses.rows(); i+=PlannerParams::SKIP_COLLISION_CHECK){
                indices.push_back(i);
            }

            Eigen::MatrixXd selected_poses(indices.size(), path.poses.cols());
            for (size_t i = 0; i < indices.size(); ++i) {
                selected_poses.row(i) = path.poses.row(indices[i]);
            }

            if (trailerlib_.check_trailer_collision(obses, selected_poses, kdtree)){
                // std::cout << " 1analystic_expantion check_trailer_collision find path " << std::endl;
                selected_path = path;
                return true;
            }
        } 
        // std::cout << " 1analystic_expantion 2" << std::endl;
        return false;       
    }

    bool TrailerHybridAStar::is_same_grid(HybridNode& node1, HybridNode& node2){
        if(node1.xind != node2.xind){
            return false;
        }
        if(node1.yind != node2.yind){
            return false;
        }
        if(node1.yawind != node2.yawind){
            return false;
        }
        return true;
    }

    bool TrailerHybridAStar::verify_index(HybridNode& node, Eigen::MatrixXd& obses
                            , double inityaw1, grid_search::KDTree& kdtree){
        
        if ((node.xind - config_.minx) >= config_.xw){
            return false;
        }else if((node.xind - config_.minx) <= 0){
            return false;
        }
        if ((node.yind - config_.miny) >= config_.yw){
            return false;
        }else if((node.yind - config_.miny) <= 0){
            return false;
        }

         // check collisiton
        Eigen::VectorXd steps = node.poses.col(4).array()*PlannerParams::MOTION_RESOLUTION;
        Eigen::VectorXd yaws1;
        bool can_estimate = trailerlib_.calc_trailer_yaw_from_xyyaw(node.poses, inityaw1, steps, yaws1);
        // if(!can_estimate){
        //         return false;
        // }
        std::vector<int> indices;
        for(int i=0; i<node.poses.rows(); i+=PlannerParams::SKIP_COLLISION_CHECK){
            indices.push_back(i);
        }

        Eigen::MatrixXd selected_poses(indices.size(), node.poses.cols());
        for (size_t i = 0; i < indices.size(); ++i) {
            selected_poses.row(i) = node.poses.row(indices[i]);
        }

        if (!trailerlib_.check_trailer_collision(obses,  selected_poses, kdtree)){
            return false;
        }
        return true;
    }


    HybridNode TrailerHybridAStar::calc_next_node(HybridNode& current, int c_id, double u, double d){

        double arc_l = PlannerParams::XY_GRID_RESOLUTION*1.5;
        int nlist = std::floor(arc_l/PlannerParams::MOTION_RESOLUTION) + 1;
        Eigen::MatrixXd poses(nlist, 5);
        Eigen::VectorXd current_pose = current.poses.row(current.poses.rows()-1);
        
        // std::cout<< " current pose " << current.poses << " nlist " << nlist << std::endl;
        if( nlist >0 ){
            poses.row(0)[0] = current_pose.x() + d * PlannerParams::MOTION_RESOLUTION*cos(current_pose[2]);
            poses.row(0)[1] = current_pose.y() + d * PlannerParams::MOTION_RESOLUTION*sin(current_pose[2]);
            poses.row(0)[2] = math_utility::pi_to_pi(current_pose[2] 
                                + d*PlannerParams::MOTION_RESOLUTION/PlannerParams::WB * tan(u));
            poses.row(0)[3] = math_utility::pi_to_pi(current_pose[3] 
                        + d*PlannerParams::MOTION_RESOLUTION/PlannerParams::LT*sin(current_pose[2]-current_pose[3]));
        }

        for(int i=0; i<(nlist-1); i++){
            poses.row(i+1)[0] = poses.row(i)[0] 
                                + d * PlannerParams::MOTION_RESOLUTION*cos(poses.row(i)[2]);
            poses.row(i+1)[1] = poses.row(i)[1] 
                                + d * PlannerParams::MOTION_RESOLUTION*sin(poses.row(i)[2]);
            poses.row(i+1)[2] = math_utility::pi_to_pi(poses.row(i)[2] 
                                + d*PlannerParams::MOTION_RESOLUTION/PlannerParams::WB * tan(u));  
            poses.row(i+1)[3] = math_utility::pi_to_pi(poses.row(i)[3] 
                                + d*PlannerParams::MOTION_RESOLUTION/PlannerParams::LT*sin(poses.row(i)[2]-poses.row(i)[3]));
        }

        Eigen::VectorXd updated_pose = poses.row(poses.rows()-1);
        // std::cout<< " cpose " << poses << std::endl;
        int xind = std::round(updated_pose[0]/config_.xyreso);
        int yind = std::round(updated_pose[1]/config_.xyreso);
        int yawind = std::round(updated_pose[2]/config_.yawreso);
        
        // std::cout<< " xind " << xind << " yind " << yind << " yawind " << yawind << std::endl;

        double addedcost = 0.0;
        double direction = 1.0;
        // std::cout<< " =========================cost=========================== " << std::endl;
        if (d > 0){
            direction = 1.0;
            addedcost += std::abs(arc_l);
        }else{
            direction = 0.0;
            addedcost += std::abs(arc_l) * PlannerParams::BACK_COST;
        }
        // std::cout<< " =========== addedcost1 " << addedcost << std::endl;
        // swich back penalty
        if (direction != current.direction){ // switch back penalty
            addedcost += PlannerParams::SB_COST;
        }
        // std::cout<< " =========== addedcost2 " << addedcost << std::endl;
        // steer penalyty
        addedcost += PlannerParams::STEER_COST*std::abs(u);
        // std::cout<< " =========== addedcost3 " << addedcost << std::endl;
        // steer change penalty
        addedcost += PlannerParams::STEER_CHANGE_COST*std::abs(current.steer - u);
        // std::cout<< " =========== addedcost4 " << addedcost << std::endl;
        double total_angle_diff = (poses.col(2) - poses.col(3))
                    .unaryExpr([](double angle) { return std::abs(math_utility::pi_to_pi(angle)); })
                    .sum();
        // jacknif cost
        addedcost += PlannerParams::JACKKNIF_COST*total_angle_diff;
        // std::cout<< " =========== yawlist " <<  poses.col(3).transpose() << std::endl;
        // std::cout<< " =========== yawlist " << poses.col(2).transpose()  << " ------- " << poses.col(3).transpose() << std::endl;
        // std::cout<< " =========== addedcost5 " << addedcost  << " total_angle_diff " << total_angle_diff << std::endl;
        double cost = current.cost + addedcost;
        // std::cout<< " =========== addedcost6 " << addedcost << std::endl;
        poses.col(4).setConstant(direction);
        HybridNode next_node(xind, yind, yawind, direction, poses, u, cost, c_id);
        return next_node;
    }

    
    

    
}


