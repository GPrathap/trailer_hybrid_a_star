#include "lib_cpp_hybrid_a_star/grid_a_star.hpp"


namespace grid_search
{

    GridAStar::GridAStar(){

    }


    void GridAStar::calc_obstacle_map(Eigen::MatrixXd& obses, double& reso, double& vr){

        // Eigen::VectorXd col_min = obses.colwise().minCoeff(); 
        // Eigen::VectorXd col_max = obses.colwise().maxCoeff();

        // minx_ = std::round(col_min[0]);
        // miny_ = std::round(col_min[1]);
        // maxx_ = std::round(col_max[0]);
        // maxy_ = std::round(col_max[1]);

        minx_ = -15;
        miny_ = -10;
        maxx_ = 15;
        maxy_ = 10;
        
        xwidth_ = maxx_ - minx_;
        ywidth_ = maxy_ - miny_;
        // std::cout<< "map min: " << col_min.transpose() << std::endl;
        // std::cout<< "map max: " << col_max.transpose() << std::endl;
        std::cout<< "xwidth_: " << xwidth_ << " " << ywidth_ << std::endl;
        obs_map.resize(xwidth_, ywidth_);
        obs_map.setZero();

        PointCloud cloud;
        cloud.points = obses;

        // std::cout<< " ============= "<< obses.rows()*obses.cols() << std::endl;
        // std::cout<< " ============= "<< obses << std::endl;

        KDTree kd_tree(2 , cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
        kd_tree.buildIndex();
        const size_t num_results = 1;
        std::vector<size_t> ret_indexes(num_results); 
        std::vector<double> out_dist_sqr(num_results);
        nanoflann::KNNResultSet<double> resultSet(num_results);
        for(int ix=0; ix<xwidth_; ix++ ){
            int x = ix + minx_;
            for(int iy=0; iy<ywidth_; iy++ ){
                int y = iy + miny_;
                const double query_pt[2] = {x*1.0, y*1.0};
                resultSet.init(ret_indexes.data(), out_dist_sqr.data());
                kd_tree.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParameters(10));
                double dis = std::sqrt(out_dist_sqr[0]);
                if (dis <= vr/reso){
                    // obs_map(ix,iy) = 1;
                } 
            }
        }
        // std::cout << "Obs map: " << obs_map << std::endl;
    }

    int GridAStar::calc_index(const Node& node){
        return (int)((node.pose.y() - miny_)*xwidth_ + (node.pose.x() - minx_));
    }

    double GridAStar::calc_cost(const Node& n, const Node& ngoal) {
        return (n.cost + h(n.pose.x() - ngoal.pose.x(), n.pose.y() - ngoal.pose.y()));
    }

    double GridAStar::h(int x, int y){
        return sqrt( pow(x,2) + pow(y,2));
    }

    Eigen::MatrixXd GridAStar::get_motion_model(){
        Eigen::MatrixXd motion(8, 3);
        motion<< 1, 0, 1, 0, 1, 1, -1, 0, 1, 0, -1, 1, -1, -1, sqrt(2), -1, 1, sqrt(2), 1, -1, sqrt(2), 1, 1, sqrt(2);
        return motion;
    }

    void GridAStar::calc_dist_policy(Eigen::Vector2d g
                        , Eigen::MatrixXd obses, double reso, double vr
                        , Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& pmap){
        
        Eigen::Vector2i goal_pose (g.x()/reso, g.y()/reso);
        Node ngoal(goal_pose, 0.0, -1);  
        for(int i=0; i<obses.rows(); i++){
            obses.row(i) = obses.row(i)/reso;
        }   
        calc_obstacle_map(obses, reso, vr);

        std::unordered_map<int, Node> open;
        std::unordered_map<int, Node> closed;
        open[calc_index(ngoal)] = ngoal;

        Eigen::MatrixXd motion = get_motion_model();
        int nmotion = motion.rows();
        // std::priority_queue<int> pq;
        std::priority_queue<CostNode, std::vector<CostNode>, CompareNode> pq;
        CostNode cost_goal(ngoal.cost, calc_index(ngoal));
        pq.push(cost_goal);
        int counter = 0;
        while( true ){
            if(open.empty()){
                std::cerr << " Finish Search " << std::endl;
                break;
            }
            CostNode next_node = pq.top();
            pq.pop();
            int c_id = next_node.id;
            Node current = open[c_id];
            // if( (current->pose.x() == ngoal.pose.x())&& (current->pose.y() == ngoal.pose.y())){
            //     std::cout << " Goal!! " << std::endl;
            //     closed[c_id] = current;
            //     break;
            // }
            // std::cout<< "----------------------------------------------" <<  std::endl;
            // std::cout<< "current->pind " <<  current.pose.transpose() <<  std::endl;
            open.erase(c_id);
            closed[c_id] = current;
            for(int i=0; i<nmotion; i++){
                // std::cout<< "current pose " <<  current->pose.transpose() <<  std::endl;
                Node node(Eigen::Vector2i(current.pose.x() + motion.row(i)[0], current.pose.y() + motion.row(i)[1])
                                , current.cost + motion.row(i)[2], c_id);
                
                // std::cout<< "node.pind " <<  node.pind <<  std::endl;
                if(!verify_node(node)){
                    // std::cout<< " not verieid " << node.pose.transpose() << std::endl;
                    continue;
                }
                int node_ind = calc_index(node);
                // std::cout<< "node pose " <<  node.pose.transpose() << " " << c_id << " " << node_ind <<  std::endl;
                if (closed.find(node_ind) != closed.end()) {
                    continue;
                }
                if (open.find(node_ind) != open.end()) {
                    if(open[node_ind].cost > node.cost){
                        open[node_ind].cost = node.cost;
                        open[node_ind].pind = c_id;
                        // std::cout<< "open node  " <<  node.cost << " " << c_id << " " << node.pose.transpose()  <<  std::endl;
                        // return;
                    }
                }else{
                    open[node_ind] = node;
                    CostNode cost_node(node.cost, calc_index(node));
                    // std::cout<< " close node " << node.pind << " cost " << node.cost << " " << node.pose.transpose() << std::endl; 
                    pq.push(cost_node);
                }
            }
            // counter++;
            // if(counter == 5){
            //     break;
            // }
        }

        std::cout<< "  grid search length(closed) " << closed.size() << std::endl;
        calc_policy_map(closed, pmap);
    }


    Node GridAStar::search_min_cost_node(const std::unordered_map<int, Node*> open, const Node ngoal){
        double mcost = std::numeric_limits<double>::infinity();
        Node mnode = ngoal;
        for (auto it = open.begin(); it != open.end(); ++it) {
            double cost = it->second->cost + h(it->second->pose.x() - ngoal.pose.x(), it->second->pose.y() - ngoal.pose.y());
            if(mcost > cost){
                mnode = *(it->second);
                mcost = cost;
            }
        }
        return mnode;
    }

    bool GridAStar::verify_node(Node node){
        // std::cout<< minx_ << " " << xwidth_ << miny_ << " " << ywidth_ << std::endl;

        if ((node.pose.x() - minx_) >= xwidth_){
            return false;
        } else if (( node.pose.x() - minx_) < 0 ){
            return false;
        } 
        if ((node.pose.y() - miny_) >= ywidth_){
            return false;
        } else if (( node.pose.y() - miny_) < 0 ){
            return false;
        }       
        // std::cout<< obs_map << std::endl;
        if (obs_map(node.pose.x()-minx_, node.pose.y()-miny_) > 0){
            return false;
        }
        return true;                    
    }

    void GridAStar::calc_policy_map(std::unordered_map<int, Node>& closed
                                                ,  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& pmap){
        pmap.resize(xwidth_, ywidth_);
        // pmap << 27.24264,26.24264,25.24264,24.24264,23.24264,22.24264,21.24264,20.24264,19.24264,18.24264,17.24264,16.24264,15.82843,16.24264,16.65685,17.07107,17.48528,17.89949,18.89949,std::numeric_limits<double>::infinity(),
        //         27.65685,26.65685,25.65685,24.65685,23.65685,22.65685,21.65685,20.65685,19.65685,18.65685,17.65685,std::numeric_limits<double>::infinity(),14.82843,15.24264,15.65685,16.07107,16.48528,17.48528,18.48528,std::numeric_limits<double>::infinity(),
        //         28.07107,27.07107,26.07107,25.07107,24.07107,23.07107,22.07107,21.07107,20.07107,19.07107,18.65685,std::numeric_limits<double>::infinity(),13.82843,14.24264,14.65685,15.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),18.89949,std::numeric_limits<double>::infinity(),
        //         28.48528,27.48528,26.48528,25.48528,24.48528,23.48528,22.48528,21.48528,20.48528,20.07107,19.65685,std::numeric_limits<double>::infinity(),12.82843,13.24264,13.65685,14.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),19.89949,std::numeric_limits<double>::infinity(),
        //         28.89949,27.89949,26.89949,25.89949,24.89949,23.89949,22.89949,21.89949,21.48528,21.07107,20.65685,std::numeric_limits<double>::infinity(),11.82843,12.24264,12.65685,13.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),20.89949,std::numeric_limits<double>::infinity(),
        //         29.31371,28.31371,27.31371,26.31371,25.31371,24.31371,23.31371,22.89949,22.48528,22.07107,21.65685,std::numeric_limits<double>::infinity(),10.82843,11.24264,11.65685,12.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),21.89949,std::numeric_limits<double>::infinity(),
        //         29.72792,28.72792,27.72792,26.72792,25.72792,24.72792,24.31371,23.89949,23.48528,23.07107,22.65685,std::numeric_limits<double>::infinity(),9.82843,10.24264,10.65685,11.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),22.89949,std::numeric_limits<double>::infinity(),
        //         30.14214,29.14214,28.14214,27.14214,26.14214,25.72792,25.31371,24.89949,24.48528,24.07107,23.65685,std::numeric_limits<double>::infinity(),8.82843,9.24264,9.65685,10.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),23.89949,std::numeric_limits<double>::infinity(),
        //         30.55635,29.55635,28.55635,27.55635,27.14214,26.72792,26.31371,25.89949,25.48528,25.07107,24.65685,std::numeric_limits<double>::infinity(),7.82843,8.24264,8.65685,9.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),24.89949,std::numeric_limits<double>::infinity(),
        //         30.97056,29.97056,28.97056,28.55635,28.14214,27.72792,27.31371,26.89949,26.48528,26.07107,25.65685,std::numeric_limits<double>::infinity(),6.82843,7.24264,7.65685,8.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),25.89949,std::numeric_limits<double>::infinity(),
        //         31.38478,30.38478,29.97056,29.55635,29.14214,28.72792,28.31371,27.89949,27.48528,27.07107,26.65685,std::numeric_limits<double>::infinity(),5.82843,6.24264,6.65685,7.65685,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),26.89949,std::numeric_limits<double>::infinity(),
        //         31.79899,31.38478,30.97056,30.55635,30.14214,29.72792,29.31371,28.89949,28.48528,28.07107,27.65685,std::numeric_limits<double>::infinity(),4.82843,5.24264,6.24264,7.24264,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),27.89949,std::numeric_limits<double>::infinity(),
        //         32.79899,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),3.82843,4.82843,5.82843,6.82843,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),28.89949,std::numeric_limits<double>::infinity(),
        //         33.79899,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),6.41421,5.41421,4.41421,3.41421,2.41421,1.41421,1.00000,1.41421,2.41421,3.41421,4.41421,5.41421,6.41421,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),29.89949,std::numeric_limits<double>::infinity(),
        //         34.79899,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),6.00000,5.00000,4.00000,3.00000,2.00000,1.00000,0.00000,1.00000,2.00000,3.00000,4.00000,5.00000,6.00000,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),30.89949,std::numeric_limits<double>::infinity(),
        //         33.79899,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),6.41421,5.41421,4.41421,3.41421,2.41421,1.41421,1.00000,1.41421,2.41421,3.41421,4.41421,5.41421,6.41421,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),29.89949,std::numeric_limits<double>::infinity(),
        //         32.79899,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),3.82843,4.82843,5.82843,6.82843,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),28.89949,std::numeric_limits<double>::infinity(),
        //         31.79899,31.38478,30.97056,30.55635,30.14214,29.72792,29.31371,28.89949,28.48528,28.07107,27.65685,std::numeric_limits<double>::infinity(),4.82843,5.24264,6.24264,7.24264,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),27.89949,std::numeric_limits<double>::infinity(),
        //         31.38478,30.38478,29.97056,29.55635,29.14214,28.72792,28.31371,27.89949,27.48528,27.07107,26.65685,std::numeric_limits<double>::infinity(),5.82843,6.24264,6.65685,7.65685,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),26.89949,std::numeric_limits<double>::infinity(),
        //         30.97056,29.97056,28.97056,28.55635,28.14214,27.72792,27.31371,26.89949,26.48528,26.07107,25.65685,std::numeric_limits<double>::infinity(),6.82843,7.24264,7.65685,8.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),25.89949,std::numeric_limits<double>::infinity(),
        //         30.55635,29.55635,28.55635,27.55635,27.14214,26.72792,26.31371,25.89949,25.48528,25.07107,24.65685,std::numeric_limits<double>::infinity(),7.82843,8.24264,8.65685,9.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),24.89949,std::numeric_limits<double>::infinity(),
        //         30.14214,29.14214,28.14214,27.14214,26.14214,25.72792,25.31371,24.89949,24.48528,24.07107,23.65685,std::numeric_limits<double>::infinity(),8.82843,9.24264,9.65685,10.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),23.89949,std::numeric_limits<double>::infinity(),
        //         29.72792,28.72792,27.72792,26.72792,25.72792,24.72792,24.31371,23.89949,23.48528,23.07107,22.65685,std::numeric_limits<double>::infinity(),9.82843,10.24264,10.65685,11.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),22.89949,std::numeric_limits<double>::infinity(),
        //         29.31371,28.31371,27.31371,26.31371,25.31371,24.31371,23.31371,22.89949,22.48528,22.07107,21.65685,std::numeric_limits<double>::infinity(),10.82843,11.24264,11.65685,12.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),21.89949,std::numeric_limits<double>::infinity(),
        //         28.89949,27.89949,26.89949,25.89949,24.89949,23.89949,22.89949,21.89949,21.48528,21.07107,20.65685,std::numeric_limits<double>::infinity(),11.82843,12.24264,12.65685,13.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),20.89949,std::numeric_limits<double>::infinity(),
        //         28.48528,27.48528,26.48528,25.48528,24.48528,23.48528,22.48528,21.48528,20.48528,20.07107,19.65685,std::numeric_limits<double>::infinity(),12.82843,13.24264,13.65685,14.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),19.89949,std::numeric_limits<double>::infinity(),
        //         28.07107,27.07107,26.07107,25.07107,24.07107,23.07107,22.07107,21.07107,20.07107,19.07107,18.65685,std::numeric_limits<double>::infinity(),13.82843,14.24264,14.65685,15.07107,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),18.89949,std::numeric_limits<double>::infinity(),
        //         27.65685,26.65685,25.65685,24.65685,23.65685,22.65685,21.65685,20.65685,19.65685,18.65685,17.65685,std::numeric_limits<double>::infinity(),14.82843,15.24264,15.65685,16.07107,16.48528,17.48528,18.48528,std::numeric_limits<double>::infinity(),
        //         27.24264,26.24264,25.24264,24.24264,23.24264,22.24264,21.24264,20.24264,19.24264,18.24264,17.24264,16.24264,15.82843,16.24264,16.65685,17.07107,17.48528,17.89949,18.89949,
        //         std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity();
                        
        // pmap << 17.72792206135786,17.313708498984763,16.899494936611667,16.48528137423857,16.071067811865476,15.65685424949238,15.242640687119284,14.82842712474619,14.414213562373096,14.0,14.414213562373096,14.82842712474619,15.242640687119284,15.65685424949238,16.071067811865476
        //         ,16.48528137423857,16.899494936611667,17.313708498984763,17.72792206135786,std::numeric_limits<double>::infinity(), 16.72792206135786,16.313708498984763,15.899494936611667,15.485281374238571,15.071067811865476,14.65685424949238,14.242640687119284,13.82842712474619
        //         ,13.414213562373096,13.0,13.414213562373096,13.82842712474619,14.242640687119284,14.65685424949238,15.071067811865476,15.485281374238571,15.899494936611667,16.313708498984763,16.72792206135786,std::numeric_limits<double>::infinity(), 15.727922061357859,15.313708498984763
        //         ,14.899494936611667,14.485281374238571,14.071067811865476,13.65685424949238,13.242640687119284,12.82842712474619,12.414213562373096,12.0,12.414213562373096,12.82842712474619,13.242640687119284,13.65685424949238,14.071067811865476,14.485281374238571,14.899494936611667
        //         ,15.313708498984763,15.727922061357859,std::numeric_limits<double>::infinity(), 14.727922061357859,14.313708498984763,13.899494936611667,13.485281374238571,13.071067811865476,12.65685424949238,12.242640687119284,11.82842712474619,11.414213562373096,11.0,11.414213562373096
        //         ,11.82842712474619,12.242640687119284,12.65685424949238,13.071067811865476,13.485281374238571,13.899494936611667,14.313708498984763,14.727922061357859,std::numeric_limits<double>::infinity(), 13.727922061357859,13.313708498984763,12.899494936611667,12.485281374238571
        //         ,12.071067811865476,11.65685424949238,11.242640687119284,10.82842712474619,10.414213562373096,10.0,10.414213562373096,10.82842712474619,11.242640687119284,11.65685424949238,12.071067811865476,12.485281374238571,12.899494936611667,13.313708498984763,13.727922061357859
        //         ,std::numeric_limits<double>::infinity(), 12.727922061357859,12.313708498984763,11.899494936611667,11.485281374238571,11.071067811865476,10.65685424949238,10.242640687119284,9.82842712474619,9.414213562373096,9.0,9.414213562373096,9.82842712474619,10.242640687119284
        //         ,10.65685424949238,11.071067811865476,11.485281374238571,11.899494936611667,12.313708498984763,12.727922061357859,std::numeric_limits<double>::infinity(), 12.313708498984763,11.313708498984763,10.899494936611667,10.485281374238571,10.071067811865476,9.65685424949238
        //         ,9.242640687119284,8.82842712474619,8.414213562373096,8.0,8.414213562373096,8.82842712474619,9.242640687119284,9.65685424949238,10.071067811865476,10.485281374238571,10.899494936611667,11.313708498984763,12.313708498984763,std::numeric_limits<double>::infinity()
        //         , 11.899494936611667,10.899494936611667,9.899494936611667,9.485281374238571,9.071067811865476,8.65685424949238,8.242640687119284,7.82842712474619,7.414213562373095,7.0,7.414213562373095,7.82842712474619,8.242640687119284,8.65685424949238,9.071067811865476,9.485281374238571
        //         ,9.899494936611667,10.899494936611667,11.899494936611667,std::numeric_limits<double>::infinity(), 11.485281374238571,10.485281374238571,9.485281374238571,8.485281374238571,8.071067811865476,7.65685424949238,7.242640687119285,6.82842712474619,6.414213562373095,6.0
        //         ,6.414213562373095,6.82842712474619,7.242640687119285,7.65685424949238,8.071067811865476,8.485281374238571,9.485281374238571,10.485281374238571,11.485281374238571,std::numeric_limits<double>::infinity(), 11.071067811865476,10.071067811865476,9.071067811865476
        //         ,8.071067811865476,7.0710678118654755,6.65685424949238,6.242640687119285,5.82842712474619,5.414213562373095,5.0,5.414213562373095,5.82842712474619,6.242640687119285,6.65685424949238,7.0710678118654755,8.071067811865476,9.071067811865476,10.071067811865476,11.071067811865476
        //         ,std::numeric_limits<double>::infinity(), 10.65685424949238,9.65685424949238,8.65685424949238,7.65685424949238,6.65685424949238,5.656854249492381,5.242640687119285,4.82842712474619,4.414213562373095,4.0,4.414213562373095,4.82842712474619,5.242640687119285,5.656854249492381
        //         ,6.65685424949238,7.65685424949238,8.65685424949238,9.65685424949238,10.65685424949238,std::numeric_limits<double>::infinity(), 10.242640687119284,9.242640687119284,8.242640687119284,7.242640687119285,6.242640687119285,5.242640687119285,4.242640687119286,3.82842712474619
        //         ,3.414213562373095,3.0,3.414213562373095,3.82842712474619,4.242640687119286,5.242640687119285,6.242640687119285,7.242640687119285,8.242640687119284,9.242640687119284,10.242640687119284,std::numeric_limits<double>::infinity(), 9.82842712474619,8.82842712474619,7.82842712474619
        //         ,6.82842712474619,5.82842712474619,4.82842712474619,3.82842712474619,2.8284271247461903,2.414213562373095,2.0,2.414213562373095,2.8284271247461903,3.82842712474619,4.82842712474619,5.82842712474619,6.82842712474619,7.82842712474619,8.82842712474619,9.82842712474619
        //         ,std::numeric_limits<double>::infinity(), 9.414213562373096,8.414213562373096,7.414213562373095,6.414213562373095,5.414213562373095,4.414213562373095,3.414213562373095,2.414213562373095,1.4142135623730951,1.0,1.4142135623730951,2.414213562373095,3.414213562373095
        //         ,4.414213562373095,5.414213562373095,6.414213562373095,7.414213562373095,8.414213562373096,9.414213562373096,std::numeric_limits<double>::infinity(), 9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.0,1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,std::numeric_limits<double>::infinity()
        //         , 9.414213562373096,8.414213562373096,7.414213562373095,6.414213562373095,5.414213562373095,4.414213562373095,3.414213562373095,2.414213562373095,1.4142135623730951,1.0,1.4142135623730951,2.414213562373095,3.414213562373095,4.414213562373095,5.414213562373095
        //         ,6.414213562373095,7.414213562373095,8.414213562373096,9.414213562373096,std::numeric_limits<double>::infinity(), 9.82842712474619,8.82842712474619,7.82842712474619,6.82842712474619,5.82842712474619,4.82842712474619,3.82842712474619,2.8284271247461903,2.414213562373095
        //         ,2.0,2.414213562373095,2.8284271247461903,3.82842712474619,4.82842712474619,5.82842712474619,6.82842712474619,7.82842712474619,8.82842712474619,9.82842712474619,std::numeric_limits<double>::infinity(), 10.242640687119284,9.242640687119284,8.242640687119284
        //         ,7.242640687119285,6.242640687119285,5.242640687119285,4.242640687119286,3.82842712474619,3.414213562373095,3.0,3.414213562373095,3.82842712474619,4.242640687119286,5.242640687119285,6.242640687119285,7.242640687119285,8.242640687119284,9.242640687119284
        //         ,10.242640687119284,std::numeric_limits<double>::infinity(), 10.65685424949238,9.65685424949238,8.65685424949238,7.65685424949238,6.65685424949238,5.656854249492381,5.242640687119285,4.82842712474619,4.414213562373095,4.0,4.414213562373095,4.82842712474619
        //         ,5.242640687119285,5.656854249492381,6.65685424949238,7.65685424949238,8.65685424949238,9.65685424949238,10.65685424949238,std::numeric_limits<double>::infinity(), 11.071067811865476,10.071067811865476,9.071067811865476,8.071067811865476,7.0710678118654755
        //         ,6.65685424949238,6.242640687119285,5.82842712474619,5.414213562373095,5.0,5.414213562373095,5.82842712474619,6.242640687119285,6.65685424949238,7.0710678118654755,8.071067811865476,9.071067811865476,10.071067811865476,11.071067811865476
        //         ,std::numeric_limits<double>::infinity(), 11.485281374238571,10.485281374238571,9.485281374238571,8.485281374238571,8.071067811865476,7.65685424949238,7.242640687119285,6.82842712474619,6.414213562373095,6.0,6.414213562373095,6.82842712474619,7.242640687119285
        //         ,7.65685424949238,8.071067811865476,8.485281374238571,9.485281374238571,10.485281374238571,11.485281374238571,std::numeric_limits<double>::infinity(), 11.899494936611667,10.899494936611667,9.899494936611667,9.485281374238571,9.071067811865476
        //         ,8.65685424949238,8.242640687119284,7.82842712474619,7.414213562373095,7.0,7.414213562373095,7.82842712474619,8.242640687119284,8.65685424949238,9.071067811865476,9.485281374238571,9.899494936611667,10.899494936611667,11.899494936611667
        //         ,std::numeric_limits<double>::infinity(), 12.313708498984763,11.313708498984763,10.899494936611667,10.485281374238571,10.071067811865476,9.65685424949238,9.242640687119284,8.82842712474619,8.414213562373096,8.0,8.414213562373096,8.82842712474619,9.242640687119284
        //         ,9.65685424949238,10.071067811865476,10.485281374238571,10.899494936611667,11.313708498984763,12.313708498984763,std::numeric_limits<double>::infinity(), 12.727922061357859,12.313708498984763,11.899494936611667,11.485281374238571,11.071067811865476
        //         ,10.65685424949238,10.242640687119284,9.82842712474619,9.414213562373096,9.0,9.414213562373096,9.82842712474619,10.242640687119284,10.65685424949238,11.071067811865476,11.485281374238571,11.899494936611667,12.313708498984763,12.727922061357859
        //         ,std::numeric_limits<double>::infinity(), 13.727922061357859,13.313708498984763,12.899494936611667,12.485281374238571,12.071067811865476,11.65685424949238,11.242640687119284,10.82842712474619,10.414213562373096,10.0,10.414213562373096,10.82842712474619
        //         ,11.242640687119284,11.65685424949238,12.071067811865476,12.485281374238571,12.899494936611667,13.313708498984763,13.727922061357859,std::numeric_limits<double>::infinity(), 14.727922061357859,14.313708498984763,13.899494936611667,13.485281374238571
        //         ,13.071067811865476,12.65685424949238,12.242640687119284,11.82842712474619,11.414213562373096,11.0,11.414213562373096,11.82842712474619,12.242640687119284,12.65685424949238,13.071067811865476,13.485281374238571,13.899494936611667,14.313708498984763
        //         ,14.727922061357859,std::numeric_limits<double>::infinity(), 15.727922061357859,15.313708498984763,14.899494936611667,14.485281374238571,14.071067811865476,13.65685424949238,13.242640687119284,12.82842712474619,12.414213562373096,12.0,12.414213562373096
        //         ,12.82842712474619,13.242640687119284,13.65685424949238,14.071067811865476,14.485281374238571,14.899494936611667,15.313708498984763,15.727922061357859,std::numeric_limits<double>::infinity(), 16.72792206135786,16.313708498984763,15.899494936611667
        //         ,15.485281374238571,15.071067811865476,14.65685424949238,14.242640687119284,13.82842712474619,13.414213562373096,13.0,13.414213562373096,13.82842712474619,14.242640687119284,14.65685424949238,15.071067811865476,15.485281374238571,15.899494936611667
        //         ,16.313708498984763,16.72792206135786,std::numeric_limits<double>::infinity(), 17.72792206135786,17.313708498984763,16.899494936611667,16.48528137423857,16.071067811865476,15.65685424949238,15.242640687119284,14.82842712474619,14.414213562373096,14.0
        //         ,14.414213562373096,14.82842712474619,15.242640687119284,15.65685424949238,16.071067811865476,16.48528137423857,16.899494936611667,17.313708498984763,17.72792206135786,std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()
        //         ,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()
        //         ,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()
        //         ,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()
        //         ,std::numeric_limits<double>::infinity();


        pmap.setConstant(std::numeric_limits<double>::infinity());
        for (auto it = closed.begin(); it != closed.end(); ++it) {
            pmap(it->second.pose.x()-minx_, it->second.pose.y()-miny_) = it->second.cost;
        }
    }

    void GridAStar::calc_astar_path(Eigen::Vector2d s, Eigen::Vector2d g
                        , Eigen::MatrixXd obses, double reso, double vr){
        
        Eigen::Vector2i goal_pose (g.x()/reso, g.y()/reso);
        Eigen::Vector2i start_pose (s.x()/reso, s.y()/reso);
        Node ngoal(goal_pose, 0.0, -1);  
        Node nstart(start_pose, 0.0, -1);  

        for(int i=0; i<obses.rows(); i++){
            obses.row(i) = obses.row(i)/reso;
        }   
        calc_obstacle_map(obses, reso, vr);

        std::unordered_map<int, Node> open;
        std::unordered_map<int, Node> closed;

        open[calc_index(nstart)] = nstart;

        Eigen::MatrixXd motion = get_motion_model();
        int nmotion = motion.rows();
        std::priority_queue<CostNode, std::vector<CostNode>, CompareNode> pq;
        CostNode cost_start(calc_cost(nstart, ngoal), calc_index(nstart));
        pq.push(cost_start);

        // std::cout<< "start " << start_pose.transpose() << " goal " << goal_pose.transpose() << std::endl; 
        while( true ){
            if(open.empty()){
                std::cerr << " Error: No open set " << std::endl;
                break;
            }

            if(pq.empty()){
                std::cerr << " Can not find a path " << std::endl;
                break;
            }
            CostNode next_node = pq.top();
            pq.pop();
            int c_id = next_node.id;
            Node current = open[c_id];
            // std::cout<< "c_id " << c_id << " current " << current.pose.transpose() << std::endl; 

            if( (current.pose.x() == ngoal.pose.x())&& (current.pose.y() == ngoal.pose.y())){
                std::cout << " Goal!! " << std::endl;
                closed[c_id] = current;
                break;
            }

            open.erase(c_id);
            closed[c_id] = current;

            for(int i=0; i<nmotion; i++){
                Node node(Eigen::Vector2i(current.pose.x() + motion.row(i)[0], current.pose.y() + motion.row(i)[1])
                                , current.cost + motion.row(i)[2], c_id);
                if(!verify_node(node)){
                    continue;
                }
                int node_ind = calc_index(node);
                // std::cout<< "node_ind " << node_ind << " , " << node.pose.transpose() << std::endl; 
                if (closed.find(node_ind) != closed.end()) {
                    continue;
                }
                if (open.find(node_ind) != open.end()) {
                    if(open[node_ind].cost > node.cost){
                        open[node_ind].cost = node.cost;
                        open[node_ind].pind = c_id;
                        // std::cout<< "open cost " << node.cost << " , " << c_id << std::endl; 
                    }
                }else{
                    open[node_ind] = node;
                    CostNode cost_node(calc_cost(node, ngoal), calc_index(node));
                    pq.push(cost_node);
                    // std::cout<< "enqueue " <<cost_node.cost << " , " << cost_node.id << std::endl; 
                }
            }
        }
        get_final_path(closed, ngoal, nstart, reso);
    }

    void GridAStar::get_final_path(std::unordered_map<int, Node>& closed, const Node ngoal
                            , const Node nstart, const double reso){

        std::vector<Eigen::Vector2i> path;
        Eigen::Vector2i r = ngoal.pose;
        path.push_back(r);
        int nid = calc_index(ngoal);
        if(closed.empty()){
            std::cerr << " Can not find a path " << std::endl;
            return;
        }
        while(true){
            if (closed.find(nid) == closed.end()) {
                std::cerr << " Can not find the full path " << std::endl;
                break;
            }
            Node next = closed[nid];
            path.push_back(next.pose);
            nid = next.pind;
            if((path.back().x() == nstart.pose.x()) & (path.back().y() == nstart.pose.y())){
                std::cout<< "find grid astar path " << std::endl;
                break;
            }
        }
        
        final_path_.clear();
        for (auto it = path.rbegin(); it != path.rend(); ++it) {
            Eigen::Vector2d next_pose(it->x()*reso, it->y()*reso);
            final_path_.push_back(next_pose);
        }
    }
}

// int main(int argc, char *argv[]){
    
//     Eigen::Vector2d s(10.0, 10.0);
//     Eigen::Vector2d g(50.0, 50.0);

//     Eigen::MatrixXd obss(320, 2);
//     int obs_index = 0;
//     for(int i=0; i< 60; i++){
//         obss.row(obs_index) << i*1.0, 0.0;
//         obs_index++;
//     }
//     for(int i=0; i< 60; i++){
//         obss.row(obs_index) << 60.0, i*1.0;
//         obs_index++;
//     }
//     for(int i=0; i< 60; i++){
//         obss.row(obs_index) << i*1.0, 60.0;
//         obs_index++;
//     }
//     for(int i=0; i< 60; i++){
//         obss.row(obs_index) << 0.0, i*1.0;
//         obs_index++;
//     }
//     for(int i=0; i< 40; i++){
//         obss.row(obs_index) << 20.0, i*1.0;
//         obs_index++;
//     }
//     for(int i=0; i< 40; i++){
//         obss.row(obs_index) << 40.0, 60.0-i*1.0;
//         obs_index++;
//     }

//     double VEHICLE_RADIUS = 5.0;
//     double GRID_RESOLUTION = 1.0;

//     grid_search::GridAStar grid_a_star;
//     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> pmap;
//     // grid_a_star.calc_dist_policy(s, g, obss, GRID_RESOLUTION, VEHICLE_RADIUS, pmap);
//     grid_a_star.calc_astar_path(s, g, obss, GRID_RESOLUTION, VEHICLE_RADIUS);
//     std::vector<Eigen::Vector2d> path_est = grid_a_star.getPath();
//     std::vector<double> path_x, path_y;
//     for(int i=0; i< path_est.size(); i++){
//         path_x.push_back(path_est[i][0]);
//         path_y.push_back(path_est[i][1]);
//     }

//     plt::figure();
//     std::vector<double> ox, oy;
//     for(int i=0; i< obss.rows(); i++){
//         ox.push_back(obss.row(i)[0]);
//         oy.push_back(obss.row(i)[1]);
//     }
//     plt::plot(ox, oy, ".r");
//     plt::plot(path_x, path_y, ".y");
//     plt::plot({s.x()}, {s.y()}, "bo"); // Start point
//     plt::plot({g.x()}, {g.y()}, "go");    // End point
//     plt::show();

// }
