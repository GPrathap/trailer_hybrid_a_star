#include <casadi/casadi.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include "matplotlibcpp.h"
#include <nlohmann/json.hpp>
#include "Eigen/Dense"

using namespace casadi;
using namespace std;
using json = nlohmann::json;

namespace plt = matplotlibcpp;


class TrajectoryOptimizer {
private:
    int n; // Number of points
    DM data;
    // DM A_p;
    // DM b_p;
    double curvature_weight;
    double curvature_rate_weight;
    double distance_weight;
    double max_steering_angle;
    double v_max;
    double delta_t;
    double wheelbase;
    vector<double> x_list, y_list, theta_list, d_list, v_list;
    vector<Eigen::Vector2d> selected_points;
    Eigen::Vector2d m_lin;
    double b_middle;
    Eigen::MatrixXd Ab_;

    vector<double> x0_guess;

    vector<MX> g;           // Constraints
    vector<double> lb_g;    // Lower bounds
    vector<double> ub_g;    // Upper bounds

public:
    TrajectoryOptimizer(): curvature_weight(1.0), curvature_rate_weight(1.0), distance_weight(1.0), 
          max_steering_angle(M_PI / 6.0), v_max(1.5), delta_t(0.1), wheelbase(1.0) {
    }

    void set_ref_trajectory(const DM& data){
        n = data.size1();
        x_list.clear();
        y_list.clear();
        theta_list.clear();
        for (int i = 0; i < n; ++i) {
            x_list.push_back(data(i, 0).scalar());
            y_list.push_back(data(i, 1).scalar());
            theta_list.push_back(data(i, 2).scalar());
            d_list.push_back(0.1);
        }
    }

    void set_ref_trajectory(const Eigen::MatrixXd& data){
        n = data.rows();
        x_list.clear();
        y_list.clear();
        theta_list.clear();
        d_list.clear();
        v_list.clear();
        
        for (int i = 0; i < n; i++) {
            x_list.push_back(data.row(i)[0]);
            y_list.push_back(data.row(i)[1]);
            theta_list.push_back(data.row(i)[2]);
            d_list.push_back(data.row(i)[3]);
            v_list.push_back(v_max);
        }
        n = d_list.size();
        std::vector<double> delta(n, 0.0);
        x0_guess.clear();
        x0_guess.insert(x0_guess.end(), x_list.begin(), x_list.end()-1);
        x0_guess.insert(x0_guess.end(), y_list.begin(), y_list.end()-1);
        x0_guess.insert(x0_guess.end(), theta_list.begin(), theta_list.end()-1);
        x0_guess.insert(x0_guess.end(), v_list.begin(), v_list.end()-1);
        x0_guess.insert(x0_guess.end(), delta.begin(), delta.end()-1);
        x0_guess.insert(x0_guess.end(), d_list.begin(), d_list.end()-2);
    }

    void set_boundary_constraints(const Eigen::MatrixXd& Ab){

        Ab_ = Ab;
        // A_p = DM::horzcat({{-0.0645709, 0.0645709}, {-0.997913, 0.997913}});
        // b_p = DM::vertcat({0.737278, 0.632449});

        // DM lin1 = A_p(casadi::Slice(0, 1), casadi::Slice());
        // DM lin2 = A_p(casadi::Slice(1, 2), casadi::Slice());


        // std::cout<< "=======================================ab================================" << std::endl;
        // std::cout<< Ab_ << std::endl;
        Eigen::Vector2d lin1 = Ab_.row(0).head(2);
        double b1 = Ab_.row(0)[2];
        Eigen::Vector2d lin2 = Ab_.row(1).head(2);
        double b2 = Ab_.row(1)[2];
        // std::cout<<  " lin1 " << lin1.transpose() << std::endl;
        // std::cout<<  " lin2 " << lin2.transpose() << std::endl;

        // Eigen::Vector2d m_lin;
        // double b_middle;
        if (static_cast<double>((lin1 + lin2).norm()) < 1e-1) {  // Tolerance to check if A1 == -A2
            // Compute the middle line
            m_lin = lin1;  
            b_middle = static_cast<double>((b1 - b2) / 2.0);
            std::cout << "The lines are  parallel" << std::endl;
        } else {
            
            m_lin = (lin1 + lin2).array()/2.0;
            b_middle = static_cast<double>((b1 + b2) / 2.0);
        }
        std::cout<< "=======================================ab================================" << std::endl;
        std::cout<< m_lin.transpose() << " b " << b_middle << std::endl;


        // std::cout<< "=======================================adb================================" << std::endl;
        // std::cout<< A_p << std::endl;
        // std::cout<<  " lin1g " << lin1 << std::endl;
        // std::cout<<  " lin2g " << lin2 << std::endl;
        

        // if (static_cast<double>(norm_2(lin1 + lin2)) < 1e-3) {  // Tolerance to check if A1 == -A2
        //     // Compute the middle line
        //     m_lin = lin1;  
        //     b_middle = static_cast<double>((b_p(0) - b_p(1)) / 2.0);
        // } else {
        //     std::cout << "The lines are not parallel, middle line calculation is not applicable."
        //             << std::endl;
        //     m_lin = (lin1 + lin2)/2.0;
        //     b_middle = static_cast<double>((b_p(0) + b_p(1)) / 2.0);
        // }

        std::cout<< m_lin << " b " << b_middle << std::endl;
        // Line parameters: Ax + b = 0
        Eigen::Vector2d A = m_lin; // Normal vector
        double b = -b_middle; // Offset

        // Compute direction vector perpendicular to A
        Eigen::Vector2d direction(-A(1), A(0)); // Rotate A by 90 degrees
        direction.normalize();           // Normalize the direction vector

        // Start point and step size
        Eigen::Vector2d target_point(x0_guess[0], x0_guess[n+1]); // Starting point on the line
        double step_size = 0.1;         // Step size for interpolation
        int max_points = n;           // Maximum number of points to generate
        selected_points.clear();
        
        Eigen::RowVectorXd A_(2);
        A_ << A(0), A(1);
        Eigen::Vector2d  start_point = closestPointToLine(A_, b_middle, target_point);
        // Container to store the points
        selected_points.push_back(start_point); // Add the starting point
        // Generate points along the line
        Eigen::Vector2d current_point = start_point;
        for (int i = 0; i < max_points; ++i) {
            // Compute the next point
            Eigen::Vector2d next_point = current_point + step_size * direction;
            // Check the constraint: Ax + b < 0
            // if (A.dot(next_point) + b >= 0) {
            //     std::cout << " Stop if the constraint is violated" << std::endl;
            //     break; // 
            // }
            // Add the next point to the list
            selected_points.push_back(next_point);
            // Update the current point
            current_point = next_point;
        }
    }

    Eigen::VectorXd closestPointToLine(const Eigen::RowVectorXd& A, double b, const Eigen::VectorXd& p) {
        // Ensure A and p are compatible
        assert(A.cols() == p.size() && "Dimension mismatch between A and p");

        // Compute the projection scalar
        double scale = (A.dot(p) - b) / A.squaredNorm();

        // Compute the closest point
        Eigen::VectorXd q = p - A.transpose() * scale;

        return q;
    }

    void add_initial_constraints(const MX& x, const MX& y, const MX& theta) {
        g.push_back(x(0) - x_list[0]);
        lb_g.push_back(0);
        ub_g.push_back(0);

        g.push_back(y(0) - y_list[0]);
        lb_g.push_back(0);
        ub_g.push_back(0);

        g.push_back(theta(0) - theta_list[0]);
        lb_g.push_back(0);
        ub_g.push_back(0);
    }

    void add_kinematic_constraints(const MX& x, const MX& y, const MX& theta, const MX& v, const MX& delta) {
        for (int i = 0; i < n - 2; ++i) {
            g.push_back(x(i + 1) - (x(i) + delta_t * v(i) * cos(theta(i))));
            lb_g.push_back(0.0);
            ub_g.push_back(0.0);

            g.push_back(y(i + 1) - (y(i) + delta_t * v(i) * sin(theta(i))));
            lb_g.push_back(0.0);
            ub_g.push_back(0.0);

            g.push_back(theta(i + 1) - (theta(i) + delta_t * v(i) * (tan(delta(i)) / wheelbase)));
            lb_g.push_back(0.0);
            ub_g.push_back(0.0);
        }
    }

    void add_velocity_and_steering_constraints(const MX& v, const MX& delta) {
        for (int i = 0; i < n - 1; ++i) {
            g.push_back(v(i));
            lb_g.push_back(0.0);
            ub_g.push_back(v_max);

            g.push_back(delta(i));
            lb_g.push_back(-max_steering_angle);
            ub_g.push_back(max_steering_angle);
        }
    }

    void add_boundary_constraints(const MX& x, const MX& y) {
        for (int i = 0; i < n - 1; ++i) {
            for (int j = 0; j < 2; ++j) {
                g.push_back(x(i) * Ab_(j, 0) + y(i) * Ab_(j, 1));
                lb_g.push_back(-casadi::inf);
                ub_g.push_back(Ab_(j, 2)); // Cast to double
            }
        }
    }

    void add_distance_constraints(const MX& x, const MX& y, const MX& d) {
        for (int i = 0; i < n - 2; ++i) {
            g.push_back(pow(d(i), 2) - (pow(x(i + 1) - x(i), 2) + pow(y(i + 1) - y(i), 2)));
            lb_g.push_back(0.0);
            ub_g.push_back(0.0);
        }
    }

    std::map<string, DM> optimize() {
        // Variables
        MX x = MX::sym("x", n - 1);
        MX y = MX::sym("y", n - 1);
        MX theta = MX::sym("theta", n - 1);
        MX v = MX::sym("v", n - 1);
        MX delta = MX::sym("delta", n - 1);
        MX d = MX::sym("d", n - 2);

        // Add constraints
        add_initial_constraints(x, y, theta);
        add_kinematic_constraints(x, y, theta, v, delta);
        add_velocity_and_steering_constraints(v, delta);
        add_boundary_constraints(x, y);
        // add_distance_constraints(x, y, d);

        // Objective function
        MX objective = 0.0;
        // for (int i = 0; i < n - 3; ++i) {
        //     MX dds_x = x(i + 2) - 2 * x(i + 1) + x(i);
        //     MX dds_y = y(i + 2) - 2 * y(i + 1) + y(i);
        //     objective += curvature_weight * (dds_x * dds_x + dds_y * dds_y);
        // }
        // for (int i = 0; i < n - 4; ++i) {
        //     MX ddds_x = x(i + 3) - 3 * x(i + 2) + 3 * x(i + 1) - x(i);
        //     MX ddds_y = y(i + 3) - 3 * y(i + 2) + 3 * y(i + 1) - y(i);
        //     objective += curvature_rate_weight * (ddds_x * ddds_x + ddds_y * ddds_y);
        // }

        

        // std::cout<< m_lin << std::endl;
        // std::cout<< lin1 << std::endl;
        // std::cout<< A_p(0, 0) << "," << A_p(0, 1) << std::endl;
        // // std::cout<< A_p << std::endl;
        // std::cout<< lin2 << std::endl;
        // std::cout<< A_p(1, 0) << "," << A_p(1, 1) << std::endl;
    
        for (int i = 0; i < n - 1; ++i) {
            // auto distance = (static_cast<double>(m_lin(0).scalar())*x(i) + static_cast<double>(m_lin(1).scalar())*y(i) 
            //             -  b_middle) / std::sqrt(std::pow(static_cast<double>(m_lin(0).scalar()), 2) 
            //             + std::pow(static_cast<double>(m_lin(1).scalar()), 2));
            // objective += distance*distance;
            objective += (x(i)-selected_points[i][0])*(x(i)-selected_points[i][0]);
            objective += (y(i)-selected_points[i][1])*(y(i)-selected_points[i][1]);
        }

        

        // for (int i = 0; i < n - 2; ++i) {
        //     MX distance_sq = MX::pow(x(i + 1) - x(i), 2) + MX::pow((y(i + 1) - y(i)), 2);
        //     objective += distance_sq + distance_weight * d(i);
        // }

        // NLP solver
        MX nlp_vars = vertcat(x, y, theta, v, delta, d);
        MXDict nlp = {{"x", nlp_vars}, {"f", objective}, {"g", vertcat(g)}};
        Dict opts;
        opts["ipopt.print_level"] = 0;
        opts["ipopt.tol"] = 1e-1;
        opts["ipopt.hessian_approximation"] = "limited-memory";
        Function solver = nlpsol("nlpsol", "ipopt", nlp, opts);
        std::map<string, DM> arg = {{"x0", x0_guess}, {"lbg", lb_g}, {"ubg", ub_g}};
        return solver(arg);
    }

    void get_feasible_trajectory(){

        std::map<string, DM> solution = optimize();
        // // Extract optimized values
        vector<double> x_opt(solution["x"](Slice(0, n - 1)).nonzeros());
        vector<double> y_opt(solution["x"](Slice(n - 1, 2 * (n - 1))).nonzeros());
        vector<double> theta_opt(solution["x"](Slice(2 * (n - 1), 3 * (n - 1))).nonzeros());
    
        // std::cout << "Optimized Path: \n";
        // for (int i = 0; i < x_opt.size(); ++i) {
        //     std::cout << "X: " << x_opt[i] << ", Y: " << y_opt[i] << std::endl;
        // }

        // std::vector<double> x_values, y_values;
        // for (int i = 0; i < data.size1(); ++i) {
        //     x_values.push_back(static_cast<double>(data(i, 0)));  // First column
        //     y_values.push_back(static_cast<double>(data(i, 1)));  // Second column
        // }
        // Generate x values from -5 to 5, with 400 points
        int num_points = 400;
        std::vector<double> x_vals(num_points);
        for (int i = 0; i < num_points; ++i) {
            x_vals[i] = -5.0 + i * (10.0 / (num_points - 1));  // from -5 to 5
        }

        // Calculate y values for the lines based on A_p and b_p
        std::vector<double> y1_vals(num_points);
        std::vector<double> y2_vals(num_points);
        std::vector<double> y3_vals(num_points);


        for (int i = 0; i < num_points; ++i) {
            y1_vals[i] = (static_cast<double>(Ab_(0, 2)) -  static_cast<double>(Ab_(0, 0)) * x_vals[i]) /  static_cast<double>(Ab_(0, 1));
            y2_vals[i] = (static_cast<double>(Ab_(1, 2)) -  static_cast<double>(Ab_(1, 0)) * x_vals[i]) /  static_cast<double>(Ab_(1, 1));
            y3_vals[i] = (b_middle -  static_cast<double>(m_lin(0)) * x_vals[i]) /  static_cast<double>(m_lin(1));
        }

        std::vector<double> middle_x, middle_y;
        for(auto p : selected_points){
            middle_x.push_back(p.x());
            middle_y.push_back(p.y());
        }

        // Plot the lines with labels
        plt::figure();
        // plt::rcparam("font.size", 25);
        // plt::plot(x_vals, y1_vals, "r-", {{"label", "y = x^2"}});
        plt::named_plot("line boundary constraints ", x_vals, y1_vals, "r-");
        plt::named_plot("line boundary constraints ", x_vals, y2_vals, "b-");
        plt::named_plot("center line", x_vals, y3_vals, "k-");
        plt::named_plot("initial trajectory", x_list, y_list, "bo--");
        plt::named_plot("refined trajectory", x_opt, y_opt, "ro-");
        plt::named_plot("reference trajectory", middle_x, middle_y, "g");
        // plt::plot(x_vals, y2_vals, "b-");
        // plt::plot(x_vals, y3_vals, "k-");
        // plt::plot(x_list, y_list, "bo--");  // Blue line with circle markers
        // plt::plot(x_opt, y_opt, "ro-");
        // plt::plot(middle_x, middle_y, "g");
        // plt::legend();   
        // plt::legend();
        // plt::title("Trajectory Optimization");
        plt::xlabel("X (m)");
        plt::ylabel("Y (m)");
        plt::show();

    }
};




int main(){

    // DM data = DM({
    // {0.0714334, 0.00638967, 1.66001},
    // {0.171036, 0.015299, 1.66001},
    // {0.270638, 0.0242084, 1.66001},
    // {0.37024, 0.0331177, 1.66001},
    // {0.469843, 0.0420271, 1.66001},
    // {0.569445, 0.0509365, 1.66001},
    // {0.669047, 0.0598458, 1.66001},
    // {0.76865, 0.0687552, 1.66001},
    // {0.868252, 0.0776645, 1.66001},
    // {0.967854, 0.0865739, 1.66001},
    // {1.06746, 0.0954833, 1.66001},
    // {1.16706, 0.104393, 1.66001},
    // {1.26666, 0.113302, 1.66001},
    // {1.36626, 0.122211, 1.66001},
    // {1.46587, 0.131121, 1.66001},
    // {1.56547, 0.14003, 1.66001},
    // {1.66507, 0.148939, 1.66001},
    // {1.76467, 0.157849, 1.66001},
    // {1.86428, 0.166758, 1.66001},
    // {1.96388, 0.175668, 1.66001},
    // {2.06348, 0.184577, 1.66001},
    // {2.16308, 3.193486, 1.66001},
    // {2.26268, 0.202396, 1.66001},
    // {2.36229, 0.211305, 1.66001},
    // {2.46189, 3.220214, 1.66001},
    // {2.56149, 0.229124, 1.66001},
    // {2.66109, 0.238033, 1.66001},
    // {2.7607, 0.246942, 1.66001},
    // {2.8603, 0.255852, 1.66001},
    // {2.9599, 0.264761, 1.66001},
    // {3.0595, 0.27367, 1.66001},
    // {3.15911, 3.28258, 1.66001},
    // {3.25871, 0.291489, 1.66001},
    // {3.35831, 0.300399, 1.66001},
    // {3.45791, 0.309308, 1.66001},
    // {3.55751, 3.318217, 1.66001},
    // {3.65712, 0.327127, 1.66001},
    // {3.75672, 0.336036, 1.66001},
    // {3.85632, 0.444945, 1.66001},
    // {3.95592, 0.453855, 1.66001},
    // {4.05553, 0.462764, 1.66001},
    // {4.15513, 0.471673, 1.66001}});


    Eigen::MatrixXd data(40, 4);
    data << 0.0806606,0.0445364,0.504489,0,
        0.168203,0.0928725,0.504489,0.1,
        0.255745,0.141209,0.504489,0.1,
        0.343287,0.189545,0.504489,0.1,
        0.430829,0.237881,0.504489,0.1,
        0.518371,0.286217,0.504489,0.1,
        0.605913,0.334553,0.504489,0.1,
        0.693456,0.382889,0.504489,0.1,
        0.780998,0.431225,0.504489,0.1,
        0.86854,0.479561,0.504489,0.1,
        0.956082,0.527897,0.504489,0.1,
        1.04362,0.576233,0.504489,0.1,
        1.13117,0.624569,0.504489,0.1,
        1.21871,0.672905,0.504489,0.1,
        1.30625,0.721241,0.504489,0.1,
        1.39379,0.769577,0.504489,0.1,
        1.48133,0.817913,0.504489,0.1,
        1.56888,0.866249,0.504489,0.1,
        1.65642,0.914585,0.504489,0.1,
        1.74396,0.962921,0.504489,0.1,
        1.8315,1.01126,0.504489,0.1,
        1.91905,1.05959,0.504489,0.1,
        2.00659,1.10793,0.504489,0.1,
        2.09413,1.15627,0.504489,0.1,
        2.18167,1.2046,0.504489,0.1,
        2.26921,1.25294,0.504489,0.1,
        2.35676,1.30127,0.504489,0.1,
        2.4443,1.34961,0.504489,0.1,
        2.53184,1.39795,0.504489,0.1,
        2.61938,1.44628,0.504489,0.1,
        2.70692,1.49462,0.504489,0.1,
        2.79447,1.54295,0.504489,0.1,
        2.88201,1.59129,0.504489,0.1,
        2.96955,1.63963,0.504489,0.1,
        3.05709,1.68796,0.504489,0.1,
        3.14464,1.7363,0.504489,0.1,
        3.23218,1.78463,0.504489,0.1,
        3.31972,1.83297,0.504489,0.1,
        3.40726,1.88131,0.504489,0.1,
        3.4948,1.92964,0.504489,0.1;

    DM A_p = DM::horzcat({{0.0213858,  -0.999771}, {-0.0213858,   0.999771}});
    DM b_p = DM::vertcat({0.286652, 0.912916});

    Eigen::MatrixXd Ab(2, 3);
    Ab << -0.0645709, -0.997913, 0.737278, 0.0645709, 0.997913, 0.632449 ;


//         Eigen::MatrixXd data(46, 4);
//     data << 0.0715195,0.0292521,0.388248,0,
// 0.164077,0.0671089,0.388248,0.1,
// 0.256634,0.104966,0.388248,0.1,
// 0.349192,0.142822,0.388248,0.1,
// 0.441749,0.180679,0.388248,0.1,
// 0.534306,0.218536,0.388248,0.1,
// 0.626864,0.256393,0.388248,0.1,
// 0.719421,0.294249,0.388248,0.1,
// 0.811978,0.332106,0.388248,0.1,
// 0.904536,0.369963,0.388248,0.1,
// 0.997093,0.40782,0.388248,0.1,
// 1.08965,0.445677,0.388248,0.1,
// 1.18221,0.483533,0.388248,0.1,
// 1.27477,0.52139,0.388248,0.1,
// 1.36732,0.559247,0.388248,0.1,
// 1.45988,0.597104,0.388248,0.1,
// 1.55244,0.63496,0.388248,0.1,
// 1.64499,0.672817,0.388248,0.1,
// 1.73755,0.710674,0.388248,0.1,
// 1.83011,0.748531,0.388248,0.1,
// 1.92267,0.786387,0.388248,0.1,
// 2.01522,0.824244,0.388248,0.1,
// 2.10778,0.862101,0.388248,0.1,
// 2.20034,0.899958,0.388248,0.1,
// 2.2929,0.937814,0.388248,0.1,
// 2.38545,0.975671,0.388248,0.1,
// 2.47801,1.01353,0.388248,0.1,
// 2.57057,1.05138,0.388248,0.1,
// 2.66313,1.08924,0.388248,0.1,
// 2.75568,1.1271,0.388248,0.1,
// 2.84824,1.16496,0.388248,0.1,
// 2.9408,1.20281,0.388248,0.1,
// 3.03336,1.24067,0.388248,0.1,
// 3.12591,1.27853,0.388248,0.1,
// 3.21847,1.31638,0.388248,0.1,
// 3.31103,1.35424,0.388248,0.1,
// 3.40358,1.3921,0.388248,0.1,
// 3.49614,1.42995,0.388248,0.1,
// 3.5887,1.46781,0.388248,0.1,
// 3.68126,1.50567,0.388248,0.1,
// 3.77381,1.54352,0.388248,0.1,
// 3.86637,1.58138,0.388248,0.1,
// 3.95893,1.61924,0.388248,0.1,
// 4.05149,1.65709,0.388248,0.1,
// 4.14404,1.69495,0.388248,0.1,
// 4.2366,1.73281,0.388248,0.1;

//     Eigen::MatrixXd Ab(2, 3);
//     Ab <<  -0.000866591,           -1,     0.466127, 0.0043336 ,    0.999991 ,    0.858569;  


    TrajectoryOptimizer traj_opt;
    // traj_opt.set_ref_trajectory(data);
    traj_opt.set_ref_trajectory(data);
    traj_opt.set_boundary_constraints(Ab);
    traj_opt.get_feasible_trajectory();

    // vector<double> x_list, y_list, theta_list;
    // int n = data.size1();

    // std::cout<< " n " << n << " cols " << data.size2() << std::endl;

    // for (int i = 0; i < n; ++i) {
    //     x_list.push_back(data(i, 0).scalar());
    //     y_list.push_back(data(i, 1).scalar());
    //     theta_list.push_back(data(i, 2).scalar());  // Sample orientation values
    //     std::cout<< " " << data(i, 0).scalar() << " " << data(i, 1).scalar() << " " << data(i, 2).scalar() << std::endl;
    // }

    // double curvature_weight = 1.0;
    // double curvature_rate_weight = 1.0;
    // double distance_weight = 1.0;
    // double max_steering_angle = 0.6;  // 30 degrees in radians
    // double v_max = 0.5;
    // double delta_t = 0.1;
    // double wheelbase = 1.0;

    // Initial fixed pose
    // double x0_val = x_list[0];
    // double y0_val = y_list[0];
    // double theta0_val = theta_list[0];



}
