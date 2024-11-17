    #include <casadi/casadi.hpp>
    #include <vector>
    #include <iostream>
    #include <cmath>
    #include "matplotlibcpp.h"
    #include <nlohmann/json.hpp>

    using namespace casadi;
    using namespace std;
    namespace plt = matplotlibcpp;

    using json = nlohmann::json;

    int main() {
        // Parameters and fixed initial pose
        DM data = DM({
        {0.0714334, 0.00638967, 1.66001},
        {0.171036, 0.015299, 1.66001},
        {0.270638, 0.0242084, 1.66001},
        {0.37024, 0.0331177, 1.66001},
        {0.469843, 0.0420271, 1.66001},
        {0.569445, 0.0509365, 1.66001},
        {0.669047, 0.0598458, 1.66001},
        {0.76865, 0.0687552, 1.66001},
        {0.868252, 0.0776645, 1.66001},
        {0.967854, 0.0865739, 1.66001},
        {1.06746, 0.0954833, 1.66001},
        {1.16706, 0.104393, 1.66001},
        {1.26666, 0.113302, 1.66001},
        {1.36626, 0.122211, 1.66001},
        {1.46587, 0.131121, 1.66001},
        {1.56547, 0.14003, 1.66001},
        {1.66507, 0.148939, 1.66001},
        {1.76467, 0.157849, 1.66001},
        {1.86428, 0.166758, 1.66001},
        {1.96388, 0.175668, 1.66001},
        {2.06348, 0.184577, 1.66001},
        {2.16308, 3.193486, 1.66001},
        {2.26268, 0.202396, 1.66001},
        {2.36229, 0.211305, 1.66001},
        {2.46189, 3.220214, 1.66001},
        {2.56149, 0.229124, 1.66001},
        {2.66109, 0.238033, 1.66001},
        {2.7607, 0.246942, 1.66001},
        {2.8603, 0.255852, 1.66001},
        {2.9599, 0.264761, 1.66001},
        {3.0595, 0.27367, 1.66001},
        {3.15911, 3.28258, 1.66001},
        {3.25871, 0.291489, 1.66001},
        {3.35831, 0.300399, 1.66001},
        {3.45791, 0.309308, 1.66001},
        {3.55751, 3.318217, 1.66001},
        {3.65712, 0.327127, 1.66001},
        {3.75672, 0.336036, 1.66001},
        {3.85632, 0.444945, 1.66001},
        {3.95592, 0.453855, 1.66001},
        {4.05553, 0.462764, 1.66001},
        {4.15513, 0.471673, 1.66001}});

        DM A_p = DM::horzcat({{-0.0645709, 0.0645709}, {-0.997913, 0.997913}});
        DM b_p = DM::vertcat({0.737278, 0.632449});

        vector<double> x_list, y_list, theta_list;
        int n = data.size1();

        std::cout<< " n " << n << " cols " << data.size2() << std::endl;

        for (int i = 0; i < n; ++i) {
            x_list.push_back(data(i, 0).scalar());
            y_list.push_back(data(i, 1).scalar());
            theta_list.push_back(data(i, 2).scalar());  // Sample orientation values
            std::cout<< " " << data(i, 0).scalar() << " " << data(i, 1).scalar() << " " << data(i, 2).scalar() << std::endl;
        }

        double curvature_weight = 1.0;
        double curvature_rate_weight = 1.0;
        double distance_weight = 1.0;
        double v_max = 1.0;
        double delta_t = 0.1;
        double wheelbase = 1.0;

        // Initial fixed pose
        double x0_val = x_list[0];
        double y0_val = y_list[0];
        double theta0_val = theta_list[0];

        // Define optimization variables for remaining poses, velocities, and steering angles
        MX x = MX::sym("x", n - 1);
        MX y = MX::sym("y", n - 1);
        MX theta = MX::sym("theta", n - 1);
        MX v = MX::sym("v", n - 1);
        MX delta = MX::sym("delta", n - 1);

        // Objective function
        MX objective = 0.0;
        
        // for (int i = 0; i < n - 3; ++i) {
        //     MX dds_x = x(i+2) - 2 * x(i+1) + x(i);
        //     MX dds_y = y(i+2) - 2 * y(i+1) + y(i);
        //     objective += curvature_weight * (dds_x*dds_x + dds_y*dds_y);
        // }
        
        // for (int i = 0; i < n - 4; ++i) {
        //     MX ddds_x = x(i+3) - 3 * x(i+2) + 3 * x(i+1) - x(i);
        //     MX ddds_y = y(i+3) - 3 * y(i+2) + 3 * y(i+1) - y(i);
        //     objective += curvature_rate_weight * (ddds_x*ddds_x + ddds_y*ddds_y);
        // }

        // Penalty for distance between consecutive points
        // for (int i = 0; i < n - 2; ++i) {
        //     MX distance_sq = MX::pow(x(i+1) - x(i), 2) + MX::pow((y(i+1) - y(i)), 2);
        //     std::cout<< distance_sq << std::endl;
        //     objective += distance_sq;
        // }

        // Constraints
        vector<MX> g;
        vector<double> lb_g, ub_g;

        // Kinematic constraints for car-like dynamics
        for (int i = 0; i < n - 2; ++i) {
            g.push_back(x(i+1) - (x(i) + delta_t * v(i) * cos(theta(i))));
            lb_g.push_back(0.0);
            ub_g.push_back(0.0);

            g.push_back(y(i+1) - (y(i) + delta_t * v(i) * sin(theta(i))));
            lb_g.push_back(0.0);
            ub_g.push_back(0.0);

            g.push_back(theta(i+1) - (theta(i) + delta_t * v(i) * (tan(delta(i)) / wheelbase)));
            lb_g.push_back(0.0);
            ub_g.push_back(0.0);
        }

        // Velocity and steering angle constraints
        double max_steering_angle = M_PI / 6.0;  // 30 degrees in radians
        // double max_steering_angle = 0.1;  // 30 degrees in radians
        std::cout<< " max_steering_angle " << max_steering_angle  << std::endl;
        for (int i = 0; i < n - 1; ++i) {
            g.push_back(v(i));
            lb_g.push_back(0.0);
            ub_g.push_back(v_max);

            g.push_back(delta(i));
            lb_g.push_back(-max_steering_angle);
            ub_g.push_back(max_steering_angle);
        }

        // // Boundary constraints
        // for (int i = 0; i < n - 1; ++i) {
        //     for (int j = 0; j < 2; ++j) {
        //         g.push_back(x(i) * A_p(j, 0) + y(i) * A_p(j, 1));
        //         lb_g.push_back(-casadi::inf);
        //         ub_g.push_back(static_cast<double>(b_p(j).scalar())); // Cast to double
        //     }
        // }



        // Define NLP problem without the fixed initial pose in the variables
        MX nlp_vars = vertcat(x, y, theta, v, delta);

        // std::cout<< " " << vars.size1() << " " << vars.size2() << std::endl;
        // Set up IPOPT solver
        Dict opts;
        opts["ipopt.print_level"] = 5;
        opts["print_time"] = false;
        opts["ipopt.tol"] = 1e-6;
        opts["ipopt.max_iter"] = 1000;
        // opts["ipopt.constr_viol_tol"] = 1e-4;
        
        auto g_cons = vertcat(g);
        MXDict nlp = {{"x", nlp_vars}, {"f", objective}, {"g", g_cons}};
        Function solver = nlpsol("nlpsol", "ipopt", nlp, opts);

        // Initial guess (excluding the fixed initial pose)
        vector<double> x0_guess;
        for (int i = 1; i < n; ++i){
            x0_guess.push_back(x_list[i]);
        } 
        for (int i = 1; i < n; ++i) {
            x0_guess.push_back(y_list[i]);
        }
        for (int i = 1; i < n; ++i){
            x0_guess.push_back(theta_list[i]);
        } 
        for (int i = 0; i < n - 1; ++i) {
            x0_guess.push_back(v_max);
        }
        for (int i = 0; i < n - 1; ++i) {
            x0_guess.push_back(0.0);
        }

        vector<double> lbx, ubx;

        std::map<std::string, DM> arg;
        arg["lbg"] = lb_g;
        arg["ubg"] = ub_g;
        arg["lbx"] = -DM::inf(); // Lower bounds
        arg["ubx"] = DM::inf();  // Upper bounds
        arg["x0"] = x0_guess;
        // Solve the problem
        std::map<string, DM> solution = solver(arg);

        // // Extract optimized values
        vector<double> x_opt(solution["x"](Slice(0, n - 1)).nonzeros());
        vector<double> y_opt(solution["x"](Slice(n - 1, 2 * (n - 1))).nonzeros());
        vector<double> theta_opt(solution["x"](Slice(2 * (n - 1), 3 * (n - 1))).nonzeros());
        vector<double> v_opt(solution["x"](Slice(3 * (n - 1), 4 * (n - 1))).nonzeros());
        vector<double> delta_opt(solution["x"](Slice(4 * (n - 1), 5 * (n - 1))).nonzeros());

        cout << "Optimized Path: \n";
        for (int i = 0; i < x_opt.size(); ++i) {
            cout << "X: " << x_opt[i] << ", Y: " << y_opt[i] << endl;
        }

        std::vector<double> x_values, y_values;
        for (int i = 0; i < data.size1(); ++i) {
            x_values.push_back(static_cast<double>(data(i, 0)));  // First column
            y_values.push_back(static_cast<double>(data(i, 1)));  // Second column
        }

        // Generate x values from -5 to 5, with 400 points
        int num_points = 400;
        std::vector<double> x_vals(num_points);
        for (int i = 0; i < num_points; ++i) {
            x_vals[i] = -5.0 + i * (10.0 / (num_points - 1));  // from -5 to 5
        }

        // Calculate y values for the lines based on A_p and b_p
        std::vector<double> y1_vals(num_points);
        std::vector<double> y2_vals(num_points);

        for (int i = 0; i < num_points; ++i) {
            y1_vals[i] = (static_cast<double>(b_p(0)) -  static_cast<double>(A_p(0, 0)) * x_vals[i]) /  static_cast<double>(A_p(0, 1));
            y2_vals[i] = (static_cast<double>(b_p(1)) -  static_cast<double>(A_p(1, 0)) * x_vals[i]) /  static_cast<double>(A_p(1, 1));
        }

        // Plot the lines with labels
        plt::figure();
        plt::plot(x_vals, y1_vals, "r-");
        plt::plot(x_vals, y2_vals, "b-");
        plt::plot(x_values, y_values, "bo--");  // Blue line with circle markers
        plt::plot(x_opt, y_opt, "ro-");
        plt::legend();
        plt::title("Trajectory Optimization");
        plt::xlabel("X position");
        plt::ylabel("Y position");
        plt::show();

        return 0;
    }
