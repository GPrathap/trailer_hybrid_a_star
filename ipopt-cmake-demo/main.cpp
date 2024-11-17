#include <casadi/casadi.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include "matplotlibcpp.h"
#include <nlohmann/json.hpp>

using namespace casadi;
using namespace std;
using json = nlohmann::json;

namespace plt = matplotlibcpp;


class TrajectoryOptimizer {
private:
    int n; // Number of points
    DM data;
    DM A_p;
    DM b_p;
    double curvature_weight;
    double curvature_rate_weight;
    double distance_weight;
    double max_steering_angle;
    double v_max;
    double delta_t;
    double wheelbase;
    vector<double> x_list, y_list, theta_list;

    vector<MX> g;           // Constraints
    vector<double> lb_g;    // Lower bounds
    vector<double> ub_g;    // Upper bounds

public:
    TrajectoryOptimizer(const DM& data) 
        : data(data), curvature_weight(1.0), curvature_rate_weight(1.0), distance_weight(1.0), 
          max_steering_angle(M_PI / 6.0), v_max(0.5), delta_t(0.1), wheelbase(1.0) {
        n = data.size1();
        A_p = DM::horzcat({{-0.0645709, 0.0645709}, {-0.997913, 0.997913}});
        b_p = DM::vertcat({0.737278, 0.632449});
        for (int i = 0; i < n; ++i) {
            x_list.push_back(data(i, 0).scalar());
            y_list.push_back(data(i, 1).scalar());
            theta_list.push_back(data(i, 2).scalar());
        }
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
                g.push_back(x(i) * A_p(j, 0) + y(i) * A_p(j, 1));
                lb_g.push_back(-casadi::inf);
                ub_g.push_back(static_cast<double>(b_p(j).scalar())); // Cast to double
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
        add_distance_constraints(x, y, d);

        // Objective function
        MX objective = 0.0;
        for (int i = 0; i < n - 3; ++i) {
            MX dds_x = x(i + 2) - 2 * x(i + 1) + x(i);
            MX dds_y = y(i + 2) - 2 * y(i + 1) + y(i);
            objective += curvature_weight * (dds_x * dds_x + dds_y * dds_y);
        }
        for (int i = 0; i < n - 4; ++i) {
            MX ddds_x = x(i + 3) - 3 * x(i + 2) + 3 * x(i + 1) - x(i);
            MX ddds_y = y(i + 3) - 3 * y(i + 2) + 3 * y(i + 1) - y(i);
            objective += curvature_rate_weight * (ddds_x * ddds_x + ddds_y * ddds_y);
        }
        for (int i = 0; i < n - 2; ++i) {
            MX distance_sq = MX::pow(x(i + 1) - x(i), 2) + MX::pow((y(i + 1) - y(i)), 2);
            objective += distance_sq + distance_weight * d(i);
        }

        // NLP solver
        MX nlp_vars = vertcat(x, y, theta, v, delta, d);
        MXDict nlp = {{"x", nlp_vars}, {"f", objective}, {"g", vertcat(g)}};
        Dict opts;
        opts["ipopt.print_level"] = 5;
        opts["ipopt.tol"] = 1e-6;
        Function solver = nlpsol("nlpsol", "ipopt", nlp, opts);

        vector<double> x0_guess(nlp_vars.size1(), 0.0);
        std::map<string, DM> arg = {{"x0", x0_guess}, {"lbg", lb_g}, {"ubg", ub_g}};
        return solver(arg);
    }
};




int main(){

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
    double max_steering_angle = 0.6;  // 30 degrees in radians
    double v_max = 0.5;
    double delta_t = 0.1;
    double wheelbase = 1.0;

    // Initial fixed pose
    double x0_val = x_list[0];
    double y0_val = y_list[0];
    double theta0_val = theta_list[0];



}
