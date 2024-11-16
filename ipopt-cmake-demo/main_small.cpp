    #include <casadi/casadi.hpp>
    #include <iostream>
    #include <fstream>
    #include <vector>

    using namespace casadi;

    int main() {
        // Parameters and fixed initial pose
        std::vector<std::vector<double>> data = {
            {0.0714334, 0.00638967, 1.66001}, {0.171036, 0.015299, 1.66001},
            {0.270638, 0.0242084, 1.66001}, {0.37024, 0.0331177, 1.66001},
            {0.469843, 0.0420271, 1.66001}, {0.569445, 0.0509365, 1.66001},
            {0.669047, 0.0598458, 1.66001}, {0.76865, 0.0687552, 1.66001},
            {0.868252, 0.0776645, 1.66001}, {0.967854, 0.0865739, 1.66001},
            {1.06746, 0.0954833, 1.66001}, {1.16706, 0.104393, 1.66001},
            {1.26666, 0.113302, 1.66001}, {1.36626, 0.122211, 1.66001},
            {1.46587, 0.131121, 1.66001}, {1.56547, 0.14003, 1.66001},
            {1.66507, 0.148939, 1.66001}, {1.76467, 0.157849, 1.66001},
            {1.86428, 0.166758, 1.66001}, {1.96388, 0.175668, 1.66001},
            {2.06348, 0.184577, 1.66001}, {2.16308, 3.193486, 1.66001},
            {2.26268, 0.202396, 1.66001}, {2.36229, 0.211305, 1.66001},
            {2.46189, 3.220214, 1.66001}, {2.56149, 0.229124, 1.66001},
            {2.66109, 0.238033, 1.66001}, {2.7607, 0.246942, 1.66001},
            {2.8603, 0.255852, 1.66001}, {2.9599, 0.264761, 1.66001},
            {3.0595, 0.27367, 1.66001}, {3.15911, 3.28258, 1.66001},
            {3.25871, 0.291489, 1.66001}, {3.35831, 0.300399, 1.66001},
            {3.45791, 0.309308, 1.66001}, {3.55751, 3.318217, 1.66001},
            {3.65712, 0.327127, 1.66001}, {3.75672, 0.336036, 1.66001},
            {3.85632, 0.444945, 1.66001}, {3.95592, 0.453855, 1.66001},
            {4.05553, 0.462764, 1.66001}, {4.15513, 0.471673, 1.66001}
        };

        // DM A_p = DM::diag(DM{-0.0645709, 0.0645709}) * DM{{1, 0}, {0, 1}};
        // DM b_p = DM{0.737278, 0.632449};

        std::vector<double> x_list, y_list, theta_list;
        for (const auto& row : data) {
            x_list.push_back(row[0]);
            y_list.push_back(row[1]);
            theta_list.push_back(row[2]);
        }

        int n = x_list.size();
        double delta_t = 0.1, distance_weight = 1.0, v_max = 1.0, wheelbase = 1.0;

        // Initial pose
        double x0_val = x_list[0], y0_val = y_list[0];

        // Define CasADi variables
        MX x = MX::sym("x", n - 1);
        MX y = MX::sym("y", n - 1);
        MX theta = MX::sym("theta", n - 1);
        MX v = MX::sym("v", n - 1);
        MX delta = MX::sym("delta", n - 1);

        // // Objective function
        MX objective = 0;
        for (int i = 0; i < n - 2; ++i) {
            MX distance_sq = MX::pow(x(i + 1) - x(i), 2) + MX::pow(y(i + 1) - y(i), 2);
            objective += distance_weight * distance_sq;
        }

        // Constraints (empty in this example)
        std::vector<MX> g;
        std::vector<double> lb_g, ub_g;

        // NLP
        // MX vars = MX::vertcat({x, y, theta, v, delta});
        MX vars = vertcat(x, y, theta, v, delta);

        std::cout<< " " << vars.size1() << " " << vars.size2() << std::endl;
        // MX constraints = MX::vertcat(g);
        std::map<std::string, MX> nlp = {{"x", vars}, {"f", objective}};

        Dict opts;
        opts["ipopt.print_level"] = 5;
        opts["ipopt.tol"] = 1e-6;
        opts["ipopt.max_iter"] = 1000;

        Function solver = nlpsol("solver", "ipopt", nlp, opts);

        // Initial guess
        std::vector<double> x0_guess;
        x0_guess.insert(x0_guess.end(), x_list.begin() + 1, x_list.end());
        x0_guess.insert(x0_guess.end(), y_list.begin() + 1, y_list.end());
        x0_guess.insert(x0_guess.end(), theta_list.begin() + 1, theta_list.end());
        x0_guess.insert(x0_guess.end(), n - 1, v_max);
        x0_guess.insert(x0_guess.end(), n - 1, 0);

        std::map<std::string, DM> arg;
        arg["x0"] = x0_guess;

        DMDict solution = solver(arg);

        // // Extract optimized values
        DM x_opt = DM::vertcat({x0_val, solution["x"](Slice(0, n - 1))});
        DM y_opt = DM::vertcat({y0_val, solution["x"](Slice(n - 1, 2 * (n - 1)))});

        std::cout << "Optimized Path:\n";
        std::cout << "x: " << x_opt << "\n";
        std::cout << "y: " << y_opt << "\n";

        return 0;
    }
