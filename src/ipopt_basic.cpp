#include <IpIpoptApplication.hpp>
#include <IpTNLP.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using Eigen::VectorXd;
using std::vector;

class TrajectoryOptimization : public Ipopt::TNLP {
public:
    size_t n;                  // Number of waypoints
    double curvature_weight;
    double curvature_rate_weight;
    double deviation_weight;
    double delta_t;
    double v_max;
    VectorXd x_init, y_init, d_init, theta_list;
    vector<double> x_opt, y_opt, d_opt;

    TrajectoryOptimization(
        size_t n_,
        double curvature_weight_,
        double curvature_rate_weight_,
        double deviation_weight_,
        double delta_t_,
        double v_max_,
        VectorXd x_init_,
        VectorXd y_init_,
        VectorXd d_init_,
        VectorXd theta_list_)
        : n(n_), curvature_weight(curvature_weight_), curvature_rate_weight(curvature_rate_weight_),
          deviation_weight(deviation_weight_), delta_t(delta_t_), v_max(v_max_),
          x_init(std::move(x_init_)), y_init(std::move(y_init_)), d_init(std::move(d_init_)),
          theta_list(std::move(theta_list_))
    {}

    bool get_nlp_info(
        Ipopt::Index& n_var, Ipopt::Index& n_cons, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
        Ipopt::TNLP::IndexStyleEnum& index_style) override
    {
        n_var = 3 * n;
        n_cons = 4 * n - 2;
        nnz_jac_g = 3 * n; // Assuming sparse Jacobian for constraints
        nnz_h_lag = n * n; // Approximate Hessian nonzeros
        index_style = Ipopt::TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n_var, double* x_l, double* x_u,
        Ipopt::Index n_cons, double* g_l, double* g_u) override
    {
        for (size_t i = 0; i < 3 * n; ++i) {
            x_l[i] = -1e19;
            x_u[i] = 1e19;
        }

        for (size_t i = 0; i < 4 * n - 2; ++i) {
            g_l[i] = 0.0;
            g_u[i] = 0.0;
        }

        for (size_t i = 0; i < n - 1; ++i) {
            g_l[3 * n + i] = -v_max;
            g_u[3 * n + i] = v_max;
            g_l[3 * n + (n - 1) + i] = -v_max;
            g_u[3 * n + (n - 1) + i] = v_max;
        }
        return true;
    }

    bool get_starting_point(
        Ipopt::Index n_var, bool init_x, double* x,
        bool init_z, double* z_L, double* z_U,
        Ipopt::Index n_cons, bool init_lambda, double* lambda) override
    {
        for (size_t i = 0; i < n; ++i) {
            x[i] = x_init(i);
            x[n + i] = y_init(i);
            x[2 * n + i] = d_init(i);
        }
        return true;
    }

    bool eval_f(
        Ipopt::Index n_var, const double* x, bool new_x, double& obj_value) override
    {
        obj_value = 0.0;
        for (size_t i = 0; i < n - 2; ++i) {
            obj_value += curvature_weight * std::pow(x[i + 2] - 2 * x[i + 1] + x[i], 2);
        }

        for (size_t i = 0; i < n - 3; ++i) {
            obj_value += curvature_rate_weight * std::pow(x[i + 3] - 3 * x[i + 2] + 3 * x[i + 1] - x[i], 2);
        }

        for (size_t i = 0; i < n; ++i) {
            obj_value += deviation_weight * std::pow(x[i + 2 * n], 2);
        }

        return true;
    }

    bool eval_grad_f(
        Ipopt::Index n_var, const double* x, bool new_x, double* grad_f) override
    {
        for (size_t i = 0; i < n_var; ++i) {
            grad_f[i] = 0.0;  // Implement appropriate gradient calculation here
        }
        return true;
    }

    bool eval_g(
        Ipopt::Index n_var, const double* x, bool new_x,
        Ipopt::Index n_cons, double* g) override
    {
        for (size_t i = 0; i < n; ++i) {
            double theta = theta_list(i) + M_PI / 2;
            g[i] = x[i] - x_init(i) - x[2 * n + i] * std::cos(theta);
            g[n + i] = x[n + i] - y_init(i) - x[2 * n + i] * std::sin(theta);
            g[2 * n + i] = x[2 * n + i];
        }

        for (size_t i = 0; i < n - 1; ++i) {
            g[3 * n + i] = (x[i + 1] - x[i]) / delta_t;
            g[3 * n + (n - 1) + i] = (x[n + i + 1] - x[n + i]) / delta_t;
        }
        return true;
    }

    bool eval_jac_g(
        Ipopt::Index n_var, const double* x, bool new_x,
        Ipopt::Index n_cons, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, double* values) override
    {
        // Define Jacobian structure and values here
        return true;
    }

    bool eval_h(
        Ipopt::Index n_var, const double* x, bool new_x, double obj_factor,
        Ipopt::Index n_cons, const double* lambda, bool new_lambda,
        Ipopt::Index nele_hess, Ipopt::Index* iRow, Ipopt::Index* jCol, double* values) override
    {
        // Define Hessian structure and values here
        return true;
    }

    void finalize_solution(
        Ipopt::SolverReturn status, Ipopt::Index n_var, const double* x,
        const double* z_L, const double* z_U, Ipopt::Index n_cons, const double* g,
        const double* lambda, double obj_value, const Ipopt::IpoptData* ip_data,
        Ipopt::IpoptCalculatedQuantities* ip_cq) override
    {
        x_opt = std::vector<double>(x, x + n);
        y_opt = std::vector<double>(x + n, x + 2 * n);
        d_opt = std::vector<double>(x + 2 * n, x + 3 * n);
        plot_trajectory();
    }

    void plot_trajectory() {
        std::vector<double> x_init_v(x_init.data(), x_init.data() + x_init.size());
        std::vector<double> y_init_v(y_init.data(), y_init.data() + y_init.size());

        plt::plot(x_init_v, y_init_v, "r--");
        plt::plot(x_opt, y_opt, "b-");
        plt::legend();
        plt::title("Trajectory Optimization");
        plt::xlabel("X position");
        plt::ylabel("Y position");
        plt::show();
    }
};

int main() {
    size_t n = 45;
    double curvature_weight = 1.0;
    double curvature_rate_weight = 1.0;
    double deviation_weight = 1.0;
    double delta_t = 1.0;
    double v_max = 1.0;

    VectorXd x_init(n), y_init(n), d_init(n), theta_list(n);
    x_init.setRandom();
    y_init.setRandom();
    d_init.setZero();
    theta_list.setRandom();

    auto trajectory_opt = std::make_shared<TrajectoryOptimization>(
        n, curvature_weight, curvature_rate_weight, deviation_weight, delta_t, v_max,
        x_init, y_init, d_init, theta_list);

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    app->Initialize();
    app->OptimizeTNLP(trajectory_opt);

    return 0;
}
