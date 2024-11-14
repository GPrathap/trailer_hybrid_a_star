#include <cassert>
#include <iostream>
#include "ipopt_interface_traj.hpp"
#include "IpTypes.hpp"

TrajectoryTracking::TrajectoryTracking()
{

}

//destructor
TrajectoryTracking::~TrajectoryTracking()
{

}

void TrajectoryTracking::init(size_t n_,
        double curvature_weight_,
        double curvature_rate_weight_,
        double deviation_weight_,
        double delta_t_,
        double v_max_,
        VectorXd x_init_,
        VectorXd y_init_,
        VectorXd d_init_,
        VectorXd theta_list_){
          
        init_n = n_;
        curvature_weight = curvature_weight_;
        curvature_rate_weight = curvature_rate_weight_;
        deviation_weight = deviation_weight_;
        delta_t = delta_t_;
        v_max = v_max_;
        x_init = x_init_;
        y_init = y_init_;
        d_init = d_init_;
        theta_list = theta_list_;
}

// returns the size of the problem
bool TrajectoryTracking::get_nlp_info(Ipopt::Index& n_var, 
                               Ipopt::Index& n_cons, 
                               Ipopt::Index& nnz_jac_g,
                               Ipopt::Index& nnz_h_lag, 
                               Ipopt::TNLP::IndexStyleEnum& index_style){

  // // The problem described in HS071_NLP.hpp has 4 variables, x[0] through x[3]
  // n_var = 4;
  
  // // one equality constraint and one inequality constraint
  // n_cons = 2;
  
  // // in this example the jacobian is dense and contains 8 nonzeros
  // nnz_jac_g = 8;
  
  // // the hessian is also dense and has 16 total nonzeros, but we
  // // only need the lower left corner (since it is symmetric)
  // nnz_h_lag = 10;
  
  // // use the C style indexing (0-based)
  // index_style = Ipopt::TNLP::C_STYLE;
  std::cout<< " get_nlp_info -----> " << init_n << std::endl;
  n_var = 3 * init_n;
  n_cons = 4 * init_n - 2;
  nnz_jac_g = 3 * init_n; // Assuming sparse Jacobian for constraints
  nnz_h_lag = 0; // Approximate Hessian nonzeros
  index_style = Ipopt::TNLP::C_STYLE;
  std::cout<< " get_nlp_info -----> x " << n_var << " n_cons " << n_cons << " nnz_jac_g " << nnz_jac_g << " nnz_h_lag " << nnz_h_lag << std::endl;
  return true;
}


// returns the variable bounds
bool TrajectoryTracking::get_bounds_info(Ipopt::Index n_var, 
                                  Ipopt::Number* x_l, 
                                  Ipopt::Number* x_u,
                                  Ipopt::Index n_cons, 
                                  Ipopt::Number* g_l, 
                                  Ipopt::Number* g_u)
{
  std::cout<< " get_bounds_info -----> " << n_var  << " " << n_cons << " " << 4 * init_n - 2 << std::endl;
  
  for (Ipopt::Index i = 0; i < n_var; ++i) {
            x_l[i] = -1000.0;
            x_u[i] = 10000.0;
            // std::cout<< "  x_u " << i << std::endl;
  }

  for (Ipopt::Index i = 0; i < n_cons; ++i) {
      g_l[i] = 0.0;
      g_u[i] = 0.0;
  }

  for (size_t i = 0; i < init_n - 1; ++i) {
      g_l[3 * init_n + i] = -v_max;
      g_u[3 * init_n + i] = v_max;
      g_l[3 * init_n + (init_n - 1) + i] = -v_max;
      g_u[3 * init_n + (init_n - 1) + i] = v_max;
      std::cout<< i << std::endl;
  }
  std::cout<< " get_bounds_info -----> " << std::endl;
  
  return true;
}



// returns the initial point for the problem
bool TrajectoryTracking::get_starting_point(Ipopt::Index n_var, 
                                     bool init_x, 
                                     Ipopt::Number* x,
                                     bool init_z, 
                                     Ipopt::Number* z_L, 
                                     Ipopt::Number* z_U,
                                     Ipopt::Index n_cons, 
                                     bool init_lambda,
                                     Ipopt::Number* lambda)
{
  // // Here, we assume we only have starting values for x, if you code
  // // your own NLP, you can provide starting values for the dual variables
  // // if you wish
  // assert(init_x == true);
  // assert(init_z == false);
  // assert(init_lambda == false);
  
  // // initialize to the given starting point
  // x[0] = 1.0;
  // x[1] = 5.0;
  // x[2] = 5.0;
  // x[3] = 1.0;
  std::cout<< " get_starting_point -----> " << std::endl;
  // return true;
  for (size_t i = 0; i < init_n; ++i) {
            x[i] = x_init(i);
            x[init_n + i] = y_init(i);
            x[2 * init_n + i] = d_init(i);
        }
  std::cout<< " get_starting_point -----> " << std::endl;
  return true;
}

// returns the value of the objective function
bool TrajectoryTracking::eval_f(Ipopt::Index n_var, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value){
  // assert(n == 4);
  // obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
  // return true;
  std::cout<< " eval_f -----> " << std::endl;
  obj_value = 0.0;
  for (size_t i = 0; i < init_n - 2; ++i) {
      obj_value += curvature_weight * std::pow(x[i + 2] - 2 * x[i + 1] + x[i], 2);
  }

  for (size_t i = 0; i < init_n - 3; ++i) {
      obj_value += curvature_rate_weight * std::pow(x[i + 3] - 3 * x[i + 2] + 3 * x[i + 1] - x[i], 2);
  }

  for (size_t i = 0; i < init_n; ++i) {
      obj_value += deviation_weight * std::pow(x[i + 2 * init_n], 2);
  }
  std::cout<< " eval_f -----> " << std::endl;
  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool TrajectoryTracking::eval_grad_f(Ipopt::Index n_var, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
  // assert(n == 4);
  
  // grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
  // grad_f[1] = x[0] * x[3];
  // grad_f[2] = x[0] * x[3] + 1;
  // grad_f[3] = x[0] * (x[0] + x[1] + x[2]);
  std::cout<< " eval_grad_f -----> " << std::endl;
  // return true;
  // for (size_t i = 0; i < n_var; ++i) {
  //     grad_f[i] = 0.0;  // Implement appropriate gradient calculation here
  // }
  std::cout<< " eval_grad_f -----> " << std::endl;
  return true;

}


// return the value of the constraints: g(x)
bool TrajectoryTracking::eval_g(Ipopt::Index n_var, 
                         const Ipopt::Number* x, 
                         bool new_x, 
                         Ipopt::Index n_cons, 
                         Ipopt::Number* g)
{
  // assert(n == 4);
  // assert(m == 2);
  
  // g[0] = x[0] * x[1] * x[2] * x[3];
  // g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
  std::cout<< " eval_g -----> " << std::endl;
  // return true;
   for (size_t i = 0; i < init_n; ++i) {
        double theta = theta_list(i) + M_PI / 2;
        g[i] = x[i] - x_init(i) - x[2 * init_n + i] * std::cos(theta);
        g[init_n + i] = x[init_n + i] - y_init(i) - x[2 * init_n + i] * std::sin(theta);
        g[2 * init_n + i] = x[2 * init_n + i];
    }

    for (size_t i = 0; i < init_n - 1; ++i) {
        g[3 * init_n + i] = (x[i + 1] - x[i]) / delta_t;
        g[3 * init_n + (init_n - 1) + i] = (x[init_n + i + 1] - x[init_n + i]) / delta_t;
    }
    std::cout<< " eval_g -----> " << std::endl;
    return true;
}

// return the structure or values of the jacobian
bool TrajectoryTracking::eval_jac_g(Ipopt::Index n, 
                             const Ipopt::Number* x, bool new_x,
                             Ipopt::Index m, 
                             Ipopt::Index nele_jac, 
                             Ipopt::Index* iRow, 
                             Ipopt::Index *jCol,
                             Ipopt::Number* values)
{
  
  // for (Ipopt::Index i = 0; i < n; ++i) {
  //       values[i] = 0.0;
  // }
  std::cout<< " eval_jac_g -----> " << n  << " " << init_n << std::endl;
  return true;
}

//return the structure or values of the hessian
bool TrajectoryTracking::eval_h(Ipopt::Index n, 
                         const Ipopt::Number* x, 
                         bool new_x,
                         Ipopt::Number obj_factor, 
                         Ipopt::Index m, 
                         const Ipopt::Number* lambda,
                         bool new_lambda, 
                         Ipopt::Index nele_hess, 
                         Ipopt::Index* iRow,
                         Ipopt::Index* jCol, 
                         Ipopt::Number* values)
{
  // if (values == NULL) {
  //   // return the structure. This is a symmetric matrix, fill the lower left
  //   // triangle only.
    
  //   // the hessian for this problem is actually dense
  //   Ipopt::Index idx=0;
  //   for (Ipopt::Index row = 0; row < 4; row++) {
  //     for (Ipopt::Index col = 0; col <= row; col++) {
  //       iRow[idx] = row;
  //       jCol[idx] = col;
  //       idx++;
  //     }
  //   }
    
  //   assert(idx == nele_hess);
  // }
  // else {
  //   // return the values. This is a symmetric matrix, fill the lower left
  //   // triangle only
    
  //   // fill the objective portion
  //   values[0] = obj_factor * (2*x[3]); // 0,0
    
  //   values[1] = obj_factor * (x[3]);   // 1,0
  //   values[2] = 0.;                    // 1,1
    
  //   values[3] = obj_factor * (x[3]);   // 2,0
  //   values[4] = 0.;                    // 2,1
  //   values[5] = 0.;                    // 2,2
    
  //   values[6] = obj_factor * (2*x[0] + x[1] + x[2]); // 3,0
  //   values[7] = obj_factor * (x[0]);                 // 3,1
  //   values[8] = obj_factor * (x[0]);                 // 3,2
  //   values[9] = 0.;                                  // 3,3
    
    
  //   // add the portion for the first constraint
  //   values[1] += lambda[0] * (x[2] * x[3]); // 1,0
    
  //   values[3] += lambda[0] * (x[1] * x[3]); // 2,0
  //   values[4] += lambda[0] * (x[0] * x[3]); // 2,1
    
  //   values[6] += lambda[0] * (x[1] * x[2]); // 3,0
  //   values[7] += lambda[0] * (x[0] * x[2]); // 3,1
  //   values[8] += lambda[0] * (x[0] * x[1]); // 3,2
    
  //   // add the portion for the second constraint
  //   values[0] += lambda[1] * 2; // 0,0
    
  //   values[2] += lambda[1] * 2; // 1,1
    
  //   values[5] += lambda[1] * 2; // 2,2
    
  //   values[9] += lambda[1] * 2; // 3,3
  // }
  std::cout<< " eval_h -----> " << n << std::endl;
  return true;
}


void TrajectoryTracking::finalize_solution(Ipopt::SolverReturn status,
                                    Ipopt::Index n, 
                                    const Ipopt::Number* x, 
                                    const Ipopt::Number* z_L, 
                                    const Ipopt::Number* z_U,
                                    Ipopt::Index m,
                                    const Ipopt::Number* g, 
                                    const Ipopt::Number* lambda,
                                    Ipopt::Number obj_value,
                                    const Ipopt::IpoptData* ip_data,
                                    Ipopt::IpoptCalculatedQuantities* ip_cq)
{
  // // here is where we would store the solution to variables, or write to a file, etc
  // // so we could use the solution.
  
  // // For this example, we write the solution to the console
  // std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
  // for (Ipopt::Index i=0; i<n; i++) {
  //   std::cout << "x[" << i << "] = " << x[i] << std::endl;
  // }
  
  // std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
  // for (Ipopt::Index i=0; i<n; i++) {
  //   std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  // }
  // for (Ipopt::Index i=0; i<n; i++) {
  //   std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  // }
  
  // std::cout << std::endl << std::endl << "Objective value" << std::endl;
  // std::cout << "f(x*) = " << obj_value << std::endl;
  
  // std::cout << std::endl << "Final value of the constraints:" << std::endl;
  // for (Ipopt::Index i=0; i<m ;i++) {
  //   std::cout << "g(" << i << ") = " << g[i] << std::endl;
  // }
  std::cout<< " finalize_solution -----> " << std::endl;
  x_opt = std::vector<double>(x, x + n);
  y_opt = std::vector<double>(x + n, x + 2 * n);
  d_opt = std::vector<double>(x + 2 * n, x + 3 * n);
  std::cout<< " finalize_solution -----> " << std::endl;
  // plot_trajectory();

}

void TrajectoryTracking::plot_trajectory() {
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
