// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10
// Modified: brian paden Aug-2017

#include <IpIpoptApplication.hpp>
#include <iostream>

#include "ipopt_interface_traj.hpp"
// #include "ipopt_interface.hpp"


// using namespace Ipopt;

int main(int argv, char* argc[])
{
  // Create a new instance of your nlp
  //  (use a SmartPtr, not raw)
  TrajectoryTracking* traj_tracker = new TrajectoryTracking();
  size_t n = 42;
  double curvature_weight = 1.0;
  double curvature_rate_weight = 1.0;
  double deviation_weight = 1.0;
  double delta_t = 1.0;
  double v_max = 1.0;

  Eigen::MatrixXd data(n, 3);
  data << 0.0714334, 0.00638967, 1.66001
    ,0.171036, 0.015299, 1.66001
    ,0.270638, 0.0242084, 1.66001
    ,0.37024, 0.0331177, 1.66001
    ,0.469843, 0.0420271, 1.66001
    ,0.569445, 0.0509365, 1.66001
    ,0.669047, 0.0598458, 1.66001
    ,0.76865, 0.0687552, 1.66001
    ,0.868252, 0.0776645, 1.66001
    ,0.967854, 0.0865739, 1.66001
    ,1.06746, 0.0954833, 1.66001
    ,1.16706, 0.104393, 1.66001
    ,1.26666, 0.113302, 1.66001
    ,1.36626, 0.122211, 1.66001
    ,1.46587, 0.131121, 1.66001
    ,1.56547, 0.14003, 1.66001
    ,1.66507, 0.148939, 1.66001
    ,1.76467, 0.157849, 1.66001
    ,1.86428, 0.166758, 1.66001
    ,1.96388, 0.175668, 1.66001
    ,2.06348, 0.184577, 1.66001
    ,2.16308, 0.193486, 1.66001
    ,2.26268, 0.202396, 1.66001
    ,2.36229, 0.211305, 1.66001
    ,2.46189, 0.220214, 1.66001
    ,2.56149, 0.229124, 1.66001
    ,2.66109, 0.238033, 1.66001
    ,2.7607, 0.246942, 1.66001
    ,2.8603, 0.255852, 1.66001
    ,2.9599, 0.264761, 1.66001
    ,3.0595, 0.27367, 1.66001
    ,3.15911, 0.28258, 1.66001
    ,3.25871, 0.291489, 1.66001
    ,3.35831, 0.300399, 1.66001
    ,3.45791, 0.309308, 1.66001
    ,3.55751, 0.318217, 1.66001
    ,3.65712, 0.327127, 1.66001
    ,3.75672, 0.336036, 1.66001
    ,3.85632, 0.344945, 1.66001
    ,3.95592, 0.353855, 1.66001
    ,4.05553, 0.362764, 1.66001
    ,4.15513, 0.371673, 1.66001;

  VectorXd x_init(n), y_init(n), d_init(n), theta_list(n);
  x_init = data.col(0);
  y_init = data.col(1);
  d_init.setZero();
  theta_list = data.col(2);

  traj_tracker->init(n, curvature_weight, curvature_rate_weight, deviation_weight, delta_t, v_max,
        x_init, y_init, d_init, theta_list);
        
  Ipopt::SmartPtr<Ipopt::TNLP> mynlp = traj_tracker;


  
  
  // Create a new instance of IpoptApplication
  //  (use a SmartPtr, not raw)
  // We are using the factory, since this allows us to compile this
  // example with an Ipopt Windows DLL
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
  // app->RethrowNonIpoptException(true);
  
  // Change some options
  // Note: The following choices are only examples, they might not be
  //       suitable for your optimization problem.
  app->Options()->SetNumericValue("tol", 1e-7);
  app->Options()->SetIntegerValue("print_level", 5);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt_output.txt"); // Output to file
  app->Options()->SetStringValue("check_derivatives_for_naninf", "yes"); // Check derivatives for NaN/Inf
  app->Options()->SetStringValue("print_user_options", "yes"); // Print user options (problem-specific)
  app->Options()->SetStringValue("nlp_file", "problem_formulation.nlp"); // Save NLP formulation

  // app->Options()->SetStringValue("warm_start_init_point", "true");
  // The following overwrites the default name (ipopt.opt) of the
  // options file
  // app->Options()->SetStringValue("option_file_name", "hs071.opt");
  
  // Initialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Ipopt::Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    return (int) status;
  }
  
  // // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp);
  
  if (status == Ipopt::Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
  }
  else {
    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
  }
  

  // return (int) status;
}
