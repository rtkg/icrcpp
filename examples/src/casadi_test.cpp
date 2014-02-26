#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <symbolic/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <symbolic/stl_vector_tools.hpp>

using namespace CasADi;
using namespace std;
/**
 *  Example program demonstrating parametric NLPs in CasADi
 *  Note that there is currently no support for parametric sensitivities via this feature (although it would make a lot of sense).
 *  For parametric sensitivities, see the parametric_sensitivities.cpp example which calculates sensitivitied via the sIPOPT extension
 *  to IPOPT.
 * 
 *  Joel Andersson, K.U. Leuven 2012
 */

int main(){
    
  /** Test problem (Ganesh & Biegler, A reduced Hessian strategy for sensitivity analysis of optimal flowsheets, AIChE 33, 1987, pp. 282-296)
   * 
   *    min     x1^2 + x2^2 + x3^2
   *    s.t.    6*x1 + 3&x2 + 2*x3 - pi = 0
   *            p2*x1 + x2 - x3 - 1 = 0
   *            x1, x2, x3 >= 0
   * 
   */
  
  // Optimization variables
  SXMatrix x = ssym("x",3);
SXMatrix y = ssym("y",3);
DMatrix D(2,2);
 D(0,0)=0; D(0,1)=1; D(1,0)=10; D(1,1)=11;

 std::vector<int> ind; ind.push_back(0); ind.push_back(1);
 D(ind,1)=77;

 std::cout<<"D: "<<std::endl<<D<<std::endl;

// SXMatrix fun=x.mul(y);
 return 0;

//std::cout<<"FUN: "<<std::endl<<fun<<std::endl;

  // Parameters
  SXMatrix p = ssym("p",2);
  
  // Objective
  SXMatrix f = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
  
  // Constraints
  SXMatrix g = vertcat(
       6*x[0] + 3*x[1] + 2*x[2] - p[0],
    p[1]*x[0] +   x[1] -   x[2] -    1
  );
  
  std::cout<<"CONSTRAINTS: ..."<<std::endl<<g<<std::endl;

  // Infinity
  double inf = numeric_limits<double>::infinity();
  
  // Initial guess and bounds for the optimization variables
  double x0_[]  = {0.15, 0.15, 0.00};
  double lbx_[] = {0.00, 0.00, 0.00};
  double ubx_[] = { inf,  inf,  inf};
  vector<double> x0(x0_,x0_+3);
  vector<double> lbx(lbx_,lbx_+3);
  vector<double> ubx(ubx_,ubx_+3);
  
  // Nonlinear bounds
  double lbg_[] = {0.00, 0.00};
  double ubg_[] = {0.00, 0.00};
  vector<double> lbg(lbg_,lbg_+2);
  vector<double> ubg(ubg_,ubg_+2);
  
  // Original parameter values
  double p0_[]  = {5.00,1.00};
  vector<double> p0(p0_,p0_+2);

  // NLP
  SXFunction nlp(nlpIn("x",x,"p",p),nlpOut("f",f,"g",g));

  std::cout<<"NLP: "<<std::endl<<nlp<<std::endl;

  // Create NLP solver
  IpoptSolver solver(nlp);
  Dictionary nlp_solver_options;
  nlp_solver_options["print_level"] = 0;

  // nlp_solver.setOption("stabilized_qp_solver",QPStabilizer::creator);
  // Dictionary stabilized_qp_solver_options;
  // stabilized_qp_solver_options["qp_solver"] = NLPQPSolver::creator;
  // Dictionary qp_solver_options;
  // qp_solver_options["nlp_solver"]= IpoptSolver::creator;
  // Dictionary nlp_solver_options;
  // nlp_solver_options["print_level"] = 0;
  // nlp_solver_options["print_time"] = 0;
  // nlp_solver_options["tol"] = 1e-16;
  // nlp_solver_options["constr_viol_tol"] = 1e-16;
  // nlp_solver_options["dual_inf_tol"] = 1e-16;
  // nlp_solver_options["compl_inf_tol"] = 1e-16;
  // qp_solver_options["nlp_solver_options"] = nlp_solver_options;
  // stabilized_qp_solver_options["qp_solver_options"] = qp_solver_options;

   nlp_solver.setOption("nlp_solver_options",nlp_solver_options);


  solver.init();
  
  // Solve NLP
  solver.setInput( x0, "x0");
  solver.setInput( p0, "p");
  solver.setInput(lbx, "lbx");
  solver.setInput(ubx, "ubx");
  solver.setInput(lbg, "lbg");
  solver.setInput(ubg, "ubg");

  solver.evaluate();
  
  // // Print the solution
  // cout << "-----" << endl;
  // cout << "Optimal solution for p = " << solver.input("p").getDescription() << ":" << endl;
  // cout << setw(30) << "Objective: " << solver.output("f").getDescription() << endl;
  // cout << setw(30) << "Primal solution: " << solver.output("x").getDescription() << endl;
  // cout << setw(30) << "Dual solution (x): " << solver.output("lam_x").getDescription() << endl;
  // cout << setw(30) << "Dual solution (g): " << solver.output("lam_g").getDescription() << endl;
  
  // // Change the parameter and resolve
  // p0[0] = 4.5;
  // solver.setInput( p0, "p");
  // solver.evaluate();
  
  // Print the new solution
  cout << "-----" << endl;
  cout << "Optimal solution for p = " << solver.input("p").getDescription() << ":" << endl;
  cout << setw(30) << "Objective: " << solver.output("f").getDescription() << endl;
  cout << setw(30) << "Primal solution: " << solver.output("x").getDescription() << endl;
  cout << setw(30) << "Dual solution (x): " << solver.output("lam_x").getDescription() << endl;
  cout << setw(30) << "Dual solution (g): " << solver.output("lam_g").getDescription() << endl;
  
  return 0;
}
