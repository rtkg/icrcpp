#include "../../include/icr.h"
#include <iostream>
#include <sys/time.h>
#include <time.h>
#include <symbolic/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <vector>
#include <iomanip>

#define KAPPA 100000

using namespace ICR;
using namespace CasADi;

FX generateCodeAndCompile(FX fcn, const std::string& name){
  std::cout << "Generating code for " << name << std::endl;

  // // Convert to an SXFunction (may or may not improve efficiency)
  // if( is_a<MXFunction>(fcn)){
  //   fcn = SXFunction(shared_cast<MXFunction>(fcn));
  //   fcn.init();
  // }

  // Generate C code
  fcn.generateCode(name + ".c");

  // Compilation command
  string compile_command = "gcc -fPIC -shared -O3 " + name + ".c -o " + name + ".so";

  // Compile the c-code
  int flag = system(compile_command.c_str());
  casadi_assert_message(flag==0, "Compilation failed");

  // Load the generated function for evaluation
  ExternalFunction fcn_e("./" + name + ".so");
  return fcn_e;
}

//===========GLOBAL===========================
uint K=2;
uint T=4;

DiscreteWrenchSpace* tws_;
 void computeInitialSolution(CasADi::DMatrix& initial_solution, const std::vector<std::vector<Eigen::MatrixXd* > >& Wh_i_list)
{
  //wrenches in Ph_list and Wh_i_list are stored as column - matrices

  uint H=Wh_i_list.size();
  uint nV=0; //number of decision variables
  uint nC=0; //number of constraints
  uint K=Wh_i_list[0][0]->rows(); //dimension
  uint L=Wh_i_list[0][0]->cols(); //number of wrenches in a wrench cone

  uint T=tws_->getNumWrenches();

  for (uint h=0; h<H;h++)
    {
      nV+=Wh_i_list[h].size()*(L+1)+2;
      nC+=2*Wh_i_list[h].size()+T;
    }


  //Extract the task wrenches and put them in a CasADi::DMatrix
  CasADi::DMatrix TW(K,T);
  uint count=0;
  for (int j=0; j<T; j++)
    for (int i=0; i<K; i++)
      {
	TW(i,j)=tws_->getWrenches().get()[count];
	count++;
      }
  
  //decision variable vector
  CasADi::SXMatrix x = CasADi::ssym("x",nV); 
  CasADi::DMatrix x0=CasADi::DMatrix::zeros(nV,1);

  //Objective
  CasADi::SXMatrix f(1,1); 

  //Constraints
  CasADi::SXMatrix g=CasADi::SXMatrix::zeros(nC,1); 
  CasADi::DMatrix ub_g=CasADi::DMatrix::inf(nC,1);
  CasADi::DMatrix lb_g=CasADi::DMatrix::inf(nC,1)*(-1);

  //Bounds
  CasADi::DMatrix ub_x=CasADi::DMatrix::inf(nV,1);
  CasADi::DMatrix lb_x=CasADi::DMatrix::inf(nV,1)*(-1);

  //Index Sets for the decision variables
  double eh_ind;
  double deh_ind;
  std::vector<int> lmbdh_ind;
  std::vector<std::vector<int> > lmbdh_i_ind;
  std::vector<int> zh_i_ind;

  //Index Sets for the constraints
  std::vector<int> tc_ind(T);
  std::vector<int> eqc_ind;
  std::vector<int> sc_ind;

  uint h_start_x_ind=0;
  uint h_start_c_ind=0;

//HAAAAAACKKK REPLACE WITH READING QHULL NORMALS FROM THE G>WS
  std::vector<double*> hps;
  double n[2]={0.6402,-0.7682};
  hps.push_back(n);
  double n2[2]={-0.9950,-0.0995};
  hps.push_back(n2);
  double n3[2]={0.2169,0.9762};
  hps.push_back(n3);
  double n4[2]={0.9829, -0.1843};
  hps.push_back(n4);
  //  facetT* curr_f=grasp_->getGWS()->conv_hull_.beginFacet().getFacetT();

  for (uint h=0; h<H;h++)
    {
   //================FORM INDEX SETS========================================
     uint Ih=Wh_i_list[h].size();
      uint nVh=2+Ih*(L+1);

      //Build the h-th index lists for nh, eh, lbmd_h_i and zh_i for eache Hyperplane h
  lmbdh_ind.resize(Ih*L);
      lmbdh_i_ind.resize(Ih);
      zh_i_ind.resize(Ih);

   
      eh_ind=h_start_x_ind;
      deh_ind=h_start_x_ind+1;

      for (int i=0; i<Ih; i++)
	{
	  lmbdh_i_ind[i].resize(L);
	  for(int l=0; l<L; l++)
	    lmbdh_i_ind[i][l]=L*i+l+h_start_x_ind+2;

	  zh_i_ind[i]=h_start_x_ind+2+L*Ih+i;
	}

  for (uint i=0; i<Ih*L; i++)
	lmbdh_ind[i]=h_start_x_ind+2+i;

      //Build the h-th index lists for the constriants for each Hyperplane h
      uint nCh=2*Ih+T;
      eqc_ind.resize(Ih);
      sc_ind.resize(Ih);

      for (uint t=0; t<T;t++)
	tc_ind[t]=h_start_c_ind+t;

      for (uint i=0; i<Ih;i++)
	{
	  eqc_ind[i]=h_start_c_ind+T+i;
	  sc_ind[i]=h_start_c_ind+T+Ih+i;
	}

      h_start_x_ind+=nVh;
   h_start_c_ind+=nCh;

      //================OBJECTIVE========================================
   f=f-x(deh_ind)-KAPPA*mul(x(zh_i_ind).trans(),CasADi::SXMatrix::ones(Ih,1));

      //==================== BOUNDS ========================================
      //de's have to be positive
   lb_x(deh_ind)=CasADi::DMatrix::zeros(1,1);

      //Slack variables have to be negative
   ub_x(zh_i_ind)=CasADi::DMatrix::zeros(Ih,1);

      //Lambdas have to be  positive
      for (int i=0; i<Ih;i++)
      	lb_x(lmbdh_i_ind[i])=CasADi::DMatrix::zeros(L,1);

  //====================HARD CONSTRAINTS ========================================

 //Constraints involving the task wrenches
      //get the h-th normal of the prototype grasp's GWS
      CasADi::DMatrix ah=CasADi::DMatrix::zeros(K,1);
      for (int k=0; k<K;k++)
	{
	  //HAAAAAACKKK REPLACE WITH READING QHULL NORMALS FROM THE G>WS
	  //        ah(k)=curr_f->normal[k]*(-1);
	  ah(k)=hps[h][k];

	}

      g(tc_ind) = mul(x(eh_ind),-CasADi::SXMatrix::ones(T,1));
      ub_g(tc_ind) = mul(TW.trans(),ah);

     //Lambdas have to sum up to one
      //Build the Equality matrix
      std::vector<int> col_ind(L);
      CasADi::SXMatrix E=CasADi::SXMatrix::zeros(Ih,Ih*L);
      for (uint i=0; i<Ih;i++)
     	{
     	  for (uint j=0; j<L;j++)
     	    col_ind[j]=i*L+j;
	     
     	  E(i,col_ind)=CasADi::SXMatrix::ones(1,L);
     	}
      //set the Equality constraints
      g(eqc_ind)=mul(E,x(lmbdh_ind));

      ub_g(eqc_ind)=CasADi::DMatrix::ones(Ih,1);
      lb_g(eqc_ind)=CasADi::DMatrix::ones(Ih,1);

      //==================== SOFT CONSTRAINTS ======================================
      for (uint i=0; i<Ih; i++)
	{
	  //Extract the i-th Wrench cone and store it in a CasADi::DMatrix
	  CasADi::SXMatrix Whi(Wh_i_list[h][i]->rows(),Wh_i_list[h][i]->cols());
	  for (uint r=0; r<Wh_i_list[h][i]->rows(); r++)
	    for (uint c=0;c<Wh_i_list[h][i]->cols() ; c++)
	      Whi(r,c)=(*Wh_i_list[h][i])(r,c);


  g(sc_ind[i])=mul(CasADi::SXMatrix(ah.trans()), mul(Whi,x(lmbdh_i_ind[i]))) +x(eh_ind)+x(deh_ind)+x(zh_i_ind[i]);
  ub_g(sc_ind[i])=CasADi::DMatrix::zeros(1,1);
	}

      //      curr_f=curr_f->next;
    }		  

  //  std::cout<<"f: "<<std::endl<<f<<std::endl;
  //  std::cout<<"g: "<<std::endl<<g<<std::endl;
  //  std::cout<<"ub_g: "<<std::endl<<ub_g<<std::endl;
  //  std::cout<<"lb_g: "<<std::endl<<lb_g<<std::endl;
  //  std::cout<<"ub_x: "<<std::endl<<ub_x<<std::endl;
  //  std::cout<<"lb_x: "<<std::endl<<lb_x<<std::endl;

  CasADi::SXFunction lp(nlpIn("x",x),nlpOut("f",f,"g",g));
  CasADi::IpoptSolver lp_solver(lp);

  lp_solver.setOption("print_level",0);
  lp_solver.setOption("linear_solver","ma97");
  lp_solver.setOption("mu_strategy","adaptive");
  lp_solver.init();

  lp_solver.setInput( x0, "x0");
  lp_solver.setInput(lb_x, "lbx");
  lp_solver.setInput(ub_x, "ubx");
  lp_solver.setInput(lb_g, "lbg");
  lp_solver.setInput(ub_g, "ubg");

  lp_solver.evaluate();

  std::cout << std::setw(30) << "Initial Objective: " << lp_solver.output("f").getDescription() << std::endl;
  std::cout << std::setw(30) << "Initial Primal solution: " << lp_solver.output("x").getDescription() << std::endl;


  //Extract the initial vectors from the solution
  initial_solution=CasADi::DMatrix();
  //curr_f=grasp_->getGWS()->conv_hull_.beginFacet().getFacetT();
      h_start_x_ind=0;
  CasADi::DMatrix solution=lp_solver.output("x");
  for (uint h=0; h<H; h++)
    {
      uint Ih=Wh_i_list[h].size();
      uint nVh=2+Ih*(L+1);
    
      eh_ind=h_start_x_ind;
      deh_ind=h_start_x_ind+1;
      std::vector<int> lmbdh_zh_ind(Ih*(L+1));
      for (int j=0; j<lmbdh_zh_ind.size(); j++)
	lmbdh_zh_ind[j]=h_start_x_ind+2+j;
  
      //get the h-th normal
      CasADi::DMatrix ah=CasADi::DMatrix::zeros(K,1);
      for (int k=0; k<K;k++)
	{
	  //HAAAAAACKKK REPLACE WITH READING QHULL NORMALS FROM THE G>WS
	  //        ah(k)=curr_f->normal[k]*(-1);
	  ah(k)=hps[h][k];
	}


   double deh=solution(deh_ind).toScalar();
    if(deh > 1e-5)
      {
      initial_solution=CasADi::vertcat(initial_solution, ah/deh); 
      initial_solution=CasADi::vertcat(initial_solution, solution(eh_ind)/deh); 
      }
    else
      {
	initial_solution=CasADi::vertcat(initial_solution, ah*(1e5)); 
      initial_solution=CasADi::vertcat(initial_solution, solution(eh_ind)*1e5); 
     }
      initial_solution=CasADi::vertcat(initial_solution, solution(lmbdh_zh_ind));



      h_start_x_ind+=nVh;
      // curr_f=curr_f->next;
    }

  std::cout<<"initial solution:"<<std::endl<<initial_solution<<std::endl;
 
}
  //-------------------------------------------------------------------

void createNLPSolver(const std::vector<Eigen::MatrixXd*>& Ph_list, const std::vector<std::vector<Eigen::MatrixXd* > >& Wh_i_list,IpoptSolver& nlp_solver)
{
  //wrenches in Ph_list and Wh_i_list are stored as column - matrices

  uint H=Ph_list.size();
  uint nV=0; //number of decision variables
  uint nC=0; //number of constraints
  uint K=Wh_i_list[0][0]->rows(); //dimension
  uint L=Wh_i_list[0][0]->cols(); //number of wrenches in a wrench cone
  uint T=tws_->getNumWrenches();

  for (uint h=0; h<H;h++)
    {
      nV+=Wh_i_list[h].size()*(L+1)+K+1;
      nC+=2*Wh_i_list[h].size()+T+Ph_list[h]->cols();
    }


  //Extract the task wrenches and put them in a SXMatrix
  SXMatrix TW(K,T);
  uint count=0;
  for (int j=0; j<T; j++)
    for (int i=0; i<K; i++)
      {
	TW(i,j)=tws_->getWrenches().get()[count];
	count++;
      }

  //decision variable vector & initial solution
  SXMatrix x = ssym("x",nV); 
  DMatrix x0(nV,1,0);
  computeInitialSolution(x0,Wh_i_list);

  //Objective
  SXMatrix f(1,1); 

  //Constraints
  SXMatrix g(nC,1); 
  DMatrix ub_g=DMatrix::inf(nC,1);
  DMatrix lb_g=DMatrix::inf(nC,1)*(-1);

  //Bounds
  DMatrix ub_x=DMatrix::inf(nV,1);
  DMatrix lb_x=DMatrix::inf(nV,1)*(-1);

  //Index Sets for the decision variables
  std::vector<int> nh_ind(K);
  double eh_ind;
  std::vector<int> lmbdh_ind;
  std::vector<std::vector<int> > lmbdh_i_ind;
  std::vector<int> zh_i_ind;

  //Index Sets for the constraints
  std::vector<int> phc_ind;
  std::vector<int> tc_ind(T);
  std::vector<int> eqc_ind;
  std::vector<int> sc_ind;

  uint h_start_x_ind=0;
  uint h_start_c_ind=0;
  for (uint h=0; h<H;h++)
    {
      //===============FORM INDEX SETS========================================
      uint Ih=Wh_i_list[h].size();
      uint n_Ph=Ph_list[h]->cols();

      //Build the h-th index lists for nh, eh, lbmd_h_i and zh_i for each Hyperplane h
      uint nVh=K+1+Ih*(L+1);
 lmbdh_ind.resize(Ih*L);
      lmbdh_i_ind.resize(Ih);
      zh_i_ind.resize(Ih);

      for (uint j=0; j<K; j++)
	nh_ind[j]=h_start_x_ind+j;

      eh_ind=h_start_x_ind+K;

      for (uint i=0; i<Ih; i++)
	{
	  lmbdh_i_ind[i].resize(L);
	  for(uint l=0; l<L; l++)
	    lmbdh_i_ind[i][l]=L*i+l+h_start_x_ind+K+1;

	  zh_i_ind[i]=h_start_x_ind+K+1+L*Ih+i;
	}

  for (uint i=0; i<Ih*L; i++)
	lmbdh_ind[i]=h_start_x_ind+K+1+i;

      //Build the h-th index lists for the constriants for each Hyperplane h
      uint nCh=2*Ih+T+n_Ph;  
      phc_ind.resize(n_Ph);
      eqc_ind.resize(Ih);
      sc_ind.resize(Ih);

      for (uint p=0; p<n_Ph;p++)
	phc_ind[p]=h_start_c_ind+p;

      for (uint t=0; t<T;t++)
	tc_ind[t]=h_start_c_ind+n_Ph+t;

      for (uint i=0; i<Ih;i++)
	{
	  eqc_ind[i]=h_start_c_ind+n_Ph+T+i;
	  sc_ind[i]=h_start_c_ind+n_Ph+T+Ih+i;
	}

      h_start_x_ind+=nVh;
      h_start_c_ind+=nCh;

      //================OBJECTIVE========================================
      f=f+0.5*mul(x(nh_ind).trans(),x(nh_ind))-KAPPA*mul(x(zh_i_ind).trans(),SXMatrix(Ih,1,1));
 
      //================HARD CONSTRAINTS========================================
      //Extract the wrenches of the prototype grasp and store them in a SXMatrix
      SXMatrix Ph(Ph_list[h]->rows(),Ph_list[h]->cols());
      for (uint i=0; i<Ph_list[h]->rows(); i++)
	for (uint j=0; j<Ph_list[h]->cols(); j++)
	  Ph(i,j)=(*Ph_list[h])(i,j);

      //Constraints involving the prototype grasps's wrenches Ph_i
      g(phc_ind)=mul(Ph.trans(),x(nh_ind))+mul(x(eh_ind),SXMatrix(K,1,1));
      ub_g(phc_ind)=DMatrix(n_Ph,1,-1);


      // for (uint j=0; j<phc_ind.size();j++)
      // 	std::cout<<phc_ind[j]<<" ";

      //Constraints involving the task wrenches
      g(tc_ind) = -mul(TW.trans(),x(nh_ind))-mul(x(eh_ind),SXMatrix(T,1,1));
      ub_g(tc_ind)=DMatrix(T,1,0); //0 for close fit, -1 for middle fit
             

      //Lambdas have to sum up to one
      //Build the Equality matrix
      std::vector<int> col_ind(L);
      CasADi::SXMatrix E(Ih,Ih*L,0);
      for (uint i=0; i<Ih;i++)
	{
	  for (uint j=0; j<L;j++)
	    col_ind[j]=i*L+j;
	     
	  E(i,col_ind)=CasADi::SXMatrix(1,L,1);
	}
      //set the Equality constraints
      g(eqc_ind)=mul(E,x(lmbdh_ind));
      ub_g(eqc_ind)=DMatrix(Ih,1,1);
      lb_g(eqc_ind)=DMatrix(Ih,1,1);

      //=================SOFT CONSTRAINTS========================================
      for (uint i=0; i<Ih; i++)
	{
	  //Extract the i-th Wrench cone and store it in a SXMatrix
	  SXMatrix Whi(Wh_i_list[h][i]->rows(),Wh_i_list[h][i]->cols());
	  for (uint r=0; r<Wh_i_list[h][i]->rows(); r++)
	    for (uint c=0;c<Wh_i_list[h][i]->cols() ; c++)
	      Whi(r,c)=(*Wh_i_list[h][i])(r,c);

	  g(sc_ind[i])=mul(x(nh_ind).trans(),mul(Whi,x(lmbdh_i_ind[i])))+x(eh_ind)+x(zh_i_ind[i]);
	  ub_g(sc_ind[i])=DMatrix(1,1,-1);

	}

      //==================== BOUNDS ========================================
      //Slack variables have to be negative
      ub_x(zh_i_ind)=DMatrix(Ih,1,0);

      //Lambdas have to be  positive
      for (uint i=0; i<Ih;i++)
	lb_x(lmbdh_i_ind[i])=DMatrix(L,1,0);

    }


  //  std::cout<<"f: "<<std::endl<<f<<std::endl;
  //  std::cout<<"g: "<<std::endl<<g<<std::endl;
  // std::cout<<"ub_g: "<<std::endl<<ub_g<<std::endl;
  // std::cout<<"lb_g: "<<std::endl<<lb_g<<std::endl;
  // std::cout<<"ub_x: "<<std::endl<<ub_x<<std::endl;
  // std::cout<<"lb_x: "<<std::endl<<lb_x<<std::endl;

  // Convert MXFunction to SXFunction before code generation (may or may not improve efficiency)
  bool expand = true;


 // NLP function
  FX nlp = SXFunction(nlpIn("x",x),nlpOut("f",f,"g",g));
  nlp.init();




//  // NLP function
// SXFunction nlp(nlpIn("x",x),nlpOut("f",f,"g",g));
//   nlp.init();

  // Gradient of the objective
  FX grad_f = nlp.gradient("x","f");
  grad_f.init();

  // Jacobian of the constraints
  FX jac_g = nlp.jacobian("x","g");
  jac_g.init();

  // Hessian of the lagrangian
  FX grad_lag = nlp.derivative(0,1);
  FX hess_lag = grad_lag.jacobian(NL_X,NL_NUM_OUT+NL_X,false,true);
  hess_lag.init();

  // Codegen and compile
  nlp = generateCodeAndCompile(nlp,"nlp");
  grad_f = generateCodeAndCompile(grad_f,"grad_f");
  jac_g = generateCodeAndCompile(jac_g,"jac_g");
  hess_lag = generateCodeAndCompile(hess_lag,"hess_lag");

  // Create an NLP solver passing derivative information
nlp_solver=IpoptSolver(nlp);
  nlp_solver.setOption("grad_f",grad_f);
  nlp_solver.setOption("jac_g",jac_g);
  nlp_solver.setOption("hess_lag",hess_lag);
  nlp_solver.init();

  nlp_solver.setInput( x0, "x0");
  nlp_solver.setInput(lb_x, "lbx");
  nlp_solver.setInput(ub_x, "ubx");
  nlp_solver.setInput(lb_g, "lbg");
  nlp_solver.setInput(ub_g, "ubg");
  // nlp_solver.setOption("print_level",0);



  nlp_solver.evaluate();

  //OOOOLDDDD

  // SXFunction nlp(nlpIn("x",x),nlpOut("f",f,"g",g));
  // nlp_solver=IpoptSolver(nlp);
  //  nlp_solver.setOption("print_level",0);
  // nlp_solver.init();

  // nlp_solver.setInput( x0, "x0");
  // nlp_solver.setInput(lb_x, "lbx");
  // nlp_solver.setInput(ub_x, "ubx");
  // nlp_solver.setInput(lb_g, "lbg");
  // nlp_solver.setInput(ub_g, "ubg");

  // nlp_solver.evaluate();

  //  std::cout << std::setw(30) << "Objective: " << nlp_solver.output("f").getDescription() << std::endl;
  // std::cout << std::setw(30) << "Primal solution: " << nlp_solver.output("x").getDescription() << std::endl; 
}
//-------------------------------------------------------------------------------------------------
void computeWarpedHyperplanes()
{
  std::vector<Eigen::MatrixXd*> Ph_list;
  Eigen::MatrixXd mat(2,2);
  mat(0,0)=0.7; mat(1,0)=2; mat(0,1)=-0.5; mat(1,1)=1;
  Ph_list.push_back(new Eigen::MatrixXd(mat));
  // mat(0,0)=0.7; mat(1,0)=2; mat(0,1)=1; mat(1,1)=-1;
  // Ph_list.push_back(new Eigen::MatrixXd(mat));
  // mat(0,0)=1; mat(1,0)=-1; mat(0,1)=-0.8; mat(1,1)=-0.6; 
  // Ph_list.push_back(new Eigen::MatrixXd(mat));
  // mat(0,0)=-0.8; mat(1,0)=-0.6; mat(0,1)=-0.5; mat(1,1)=1;
  // Ph_list.push_back(new Eigen::MatrixXd(mat));


  std::vector<std::vector<Eigen::MatrixXd* > > Wh_i_list(Ph_list.size());
  mat(0,0)=-0.6; mat(1,0)=0.6; mat(0,1)=-0.2; mat(1,1)=1;
  Wh_i_list[0].push_back(new Eigen::MatrixXd(mat));
  // mat(0,0)=0.2; mat(1,0)=1.8; mat(0,1)=0.8; mat(1,1)=1.5;
  // Wh_i_list[0].push_back(new Eigen::MatrixXd(mat));
  // mat(0,0)=0; mat(1,0)=1.1; mat(0,1)=1; mat(1,1)=1;
  // Wh_i_list[0].push_back(new Eigen::MatrixXd(mat));

  // mat(0,0)=0.2; mat(1,0)=1.8; mat(0,1)=0.8; mat(1,1)=1.5;
  // Wh_i_list[1].push_back(new Eigen::MatrixXd(mat));
  // mat(0,0)=0.5; mat(1,0)=-0.55; mat(0,1)=1.2; mat(1,1)=0.1;
  // Wh_i_list[1].push_back(new Eigen::MatrixXd(mat));

  // mat(0,0)=0.5; mat(1,0)=-0.55; mat(0,1)=1.2; mat(1,1)=0.1;
  // Wh_i_list[2].push_back(new Eigen::MatrixXd(mat));
  // mat(0,0)=-1; mat(1,0)=-0.6; mat(0,1)=-0.4; mat(1,1)=-0.6;
  // Wh_i_list[2].push_back(new Eigen::MatrixXd(mat));

  // mat(0,0)=-1; mat(1,0)=-0.6; mat(0,1)=-0.4; mat(1,1)=-0.6;
  // Wh_i_list[3].push_back(new Eigen::MatrixXd(mat));
  // mat(0,0)=-0.6; mat(1,0)=0.6; mat(0,1)=-0.2; mat(1,1)=1;
  // Wh_i_list[3].push_back(new Eigen::MatrixXd(mat));

  IpoptSolver nlp_solver;
  createNLPSolver(Ph_list,Wh_i_list,nlp_solver);

  for (int i=0; i<Ph_list.size(); i++)
    delete Ph_list[i];

  for (int i=0; i<Wh_i_list.size(); i++)
    for (int j=0; j<Wh_i_list[i].size(); j++)
      delete Wh_i_list[i][j];
}
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
int main()
{ 

  // ==============TWS===============================================
  Eigen::MatrixXd tmp(T,K);
  tmp<<0.5, 0.4, -0.1, 0.6, -0.4, -0.4, 0.1, -0.3;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tw_e=tmp.transpose();
  double* a;
  eigenMatrixToDoubleArray(tw_e,a);
  tws_=new DiscreteWrenchSpace(K,SharedDoublePtr(a),T);

  computeWarpedHyperplanes();

 
  // orgQhull::Qhull conv_hull_;
  //  conv_hull_.runQhull("", dim,num_t_w,a ,"Qx Qt");

  //  conv_hull_.defineVertexNeighborFacets();
  //  double area_=conv_hull_.area();
  //  double volume_=conv_hull_.volume();
  //  int  num_vtx_=conv_hull_.vertexCount();
  //  int  num_facets_=conv_hull_.facetCount();

  //  facetT* curr_f=conv_hull_.beginFacet().getFacetT();
  //  double r_oc_insphere_=-(curr_f->offset);


  //  for(int i=0;i< num_facets_;i++)
  //    {
  //      curr_f->id=i; //Replaces the Qhull runtime indexing with indices 0 - num_facets_
  //      r_oc_insphere_ = (-(curr_f->offset) < r_oc_insphere_) ? -(curr_f->offset) : r_oc_insphere_;
  //      std::cout<<"normal: "<<-Eigen::Map<Eigen::Matrix<double,1,2> >(curr_f->normal)<<" offset: "<<-(curr_f->offset)  <<std::endl;
    
   
  //      curr_f=curr_f->next;
  //    }
  
  //  std::cout<<"num_facets: "<<num_facets_<<std::endl;
  //  std::cout<<"r_oc_insphere: "<<r_oc_insphere_<<std::endl;


  // tw_a.reset(a);

  // DiscreteWrenchSpace tws(dim,tw_a,num_t_w);
  // tws.computeConvexHull(); 

  // std::cout<<tws<<std::endl;

  delete tws_;

  return 0;
}
