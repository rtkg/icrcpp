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

//===========GLOBAL===========================
uint K=2;
uint T=4;

DiscreteWrenchSpace* tws_;

//-------------------------------------------------------------------------------------------------
void computeInitialSolution(DMatrix& solution, const std::vector<std::vector<Eigen::MatrixXd* > >& Wh_i_list)
{
  //wrenches in Ph_list and Wh_i_list are stored as column - matrices

  uint H=Wh_i_list.size();
  uint nV=0; //number of decision variables
  uint L=2;
  for (uint i=0; i<H;i++)
    nV+=Wh_i_list[i].size()*(L+1)+2;

  //Extract the task wrenches and put them in a DMatrix
  DMatrix TW(K,tws_->getNumWrenches());
  uint count=0;
  for (int j=0; j<tws_->getNumWrenches(); j++)
    for (int i=0; i<K; i++)
      {
	TW(i,j)=tws_->getWrenches().get()[count];
	count++;
      }

  SXMatrix x = ssym("x",nV); //decision variable vector
  SXMatrix f(1,1); //Objective
  DMatrix x0(nV,1,0);
  double eh_ind;
  double deh_ind;
  std::vector<std::vector<int> > lmbdh_i_ind;
  std::vector<int> zh_i_ind;
  SXMatrix g,  g_T,g_lmbde,g_s;
  DMatrix  lb_g, ub_g, lb_g_T, ub_g_T, lb_g_lmbde, ub_g_lmbde,lb_g_s,ub_g_s;
  DMatrix ub_x(nV,1,std::numeric_limits<double>::infinity());
  DMatrix lb_x(nV,1,std::numeric_limits<double>::infinity()*(-1));
  for (uint h=0; h<H;h++)
    {
 
      //Build the h-th index lists for nh, eh, lbmd_h_i and zh_i for eache Hyperplane h
      uint Ih=Wh_i_list[h].size();
      uint nVh=2+Ih*(L+1);
      lmbdh_i_ind.resize(Ih);
      zh_i_ind.resize(Ih);

   
      eh_ind=h*nVh;
      deh_ind=h*nVh+1;

      for (int i=0; i<Ih; i++)
	{
	  lmbdh_i_ind[i].resize(L);
	  for(int l=0; l<L; l++)
	    lmbdh_i_ind[i][l]=L*i+l+h*nVh+2;

	  zh_i_ind[i]=h*nVh+2+L*Ih+i;
	}

       //================OBJECTIVE========================================
      f=f+x(eh_ind)-KAPPA*x(zh_i_ind).trans().mul(DMatrix(Ih,1,1));

      //==================== BOUNDS ========================================
      //de's have to be positive
      lb_x(deh_ind)=DMatrix(1,1,0);

      //Slack variables have to be negative
      ub_x(zh_i_ind)=DMatrix(Ih,1,0);

      //Lambdas have to be  positive
      for (int i=0; i<Ih;i++)
      	lb_x(lmbdh_i_ind[i])=DMatrix(L,1,0);

           //==================== EQUALITY CONSTRAINTS ========================================
          //Lambda_i have to sum up to one
      for (int i=0; i<Ih;i++)
	{
	  g_lmbde=vertcat(g_lmbde,(x(lmbdh_i_ind[i]).trans()).mul(DMatrix(L,1,1)));
	  ub_g_lmbde=vertcat(ub_g_lmbde,DMatrix(1,1,1));
	  lb_g_lmbde=vertcat(lb_g_lmbde,DMatrix(1,1,1));
	}
  
           //==================== INEQUALITY CONSTRAINTS ======================================
   //Constraints involving the task wrenches


      // g_T = vertcat(g_T, (-1)*(x(nh_ind).trans().mul(TW)).trans()-x(eh_ind).mul(DMatrix(T,1,1)));
      // ub_g_T = vertcat(ub_g_T,DMatrix(T,1,0));
      // lb_g_T =vertcat(lb_g_T, DMatrix(T,1,std::numeric_limits<double>::infinity()*(-1)));

      //Soft constraints involving the slack variables

    }

  g=vertcat(g_T,g_lmbde); g=vertcat(g,g_s);
  ub_g=vertcat(ub_g_T,ub_g_lmbde); ub_g=vertcat(ub_g,ub_g_s);
  lb_g=vertcat(lb_g_T,lb_g_lmbde); lb_g=vertcat(lb_g,lb_g_s);


  // std::cout<<"f: "<<std::endl<<f<<std::endl;
  std::cout<<"g: "<<std::endl<<g<<std::endl;
  std::cout<<"ub_g: "<<std::endl<<ub_g<<std::endl;
  std::cout<<"lb_g: "<<std::endl<<lb_g<<std::endl;
   // std::cout<<"ub_x: "<<std::endl<<ub_x<<std::endl;
   // std::cout<<"lb_x: "<<std::endl<<lb_x<<std::endl;

  // SXFunction nlp(nlpIn("x",x),nlpOut("f",f,"g",g));
  // nlp_solver=IpoptSolver(nlp);
  // nlp_solver.setOption("print_level",0);


  // nlp_solver.init();

  // nlp_solver.setInput( x0, "x0");
  // nlp_solver.setInput(lb_x, "lbx");
  // nlp_solver.setInput(ub_x, "ubx");
  // nlp_solver.setInput(lb_g, "lbg");
  // nlp_solver.setInput(ub_g, "ubg");

  // nlp_solver.evaluate();

  // std::cout << std::setw(30) << "Objective: " << nlp_solver.output("f").getDescription() << std::endl;
  // std::cout << std::setw(30) << "Primal solution: " << nlp_solver.output("x").getDescription() << std::endl;

}
//-------------------------------------------------------------------------------------------------
void createNLPSolver(const std::vector<Eigen::MatrixXd*>& Ph_list, const std::vector<std::vector<Eigen::MatrixXd* > >& Wh_i_list,IpoptSolver& nlp_solver)
{
  //wrenches in Ph_list and Wh_i_list are stored as column - matrices

  uint H=Ph_list.size();
  uint nV=0; //number of decision variables
  uint L=2;
  for (uint i=0; i<H;i++)
    nV+=Wh_i_list[i].size()*(L+1)+K+1;

  //Extract the task wrenches and put them in a DMatrix
  DMatrix TW(K,tws_->getNumWrenches());
  uint count=0;
  for (int j=0; j<tws_->getNumWrenches(); j++)
    for (int i=0; i<K; i++)
      {
	TW(i,j)=tws_->getWrenches().get()[count];
	count++;
      }

  SXMatrix x = ssym("x",nV); //decision variable vector
  SXMatrix f(1,1); //Objective
  DMatrix x0;
  computeInitialSolution(x0,Wh_i_list);
  std::vector<int> nh_ind(K);
  double eh_ind;
  std::vector<std::vector<int> > lmbdh_i_ind;
  std::vector<int> zh_i_ind;
  SXMatrix g, g_P, g_T,g_lmbde,g_s;
  DMatrix  lb_g, ub_g, lb_g_P, lb_g_T, ub_g_P, ub_g_T, lb_g_lmbde, ub_g_lmbde,lb_g_s,ub_g_s;
  DMatrix ub_x(nV,1,std::numeric_limits<double>::infinity());
  DMatrix lb_x(nV,1,std::numeric_limits<double>::infinity()*(-1));
  for (uint h=0; h<H;h++)
    {
 
      //Build the h-th index lists for nh, eh, lbmd_h_i and zh_i for eache Hyperplane h
      uint Ih=Wh_i_list[h].size();
      uint nVh=K+1+Ih*(L+1);
      lmbdh_i_ind.resize(Ih);
      zh_i_ind.resize(Ih);

      for (int j=0; j<K; j++)
        nh_ind[j]=h*nVh+j;

      eh_ind=h*nVh+K;

      for (int i=0; i<Ih; i++)
	{
	  lmbdh_i_ind[i].resize(L);
	  for(int l=0; l<L; l++)
	    lmbdh_i_ind[i][l]=L*i+l+h*nVh+K+1;

	  zh_i_ind[i]=h*nVh+K+1+L*Ih+i;
	}
      //================OBJECTIVE========================================
      f=f+0.5*x(nh_ind).trans().mul(x(nh_ind))-KAPPA*x(zh_i_ind).trans().mul(DMatrix(Ih,1,1));

      //================HARD CONSTRAINTS========================================
      //Extract the wrenches of the prototype grasp and store them in a DMatrix
      DMatrix Ph(Ph_list[h]->rows(),Ph_list[h]->cols());
      for (int i=0; i<Ph_list[h]->rows(); i++)
	for (int j=0; j<Ph_list[h]->cols(); j++)
	  Ph(i,j)=(*Ph_list[h])(i,j);

      //Constraints involving the prototype grasps's wrenches
      g_P = vertcat(g_P, (x(nh_ind).trans().mul(Ph)).trans()+x(eh_ind).mul(DMatrix(K,1,1)));
      ub_g_P = vertcat(ub_g_P,DMatrix(Ph_list[h]->cols(),1,-1));
      lb_g_P =vertcat(lb_g_P, DMatrix(Ph_list[h]->cols(),1,std::numeric_limits<double>::infinity()*(-1)));

      //Constraints involving the task wrenches
      g_T = vertcat(g_T, (-1)*(x(nh_ind).trans().mul(TW)).trans()-x(eh_ind).mul(DMatrix(T,1,1)));
      ub_g_T = vertcat(ub_g_T,DMatrix(T,1,0));
      lb_g_T =vertcat(lb_g_T, DMatrix(T,1,std::numeric_limits<double>::infinity()*(-1)));


      //Lambda_i have to sum up to one
      for (int i=0; i<Ih;i++)
	{
	  g_lmbde=vertcat(g_lmbde,(x(lmbdh_i_ind[i]).trans()).mul(DMatrix(L,1,1)));
	  ub_g_lmbde=vertcat(ub_g_lmbde,DMatrix(1,1,1));
	  lb_g_lmbde=vertcat(lb_g_lmbde,DMatrix(1,1,1));
	}

      //==================== BOUNDS ========================================
      //Slack variables have to be negative
      ub_x(zh_i_ind)=DMatrix(Ih,1,0);

      //Lambdas have to be  positive
      for (int i=0; i<Ih;i++)
	lb_x(lmbdh_i_ind[i])=DMatrix(L,1,0);


      //=================SOFT CONSTRAINTS========================================
      for (int i=0; i<Ih; i++)
	{
	  //Extract the i-th Wrench cone and store it in a DMatrix
	  DMatrix Whi(Wh_i_list[h][i]->rows(),Wh_i_list[h][i]->cols());
	  for (int r=0; r<Wh_i_list[h][i]->rows(); r++)
	    for (int c=0;c<Wh_i_list[h][i]->cols() ; c++)
	      Whi(r,c)=(*Wh_i_list[h][i])(r,c);

	  g_s=vertcat(g_s, (x(lmbdh_i_ind[i]).trans().mul((x(nh_ind).trans().mul(Whi)).trans())).trans()+x(eh_ind)+x(zh_i_ind[i]));
          ub_g_s=vertcat(ub_g_s,DMatrix(1,1,-1));
          lb_g_s=vertcat(lb_g_s,DMatrix(1,1,std::numeric_limits<double>::infinity()*(-1)));
	}

    }

  g=vertcat(g_P,g_T); g=vertcat(g,g_lmbde); g=vertcat(g,g_s);
  ub_g=vertcat(ub_g_P,ub_g_T);  ub_g=vertcat(ub_g,ub_g_lmbde); ub_g=vertcat(ub_g,ub_g_s);
  lb_g=vertcat(lb_g_P,lb_g_T);  lb_g=vertcat(lb_g,lb_g_lmbde); lb_g=vertcat(lb_g,lb_g_s);

  // std::cout<<"f: "<<std::endl<<f<<std::endl;
  // std::cout<<"g: "<<std::endl<<g<<std::endl;
  // std::cout<<"ub_g: "<<std::endl<<ub_g<<std::endl;
  // std::cout<<"lb_g: "<<std::endl<<lb_g<<std::endl;
  // std::cout<<"ub_x: "<<std::endl<<ub_x<<std::endl;
  // std::cout<<"lb_x: "<<std::endl<<lb_x<<std::endl;

  SXFunction nlp(nlpIn("x",x),nlpOut("f",f,"g",g));
  nlp_solver=IpoptSolver(nlp);
  nlp_solver.setOption("print_level",0);


  nlp_solver.init();

  nlp_solver.setInput( x0, "x0");
  nlp_solver.setInput(lb_x, "lbx");
  nlp_solver.setInput(ub_x, "ubx");
  nlp_solver.setInput(lb_g, "lbg");
  nlp_solver.setInput(ub_g, "ubg");

  nlp_solver.evaluate();

  std::cout << std::setw(30) << "Objective: " << nlp_solver.output("f").getDescription() << std::endl;
  std::cout << std::setw(30) << "Primal solution: " << nlp_solver.output("x").getDescription() << std::endl;

}
//-------------------------------------------------------------------------------------------------
void computeWarpedHyperplanes()
{
  std::vector<Eigen::MatrixXd*> Ph_list;
  Eigen::MatrixXd mat(2,2);
  mat(0,0)=0.7; mat(1,0)=2; mat(0,1)=-0.5; mat(1,1)=1;
  Ph_list.push_back(new Eigen::MatrixXd(mat));
  mat(0,0)=0.7; mat(1,0)=2; mat(0,1)=1; mat(1,1)=-1;
  Ph_list.push_back(new Eigen::MatrixXd(mat));
  mat(0,0)=1; mat(1,0)=-1; mat(0,1)=-0.8; mat(1,1)=-0.6; 
  Ph_list.push_back(new Eigen::MatrixXd(mat));
  mat(0,0)=-0.8; mat(1,0)=-0.6; mat(0,1)=-0.5; mat(1,1)=1;
  Ph_list.push_back(new Eigen::MatrixXd(mat));


  std::vector<std::vector<Eigen::MatrixXd* > > Wh_i_list(Ph_list.size());
  mat(0,0)=-0.6; mat(1,0)=0.6; mat(0,1)=-0.2; mat(1,1)=1;
  Wh_i_list[0].push_back(new Eigen::MatrixXd(mat));
  mat(0,0)=0.2; mat(1,0)=1.8; mat(0,1)=0.8; mat(1,1)=1.5;
  Wh_i_list[0].push_back(new Eigen::MatrixXd(mat));
  //  mat(0,0)=0; mat(1,0)=1.1; mat(0,1)=1; mat(1,1)=1;
  //  Wh_i_list[0].push_back(new Eigen::MatrixXd(mat));

  mat(0,0)=0.2; mat(1,0)=1.8; mat(0,1)=0.8; mat(1,1)=1.5;
  Wh_i_list[1].push_back(new Eigen::MatrixXd(mat));
  mat(0,0)=0.5; mat(1,0)=-0.55; mat(0,1)=1.2; mat(1,1)=0.1;
  Wh_i_list[1].push_back(new Eigen::MatrixXd(mat));

  mat(0,0)=0.5; mat(1,0)=-0.55; mat(0,1)=1.2; mat(1,1)=0.1;
  Wh_i_list[2].push_back(new Eigen::MatrixXd(mat));
  mat(0,0)=-1; mat(1,0)=-0.6; mat(0,1)=-0.4; mat(1,1)=-0.6;
  Wh_i_list[2].push_back(new Eigen::MatrixXd(mat));

  mat(0,0)=-1; mat(1,0)=-0.6; mat(0,1)=-0.4; mat(1,1)=-0.6;
  Wh_i_list[3].push_back(new Eigen::MatrixXd(mat));
  mat(0,0)=-0.6; mat(1,0)=0.6; mat(0,1)=-0.2; mat(1,1)=1;
  Wh_i_list[3].push_back(new Eigen::MatrixXd(mat));

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
