#include "../../include/icr.h"
#include <iostream>
#include <sys/time.h>
#include <time.h>
#include <symbolic/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <vector>

using namespace ICR;
using namespace CasADi;

//===========GLOBAL===========================
uint K=2;
uint T=4;

DiscreteWrenchSpace* tws_;

//-------------------------------------------------------------------------------------------------
void createNLPSolver(const std::vector<Eigen::MatrixXd*>& Ph_list, const std::vector<std::vector<Eigen::MatrixXd* > >& Wh_n_list,IpoptSolver& nlp_solver)
{
  uint H=Ph_list.size();
  uint nV=0; //number of decision variables
  uint L=2;
  for (uint i=0; i<H;i++)
    nV+=Wh_n_list[i].size()*(L+1)+K+1;

  

  SXMatrix x = ssym("x",nV); //decision variable vector
  x.print();
  std::vector<int> nh_ind(K);
  double eh_ind;
  std::vector<std::vector<int> > lmbdh_i_ind;
  std::vector<int> zh_i_ind;
 
  for (uint h=0; h<H;h++)
    {
      //Build index lists for nh, eh, lbmd_h_i and zh_i for eache Hyperplane h
      uint Ih=Wh_n_list[h].size();
      uint nVh=K+1+Ih*(L+1);
      lmbdh_i_ind.resize(Ih);
      zh_i_ind.resize(Ih);

      for (int j=0; j<K; j++)
        nh_ind[j]=h*nVh+j;

      eh_ind=h*nVh+L;

      for (int i=0; i<Ih; i++)
	{
	  lmbdh_i_ind[i].resize(L);
	  for(int l=0; l<L; l++)
	    lmbdh_i_ind[i][l]=L*(i+1)+l+h*nVh+1;

	  zh_i_ind[i]=h*nVh+K+1+L*Ih+i;
	}




    }





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


  std::vector<std::vector<Eigen::MatrixXd* > > Wh_n_list(Ph_list.size());
  mat(0,0)=-0.6; mat(1,0)=0.6; mat(0,1)=-0.2; mat(1,1)=1;
  Wh_n_list[0].push_back(new Eigen::MatrixXd(mat));
  mat(0,0)=0.2; mat(1,0)=1.8; mat(0,1)=0.8; mat(1,1)=1.5;
  Wh_n_list[0].push_back(new Eigen::MatrixXd(mat));
  //F{1}.W{3}=[0 1.1; 1 1];

  mat(0,0)=0.2; mat(1,0)=1.8; mat(0,1)=0.8; mat(1,1)=1.5;
  Wh_n_list[1].push_back(new Eigen::MatrixXd(mat));
  mat(0,0)=0.5; mat(1,0)=-0.55; mat(0,1)=1.2; mat(1,1)=0.1;
  Wh_n_list[1].push_back(new Eigen::MatrixXd(mat));

  mat(0,0)=0.5; mat(1,0)=-0.55; mat(0,1)=1.2; mat(1,1)=0.1;
  Wh_n_list[2].push_back(new Eigen::MatrixXd(mat));
  mat(0,0)=-1; mat(1,0)=-0.6; mat(0,1)=-0.4; mat(1,1)=-0.6;
  Wh_n_list[2].push_back(new Eigen::MatrixXd(mat));

  mat(0,0)=-1; mat(1,0)=-0.6; mat(0,1)=-0.4; mat(1,1)=-0.6;
  Wh_n_list[3].push_back(new Eigen::MatrixXd(mat));
  mat(0,0)=-0.6; mat(1,0)=0.6; mat(0,1)=-0.2; mat(1,1)=1;
  Wh_n_list[3].push_back(new Eigen::MatrixXd(mat));


  IpoptSolver nlp_solver;
  createNLPSolver(Ph_list,Wh_n_list,nlp_solver);


  for (int i=0; i<Ph_list.size(); i++)
    delete Ph_list[i];

  for (int i=0; i<Wh_n_list.size(); i++)
    for (int j=0; j<Wh_n_list[i].size(); j++)
      delete Wh_n_list[i][j];
}
//-------------------------------------------------------------------------------------------------
int main()
{ 
  //==============TWS===============================================
  Eigen::MatrixXd tmp(T,K);
  tmp<<0, 0, 0, 2, 1, 1, 2, 0;
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
