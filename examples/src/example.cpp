#include "../../include/icr.h"
#include <iostream>
#include <sys/time.h>
#include <time.h>

using namespace ICR;

int main()
{ 
  uint dim=2;
  uint num_t_w=4;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tmp(num_t_w,dim);
  tmp<<0, 0, 0, 2, 1, 1, 2, 0;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tw_e=tmp.transpose();

  SharedDoublePtr tw_a;
  double * a;
  eigenMatrixToDoubleArray(tw_e,a);

  std::cout<<"tw_e"<<std::endl<<tw_e<<std::endl;

  orgQhull::Qhull conv_hull_;
   conv_hull_.runQhull("", dim,num_t_w,a ,"Qx Qt");

   conv_hull_.defineVertexNeighborFacets();
   double area_=conv_hull_.area();
   double volume_=conv_hull_.volume();
   int  num_vtx_=conv_hull_.vertexCount();
   int  num_facets_=conv_hull_.facetCount();

   facetT* curr_f=conv_hull_.beginFacet().getFacetT();
   double r_oc_insphere_=-(curr_f->offset);


   for(int i=0;i< num_facets_;i++)
     {
       curr_f->id=i; //Replaces the Qhull runtime indexing with indices 0 - num_facets_
       r_oc_insphere_ = (-(curr_f->offset) < r_oc_insphere_) ? -(curr_f->offset) : r_oc_insphere_;
       std::cout<<"normal: "<<-Eigen::Map<Eigen::Matrix<double,1,2> >(curr_f->normal)<<" offset: "<<-(curr_f->offset)  <<std::endl;
   
       curr_f=curr_f->next;
     }
  
   std::cout<<"num_facets: "<<num_facets_<<std::endl;
   std::cout<<"r_oc_insphere: "<<r_oc_insphere_<<std::endl;


  // tw_a.reset(a);

  // DiscreteWrenchSpace tws(dim,tw_a,num_t_w);
  // tws.computeConvexHull(); 

  // std::cout<<tws<<std::endl;


  return 0;
}
