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
  tmp<<0.5, 0.4, -0.1, 0.6, -0.4, -0.4, 0.1, -0.3;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tw_e=tmp.transpose();

  SharedDoublePtr tw_a;
  double * a;
  eigenMatrixToDoubleArray(tw_e,a);
  tw_a.reset(a);

  DiscreteWrenchSpace tws(dim,tw_a,num_t_w);
  tws.computeConvexHull(); 

  std::cout<<tws<<std::endl;


  return 0;
}
