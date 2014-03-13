// #include "../../include/icr.h"
// #include "../../include/config.h"
// #include <iostream>
// #include <sys/time.h>
// #include <Eigen/Geometry>
#include "gurobi_c++.h"

int main()
{
  int nVar=3; 
  int nConstraints=2;
  GRBEnv env;
  GRBModel model(env);
  GRBVar* x =model.addVars(NULL,NULL,NULL,NULL,NULL,nVar);
  model.update();

  GRBLinExpr* lhsides=new GRBLinExpr[nConstraints];
  char* senses=new char[nConstraints]; std::fill_n(senses,nConstraints,GRB_LESS_EQUAL);
  double* rhsides=new double[nConstraints]; std::fill_n(rhsides,nConstraints,0); 

  double* coeff=new double[nVar]; std::fill_n(coeff,nVar,0); 
  for (uint i=0; i<nConstraints;i++)
    lhsides[i].addTerms(coeff,x,nVar);
   
  try
    {
      GRBConstr* constrs=model.addConstrs(lhsides,senses,rhsides,NULL,nConstraints);
    }
  catch(GRBException e) 
    {
      std::cout << "Error code = " << e.getErrorCode() << std::endl;
      std::cout << e.getMessage() << std::endl;
    }

  delete[] x;
  delete[] lhsides;
  delete[] senses;
  delete[] rhsides;
  delete[] coeff;

  return 0;
}
