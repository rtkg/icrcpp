#include "../../include/icr.h"
#include "../../include/config.h"
#include <iostream>
#include <sys/time.h>
#include <Eigen/Geometry>

using namespace ICR;

int main()
{
  struct timeval start, end;

  //PARAMETERS
  std::string name="cup";
  std::string save_path="/home/rkg/Data/Coding/ICR/experiments/gait/";
  uint nF=3;
  double g_scale=0.1; //scale used for gravity
  double mu=0.5;
  uint L=8;
  uint p_id=2;

  //LOAD TARGET OBJECT
  ObjectLoader obj_loader;
  obj_loader.loadObject("../models/"+name+".obj",name);
  TargetObjectPtr obj=obj_loader.getObject(); obj->computeCentroid();

  //GRASP
  VectorXui G(nF); G(0)=643; G(1)=940; G(2)=1298;
  FParamList f_parameters(nF);
  FingerParameters parameters;
  parameters.setFrictionalContact(1, L, mu);
  parameters.setWrenchIncusionTestType(Convex_Combination);
  for (uint i=0; i<nF;i++)
    f_parameters[i]=parameters;
  GraspPtr prototype_grasp(new Grasp());
  prototype_grasp->init(f_parameters,obj,G);

  //COMPOSE A NEW TWS 
  double *a=new double[6]; std::fill_n(a, 6, 0); 
  SharedDoublePtr t_wrenches(a); 
  DiscreteTaskWrenchSpacePtr tws(new DiscreteTaskWrenchSpace(6,t_wrenches,1));  // 	DiscreteTaskWrenchSpace (uint dimension, SharedDoublePtr wrenches, uint num_wrenches)

  //COMPUTE SEARCH ZONES
  SearchZonesPtr search_zones(new SearchZones(prototype_grasp));
  search_zones->setTaskWrenchSpace(tws);
  search_zones->computePrioritizedSearchZones(p_id);

  //COMPUTE ICR 
  IndependentContactRegions icr(search_zones,prototype_grasp);
  icr.computeICR(Full);
  remove((save_path+name+"_icr_1.txt").c_str());
  icr.writeICR(save_path+name+"_icr_1.txt");
  std::cout<<"icr:"<<icr<<std::endl;

  return 0;
}

