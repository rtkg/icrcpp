#include "../../include/icr.h"
#include <iostream>
#include <sys/time.h>
#include <time.h>

using namespace ICR;

int main()
{ 
  //Load a new target object 
  ObjectLoader obj_loader;
  obj_loader.loadObject("../models/Fish_5k.obj","Fish_5k");
 
  //Create a list of default finger parameters (default parameters defined in config.h) and a 
  //vector of centerpoint contact id's for the 
  FParamList f_parameters;
  FingerParameters parameters;
  parameters.setFrictionalContact (1, 6, 0.8);
  parameters.setWrenchIncusionTestType(Convex_Combination);
  parameters.setContactModelType(Single_Point); 
  parameters.setInclusionRuleType(Sphere);
  parameters.setInclusionRuleParameter(10);
  //parameters.setSoftFingerContact(1, 8, 0.8, 0.8);
  uint n_fingers=4; //number of fingers in the prototype grasp
  for (int i=0;i<n_fingers;i++)
    f_parameters.push_back(parameters);
  
  VectorXui centerpoint_ids(n_fingers);

    WrenchSpacePtr tws(new SphericalWrenchSpace(6,1e-3));
    //   centerpoint_ids=generateRandomGrasp(obj_loader.getObject(),f_parameters,tws);
    centerpoint_ids<<1929,910, 642,1554;

  //Create a prototype grasp and search zones, the parameter alpha is the scale of the largest
  //origin-centered ball contained by the Grasp wrench space of the prototype grasp


  GraspPtr prototype_grasp(new Grasp());
  prototype_grasp->init(f_parameters,obj_loader.getObject(),centerpoint_ids);
  prototype_grasp->setCenterPointIds(centerpoint_ids);
  SearchZonesPtr search_zones(new SearchZones(prototype_grasp));

  // //  Create a new Discrete Task Wrench Space
  // double *a=new double[6]; std::fill_n(a, 6, 0); 
  // SharedDoublePtr t_wrenches(a); 
  // DiscreteTaskWrenchSpacePtr tws(new DiscreteTaskWrenchSpace(6,t_wrenches,1));  // 	DiscreteTaskWrenchSpace (uint dimension, SharedDoublePtr wrenches, uint num_wrenches)
  // search_zones->setTaskWrenchSpace(tws);

   search_zones->setTaskWrenchSpace(tws);

  struct timeval start, end;
  double c_time;
  gettimeofday(&start,0);
  //search_zones->computePrioritizedSearchZones(1);
  search_zones->computeShiftedSearchZones();
  gettimeofday(&end,0);
  c_time = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
  std::cout<<"Computation time search zones: "<<c_time<<" s"<<std::endl;


  //Create and plot the Independent Contact Regions
  IndependentContactRegions icr(search_zones,prototype_grasp);

  gettimeofday(&start,0);
  icr.computeICR(BFS); //explore points for inclusion via a Breadth-First Search with the prototype grasp's centerpoints as root nodes
  gettimeofday(&end,0);
  c_time = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
  std::cout<<"Computation time icr: "<<c_time<<" s"<<std::endl;
  std::cout<<"icr: "<<icr<<std::endl;
  
  remove("/home/rkg/ros/generic/aass_icr/libicr/icrcpp/debug/icr.txt");
  icr.writeICR("/home/rkg/ros/generic/aass_icr/libicr/icrcpp/debug/icr.txt");

  return 0;
}
