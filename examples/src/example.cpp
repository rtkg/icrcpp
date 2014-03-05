#include "../../include/icr.h"
#include <iostream>
#include <sys/time.h>
#include <time.h>

using namespace ICR;

int main()
{ 
  //Load a new target object 
  ObjectLoader obj_loader;
  obj_loader.loadObject("../models/beer_can.obj","beer_can");
 
  //Create a list of default finger parameters (default parameters defined in config.h) and a 
  //vector of centerpoint contact id's for the 
  FParamList f_parameters;
  FingerParameters parameters;
  parameters.setFrictionalContact (1, 8, 0.8);
  parameters.setWrenchIncusionTestType(Convex_Combination);

  //parameters.setSoftFingerContact(1, 8, 0.8, 0.8);
  uint n_fingers=3; //number of fingers in the prototype grasp
  for (int i=0;i<n_fingers;i++)
    f_parameters.push_back(parameters);
  
  VectorXui centerpoint_ids(n_fingers);
  //centerpoint_ids=generateRandomGrasp(obj_loader.getObject(),f_parameters);
  //centerpoint_ids << 2036, 4508;  //centerpoint_ids<<346, 3131, 6168;// centerpoint_ids << 1838, 4526, 4362, 1083, 793; centerpoint_ids << 2036, 4508;
  centerpoint_ids<<346, 3131, 6168;

  //Create a prototype grasp and search zones, the parameter alpha is the scale of the largest
  //origin-centered ball contained by the Grasp wrench space of the prototype grasp
  GraspPtr prototype_grasp(new Grasp());
  prototype_grasp->init(f_parameters,obj_loader.getObject(),centerpoint_ids);
  SearchZonesPtr search_zones(new SearchZones(prototype_grasp));

  //  Create a new Discrete Task Wrench Space
  double *a=new double[6]; std::fill_n(a, 6, 0); 

  SharedDoublePtr t_wrenches(a); 
  DiscreteTaskWrenchSpacePtr tws(new DiscreteTaskWrenchSpace(6,t_wrenches,1));  // 	DiscreteTaskWrenchSpace (uint dimension, SharedDoublePtr wrenches, uint num_wrenches)
  search_zones->setTaskWrenchSpace(tws);
  //  SphericalWrenchSpacePtr tws(new SphericalWrenchSpace(6,1e-5));
  //  search_zones->setTaskWrenchSpace(tws);

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

  std::cout<<"HF: "<<icr<<std::endl;
   std::cout<<"Number of ICR points: "<<icr.getNumICRPoints()<<std::endl;


  //  parameters.setSoftFingerContact(1, 8, 0.8, 0.8);
  //   for (int i=0;i<n_fingers;i++)
  //     f_parameters[i]=parameters;

  //   prototype_grasp.reset(new Grasp());
  //   prototype_grasp->init(f_parameters,obj_loader.getObject(),centerpoint_ids);
  //   search_zones.reset(new SearchZones(prototype_grasp));
  //   search_zones->setTaskWrenchSpace(tws);
  //   search_zones->computeShiftedSearchZones();
  //   IndependentContactRegions icr2(search_zones,prototype_grasp);
  // icr2.computeICR(BFS);

  //   std::cout<<"SF: "<<icr2;
  return 0;
}
