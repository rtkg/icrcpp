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
 
  //Create a list of 5 default finger parameters (default parameters defined in config.h) and a 
  //vector of centerpoint contact id's for the 5-fingered prototype grasp
  FParamList f_parameters;
  FingerParameters parameters;
  for (int i=0;i<5;i++)
    f_parameters.push_back(parameters);
  VectorXui centerpoint_ids(5);
  centerpoint_ids << 1838, 4526, 4362, 1083, 793;

  //Create a prototype grasp and search zones, the parameter alpha is the scale of the largest
  //origin-centered ball contained by the Grasp wrench space of the prototype grasp
  GraspPtr prototype_grasp(new Grasp());
  prototype_grasp->init(f_parameters,obj_loader.getObject(),centerpoint_ids);
  //SearchZonesPtr search_zones(new SearchZones(prototype_grasp));
  //double alpha=0.5;
  //search_zones->computeShiftedSearchZones(alpha);

  //Create and plot the Independent Contact Regions
  //IndependentContactRegions icr(search_zones,prototype_grasp);
  //icr.computeICR(BFS); //explore points for inclusion via a Breadth-First Search with the prototype grasp's centerpoints as root nodes
//  std::cout<<icr;

  return 0;
}
