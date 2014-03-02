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
  uint n_fingers=3; //number of fingers in the prototype grasp
  for (int i=0;i<n_fingers;i++)
    f_parameters.push_back(parameters);
  
  VectorXui centerpoint_ids(n_fingers);
   //centerpoint_ids=generateRandomGrasp(obj_loader.getObject(),f_parameters);
   centerpoint_ids<<346, 3131, 6168;// centerpoint_ids << 1838, 4526, 4362, 1083, 793;

  //Create a prototype grasp and search zones, the parameter alpha is the scale of the largest
  //origin-centered ball contained by the Grasp wrench space of the prototype grasp
  GraspPtr prototype_grasp(new Grasp());
  prototype_grasp->init(f_parameters,obj_loader.getObject(),centerpoint_ids);
  SearchZonesPtr search_zones(new SearchZones(prototype_grasp));

  //  Create a new Discrete Task Wrench Space
  double *a;
  eigenMatrixToDoubleArray(Eigen::MatrixXd::Identity(6,6)*1.0e-4,a);
  SharedDoublePtr t_wrenches(a); 
  DiscreteTaskWrenchSpacePtr tws(new DiscreteTaskWrenchSpace(6,t_wrenches,6));  // 	DiscreteTaskWrenchSpace (uint dimension, SharedDoublePtr wrenches, uint num_wrenches)
  search_zones->setTaskWrenchSpace(tws);

  //Create a set of conditioning patches;
  std::vector< ContactRegion * > conditioning_patches(centerpoint_ids.rows());
  double inclusion_radius= 5; //radius of the sphere used for patch inclusion
  for (int n=0; n<centerpoint_ids.rows();n++)
    {
      Patch* patch(new Patch(centerpoint_ids(n), *obj_loader.getObject(),InclusionRule(inclusion_radius,Sphere,true)));
      std::cout<<(*patch)<<std::endl;
      conditioning_patches[n]=new ContactRegion();
      conditioning_patches[n]->push_back(patch);
    }

  search_zones->computeConditionedSearchZones(conditioning_patches);
  //  // //Create and plot the Independent Contact Regions
  //  // IndependentContactRegions icr(search_zones,prototype_grasp);
  //  // icr.computeICR(BFS); //explore points for inclusion via a Breadth-First Search with the prototype grasp's centerpoints as root nodes
  //  // std::cout<<icr;

  //clean up
  for (int n=0; n<centerpoint_ids.rows();n++)
    {
      for (int p=0; p<conditioning_patches[n]->size();p++)
	delete conditioning_patches[n]->at(p);

      delete conditioning_patches[n];
    }
  return 0;
}
