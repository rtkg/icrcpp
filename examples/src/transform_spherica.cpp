#include "../../include/icr.h"
#include <iostream>
#include <sys/time.h>
#include <time.h>
#include <Eigen/Geometry>

using namespace ICR;

int main()
{ 
  //Load a new target object 
  ObjectLoader obj_loader;
  obj_loader.loadObject("../models/Fish_5k.obj","Fish_5k");
  TargetObjectPtr obj=obj_loader.getObject();
 
  //Create a list of default finger parameters (default parameters defined in config.h) and a 
  //vector of centerpoint contact id's for the 
  FParamList f_parameters;
  FingerParameters parameters;
  parameters.setFrictionalContact (1, 100, 0.8);
  parameters.setWrenchIncusionTestType(Primitive);
  parameters.setContactModelType(Single_Point); 
  parameters.setInclusionRuleType(Sphere);
  parameters.setInclusionRuleParameter(10);
  //parameters.setSoftFingerContact(1, 8, 0.8, 0.8);
  ICR::uint n_fingers=4; //number of fingers in the prototype grasp
  for (int i=0;i<n_fingers;i++)
    f_parameters.push_back(parameters);
  
SphericalWrenchSpacePtr tws(new SphericalWrenchSpace(6,1e-15));

  VectorXui centerpoint_ids(n_fingers);
 centerpoint_ids<<1929,910, 642,1554;

  GraspPtr prototype_grasp(new Grasp());
  prototype_grasp->init(f_parameters,obj,centerpoint_ids);
  prototype_grasp->setCenterPointIds(centerpoint_ids);
  SearchZonesPtr search_zones(new SearchZones(prototype_grasp));
   search_zones->setTaskWrenchSpace(tws);
  search_zones->computeShiftedSearchZones();
  IndependentContactRegions icr(search_zones,prototype_grasp);
  icr.computeICR(BFS); //explore points for inclusion via a Breadth-First Search with the prototype grasp's centerpoints as root nodes
  std::cout<<"icr before transform: "<<icr<<std::endl;
  std::cout<<"Num icr before transform: "<<icr.getNumICRPoints()<<std::endl;


  Eigen::Matrix3d rot;
  rot = Eigen::AngleAxisd(0.35*M_PI, Eigen::Vector3d::UnitX())
      * Eigen::AngleAxisd(0.11*M_PI, Eigen::Vector3d::UnitY())
      * Eigen::AngleAxisd(0.52*M_PI, Eigen::Vector3d::UnitZ());
  Eigen::Vector3d trans; trans(0)=0; trans(1)=0; trans(2)=0; 
  Eigen::Affine3d T;T.setIdentity(); T.translate(trans); T.rotate(rot); //T.translate(trans);//T.rotation()=rot; T.translation()=trans;

  obj->transform(T);

 prototype_grasp.reset(new Grasp());
  prototype_grasp->init(f_parameters,obj,centerpoint_ids);
search_zones.reset(new SearchZones(prototype_grasp));
   search_zones->setTaskWrenchSpace(tws);
  search_zones->computeShiftedSearchZones();
  IndependentContactRegions icr_t(search_zones,prototype_grasp);
  icr_t.computeICR(BFS); //explore points for inclusion via a Breadth-First Search with the prototype grasp's centerpoints as root nodes
  std::cout<<"icr after transform: "<<icr_t<<std::endl;
  std::cout<<"Num icr after transform: "<<icr_t.getNumICRPoints()<<std::endl;

  return 0;
}


////BENCHMARK ON FISH_5K
//   centerpoint_ids<<1929,910, 642,1554;
//  parameters.setFrictionalContact (1, 6, 0.8);
//    WrenchSpacePtr tws(new SphericalWrenchSpace(6,1e-3));

// //PRIMITIVE INCLUSION
// Number of contact regions: 4
// Contact regions computed: true
// Centerpoint id's region 0: 1929 1951 2038 2089 2116 2055 2124 2149 2106 2151 2191 2199 2242 2211 2243 
// Centerpoint id's region 1: 910 887 929 938 990 1013 1001 1030 1073 1074 989 1048 1040 1089 1118 1047 1070 1152 1198 1222 
// Centerpoint id's region 2: 642 696 
// Centerpoint id's region 3: 1554 2496 2500 2498 2501 2493 2494 


// //CONVEX COMBINATION
// IndependentContactRegions: 
// Number of contact regions: 4
// Contact regions computed: true
// Centerpoint id's region 0: 1929 1951 2038 1857 2089 2002 2116 2055 2124 2149 2105 2135 2106 2151 2191 2199 2198 2242 2197 2201 2211 2243 
// Centerpoint id's region 1: 910 887 888 929 938 939 869 876 960 990 1013 1030 875 894 905 972 989 1001 1073 1074 1040 1089 1118 950 957 1027 1047 1048 1134 1050 1152 1198 949 987 1046 1115 1070 1102 1168 1004 1281 1222 1130 1132 1151 1197 1247 1310 1150 1263 1297 1279 1321 1333 
// Centerpoint id's region 2: 642 633 679 696 
// Centerpoint id's region 3: 1554 1630 2496 2500 2498 2501 2493 1295 1334 2494 1294 1296 1332 1234 1273 1167 1260 
