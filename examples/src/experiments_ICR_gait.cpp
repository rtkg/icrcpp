#include "../../include/icr.h"
#include "../../include/config.h"
#include <iostream>
#include <sys/time.h>
#include <Eigen/Geometry>

using namespace ICR;

//--------------------------------------------------------------------------------
void writeValues(const std::string& path, const Eigen::VectorXd& values)
{
  assert(values.size() > 0);
  FILE* f=fopen (path.c_str(),"a");
  if(!f)
    {
      std::cout<<"Error in writeValues(const std::string& path) - Couldn't open file. Exiting..."<<std::endl;
      exit(0);
    }

  for(uint id=0; id < values.size(); id++)
    fprintf(f, "%f ", values(id));

  fprintf(f,"\n");
  fclose (f);
}
// //--------------------------------------------------------------------------------
// void readGrasps(const std::string& path, Eigen::MatrixXi& G_mat)
// {
//   uint M=G_mat.rows();
//   uint N=G_mat.cols();

//   std::ifstream f;
//   f.open(path);

//   if(f.fail())
//     {
//       std::cout<<"Couldn't read file "<<path<<". Exiting..."<<std::endl;
//       exit(0);
//     }

//   for (uint i=0; i<M;i++)
//     for (uint j=0;j<N;j++)
//       {
// 	f>>G_mat(i,j);
// 	G_mat(i,j)=G_mat(i,j)-1; //revert from matalb indexing to C++ 0-based indexing
//       }
// }
//--------------------------------------------------------------------------------
int main()
{
  struct timeval start, end;

  //PARAMETERS
  std::string name="Fish_5k";
  std::string save_path="/home/rkg/Data/Coding/ICR/experiments/gait/";
  uint nF=3;
  double mu=0.8;
  uint L=8;
  VectorXui G(nF); G<<962,2263, 1110;//G<<643, 940, 1298;[963 2264 1111]
  uint p_id=0;
  double g_scale=0.1;

  //LOAD TARGET OBJECT
  ObjectLoader obj_loader;
  obj_loader.loadObject("../models/"+name+".obj",name);
  TargetObjectPtr obj=obj_loader.getObject(); obj->computeCentroid();

  //GENERATE GRASP
  FParamList f_parameters(nF);
  FingerParameters parameters; 
  parameters.setFrictionalContact(1, L, mu); 
  parameters.setWrenchIncusionTestType(Primitive);
  for(uint i=0;i<nF;i++)
    f_parameters[i]=parameters;
  GraspPtr prototype_grasp(new Grasp());
  prototype_grasp->init(f_parameters,obj,G);

//   //COMPUTE TASK
//   SharedDoublePtr t_wrenches(new double[72]); std::fill_n(t_wrenches.get(),72,0.0); 
// t_wrenches.get()[0]=1e-8; t_wrenches.get()[7]=1e-8; t_wrenches.get()[14]=1e-8; t_wrenches.get()[21]=1e-8; t_wrenches.get()[28]=1e-8; t_wrenches.get()[35]=1e-8; 
// t_wrenches.get()[36]=-1e-8; t_wrenches.get()[43]=-1e-8; t_wrenches.get()[50]=-1e-8; t_wrenches.get()[57]=-1e-8; t_wrenches.get()[64]=-1e-8; t_wrenches.get()[71]=-1e-8; 

//   WrenchSpacePtr tws(new DiscreteTaskWrenchSpace(6,t_wrenches,12)); 
//   WrenchSpacePtr tws_s(new SphericalWrenchSpace(6,1e-8*prototype_grasp->getGWS()->getOcInsphereRadius())); 

  //COMPOSE A NEW TWS OF THE 0-WRENCH AND THE GRAVITY WRENCH - gravity is acting in -y direction for the KIT models!!!
  Eigen::MatrixXd t_wrenches_e=Eigen::MatrixXd::Zero(12,1); t_wrenches_e(7)=-g_scale;
  Eigen::Vector3d g_force; g_force<<0,-g_scale,0;
#ifdef DIVIDE_OWS_BY_LAMBDA    
  double lambda=0;
  for (uint i=0; i<obj->getNumCp();i++)
    if (obj->getContactPoint(i)->getVertex()->norm() > lambda)
      lambda=obj->getContactPoint(i)->getVertex()->norm();
#else
  double lambda=1;
#endif

  t_wrenches_e.bottomRows(3)=(*obj->getCentroid()).cross(g_force)/lambda;
  //  SharedDoublePtr t_wrenches(new double[12]); double* temp=t_wrenches.get();
  double* temp=new double[12]; eigenMatrixToDoubleArray(t_wrenches_e,temp);
  SharedDoublePtr t_wrenches(temp);
  WrenchSpacePtr tws(new DiscreteTaskWrenchSpace(6,t_wrenches,2)); 


  //COMPUTE SEARCH ZONES
  SearchZonesPtr search_zones(new SearchZones(prototype_grasp));
  search_zones->setTaskWrenchSpace(tws);
  search_zones->computePrioritizedSearchZones(p_id);
  //search_zones->computeShiftedSearchZones();

  //COMPUTE ICR
  IndependentContactRegions icr(search_zones,prototype_grasp);
  icr.computeICR(BFS);

  remove((save_path+name+"_icr_1.txt").c_str()); icr.writeICR(save_path+name+"_icr_1.txt");

  return 0;
}
