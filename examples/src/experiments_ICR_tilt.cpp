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
//--------------------------------------------------------------------------------
void readGrasps(const std::string& path, Eigen::MatrixXi& G_mat)
{
  uint M=G_mat.rows();
  uint N=G_mat.cols();

  std::ifstream f;
  f.open(path);

  if(f.fail())
    {
      std::cout<<"Couldn't read file "<<path<<". Exiting..."<<std::endl;
      exit(0);
    }

  for (uint i=0; i<M;i++)
    for (uint j=0;j<N;j++)
      {
	f>>G_mat(i,j);
	G_mat(i,j)=G_mat(i,j)-1; //revert from matalb indexing to C++ 0-based indexing
      }
}
//--------------------------------------------------------------------------------
int main()
{
  struct timeval start, end;
  double c_time;

  //PARAMETERS
  std::string name="Fish_5k";
  std::string save_path="/home/rkg/Data/Coding/ICR/experiments/tilt/";
  uint nG=1000;
  uint nF=4;
  double g_scale=0.1; //scale used for gravity
  double mu=0.5;
  uint L=8;

  //LOAD GRASPS
  Eigen::MatrixXi G_mat(nG,nF);
  readGrasps((save_path+name+"_G_8_0.5.txt"),G_mat);

  //LOAD TARGET OBJECT
  ObjectLoader obj_loader;
  obj_loader.loadObject("../models/"+name+".obj",name);
  TargetObjectPtr obj=obj_loader.getObject(); obj->computeCentroid();

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

  Eigen::Vector3d c_times; 
  FParamList f_parameters(nF);
  FingerParameters parameters;
  VectorXui G(nF);
  parameters.setFrictionalContact(1, L, mu);

  for (uint i=0; i<nG;i++)
    {
      parameters.setWrenchIncusionTestType(Primitive);
      for (uint p=0;p<nF;p++)
	f_parameters[p]=parameters;

      //Current grasp indices
      for (uint n=0;n<nF;n++)
	G(n)=G_mat(i,n);

      std::cout<<"G: "<<G<<std::endl;
      //Seed the random number generator and compute a random number 0 <= c < nF
      gettimeofday(&start,0);
      srand(start.tv_usec);
      int c=rand()%nF;
      std::cout<<"c: "<<c<<std::endl;
      //COMPUTE ICR WITH PRIMITIVE WRENCH INCLUSION
      GraspPtr prototype_grasp(new Grasp());
      gettimeofday(&start,0);
      prototype_grasp->init(f_parameters,obj,G);
      SearchZonesPtr search_zones(new SearchZones(prototype_grasp));
      search_zones->setTaskWrenchSpace(tws);
      search_zones->computePrioritizedSearchZones(c);
      gettimeofday(&end,0);
      c_times(0) = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
      gettimeofday(&start,0);
      IndependentContactRegions icr_pw(search_zones,prototype_grasp);
      icr_pw.computeICR(BFS);
      gettimeofday(&end,0);
      c_times(1) = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
      icr_pw.writeICR(save_path+name+"_ICR_PW.txt");
      // std::cout<<"ICR PRIMITIVE: "<<icr_pw<<std::endl;

      //COMPUTE ICR WITH LP WRENCH INCLUSION
      for (int i=0;i<nF;i++)
	f_parameters[i].setWrenchIncusionTestType(Convex_Combination);

      prototype_grasp.reset(new Grasp());
      prototype_grasp->init(f_parameters,obj,G);
      search_zones.reset(new SearchZones(prototype_grasp));
      search_zones->setTaskWrenchSpace(tws);
      search_zones->computeShiftedSearchZones();
      gettimeofday(&start,0);
      IndependentContactRegions icr_cc(search_zones,prototype_grasp);
      icr_cc.computeICR(BFS);
      gettimeofday(&end,0);
      c_times(2) = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
      icr_cc.writeICR(save_path+name+"_ICR_CC.txt");
      // std::cout<<"ICR CONVEX COMBINATION: "<<icr_cc<<std::endl;
	 
      writeValues(save_path+name+"_CT.txt",c_times);
      std::cout<<nG-i-1 <<" checks to go..."<<std::endl;
    }

  return 0;
}
