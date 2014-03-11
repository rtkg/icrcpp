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
int main()
{
  struct timeval start, end;
  double c_time;

  //PARAMETERS
  std::string name="Fish_5k";
  std::string save_path="/home/rkg/Data/Coding/ICR/experiments/shift/";
  uint nG=1000;
  uint nF=4;
  double g_scale=0.1; //scale used for gravity
  Eigen::VectorXd mu(3); mu<<0.2,0.5,0.8;
  Eigen::VectorXd L(3); L<<6,8,10;   //consider that one wrench will be added (corresponding to the normal force)

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

  remove((save_path+name+"_P.txt").c_str());
  remove((save_path+name+"_CT.txt").c_str());
  remove((save_path+name+"_ICR_PW.txt").c_str());
  remove((save_path+name+"_ICR_CC.txt").c_str());

  uint count=0;
  Eigen::Vector3d c_times; 
  FParamList f_parameters(nF);
  FingerParameters parameters;
  VectorXui G;
  Eigen::VectorXd p_values(3,1);
  for (uint i=0;i<mu.size();i++)
    for (uint j=0;j<L.size();j++)
      for (uint n=0; n<nG;n++)
	{
	  //SET FINGER PARAMETERS
          parameters.setFrictionalContact(1, L(j), mu(i));
          p_values(0)=nF; p_values(1)=mu(i); p_values(2)=L(j);
	  writeValues(save_path+name+"_P.txt",p_values);
          parameters.setWrenchIncusionTestType(Primitive);
	  for (uint p=0;p<nF;p++)
	    f_parameters[p]=parameters;

	  //GENERATE RANDOM GRASP INDICES
	  G=generateRandomGrasp(obj,f_parameters,tws);
     
	  //COMPUTE ICR WITH PRIMITIVE WRENCH INCLUSION
	  GraspPtr prototype_grasp(new Grasp());
	  prototype_grasp->init(f_parameters,obj,G);
	  gettimeofday(&start,0);
          prototype_grasp->setCenterPointIds(G);
	  SearchZonesPtr search_zones(new SearchZones(prototype_grasp));
	  search_zones->setTaskWrenchSpace(tws);
	  search_zones->computeShiftedSearchZones();
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
	  count++;
	  std::cout<<nG*L.size()*mu.size()-count <<" checks to go..."<<std::endl;
	}

  return 0;
}
