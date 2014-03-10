#include "../include/grasp.h"
#include "../include/ows.h"
#include "../include/debug.h"
#include <iostream>
#include "assert.h"
#include <sys/time.h>

namespace ICR
{
  //------------------------------------------------------------------
  //------------------------------------------------------------------
  Grasp::Grasp() : initialized_(false) , num_fingers_(0),num_grasp_wrenches_(0){}
  //------------------------------------------------------------------
  Grasp::Grasp(Grasp const& src) :  fingers_(src.fingers_), initialized_(src.initialized_),
				    obj_(src.obj_),gws_(src.gws_),num_fingers_(src.num_fingers_),num_grasp_wrenches_(src.num_grasp_wrenches_) {}
  //------------------------------------------------------------------
  Grasp& Grasp::operator=(Grasp const& src)
  {
    if (this !=&src)
      {
	fingers_=src.fingers_;
	initialized_=src.initialized_;
	obj_=src.obj_;
	gws_=src.gws_;
	num_fingers_=src.num_fingers_;
	num_grasp_wrenches_=src.num_grasp_wrenches_;
      }
    return *this;
  }
  //------------------------------------------------------------------
  std::ostream& operator<<(std::ostream& stream, Grasp const& grasp)
  {
    stream <<'\n'<<"GRASP: "<<'\n'
	   <<"Number of fingers: " << grasp.num_fingers_ <<'\n'
	   <<"Number of grasp wrenches: " << grasp.num_grasp_wrenches_ <<'\n'
	   <<"Is initialized: "<<std::boolalpha<<grasp.initialized_<<'\n'<<'\n';

    return stream;
  }
  //------------------------------------------------------------------
  Grasp::~Grasp(){clear();}
  //------------------------------------------------------------------
  void Grasp::clear()
  {
    for(uint i=0; i < num_fingers_; i++)
      delete fingers_[i];

    fingers_.clear();
  }
  //------------------------------------------------------------------
  void Grasp::computeGWS()
  {
    SharedDoublePtr wrenches(new double[num_grasp_wrenches_*gws_->getDimension()]);
    double* wrench_it=wrenches.get();
    uint wrench_cone_cardinality=0;
    //ConstIndexListIterator curr_pt;

    for(uint i=0;i<num_fingers_;i++)
      {
	wrench_cone_cardinality=fingers_[i]->getOWS()->getLimitSurface()->getNumPrimitiveWrenches();
	// curr_pt=fingers_[i]->getCenterPointPatch()->patch_ids_.begin();
	// for(uint j=0; j < fingers_[i]->getCenterPointPatch()->patch_ids_.size(); j++)
	for(ConstIndexListIterator curr_pt = fingers_[i]->getCenterPointPatch()->patch_ids_.begin(); curr_pt != fingers_[i]->getCenterPointPatch()->patch_ids_.end();curr_pt++)
	  {
	    Eigen::Map<Matrix6Xd >(wrench_it,6,wrench_cone_cardinality)=*(fingers_[i]->getOWS()->getWrenchCone(*curr_pt)->getWrenches());
	    wrench_it+=gws_->getDimension()*wrench_cone_cardinality;
	  }
      }  

    gws_->setWrenches(wrenches,num_grasp_wrenches_);
    gws_->computeConvexHull();

  }
  //------------------------------------------------------------------
  PatchListPtr Grasp::computePatches(Finger* new_finger)
  {
    PatchListPtr patches;
    if(new_finger->getContactModel()->getModelType()==Single_Point)
      { 
	if( fingers_.size() > 0)
	  {
	    patches=fingers_[0]->getPatches();
	    return patches;
	  }

	patches=PatchListPtr(new std::vector<Patch* >);
	patches->reserve(obj_->getNumCp());
	for(uint centerpoint_id=0; centerpoint_id < obj_->getNumCp();centerpoint_id++)
	  patches->push_back(new Patch(centerpoint_id)); 
      }
    else if(new_finger->getContactModel()->getModelType()==Multi_Point)
      {
	for(uint id=0; id < fingers_.size(); id++)
	  {
	    if(dynamic_cast<MultiPointContactModel*>(fingers_[id]->getContactModel())->getInclusionRule() == dynamic_cast<MultiPointContactModel*>(new_finger->getContactModel())->getInclusionRule())
	      {
		patches=fingers_[id]->getPatches();
		return patches;
	      }
	  }
	patches=PatchListPtr(new std::vector<Patch* >);
	patches->reserve(obj_->getNumCp());
	for(uint centerpoint_id=0; centerpoint_id < obj_->getNumCp();centerpoint_id++)
	  patches->push_back(new Patch(centerpoint_id,*obj_,*dynamic_cast<MultiPointContactModel*>(new_finger->getContactModel())->getInclusionRule()));

      }
    else
      {
	std::cout<<"Error in Grasp::computePatches - Invalid contact model type. Exiting..."<<std::endl;
	exit(0);
      }
    return patches;
  }
  //------------------------------------------------------------------
  OWSPtr Grasp::computeOWS(Finger const* new_finger)
  {
    OWSPtr ows;
    for(uint id=0; id < fingers_.size(); id++)
      {
	if( *fingers_[id]->getContactModel()->getLimitSurface() == *new_finger->getContactModel()->getLimitSurface())
	  {
	    ows=fingers_[id]->getOWS();
	    return ows;
	  }
      }

    ows=OWSPtr(new OWS());
    ows->init(*obj_,*new_finger->getContactModel()->getLimitSurface());
    return ows;
  }
  //------------------------------------------------------------------
  void Grasp::addFinger(FingerParameters const& param, uint centerpoint_id)
  {

    Finger* new_finger=new Finger(param); 

    new_finger->setName(param.getName());
    new_finger->init(centerpoint_id,computePatches(new_finger),computeOWS(new_finger));
    num_grasp_wrenches_+=(new_finger->getCenterPointPatch()->patch_ids_.size())*(new_finger->getContactModel()->getLimitSurface()->getNumPrimitiveWrenches());
    fingers_.push_back(new_finger); 
  }
  //------------------------------------------------------------------
  Finger const* Grasp::getFinger(uint id) const {return fingers_.at(id);}
  //------------------------------------------------------------------
  uint Grasp::getNumFingers()const{return num_fingers_;}
  //------------------------------------------------------------------
  void Grasp::init(FParamList const& f_param_list,const TargetObjectPtr obj,VectorXui const& centerpoint_ids)
  {
    assert(f_param_list.size()>0);
    assert(f_param_list.size()==(uint)centerpoint_ids.size());
    assert((bool)obj);

    if(initialized_)//Clean up possible previous grasp
      clear();

    gws_.reset(new DiscreteWrenchSpace(6)); //Hard coded for 6D wrench space right now - should be changed to allow for 3D wrenches associated to 2D grasps ...

    num_fingers_=f_param_list.size();
    for (uint i=0;i < num_fingers_;i++)
      assert(centerpoint_ids(i) < obj->getNumCp()); //make sure all centerpoint ids are valid

    obj_=obj;//const_cast<TargetObject*>(obj);
    for (uint i=0; i< num_fingers_;i++)
      addFinger(f_param_list[i],(uint)centerpoint_ids(i));

    computeGWS();
    initialized_=true;
  }
  //------------------------------------------------------------------
  void Grasp::setCenterPointId(uint finger_id,uint centerpoint_id)
  {
    //Not tested yet!! Changing the center contact point requires updating the number of grasp
    //wrenches (the number of patch points might differ between old and new center point) and
    //recomputing the GWS
    assert(initialized_);
    assert(centerpoint_id <= obj_->getNumCp()-1);
    num_grasp_wrenches_-=(fingers_.at(finger_id)->getCenterPointPatch()->patch_ids_.size())*(fingers_.at(finger_id)->getContactModel()->getLimitSurface()->getNumPrimitiveWrenches());
    fingers_.at(finger_id)->setCenterPointId(centerpoint_id);
    num_grasp_wrenches_+=(fingers_.at(finger_id)->getCenterPointPatch()->patch_ids_.size())*(fingers_.at(finger_id)->getContactModel()->getLimitSurface()->getNumPrimitiveWrenches());
    gws_.reset(new DiscreteWrenchSpace(6)); //old gws_ needs to be deleted, since DiscreteWrenchSpace::conv_hull_.runQhull() can only be executed once in the lifetime of conv_hull_
    computeGWS();
  }
  //------------------------------------------------------------------
  void Grasp::setCenterPointIds(VectorXui const& centerpoint_ids)
  {
    //Not tested yet!!
    assert(initialized_);
    assert(centerpoint_ids.size()==(int)num_fingers_);

    for(uint finger_id=0; finger_id < num_fingers_; finger_id++)
      {
	assert(centerpoint_ids(finger_id) <= obj_->getNumCp()-1);
	num_grasp_wrenches_-=(fingers_.at(finger_id)->getCenterPointPatch()->patch_ids_.size())*(fingers_.at(finger_id)->getContactModel()->getLimitSurface()->getNumPrimitiveWrenches());
	fingers_.at(finger_id)->setCenterPointId(centerpoint_ids(finger_id));
	num_grasp_wrenches_+=(fingers_.at(finger_id)->getCenterPointPatch()->patch_ids_.size())*(fingers_.at(finger_id)->getContactModel()->getLimitSurface()->getNumPrimitiveWrenches());
      }

    gws_.reset( new DiscreteWrenchSpace(6));  //old gws_ needs to be deleted, since DiscreteWrenchSpace::conv_hull_.runQhull() can only be executed once
    computeGWS();
  }
  //------------------------------------------------------------------
  bool Grasp::isInitialized()const{return initialized_;}
  //------------------------------------------------------------------
  const TargetObjectPtr Grasp::getParentObj()const{return obj_;}
  //------------------------------------------------------------------
  const DiscreteWrenchSpacePtr Grasp::getGWS()const{return gws_;}
  //------------------------------------------------------------------
  VectorXui generateRandomGrasp(TargetObjectPtr object,const FParamList f_parameters, WrenchSpacePtr const tws)
  {
    assert(object->getNumCp() >= f_parameters.size());
    struct timeval c_time;
    uint n_fingers=f_parameters.size();
    VectorXui centerpoint_ids(n_fingers);

    while(1)
      {
	std::list<uint> l;
	for(uint i=0; i<n_fingers;i++)
	  {
	    while (1)
	      {
		gettimeofday(&c_time,0);
		srand(c_time.tv_usec);    
		bool is_unique=true;
		uint rand_val=(rand() % (object->getNumCp()));
		for (uint j=0; j<i;j++)
		  if (centerpoint_ids(j)==rand_val)
		    is_unique=false;
         
		if (is_unique)
		  {
		    centerpoint_ids(i)=rand_val;
		    break;
		  }
	      }
	  }

	Grasp grasp;
	grasp.init(f_parameters,object,centerpoint_ids);
        
	if (tws->getWrenchSpaceType()==Spherical)
          {
	    if(grasp.getGWS()->getOcInsphereRadius()  >= dynamic_cast<const SphericalWrenchSpace*>(tws.get())->getRadius())
              break;
          }
	else if (tws->getWrenchSpaceType()==Discrete)
	  {
	    //Should ascertain that the TWS contains the origin ...
	    uint nTW = dynamic_cast<const DiscreteTaskWrenchSpace*>(tws.get())->getNumWrenches(); 
	    Eigen::MatrixXd TWS(6,nTW);
	    doubleArrayToEigenMatrix(dynamic_cast<const DiscreteTaskWrenchSpace*>(tws.get())->getWrenches().get(),TWS);

            bool contained =true;
	    Eigen::MatrixXd nh(6,1); Eigen::MatrixXd eh(1,1);
	    facetT* curr_f= grasp.getGWS()->getConvexHull()->beginFacet().getFacetT();
	    for (uint h=0; h< grasp.getGWS()->getConvexHull()->facetCount();h++)
	      {
		doubleArrayToEigenMatrix(curr_f->normal,nh);
                eh(0,0)=-curr_f->offset;
                nh=-nh;

		if((TWS.transpose()*nh+eh.replicate(nTW,1)).minCoeff() <= 0)
		  {
		    contained =false;
                    break;
                  }
		curr_f=curr_f->next;          
	      }

	    if (contained)
	      break;

          }
	else
	  {
	    std::cout<<"Error in generateRandomGrasp(TargetObjectPtr object,const FParamList f_parameters,const WrenchSpace* const tws) - Invalid task wrench space type, cannot compute a valid random grasp. Exiting ..."<<std::endl;
	    exit(0);
          }
      }
    return centerpoint_ids;
  }
  //------------------------------------------------------------------
}//namespace ICR
