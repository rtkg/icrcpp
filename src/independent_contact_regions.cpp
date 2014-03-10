#include "../include/independent_contact_regions.h"
#include "../include/wrench_cone.h"
#include "../include/search_zones.h"
#include "../include/ows.h"
#include "../include/grasp.h"
#include "../include/target_object.h"
#include "../include/contact_point.h"
#include "../include/debug.h"
#include "../include/config.h"
#include "assert.h"
#include <thread>

#include <sys/time.h>
#include <time.h>

namespace ICR
{

  IndependentContactRegions::IndependentContactRegions() :  
    icr_computed_(false),
    num_contact_regions_(0) {contact_regions_.clear();}
  //--------------------------------------------------------------------
  IndependentContactRegions::IndependentContactRegions(const SearchZonesPtr search_zones,const GraspPtr grasp) : 
    search_zones_(search_zones), 
    grasp_(grasp),
    icr_computed_(false), 
    num_contact_regions_(0)
  {
    contact_regions_.clear();
    assert((bool)search_zones_);
    assert((bool)grasp_);
  }
  //--------------------------------------------------------------------
  IndependentContactRegions::IndependentContactRegions(IndependentContactRegions const& src) : 
    search_zones_(src.search_zones_), 
    grasp_(src.grasp_), 
    icr_computed_(src.icr_computed_),
    contact_regions_(src.contact_regions_),
    num_contact_regions_(src.num_contact_regions_){}
  //--------------------------------------------------------------------
  IndependentContactRegions& IndependentContactRegions::operator=(IndependentContactRegions const& src)
  {
    if (this !=&src)
      {
	search_zones_=src.search_zones_;
	grasp_=src.grasp_;
	icr_computed_=src.icr_computed_;
	contact_regions_=src.contact_regions_;
	num_contact_regions_=src.num_contact_regions_;
      }

    return *this;
  }

  //--------------------------------------------------------------------
  std::ostream& operator<<(std::ostream& stream, IndependentContactRegions const& icr)
  {
    stream <<'\n'<<"IndependentContactRegions: "<<'\n'
	   <<"Number of contact regions: " << icr.num_contact_regions_<<'\n'
	   <<"Contact regions computed: "<<std::boolalpha<<icr.icr_computed_<<'\n';
    if(icr.icr_computed_)
      {
	for(uint i=0;i<icr.num_contact_regions_;i++)
	  {
	    stream<<"Centerpoint id's region "<<i<<": ";
	    for(uint j=0; j < icr.contact_regions_[i]->size();j++)
	      {
		stream<<(*icr.contact_regions_[i])[j]->patch_ids_.front()<<" ";
	      }
	    stream<<'\n';
	  }
      }
    stream<<'\n';
    return stream;
  }
  //--------------------------------------------------------------------
  IndependentContactRegions::~IndependentContactRegions(){clear();}
  //--------------------------------------------------------------------
  void IndependentContactRegions::clear()
  {
    for(uint i=0;i<contact_regions_.size();i++)
      {
	delete contact_regions_[i];
      }
    contact_regions_.clear();
    num_contact_regions_ = 0;
    icr_computed_ = false;
  }
  //--------------------------------------------------------------------
#ifdef WITH_GUROBI
  bool IndependentContactRegions::convexCombinationSearchZoneInclusionTest(PrimitiveSearchZone* prim_sz,WrenchCone const* wc)const
  {
    try
      {
	// struct timeval start, end;
	// double c_time;

	if ((prim_sz->satisfied_wc_ids_.array() == wc->id_).any())
	  return true;

	uint L=wc->num_primitive_wrenches_;
	uint S=prim_sz->hyperplane_ids_.size();
	Eigen::MatrixXd W=wc->cone_;

	GRBEnv env;
	env.set(GRB_IntParam_OutputFlag,0);
	env.set(GRB_IntParam_Presolve,0);
	//env.set(GRB_IntParam_Method,0);

	GRBModel lp(env);// = GRBModel(env);
        double* obj = new double[L+1]; std::fill_n(obj, L, 0); obj[L]=1; //only minimize slack variable z 
	GRBVar* x =lp.addVars(NULL,NULL,obj,NULL,NULL,L+1); //LB = 0 
	lp.update();

	GRBLinExpr* lhs=new GRBLinExpr[S+1];
	double* rhs = new double[S+1];
	char* senses = new char[S+1];
	double* ones_L = new double[L]; std::fill_n(ones_L, L, 1); 
        double minus_one[1]; minus_one[0]=-1;        

	double* proj; 
	double offset; 
	Eigen::VectorXd proj_e;
	for (uint s=0;s<S; s++)
	  {
	    proj_e=search_zones_->hyperplane_normals_.row(prim_sz->hyperplane_ids_[s])*W;
	    offset=search_zones_->hyperplane_offsets_(prim_sz->hyperplane_ids_[s]);
	    eigenMatrixToDoubleArray(proj_e,proj);
           
	    //s-th inequality constraint
            lhs[s].addTerms(proj,x,L);
            lhs[s].addTerms(minus_one,(x+L),1);
            senses[s]=GRB_LESS_EQUAL;
            rhs[s]=-offset;         
      	 
	    delete[] proj;
	  }

	//equality constriant        
	lhs[S].addTerms(ones_L,x,L);        
	senses[S]=GRB_EQUAL;
	rhs[S]=1;      
    
	GRBConstr* constrs = lp.addConstrs(lhs,senses,rhs,NULL,S+1);
	lp.update();
        lp.optimize(); 

	double z=x[L].get(GRB_DoubleAttr_X); 

	delete[] x;
	delete[] constrs;
	delete[] lhs;
	delete[] rhs;
	delete[] senses;
	delete[] ones_L;
        delete[] obj;

	if (z > 1e-15)
          return false;

	prim_sz->satisfied_wc_ids_.conservativeResize(prim_sz->satisfied_wc_ids_.size()+1);
	prim_sz->satisfied_wc_ids_(prim_sz->satisfied_wc_ids_.size()-1)=wc->id_;
	return true;
      }
    catch(GRBException e) 
      {
	std::cout<<"Gurobi exception with error code "<<e.getErrorCode()<<" in IndependentContactRegions::convexCombinationSearchZoneInclusionTest(PrimitiveSearchZone* prim_sz,WrenchCone const* wc)const"<<std::endl;
	std::cout<<e.getMessage()<<std::endl;
	return false;
      }
  }
#endif
  //--------------------------------------------------------------------
  bool IndependentContactRegions::primitiveSearchZoneInclusionTest(PrimitiveSearchZone* prim_sz,WrenchCone const* wc)const
  {
    if ((prim_sz->satisfied_wc_ids_.array() == wc->id_).any())
      return true;

    bool constraints_satisfied;
    for(uint pw_count=0; pw_count < wc->num_primitive_wrenches_;pw_count++)
      {
	constraints_satisfied=true;
	for(uint hp_count=0; hp_count< prim_sz->hyperplane_ids_.size(); hp_count++)
	  {
	    if(search_zones_->hyperplane_normals_.row(prim_sz->hyperplane_ids_[hp_count]) * wc->cone_.col(pw_count)  + search_zones_->hyperplane_offsets_(prim_sz->hyperplane_ids_[hp_count]) > 0)
	      {
		constraints_satisfied=false;
		break;
	      }
	  }
	if(constraints_satisfied)
	  {
	    prim_sz->satisfied_wc_ids_.conservativeResize(prim_sz->satisfied_wc_ids_.size()+1);
	    prim_sz->satisfied_wc_ids_(prim_sz->satisfied_wc_ids_.size()-1)=wc->id_;
	    return true;
	  }
      }
    return false;
  }
  //--------------------------------------------------------------------
  bool IndependentContactRegions::searchZoneInclusionTest(uint region_id,Patch const* patch, const WrenchInclusionTestType wrench_inclusion_test_type)const
  {

    bool psz_satisfied;
    for(uint psz_id=0; psz_id < search_zones_->search_zones_[region_id]->size();psz_id++)//iterate over all primitive search zones of the queried search zone
      {

	//check if a wrench cone associated with any point in the patch satisfies the primitiveSearchZoneInclusionTest
	psz_satisfied=false;
	for(ConstIndexListIterator patch_point=patch->patch_ids_.begin(); patch_point != patch->patch_ids_.end(); patch_point++)
	  {
	    //struct timeval start, end;
	    //double c_time;
	    // gettimeofday(&start,0);

	    if(wrench_inclusion_test_type==Primitive)
              psz_satisfied=primitiveSearchZoneInclusionTest((*search_zones_->search_zones_[region_id])[psz_id] ,grasp_->getFinger(region_id)->getOWS()->getWrenchCone(*patch_point));
            else if(wrench_inclusion_test_type==Convex_Combination)
	      {
#ifdef WITH_GUROBI
		psz_satisfied=convexCombinationSearchZoneInclusionTest((*search_zones_->search_zones_[region_id])[psz_id] ,grasp_->getFinger(region_id)->getOWS()->getWrenchCone(*patch_point));
#else
		std::cout<<"Error in IndependentContactRegions::searchZoneInclusionTest(uint region_id,Patch const* patch, const WrenchInclusionTestType wrench_inclusion_test_type), icrcpp was not compiled with Gurobi-support. ICR::WrenchInclusionTestType::Convex_Combination is not supported. Exiting ..."<<std::endl;
		exit(0);
#endif
	      }
            else 
	      {
		std::cout<<"Error in IndependentContactRegions::searchZoneInclusionTest(uint region_id,Patch const* patch, const WrenchInclusionTestType wrench_inclusion_test_type) - invalid wrench inclusion test type specified! Exiting ..."<<std::endl;
		exit(0);
	      }

	    // gettimeofday(&end,0);
	    // c_time = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
	    // std::cout<<"Computation time: "<<c_time<<" s"<<std::endl;

	    if (psz_satisfied)
	      break;
	  }
	if(!psz_satisfied)
	  return false;
      }     
    return true;
  }
  //--------------------------------------------------------------------
  void IndependentContactRegions::computeContactRegionBFS(uint region_id, const WrenchInclusionTestType wrench_inclusion_test_type)
  {
    search_zones_->resetPrimitiveSearchZones(region_id); //Make sure the member ICR::PrimitiveSearchZone::satisfied_wc_ids_ is empty
    std::list<Node> nodes; //List of nodes for the breadth-first-search on the target object's vertices
    VectorXui explored_cp=VectorXui::Zero(grasp_->getParentObj()->getNumCp(),1); //set all entries to NOT_EXPLORED
    uint init_centerpoint_id=grasp_->getFinger(region_id)->getCenterPointPatch()->patch_ids_.front();
  
    contact_regions_[region_id]->reserve(grasp_->getParentObj()->getNumCp());
    nodes.push_back(Node(const_cast<ContactPoint*>(grasp_->getParentObj()->getContactPoint(init_centerpoint_id)))); //push the node corresponding to the centerpoint of the patch associated with the initial prototype grasp contact
    contact_regions_[region_id]->push_back(const_cast<Patch*>(grasp_->getFinger(region_id)->getCenterPointPatch())); 
    explored_cp(init_centerpoint_id)=EXPLORED_QUALIFIED;

    while(nodes.size() > 0)
      {   
	for(ConstIndexListIterator neighbor=nodes.front().contact_point_->getNeighborItBegin(); neighbor != nodes.front().contact_point_->getNeighborItEnd(); neighbor++)      
	  {
	    if(explored_cp(*neighbor)==NOT_EXPLORED)
	      {	     
		if(searchZoneInclusionTest(region_id,grasp_->getFinger(region_id)->getPatch(*neighbor),wrench_inclusion_test_type))
		  {
		    nodes.push_back(Node(const_cast<ContactPoint*>(grasp_->getParentObj()->getContactPoint(*neighbor))));
		    contact_regions_[region_id]->push_back(const_cast<Patch*>(grasp_->getFinger(region_id)->getPatch(*neighbor))); 
		    explored_cp(*neighbor)=EXPLORED_QUALIFIED;
		  }
		else
		  explored_cp(*neighbor)=EXPLORED_UNQUALIFIED;
	      }
	  }
	nodes.pop_front();
      }
  }
  //--------------------------------------------------------------------
  void IndependentContactRegions::computeICR(ICRType type)
  {
    assert((bool)search_zones_);
    assert((bool)grasp_);

    if(icr_computed_) //clean up possible previously computed contact regions
      clear();

    num_contact_regions_=search_zones_->num_search_zones_;
    contact_regions_.reserve(num_contact_regions_);


#ifdef MULTITHREAD_ICR_COMPUTATION
    std::vector<std::thread*> threads;
    threads.reserve(num_contact_regions_);
    for(uint region_id=0; region_id < num_contact_regions_; region_id++)
      {
	contact_regions_.push_back(new ContactRegion);
	if (type==BFS)
	  threads.push_back(new std::thread(&IndependentContactRegions::computeContactRegionBFS,this,region_id,grasp_->getFinger(region_id)->getWrenchInclusionTestType()));
	else if (type==Full)
	  threads.push_back(new std::thread(&IndependentContactRegions::computeContactRegionFull,this,region_id,grasp_->getFinger(region_id)->getWrenchInclusionTestType()));
	else
	  std::cout<<"Error: Unknown ICR type - cannot compute ICR!"<<std::endl;
      
      }
    for(uint thread_id=0; thread_id < num_contact_regions_; thread_id++)
      {
	threads[thread_id]->join();
	delete threads[thread_id];
      }
#else
    for(uint region_id=0; region_id < num_contact_regions_; region_id++)
      {
	contact_regions_.push_back(new ContactRegion);
	if (type==BFS)
	  computeContactRegionBFS(region_id,grasp_->getFinger(region_id)->getWrenchInclusionTestType());
	else if (type==Full)
	  computeContactRegionFull(region_id,grasp_->getFinger(region_id)->getWrenchInclusionTestType());
        else
	  std::cout<<"Error: Unknown ICR type - cannot compute ICR!"<<std::endl;
      }
#endif

    icr_computed_=true;
  }
  //--------------------------------------------------------------------
  bool IndependentContactRegions::writeICR(const std::string& filepath, const char* mode)const
  {
    //write the computed ICR
    if (!icr_computed_)
      {
	std::cout<<"Error in bool IndependentContactRegions::writeICR(const std::string& filepath)const - No ICR computed, cannot write to file"<<std::endl;
	return false;
      }

    FILE* f=fopen (filepath.c_str(),mode);

    if(f==NULL)
      {
	std::cout<<"Error in bool IndependentContactRegions::writeICR(const std::string& filepath)const - Could not open file "<<filepath<<". Exiting ..."<<std::endl;
	exit(0);
      }

    for (uint n=0; n<num_contact_regions_;n++)
      {
	for (uint p=0; p<contact_regions_[n]->size(); p++)
	  fprintf(f, "% d",contact_regions_[n]->at(p)->patch_ids_.front()+1); //write the centerpoint of corresponding patch, Indexing starts from 1 opposed to 0

	fprintf(f, "\n");
      }
    fclose (f);
  }
  //--------------------------------------------------------------------
  bool IndependentContactRegions::icrComputed()const{return icr_computed_;}
  //--------------------------------------------------------------------
  ContactRegion const* IndependentContactRegions::getContactRegion(uint id)const{return contact_regions_.at(id);}
  //--------------------------------------------------------------------
  uint IndependentContactRegions::getNumContactRegions()const{return num_contact_regions_;}
  //--------------------------------------------------------------------
  const SearchZonesPtr IndependentContactRegions::getSearchZones()const{return search_zones_;}
  //--------------------------------------------------------------------
  void IndependentContactRegions::setSearchZones(SearchZonesPtr sz_in)
  {
    clear();
    search_zones_ = sz_in;
  }
  //--------------------------------------------------------------------
  void IndependentContactRegions::setGrasp(GraspPtr g_in) 
  {
    clear();
    grasp_ = g_in;
  }
  //--------------------------------------------------------------------
  bool IndependentContactRegions::hasInitializedGrasp() 
  {
    return (grasp_ != NULL && grasp_->isInitialized());
  }
  //--------------------------------------------------------------------
  uint IndependentContactRegions::getNumICRPoints()const
  {
    assert(icr_computed_);
    uint n_P=0;
    for (uint n=0; n<num_contact_regions_;n++)
      for (uint p=0; p<contact_regions_[n]->size(); p++)
	n_P++;

    return n_P;
  }
  //--------------------------------------------------------------------
  void IndependentContactRegions::computeContactRegionFull(uint region_id, const WrenchInclusionTestType wrench_inclusion_test_type)
  {
    search_zones_->resetPrimitiveSearchZones(region_id); //Make sure the member ICR::PrimitiveSearchZone::satisfied_wc_ids_ is empty
    VectorXui explored_cp=VectorXui::Zero(grasp_->getParentObj()->getNumCp(),1); //set all entries to NOT_EXPLORED
    contact_regions_[region_id]->reserve(grasp_->getParentObj()->getNumCp());
 
    PatchListPtr patches=grasp_->getFinger(region_id)->getPatches();

    for (unsigned int i=0; i<patches->size(); i++)
      if(searchZoneInclusionTest(region_id,(*patches)[i],wrench_inclusion_test_type))
	{
	  contact_regions_[region_id]->push_back(const_cast<Patch*>((*patches)[i]));
          explored_cp((*patches)[i]->patch_ids_.front())=EXPLORED_QUALIFIED; //first point in the patch_id_ list is the centerpoint (hopefully ...)
	}
      else
	explored_cp((*patches)[i]->patch_ids_.front())=EXPLORED_UNQUALIFIED;
  }
  //-------------------------------------------------------------------- 
}//namespace ICR
