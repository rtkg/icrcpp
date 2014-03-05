#include "../include/search_zones.h"
#include "../include/debug.h"
#include "../include/grasp.h"
#include "assert.h"
#include <libqhullcpp/QhullVertexSet.h>
#include <gurobi_c++.h>

// #include <sys/time.h>
// #include <time.h>

namespace ICR
{
  //-------------------------------------------------------------------
  PrimitiveSearchZone::PrimitiveSearchZone(){hyperplane_ids_.clear();}
  //-------------------------------------------------------------------
  std::ostream& operator<<(std::ostream& stream, PrimitiveSearchZone const& psz)
  {
    stream <<'\n'<<"PRIMITIVE SEARCH ZONE: "<<'\n'
	   <<"Hyperplane id's: ";
       
    for(uint i=0;i<psz.hyperplane_ids_.size();i++)
      stream<<psz.hyperplane_ids_[i]<<" ";

    stream<<'\n'<<'\n';

    return stream;
  }
  //-------------------------------------------------------------------
  SearchZones::SearchZones() : num_search_zones_(0),search_zones_computed_(false)
  {
    search_zones_.clear();
  }
  //-------------------------------------------------------------------
  SearchZones::SearchZones(const GraspPtr grasp) : grasp_(grasp), num_search_zones_(0),search_zones_computed_(false)
  {
    assert(grasp_->isInitialized());
    assert(grasp_->getGWS()->containsOrigin());
    search_zones_.clear();
  }
  //-------------------------------------------------------------------
  SearchZones::SearchZones(SearchZones const& src) :  grasp_(src.grasp_), tws_(src.tws_), search_zones_(src.search_zones_),
						      num_search_zones_(src.num_search_zones_), search_zones_computed_(src.search_zones_computed_),map_vertex2finger_(src.map_vertex2finger_),
                                                      hyperplane_normals_(src.hyperplane_normals_),hyperplane_offsets_(src.hyperplane_offsets_){}
  //-------------------------------------------------------------------
  SearchZones& SearchZones::operator=(SearchZones const& src)
  {
    if (this !=&src)
      {
	grasp_=src.grasp_;
	tws_=src.tws_;
	search_zones_=src.search_zones_;
	num_search_zones_=src.num_search_zones_;
	search_zones_computed_=src.search_zones_computed_;
	map_vertex2finger_=src.map_vertex2finger_;
	hyperplane_normals_=src.hyperplane_normals_;
	hyperplane_offsets_=src.hyperplane_offsets_;
      }
    return *this;
  }
  //-------------------------------------------------------------------
  std::ostream& operator<<(std::ostream& stream,SearchZones const& sz)
  {
    stream <<'\n'<<"SEARCH ZONES: "<<'\n'
	   <<"Number of patch search zones: " << sz.num_search_zones_<<'\n'
	   <<"Search zones computed: "<<std::boolalpha<<sz.search_zones_computed_<<'\n'
	   <<"Vertex-finger map: "<<sz.map_vertex2finger_<<'\n';

    return stream;
  }
  //-------------------------------------------------------------------
  SearchZones::~SearchZones(){clear();}
  //-------------------------------------------------------------------
  SearchZone const* SearchZones::getSearchZone(uint finger_id)const{return search_zones_.at(finger_id);}
  //-------------------------------------------------------------------
  uint SearchZones::getNumSearchZones()const{return num_search_zones_;}
  //-------------------------------------------------------------------
  bool SearchZones::searchZonesComputed()const{return search_zones_computed_;}
  //-------------------------------------------------------------------
  void SearchZones::computeShiftedHyperplanes()
  {
    uint K=grasp_->getGWS()->getDimension();
    uint H=grasp_->getGWS()->num_facets_;
    hyperplane_normals_.resize(H,K);
    hyperplane_offsets_.resize(H);
    facetT* curr_f=grasp_->getGWS()->conv_hull_.beginFacet().getFacetT();
 
    if(tws_->getWrenchSpaceType() == Spherical)
      {
	double offset=tws_->getOcInsphereRadius();
	hyperplane_offsets_.fill(offset); 
     
	for(uint h=0; h < H;h++)
	  {
	    assert(offset <= (-curr_f->offset)); //Make sure that GWSinit contains the TWS and the hyperplanes are shifted inwards
	    //The direction of the normal is switched to inward-pointing in order to be consistent with the positive offset

	    //  hyperplane_normals_.row(h)=(-Eigen::Map<Eigen::Matrix<double,1,6> >(curr_f->normal));
	    hyperplane_normals_.row(h)=-Eigen::Map<Eigen::RowVectorXd>(curr_f->normal,K);
	    curr_f=curr_f->next;
	  }
      }

    else if(tws_->getWrenchSpaceType() == Discrete)
      {
	//Should assert that the TWS contains the origin ...
	Eigen::MatrixXd TW(K,dynamic_cast<DiscreteTaskWrenchSpace*>(tws_.get())->getNumWrenches());
        doubleArrayToEigenMatrix(dynamic_cast<DiscreteTaskWrenchSpace*>(tws_.get())->getWrenches().get(),TW);

      	for(uint h=0; h < H;h++)
	  {
	    Eigen::RowVectorXd normal=-Eigen::Map<Eigen::RowVectorXd>(curr_f->normal,K);
            double e = (-curr_f->offset);
	    assert(e > 0); //just in case ...

	    double e_shift=-(normal*TW).minCoeff();
	    //Should provide some tolerances for the checks below ...
	    if (e_shift < 0)
	      {
		std::cout<<"Error in SearchZones::computeShiftedHyperplanes() - the TWS does not contain the origin! Exiting ..."<<std::endl;
		exit(0);
	      }
	    if (e_shift > e)
	      {
		std::cout<<"Error in SearchZones::computeShiftedHyperplanes() - the prototype grasp's GWS does not contain the TWS! Exiting ..."<<std::endl;
		exit(0);
              }



	    hyperplane_offsets_(h)=e_shift;
	    hyperplane_normals_.row(h)=normal;
	    // std::cout<<hyperplane_normals_.row(h)<<" "<<hyperplane_offsets_(h)<<std::endl;
	    curr_f=curr_f->next;
          }

      }
    else
      std::cout<<"Error in SearchZones::computeShiftedHyperplanes() - invalid Task Wrench Space type!"<<std::endl;
  }
  //-------------------------------------------------------------------
  void SearchZones::mapFacetToFingers(const orgQhull::QhullFacet& facet,IndexList & finger_ids)const
  {
    finger_ids.clear();

    for (uint v=0; v<facet.vertices().size(); v++)
      {
  	uint v_id=(uint)facet.vertices()[v].point().id();

  	//determine to which finger the current vertex belongs
  	for(uint n=0; n < grasp_->getNumFingers(); n++)
  	  if(v_id <= map_vertex2finger_(n))
  	    {
  	      finger_ids.push_back(n);
  	      break;
  	    }
      }
  }
  //-------------------------------------------------------------------
  void SearchZones::addPrimitiveSearchZone(uint finger_id,vertexT const* curr_vtx)
  {
    setT *neighbor_facets=curr_vtx->neighbors;
    PrimitiveSearchZone* p_search_zone = new PrimitiveSearchZone();
    p_search_zone->hyperplane_ids_.reserve(qh_setsize(neighbor_facets));

    for(uint i=0; i < (uint)qh_setsize(neighbor_facets); i++)
      p_search_zone->hyperplane_ids_.push_back(((facetT*)neighbor_facets->e[i].p)->id);

    search_zones_[finger_id]->push_back(p_search_zone);
  }
  //-------------------------------------------------------------------
  void SearchZones::clear()
  {
    for(uint finger_id=0;finger_id < num_search_zones_;finger_id++)
      {
	for(uint prim_sz_id=0; prim_sz_id < search_zones_[finger_id]->size(); prim_sz_id++)
	  delete (*search_zones_[finger_id])[prim_sz_id];            

	delete search_zones_[finger_id];
      }
    search_zones_.clear();
    num_search_zones_ = 0;
    search_zones_computed_ = false;

  }
  //-------------------------------------------------------------------
  void SearchZones::initializeSearchZones()
  {
    if(search_zones_computed_) //Clear up possible previously computed search zones
      clear();

    num_search_zones_=grasp_->getNumFingers();
    search_zones_.reserve(num_search_zones_);
    map_vertex2finger_.resize(1,num_search_zones_);

    for(uint finger_id=0; finger_id < num_search_zones_; finger_id++)
      {
	search_zones_.push_back(new SearchZone());
	//reserve for the maximum number of primitive search zones for each finger (i.e all wrenches are vertexes of the gws and span a primitive search zone)
	search_zones_[finger_id]->reserve(grasp_->getFinger(finger_id)->getCenterPointPatch()->patch_ids_.size()
					  *grasp_->getFinger(finger_id)->getContactModel()->getLimitSurface()->getNumPrimitiveWrenches());

	//the elements of map_vertex2finger_ describe the relation between a vertex-index from the grasp wrench space to the finger this vertex belongs.
	//I.e. if the value of vertex id <= than map_vertex2finger_(0), then the vertex is a wrench belonging to the first finger, if vertex id <= than 
	//map_vertex2finger_(1) the vertex belongs to the second finger and so on. E.g, a 3-fingered grasp comprising patches with 2 contact points each and 10 
	//primitive wrenches in each contact point would have map_vertex2finger_=[19 39 59] 
	map_vertex2finger_(finger_id)=grasp_->getFinger(finger_id)->getCenterPointPatch()->patch_ids_.size()
	  *grasp_->getFinger(finger_id)->getContactModel()->getLimitSurface()->getNumPrimitiveWrenches()-1;
	if (finger_id > 0)
	  map_vertex2finger_(finger_id)=map_vertex2finger_(finger_id)+map_vertex2finger_(finger_id-1)+1;
      }
  }
  //--------------------------------------------------------------------
  void SearchZones::computeShiftedSearchZones()
  {
    //Make sure the TWS is valid
    assert(tws_.get() != NULL);
    assert(tws_->getDimension() == grasp_->getGWS()->getDimension());
    assert(grasp_->getGWS()->getOcInsphereRadius() >= tws_->getOcInsphereRadius());

    //Make sure the grasp is valid
    assert(grasp_->getGWS()->containsOrigin());

    initializeSearchZones();
    computeShiftedHyperplanes();
    computePrimitiveSearchZones();
 
    search_zones_computed_=true;
  }
  //--------------------------------------------------------------------
  void SearchZones::computePrioritizedSearchZones(uint finger_id)
  {
    //Make sure the TWS is valid
    assert(tws_.get() != NULL);
    assert(tws_->getDimension() == grasp_->getGWS()->getDimension());
    assert(grasp_->getGWS()->getOcInsphereRadius() >= tws_->getOcInsphereRadius());

    //Make sure the grasp is valid
    assert(grasp_->getGWS()->containsOrigin());

    initializeSearchZones();
    computePrioritizedHyperplanes(finger_id);
    computePrimitiveSearchZones();
 
    search_zones_computed_=true;
  }
  //-------------------------------------------------------------------
  void SearchZones::computePrioritizedHyperplanes(uint finger_id)
  {
    uint K=grasp_->getGWS()->getDimension();
    uint H=grasp_->getGWS()->num_facets_;
    hyperplane_normals_.resize(H,K);
    hyperplane_offsets_.resize(H);
    orgQhull::QhullFacet curr_f=grasp_->getGWS()->conv_hull_.beginFacet();
 
    if(tws_->getWrenchSpaceType() == Spherical)
      {
	std::cout<<"Error in SearchZones::computePrioritizedHyperplanes(uint finger_id) - method is not yet implemented for spherical task wrench spaces! Exiting ..."<<std::endl;
        exit(0);
      }
    else if(tws_->getWrenchSpaceType() == Discrete)
      {

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag,0);
	env.set(GRB_IntParam_Presolve,0);
       	//   env.set(GRB_IntParam_Method,2);

       	//Should assert that the TWS contains the origin ...
       	Eigen::MatrixXd TW(K,dynamic_cast<DiscreteTaskWrenchSpace*>(tws_.get())->getNumWrenches());
	doubleArrayToEigenMatrix(dynamic_cast<DiscreteTaskWrenchSpace*>(tws_.get())->getWrenches().get(),TW);

	for(uint h=0; h < H;h++)
	  {
	    Eigen::RowVectorXd normal=-Eigen::Map<Eigen::RowVectorXd>(curr_f.getFacetT()->normal,K);
	    double e = -(curr_f.getFacetT()->offset);
	    assert(e > 0); //just in case ...
	 
	    //if the finger to be prioritized does not participate in spanning the h-th facet, directly transfer hyperplane h from the GWS to the search zones and continue to the next facet
	    IndexList finger_ids;
	    mapFacetToFingers(curr_f,finger_ids);
	    IndexList::iterator it = std::find(finger_ids.begin(), finger_ids.end(), finger_id);
         
	    if (it == finger_ids.end())
	      {
	       //WITH ADDITIONAL SHIFTING
	       hyperplane_normals_.row(h)=normal;
               double e_shift=-(normal*TW).minCoeff();
               hyperplane_offsets_(h)=e_shift;
               
	       //WITHOUT ADDITIONAL SHIFTING
	       //hyperplane_offsets_(h)=e;
	    	continue;
	      } 
            
	    // std::cout<<"vertex2finger list size: "<<finger_ids.size()<<std::endl;
	    uint l=0;
	    Eigen::MatrixXd R,P;
	    for (IndexList::const_iterator it = finger_ids.begin(), end = finger_ids.end(); it != end; ++it)
	      {
		//	std::cout<<(*it)<<" ";
		Eigen::VectorXd vertex=Eigen::Map<Eigen::RowVectorXd>(curr_f.vertices()[l].getVertexT()->point,K);
		if ((*it)==finger_id)
		  {
		    //the l-th vertex of the current facet belongs to finger 'finger_id'
                    R.conservativeResize(K,R.cols()+1);
		    R.rightCols(1)=vertex;
		  }
		else
		  {
		    //the l-th vertex of the current facet doesn't belong to finger 'finger_id'
		    P.conservativeResize(K,P.cols()+1);
 		    P.rightCols(1)=vertex;
		  }

		l++;
	      }

	    //SOLVE THE TILT PROBLEM

	    GRBModel qp = GRBModel(env);
            double* lb_x = new double[K+1]; std::fill_n(lb_x, K+1,std::numeric_limits<double>::infinity()*(-1)); lb_x[K]=0; //lower bound on eh is only given for numerical reasons - solver f***s up otherwise ...
	    GRBVar* x =qp.addVars(lb_x,NULL,NULL,NULL,NULL,K+1);
	    //just some convenience pointers ...
            GRBVar* nh=x;
            GRBVar* eh=(x+K);
	    qp.update();
	 
	    //Objective
            GRBQuadExpr obj = 0.0;
            for (uint k=0; k<K;k++)
	      obj += 0.25*nh[k]*nh[k];

            qp.setObjective(obj,GRB_MINIMIZE);

	    //Equality constriants
	    GRBLinExpr* eq_c_expr=new GRBLinExpr[P.cols()];
            char* eq_c_senses = new char[P.cols()];  std::fill_n(eq_c_senses,P.cols(),GRB_EQUAL); 
            double* eq_c_rhs = new double[P.cols()]; std::fill_n(eq_c_rhs,P.cols(),0.0);           
            double** p_i = new double*[P.cols()];
	    for (uint i=0; i<P.cols();i++)
	      {
		eigenMatrixToDoubleArray(P.col(i),p_i[i]);
		eq_c_expr[i] = 0.0;
		for (uint k=0; k<K;k++)
		  eq_c_expr[i] += p_i[i][k]*nh[k];

		eq_c_expr[i]+=eh[0];             
	      }
	    GRBConstr* eq_c= qp.addConstrs(eq_c_expr,eq_c_senses,eq_c_rhs,NULL,P.cols());

	    //Inequality constriants on R
   	    GRBLinExpr* ineq_c_R_expr=new GRBLinExpr[R.cols()];
            char* ineq_c_R_senses = new char[R.cols()];  std::fill_n(ineq_c_R_senses,R.cols(),GRB_LESS_EQUAL); 
            double* ineq_c_R_rhs = new double[P.cols()]; std::fill_n(ineq_c_R_rhs,R.cols(),-1.0);           
            double** r_i = new double*[R.cols()];
   	    for (uint i=0; i<R.cols();i++)
   	      {
   		eigenMatrixToDoubleArray(R.col(i),r_i[i]);
   		ineq_c_R_expr[i] = 0.0;
   		for (uint k=0; k<K;k++)
   		  ineq_c_R_expr[i] += r_i[i][k]*nh[k];

   		ineq_c_R_expr[i]+=eh[0];             
   	      }
   	    GRBConstr* ineq_c_R= qp.addConstrs(ineq_c_R_expr,ineq_c_R_senses,ineq_c_R_rhs,NULL,R.cols());

	    //Inequality constriants on TW
   	    GRBLinExpr* ineq_c_T_expr=new GRBLinExpr[TW.cols()];
            char* ineq_c_T_senses = new char[TW.cols()];  std::fill_n(ineq_c_T_senses,TW.cols(),GRB_LESS_EQUAL); 
            double* ineq_c_T_rhs = new double[TW.cols()]; std::fill_n(ineq_c_T_rhs,TW.cols(),0.0);           
            double** t_i = new double*[TW.cols()];
   	    for (uint i=0; i<TW.cols();i++)
   	      {
   		eigenMatrixToDoubleArray(R.col(i),t_i[i]);
   		ineq_c_T_expr[i] = 0.0;
   		for (uint k=0; k<K;k++)
   		  ineq_c_T_expr[i] += t_i[i][k]*nh[k];

   		ineq_c_T_expr[i]+=eh[0];                         
   	      }
   	    GRBConstr* ineq_c_T= qp.addConstrs(ineq_c_T_expr,ineq_c_T_senses,ineq_c_T_rhs,NULL,TW.cols());

            qp.update();
	    qp.optimize();
	    int status = qp.get(GRB_IntAttr_Status);

	    if (status != GRB_OPTIMAL) 
	      {
		std::cout<<"Error in SearchZones::computePrioritizedHyperplanes(uint finger_id) - QP solution not found! Exiting ..."<<std::endl;
		exit(0);
	      }

	    double* x_opt=qp.get(GRB_DoubleAttr_X,x,K+1);

	    Eigen::MatrixXd nh_tilt(1,K);
	    doubleArrayToEigenMatrix(x_opt,nh_tilt);
	    double de=1/nh_tilt.norm(); nh_tilt.normalize();

	        hyperplane_normals_.row(h)=nh_tilt;
             	//WITH ADDITIONAL SHIFTING
               double e_shift=-(nh_tilt*TW).minCoeff();
               hyperplane_offsets_(h)=e_shift;

	       //WITHOUT ADDITIONAL SHIFTING
	       //hyperplane_offsets_(h)=x_opt[K]*de;

 
	    // std::cout<<"solution: "<<hyperplane_normals_.row(h)<<" "<<hyperplane_offsets_(h)<<std::endl;
	    // std::cout<<"x_opt: ";
	    // for (uint j=0;j<K+1; j++)
	    //   std::cout<<x_opt[j]<<" ";
	    // std::cout<<std::endl<<std::endl;

 	   


	    assert(hyperplane_offsets_(h) >= 0); //just to be sure ...           
	    curr_f=curr_f.next();

 	    //clean-up
	    for (uint i=0; i<TW.cols();i++)
	      delete[] t_i[i];
	    for (uint i=0; i<P.cols();i++)
	      delete[] p_i[i];
	    for (uint i=0; i<R.cols();i++)
	      delete[] r_i[i];

	    delete[] x; delete[] lb_x; delete[] x_opt; 
	    delete p_i; delete[] eq_c_senses; delete[] eq_c_rhs; delete[] eq_c; delete[] eq_c_expr;
	    delete r_i; delete[] ineq_c_R_senses; delete[] ineq_c_R_rhs; delete[] ineq_c_R; delete[] ineq_c_R_expr;
	    delete t_i; delete[] ineq_c_T_senses; delete[] ineq_c_T_rhs; delete[] ineq_c_T; delete[] ineq_c_T_expr; 
	  }
      }
    else
      std::cout<<"Error in SearchZones::computeShiftedHyperplanes() - invalid Task Wrench Space type!"<<std::endl;
      
  } 
  //--------------------------------------------------------------------
  void SearchZones::computePrimitiveSearchZones()
  { 
    //iterate through the vertices of the GWS and create the appropriate Primitive Search Zone for each 
    //vertex and push it back on the search zone of its corresponding finger

    orgQhull::QhullVertex curr_vtx=grasp_->getGWS()->conv_hull_.beginVertex();

    for (uint vtx_count=0; vtx_count < grasp_->getGWS()->num_vtx_;vtx_count++)
      {
	//determine to which finger the current vertex belongs and create an appropriate primitive search zone and push it on the list for the corresponding
	//finger
	for(uint finger_id=0; finger_id < num_search_zones_; finger_id++)
	  {
	    if((uint)curr_vtx.point().id() <= map_vertex2finger_(finger_id))
	      {
		addPrimitiveSearchZone(finger_id,curr_vtx.getVertexT());
		break;
	      }
	  }

	curr_vtx=curr_vtx.next();
      }
  }
  //--------------------------------------------------------------------
  Eigen::Matrix<double,Eigen::Dynamic,6> const* SearchZones::getHyperplaneNormals()const{return &hyperplane_normals_;}
  //--------------------------------------------------------------------
  Eigen::VectorXd const* SearchZones::getHyperplaneOffsets()const{return &hyperplane_offsets_;}
  //--------------------------------------------------------------------
  const GraspPtr SearchZones::getGrasp()const{return grasp_;}
  //--------------------------------------------------------------------
  WrenchSpacePtr SearchZones::getTaskWrenchSpace()const
  {
    assert(search_zones_computed_);
    return tws_;
  }
  //--------------------------------------------------------------------
  void SearchZones::resetPrimitiveSearchZones(uint sz_id)
  {
    for(uint psz_id=0; psz_id < search_zones_[sz_id]->size(); psz_id++)
      (*search_zones_[sz_id])[psz_id]->satisfied_wc_ids_.resize(0);  
  }
  //--------------------------------------------------------------------
  void SearchZones::setTaskWrenchSpace(WrenchSpacePtr tws)
  {
    //make sure the TWS is valid
    assert(tws.get() != NULL);

    tws_=tws;
  }
  //--------------------------------------------------------------------
  void SearchZones::resetSearchZones()
  {
    for(uint sz_id=0;sz_id < num_search_zones_; sz_id++)
      resetPrimitiveSearchZones(sz_id);
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
}//namespace ICR
