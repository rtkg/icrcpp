#include "../include/search_zones.h"
#include "../include/debug.h"
#include "../include/grasp.h"
#include "assert.h"
#include <libqhullcpp/QhullVertexSet.h>


#include <iomanip>
#include <sys/time.h>
#include <time.h>

namespace ICR
{
#define KAPPA 100000
 
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
  //  //-------------------------------------------------------------------
  // void SearchZones::extractPhList(std::vector<Eigen::MatrixXd*>& Ph_list)const
  // {
  //   //Ph_list stores, for each of the h facets of the prototype GWS's convex hull, the vertices spanning the h-th facet

  //   uint H=grasp_->getGWS()->getNumFacets();
  //   uint K=grasp_->getGWS()->getDimension();

  //   Ph_list.resize(H);

  //   facetT* curr_f=grasp_->getGWS()->getConvexHull()->beginFacet().getFacetT();
  //   for (uint h=0; h<H; h++)
  //     {
  // 	setT* vertices=curr_f->vertices;
  // 	uint V=(uint)qh_setsize(vertices);
  // 	Eigen::MatrixXd* Ph(new  Eigen::MatrixXd(K,V));

  // 	//put the vertices of the current facet in the columns of the matrix Ph
  // 	for(uint v=0;v<V ;v++)
  // 	  {
  // 	    coordT* p =((vertexT*)vertices->e[v].p)->point;
  //           for (uint k=0;k<K;k++) //Really shouldn't be doint that in a loop...
  // 	      (*Ph)(k,v)=p[k];
  // 	  }

  // 	Ph_list[h]=Ph;
  // 	curr_f=curr_f->next;
  //     }
  // }
  // //-------------------------------------------------------------------
  // void SearchZones::mapFacetToFinger(const orgQhull::QhullFacet& facet,IndexList & finger_ids)const
  // {
  //   finger_ids.clear();

  //   for (uint v=0; v<facet.vertices().size(); v++)
  //     {
  // 	uint v_id=(uint)facet.vertices()[v].point().id();

  // 	//determine to which finger the current vertex belongs
  // 	for(uint n=0; n < grasp_->getNumFingers(); n++)
  // 	  if(v_id <= map_vertex2finger_(n))
  // 	    {
  // 	      finger_ids.push_back(n);
  // 	      break;
  // 	    }
  //     }
  //   finger_ids.sort();
  //   finger_ids.unique();
  // }
  // //-------------------------------------------------------------------
  // void SearchZones::extractWhiList(std::vector< ContactRegion * > const & conditioning_patches,std::vector<std::vector<Eigen::MatrixXd*> >& Wh_i_list)const
  // {
  //   //Wh_i_list stores, for each of the h facets of the prototype GWS's convex hull, the wrench cones of the conditioning patches belonging to those fingers of the prototype grasp whose wrenches span the h-th facet

  //   uint H=grasp_->getGWS()->num_facets_;
  //   Wh_i_list.resize(H);
  //   IndexList finger_ids;
  //   orgQhull::QhullFacet facet=grasp_->getGWS()->conv_hull_.beginFacet();

  //   for (uint h=0; h<H; h++)
  //     {
  // 	mapFacetToFinger(facet,finger_ids);
  // 	std::vector<Eigen::MatrixXd*> Wh_i;
  // 	for (IndexList::const_iterator f_it = finger_ids.begin(), end = finger_ids.end(); f_it != end; ++f_it)
  // 	  for (uint j=0; j<conditioning_patches[*f_it]->size();j++)
  // 	    {
  // 	      Patch* patch=conditioning_patches[*f_it]->at(j);
  // 	      for (IndexList::const_iterator p_it = patch->patch_ids_.begin(), end = patch->patch_ids_.end(); p_it != end; ++p_it)
  // 		Wh_i.push_back(new Eigen::MatrixXd(*(grasp_->getFinger(*f_it)->getOWS()->getWrenchCone(*p_it)->getWrenches())));

  // 	    }

  // 	Wh_i_list[h]=Wh_i;
  // 	facet=facet.next();
  //     }
  // }
  // //-------------------------------------------------------------------
  // void SearchZones::computeConditionedHyperplanes(std::vector< ContactRegion * > const & conditioning_patches)
  // {
  //   assert(grasp_->getGWS()->convHullComputed());
  //   uint H=grasp_->getGWS()->num_facets_;
  //   std::vector<std::vector<Eigen::MatrixXd* > > Wh_i_list(H);
  //   std::vector<Eigen::MatrixXd*> Ph_list(H);
  //   CasADi::IpoptSolver nlp_solver;


  //   struct timeval start, end;
  //   double c_time;

  //   if(tws_->getWrenchSpaceType() == Spherical)
  //     {
  // 	std::cout<<"Error in SearchZones::computeConditionedHyperplanes() - Spherical Task Wrench Space based search zone computation not implemented yet!"<<std::endl;
  //     }
  //   else if(tws_->getWrenchSpaceType() == Discrete)
  //     {
  // 	std::vector<Eigen::MatrixXd*> Ph_list;
  //       extractPhList(Ph_list); //ALLOCATES MEMORY!!! Should be done with shared pointers ...
  //       extractWhiList(conditioning_patches,Wh_i_list);//ALLOCATES MEMORY!!! Should be done with shared pointers ...

  // 	createDiscreteTWSNLPSolver(Ph_list,Wh_i_list,nlp_solver);
  // 	std::cout<<"TRYING TO SOLVE..."<<std::endl;
  // 	gettimeofday(&start,0);
  // 	nlp_solver.evaluate();

  // 	// std::cout << std::setw(30) << "Objective: " << nlp_solver.output("f").getDescription() << std::endl;
  // 	// std::cout << std::setw(30) << "Primal solution: " << nlp_solver.output("x").getDescription() << std::endl;

  // 	gettimeofday(&end,0);
  // 	c_time = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
  // 	std::cout<<"Computation time SOLVER: "<<c_time<<" s"<<std::endl;
  // 	setTWSNLPHyperplanes(nlp_solver,Wh_i_list);
  //     }
  //   else
  //     {
  // 	std::cout<<"Error in SearchZones::computeConditionedHyperplanes() - invalid Task Wrench Space type!"<<std::endl;
  // 	return;
  //     }

  //   //clean up
  //   for (uint h=0; h<H; h++)
  //     {
  // 	for (uint i=0; i<Wh_i_list[h].size();i++)
  // 	  delete Wh_i_list[h][i];  

  //       delete Ph_list[h];
  //     }
  // }
  //-------------------------------------------------------------------
  void SearchZones::addShiftedPrimitiveSearchZone(uint finger_id,vertexT const* curr_vtx)
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
  // void SearchZones::computeConditionedSearchZones(std::vector< ContactRegion * > const & conditioning_patches)
  // {
  //   //Make sure the TWS is valid
  //   assert(tws_.get() != NULL);
  //   assert(tws_->getDimension() == grasp_->getGWS()->getDimension());
  //   assert(grasp_->getGWS()->getOcInsphereRadius() >= tws_->getOcInsphereRadius());

  //   //Make sure the grasp is valid
  //   assert(grasp_->getGWS()->containsOrigin());

  //   //Make sure the conditioning patches are valid
  //   assert(conditioning_patches.size()==grasp_->getNumFingers());

  //   initializeSearchZones();
  //   computeConditionedHyperplanes(conditioning_patches);
  //   computePrimitiveSearchZones();
 
  //   search_zones_computed_=true;
  // }
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
		addShiftedPrimitiveSearchZone(finger_id,curr_vtx.getVertexT());
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
