#ifndef search_zones_h___
#define search_zones_h___

#include <iostream>
#include "utilities.h"
#include "ows.h"
#include <libqhullcpp/Qhull.h>
#include <iomanip>
#include <symbolic/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <symbolic/stl_vector_tools.hpp>

namespace ICR
{
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  /*! 
   *  \brief Used in ICR::IndependentContactRegions#primitiveSearchZoneInclusionTest to determine
   *  whether a contact point associated to a wrench cone is eligible for inclusion in the independent regions
   */
  struct PrimitiveSearchZone
  {
    /*!
     * Idexes the hyperplanes in ICR::SearchZones#hyperplane_normals_ and ICR::SearchZones#hyperplane_offsets_
     */
    std::vector<uint> hyperplane_ids_;
    PrimitiveSearchZone();
    /*!
     * Helper for ICR::IndependentContactRegions#primitiveSearchZoneInclusionTest; stores the indices of already explored wrench cones;
     */
    VectorXui satisfied_wc_ids_;

    friend std::ostream& operator<<(std::ostream& stream, PrimitiveSearchZone const& psz);
  };
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  /*! 
   *  \brief Holds a std::vector of ICR::SearchZone pointers, one for each finger of the associated prototype grasp.
   */
  class SearchZones
  {
  private:

    GraspPtr grasp_;
    WrenchSpacePtr tws_;
    std::vector<SearchZone*> search_zones_;
    uint num_search_zones_;
    bool search_zones_computed_;
    RowVectorXui map_vertex2finger_;
    Eigen::Matrix<double,Eigen::Dynamic,6> hyperplane_normals_;
    Eigen::VectorXd hyperplane_offsets_;  

    void computeConditionedHyperplanes(std::vector< ContactRegion * > const & conditioning_patches);
    void computeShiftedHyperplanes();
    void initializeSearchZones();
    void addShiftedPrimitiveSearchZone(uint finger_id,vertexT const* curr_vtx);
    void resetPrimitiveSearchZones(uint sz_id);    
    void clear();
    void computePrimitiveSearchZones();
    void mapFacetToFinger(const orgQhull::QhullFacet& facet,IndexList & finger_ids)const;
    void extractPhList(std::vector<Eigen::MatrixXd*>& Ph_list)const;
    void extractWhiList(std::vector< ContactRegion * > const & conditioning_patches,std::vector<std::vector<Eigen::MatrixXd*> >& Wh_i_list)const;
    void createDiscreteTWSNLPSolver(const std::vector<Eigen::MatrixXd*>& Ph_list,const std::vector<std::vector<Eigen::MatrixXd* > >& Wh_i_list,CasADi::IpoptSolver& nlp_solver);
    void initialSolutionTWSNLP(CasADi::DMatrix& initial_solution, const std::vector<std::vector<Eigen::MatrixXd* > >& Wh_i_list);
    void setTWSNLPHyperplanes(const CasADi::IpoptSolver& nlp_solver,const std::vector<std::vector<Eigen::MatrixXd* > >& Wh_i_list);
    SearchZones();

  public:

    friend class IndependentContactRegions;

    SearchZones(const GraspPtr grasp);
    SearchZones(SearchZones const& src);
    SearchZones& operator=(SearchZones const& src);
    friend std::ostream& operator<<(std::ostream& stream,SearchZones const& sz);
    ~SearchZones();

    /*! 
     *  \brief Creates ICR::SearchZones by shifting the hyperplanes described in
     *  ICR::SearchZones#hyperplane_normals_ and ICR::SearchZones#hyperplane_offsets_ until their are
     *  tangent to a Task Wrench Space 
     */
    void computeShiftedSearchZones();
    void computeConditionedSearchZones(std::vector< ContactRegion * > const & conditioning_patches);
    const GraspPtr getGrasp()const;
    SearchZone const* getSearchZone(uint finger_id)const;
    uint getNumSearchZones()const;
    bool searchZonesComputed()const;
    void setTaskWrenchSpace(WrenchSpacePtr tws);
    /*! 
     *  \brief Empties the vectors ICR::PrimitiveSearchZone::satisfied_wc_ids_ of all primitive search
     *  zones; This is necessary when transferring search zones which were utilized to compute ICR for
     *  one object to a novel object; Note that its child function
     *  ICR::SearchZones#resetPrimitiveSearchZones is automatically called by
     *  ICR::IndependentContactRegions#computeContactRegion and thus doesn't need to be called
     *  explicitly;
     */
    void resetSearchZones();
    Eigen::Matrix<double,Eigen::Dynamic,6> const* getHyperplaneNormals()const;
    Eigen::VectorXd const*getHyperplaneOffsets()const;
    WrenchSpacePtr getTaskWrenchSpace()const;
  };
  //-------------------------------------------------------------------
  CasADi::FX generateCodeAndCompile(CasADi::FX fcn, const std::string& name);
  //-------------------------------------------------------------------
}//namespace ICR
#endif
