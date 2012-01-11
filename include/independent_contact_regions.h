#ifndef independent_contact_regions_h___
#define independent_contact_regions_h___

#include <iostream>
#include "utilities.h"

namespace ICR
{
//--------------------------------------------------------------------
//--------------------------------------------------------------------
/*! 
 *  \brief Holds shared pointers to the prototype ICR::Grasp and previously computed
 *  ICR::SearchZones; checks the contact points on the target object's surface for inclusion in the
 *  independent regions
 */
class IndependentContactRegions
{

 private:

  SearchZonesPtr search_zones_;
  GraspPtr grasp_;
  bool icr_computed_;
  std::vector<ContactRegion*> contact_regions_;
  uint num_contact_regions_;
/*! 
 *  Returns true if at least one of the primitive wrenches in wc is contained in the exterior
 *  half-space of each hyperplane defined by prim_sz; Here, an exterior half-space of a hyperplane
 *  is the half-space which does not contain the origin;
 */
  bool primitiveSearchZoneInclusionTest(PrimitiveSearchZone* prim_sz,WrenchCone const* wc)const;
  bool searchZoneInclusionTest(uint region_id,Patch const* patch)const;
  void computeContactRegion(uint region_id);
  void clear();

  IndependentContactRegions();

 public:

  IndependentContactRegions(const SearchZonesPtr search_zones,const GraspPtr grasp);
  IndependentContactRegions(IndependentContactRegions const& src);
  IndependentContactRegions& operator=(IndependentContactRegions const& src);
  friend std::ostream& operator<<(std::ostream& stream, IndependentContactRegions const& icr);
  ~IndependentContactRegions();

  void computeICR();
  bool icrComputed()const;
  ContactRegion const* getContactRegion(uint id)const;
  uint getNumContactRegions()const;
  const SearchZonesPtr getSearchZones()const;
  const GraspPtr getGrasp()const;
};
//--------------------------------------------------------------------
//--------------------------------------------------------------------
}//;namespace ICR
#endif
