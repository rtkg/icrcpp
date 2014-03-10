#ifndef icr_h___
#define icr_h___

/*! \mainpage The icrcpp library
 *
 *Author: Robert Krug robert.krug@oru.se
 *
 *The library provides computation of Independent Contact Regions (ICR) on discretized object models
 *represented by triangular meshes provided in the Wavefront .obj file format. Single- and
 *Multi-Point contact models are possible. The icrcpp library is available under the GNU General
 *Public License.
 *
 * \section concept_sec General Concept
 *
 * The central object of the library is the Grasp, which is tied to a previously loaded target
 * object. A Grasp contains a list of Finger-pointers. A Finger holds a pointer to an Object Wrench
 * Space and a pointer to a list of patches. OWS and patches are precomputed when the Finger's
 * parent grasp is initialized. A separate OWS is computed for each finger with a different contact
 * model, a separate list of patches is computed for each finger with a different inclusion
 * rule. Fingers with the same contact model/inclusion rule hold pointers to the same OWS/patch
 * list. The inclusion rule determines which vertices of the target object's mesh are eligible for
 * inclusion in a patch centered around a given center-point. In case of the single-point contact
 * model, each patch only contains its respective center-point.The class SearchZones requires a
 * force-closure prototype grasp for creation. After computing the search zones for the given
 * prototype grasp, these search zones and the grasp are used in the constructor for objects of type
 * IndependentContactRegions. All objects and functions share the namespace ICR.
 *
 *The basic program flow is as follows: 
 *- Create an ICR::ObjectLoader 
 *- Load a target object in form of a Wavefront .obj file 
 *- Create a ICR::Grasp
 *- Initialize the Grasp with a set of proper ICR::FingerParameters, center-point indices and the loaded object
 *- Create ICR::SearchZones utilizing the prototype grasp
 *- Compute the search zones given the parameters of a task wrench space
 *- Create an ICR::IndependentContactRegions object utilizing the previously created search zones and grasp
 *- Compute the independent contact regions
 *
 * See also the provided example in /icrcpp/examples/src/example.cpp.
 *
 * \section todo_sec ToDo
 *
 *
 *- In the ICR::ObjectLoader class, implement the callbacks allowing to read .obj files with texture
 *   information. Right now, program execution is stopped by an attempt to load an .obj file
 *   containing textures.
 *
 *- Maybe ICR::Grasp::computeGWS should not be called in ICR::Grasp::init, but
 *   separately. Transferring existing ICR::SearchZones to a novel object also requires creation of
 *   a new Grasp having this novel object as parent, because new Object Wrench Spaces and Patch lists are
 *   needed. However, calculating the GWS for this grasp is superfluous and could be avoided by
 *   disconnecting the GWS calculation from the initialization of the grasp and thus increase
 *   computational efficiency.
 *
 *- Generalize the ICR::InclusionRule which is responsible for generating patches on the target
 *   object's surface. Right now, the inclusion rule is governed only by the scalar
 *   ICR::InclusionRule::rule_parameter_ which gives the radius of a sphere centered at a given
 *   center-point on the target object. One possibility would be to change this rule parameter to a
 *   3x3 Matrix parametrizing an ellipsoid and add an additional rotation matrix to give the
 *   orientation of this ellipsoid w.r.t the target object. An even more sophisticated possibility
 *   would be to formulate the rule parameter as a mesh describing the respective finger tip with a
 *   given offset. Maybe also change the name to ICR::PatchRule or something like that to avoid
 *   confusion with ICR::IndependentContactRegions::searchZoneInclusionTest and
 *   ICR::IndependentContactRegions::primitiveSearchZoneInclusionTest.
 *
 *- Consistently use shared pointers for all non-elemental class members and Add clone methods to
 *   all classes to allow for deep copies
 *
 *- Make a central (global) Qhull server and remove Qhull as class member from the
 *   DiscreteWrenchSpace class since Qhull allows only one instance to be alive at a time (this would
 *   also allow to remove the DiscreteTaskWrenchSpace class which only exists for this reason)
 *
 *- Remove the wrench-space 6D-assumption which is hardcoded in places
 *
 *- Assert that a discrete TWS contains the origin in SearchZones::computeShiftedHyperplanes() 
 *
 *- Implement proper error handling opposed to the simple screen prints used right now
 *
 *- Remove the self-defined Eigen2Array and Array2Eigen functions defined in utilities.h/cpp and
 *  replace them with Eigen::Map
 *
 *
 *- Need to run qhull with the Qx option, which results in coplanar facets, otherwise there are
 *    undefined facets - should find some workaround for this to reduce the unnecessary large number
 *    of hyperplanes (edit: also the Q0 options seems to work which doesn't produce redundant planes...)
 *
 *- Implement prioritized QP search zone computation for sphericas task wrench spaces
 *
 *- Try a version of the QP-based search zone computation as one 'big' QP for all hyperplanes - should be faster
 *
 *- Also write the loaded mesh to a file if the according debug flag is set
 *
 *- Change the Discrete Wrench Space constructors to accept shared pointers instead of double*
 *
 *- Still a small memory leak in
 *    ICR::IndependentContactRegions::convexCombinationSearchZoneInclusionTest(PrimitiveSearchZone*
 *   prim_sz,WrenchCone const* wc)const - seemst to be produced by Gurobi's GRBEnv constructor
 *
 *- Add a filter method to the ICR::IndependentContactRegions class to remove duplicated points from regions 
 *
 *- Object parser always returns false ...?
 *
 *- Should make sure that the given TWS contains the origin in  VectorXui generateRandomGrasp(TargetObjectPtr object,const FParamList f_parameters, WrenchSpacePtr const tws) 
 *
 *- Should add a public method to allow dividing the OWS by lambda not only via a preprocessor flag
 *
 *- Wrench space inclusion test type should be set in the ICR::IndependentContactRegions class, not via finger parameters in the fingers
 *
 *- Pinpoint the cause of the sporadic Gurobi exceptions in IndependentContactRegions::convexCombinationSearchZoneInclusionTest(PrimitiveSearchZone* prim_sz,WrenchCone const* wc)const
 *
 *- Add a flipNormals() method to the ICR::TargetObject class to invert vertex normal directions
 *
 *- Make the Gurobi environment a member of the ICR::IndependentContactRegions class and protect it with a mutex
 */

#include "utilities.h"
#include "contact_model.h"
#include "contact_point.h"
#include "finger.h"
#include "grasp.h"
#include "independent_contact_regions.h"
#include "limit_surface.h"
#include "object_loader.h"
#include "ows.h"
#include "search_zones.h"
#include "target_object.h"
#include "wrench_cone.h"
#include "wrench_space.h"

#endif
