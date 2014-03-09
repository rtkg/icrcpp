#ifndef utilities_h___
#define utilities_h___

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <list>
#include <vector>
#include <tr1/memory>

#define EPSILON_WRENCH_CONE_ROTATION 1e-10
#define EPSILON_UNIT_NORMAL 1e-6
#define EPSILON_FORCE_CLOSURE 1e-10
#define PI 3.14159265358979323846264338327950
#define NOT_VISITED         -1
#define NOT_EXPLORED         0
#define EXPLORED_QUALIFIED   1
#define EXPLORED_UNQUALIFIED 2

namespace ICR
{
//--------------------------------------------------------------------
//--------------------------------------------------------------------

//Forward declerations
class ContactPoint;
class WrenchCone;
class OWS;
class Finger;
class PointContactModel;
class MultiPointContactModel;
class FingerParameters;
class SearchZones;
class LimitSurface;
class IndependentContactRegions;
class Grasp;
class TargetObject;
class WrenchSpace;
class SphericalWrenchSpace;
class DiscreteWrenchSpace;
class DiscreteTaskWrenchSpace;
class ObjectLoader;
struct Node;
struct PrimitiveSearchZone;
struct Patch;
struct InclusionRule;

enum ContactType {Frictionless=1, Frictional, Soft_Finger};
enum ModelType {Single_Point=1,Multi_Point};
enum RuleType {Sphere=1};
enum WrenchSpaceType {Discrete=1,Spherical};
enum ICRType {BFS=1,Full};
enum WrenchInclusionTestType {Primitive=1,Convex_Combination};

typedef unsigned int uint;
typedef Eigen::Matrix<double,6,Eigen::Dynamic> Matrix6Xd;
typedef Eigen::Matrix<uint,1, Eigen::Dynamic > RowVectorXui;
typedef Eigen::Matrix<double,6,1> Vector6d;
typedef Eigen::Matrix<uint,Eigen::Dynamic,6> MatrixX6ui;
typedef Eigen::Array<uint,Eigen::Dynamic,6> ArrayX6ui;
typedef Eigen::Array<uint,1,6> Array6ui;
typedef Eigen::Matrix<uint,Eigen::Dynamic,1> VectorXui;
typedef std::list<uint> IndexList;
typedef std::list<uint>::iterator IndexListIterator;
typedef std::list<uint>::const_iterator ConstIndexListIterator;
typedef std::vector<ContactPoint*> ContactPointList;
typedef std::vector<WrenchCone* > WrenchConeList;
typedef std::vector<PrimitiveSearchZone*> SearchZone;
typedef std::vector<FingerParameters> FParamList;
typedef std::vector<Patch*> ContactRegion;
typedef std::tr1::shared_ptr<TargetObject> TargetObjectPtr;
typedef std::tr1::shared_ptr<Grasp> GraspPtr;
typedef std::tr1::shared_ptr<SearchZones> SearchZonesPtr;
typedef std::tr1::shared_ptr<OWS> OWSPtr;
typedef std::tr1::shared_ptr<std::vector<Patch*> > PatchListPtr;
typedef std::vector<std::tr1::shared_ptr<Finger> > FingerPtrList;
typedef std::tr1::shared_ptr<IndependentContactRegions>  IndependentContactRegionsPtr;
typedef std::tr1::shared_ptr<double>  SharedDoublePtr;
typedef std::tr1::shared_ptr<WrenchSpace> WrenchSpacePtr;
typedef std::tr1::shared_ptr<DiscreteWrenchSpace> DiscreteWrenchSpacePtr;
typedef std::tr1::shared_ptr<DiscreteTaskWrenchSpace> DiscreteTaskWrenchSpacePtr;
typedef std::tr1::shared_ptr<SphericalWrenchSpace> SphericalWrenchSpacePtr;
 typedef std::tr1::shared_ptr<Eigen::Vector3d> Vector3dPtr;
//--------------------------------------------------------------------
Eigen::Matrix3d skewSymmetricMatrix(Eigen::Vector3d vector);
//--------------------------------------------------------------------
uint factorial(uint x);
//--------------------------------------------------------------------
uint dfactorial(uint x);
//--------------------------------------------------------------------
/*!
 * Column-wise conversion from an Eigen::MatrixXd to a double array. Allocates Memory for the array!!! 
 */
 void eigenMatrixToDoubleArray(Eigen::MatrixXd const & matrix, double* & array);
//--------------------------------------------------------------------
/*!
 * Column-wise conversion from a double array to an Eigen::MatrixXd
 */
 void doubleArrayToEigenMatrix( double * const array, Eigen::MatrixXd & matrix);
//--------------------------------------------------------------------

}//namespace ICR
#endif
