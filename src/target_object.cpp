#include "../include/target_object.h"
#include "../include/debug.h"
#include <assert.h>
#include <math.h>

namespace ICR
{
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  TargetObject::TargetObject() : num_cp_(0), centroid_(new Eigen::Vector3d), centroid_computed_(false){contact_points_.clear();}
  //--------------------------------------------------------------------
  TargetObject::TargetObject(std::string const& name) : name_(name), num_cp_(0), centroid_(new Eigen::Vector3d), centroid_computed_(false) {contact_points_.clear();}
  //--------------------------------------------------------------------
  TargetObject::TargetObject(TargetObject const& src) : name_(src.name_), num_cp_(src.num_cp_),
							contact_points_(src.contact_points_), centroid_(src.centroid_), centroid_computed_(src.centroid_computed_) {}
  //--------------------------------------------------------------------
  TargetObject& TargetObject::operator=(TargetObject const& src)		  
  {
    if (this !=&src)
      {
	name_=src.name_;
	num_cp_=src.num_cp_;
	contact_points_=src.contact_points_;
	centroid_=src.centroid_;
	centroid_computed_=src.centroid_computed_;
      }
    return *this;
  }
  //--------------------------------------------------------------------
  std::ostream& operator<<(std::ostream& stream, TargetObject const& obj)
  {
    stream <<'\n'<<"TARGET OBJECT: "<<'\n'
	   <<"Name: "<<obj.name_<<'\n'
	   <<"Centroid computed: "<<obj.centroid_computed_<<'\n'
	   <<"Number of contact points: "<<obj.num_cp_<<'\n'<<'\n';

    return stream;
  }
  //--------------------------------------------------------------------
  TargetObject::~TargetObject()
  {
    for (uint i=0;i<contact_points_.size();i++)
      delete contact_points_[i];
  }
  //--------------------------------------------------------------------
  std::string const TargetObject::getName() const {return name_;}
  //--------------------------------------------------------------------
  void TargetObject::setName(std::string const& name) {name_=name;}
  //--------------------------------------------------------------------
  void TargetObject::reserveCpList(uint num_cp){contact_points_.reserve(num_cp);}
  //--------------------------------------------------------------------
  // void TargetObject::findClosestIdx(double x, double y, double z) {
  //   for (uint i=0; i<num_cp_ ; i++) {
    
    
  //   }
  // }
  // //--------------------------------------------------------------------
  // void TargetObject::findClosestIdx(Eigen::Vector3d& pt) {

  // }
  //--------------------------------------------------------------------
  uint TargetObject::getNumCp() const {return num_cp_;}
  //--------------------------------------------------------------------
  void TargetObject::addContactPoint(ContactPoint const& point)
  {
    centroid_computed_=false;

    //Ensure unit normal
    assert(fabs( (*point.getVertexNormal()).norm()-1) <= EPSILON_UNIT_NORMAL);

    contact_points_.push_back(new ContactPoint(point));
    num_cp_++;
  }
  //--------------------------------------------------------------------
  void TargetObject::scaleObject(double scale)
  {
    for (uint i=0; i < contact_points_.size(); i++)
      *(contact_points_[i]->getVertex())=*(contact_points_[i]->getVertex())*scale;

    if (centroid_computed_)
      (*centroid_)=(*centroid_)*scale;
  }
  //--------------------------------------------------------------------
  ContactPoint const* TargetObject::getContactPoint(uint id) const {return contact_points_.at(id);}
  //--------------------------------------------------------------------
  void TargetObject::computeCentroid()
  {
    assert(contact_points_.size() > 0);

    centroid_->setZero();
    for (uint i=0; i<contact_points_.size(); i++)
      (*centroid_)+= (*(contact_points_[i]->getVertex()));

    (*centroid_)=(*centroid_)/contact_points_.size();

    centroid_computed_=true;
  }
  //--------------------------------------------------------------------
  Vector3dPtr TargetObject::getCentroid()const
  {
    assert(centroid_computed_);
    return centroid_;
  }
  //--------------------------------------------------------------------
  void TargetObject::transform(const Eigen::Affine3d& transform)
  {
    for (uint i=0; i<contact_points_.size(); i++)
	contact_points_[i]->transform(transform);

    if (centroid_computed_)
      (*centroid_)=transform*(*centroid_);

 }
  //--------------------------------------------------------------------
bool TargetObject::writeToFile(const std::string& points_path, const std::string& normals_path,const std::string& neighbors_path)const
  {
    remove(points_path.c_str());
    remove(normals_path.c_str());
    remove(neighbors_path.c_str());
    FILE* cp=fopen (points_path.c_str(),"a");
    FILE* vn=fopen (normals_path.c_str(),"a");
    FILE* nb=fopen (neighbors_path.c_str(),"a");
    if(!vn |!cp |!nb )
      {
	std::cout<<"Warning in TargetObject::writeToFile(const std::string& points_path, const std::string& normals_path,const std::string& neighbors_path)const - Couldn't write to file"<<std::endl;
	return false;
      }
    Eigen::Vector3d cp_vtx;
    Eigen::Vector3d cp_vtx_normal;
    for(uint id=0; id < num_cp_;id++)
      {
	cp_vtx=*getContactPoint(id)->getVertex();
	cp_vtx_normal=*getContactPoint(id)->getVertexNormal();
	fprintf(cp, "% f % f % f \n", cp_vtx(0),cp_vtx(1),cp_vtx(2));
	fprintf(vn, "% f % f % f \n", cp_vtx_normal(0),cp_vtx_normal(1),cp_vtx_normal(2));

	for (ConstIndexListIterator it= getContactPoint(id)->getNeighborItBegin(); it !=getContactPoint(id)->getNeighborItEnd(); it++)
	  fprintf(nb, "% d",*it+1);//Indexing of neighboring points starting from 1 instead of 0

	fprintf(nb, "\n");
      }
    fclose (cp); fclose (vn); fclose (nb);
    return true;
}
  //--------------------------------------------------------------------
}//namespace ICR
