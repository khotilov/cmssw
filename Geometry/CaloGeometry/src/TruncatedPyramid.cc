#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
//#include <algorithm>
#include <iostream>
//#include "assert.h"

//----------------------------------------------------------------------

TruncatedPyramid::TruncatedPyramid()
 {}

//----------------------------------------------------------------------

TruncatedPyramid::TruncatedPyramid(double dz, double theta, double phi, 
                                   double h1, double bl1, double tl1, double alpha1, 
                                   double h2, double bl2, double tl2, double alpha2)
{
  bool frontSideIsPositiveZ = (trapeziumArea(h2,tl2,bl2) < trapeziumArea(h1,tl1,bl1));
  init(dz,theta,phi,h1,bl1,tl1,alpha1,h2,bl2,tl2,alpha2,frontSideIsPositiveZ);
  // calculate the center of the front face
  int i;
  int offset;
  if (frontSideIsPositiveZ)
    offset = 4;                 // calculate the center of the front
                                // face from the second four points
  else
    offset = 0;                 // or from the first four points
  
  HepGeom::Point3D<double> position;
  for (i=0; i<4; ++i)
    position += HepGeom::Point3D<double>(corners[i + offset].x(),corners[i + offset].y(),corners[i + offset].z());

  //While waiting for *= operator
  position *= 0.25;
  setPosition(GlobalPoint(position.x(),position.y(),position.z()));
}

//----------------------------------------------------------------------

void TruncatedPyramid::init(double dz, double theta, double phi, 
			    double h1, double bl1, double tl1, double alpha1, 
			    double h2, double bl2, double tl2, double alpha2,
			    bool frontSideIsPositiveZ)
{
  corners.resize(8);
  
  double tan_alpha1 = tan(alpha1); // lower plane
  double tan_alpha2 = tan(alpha2); // upper plane

  double tan_theta_cos_phi = tan(theta) * cos(phi);
  double tan_theta_sin_phi = tan(theta) * sin(phi);

  //                       shift due to trapezoid| shift due to the                                               // approximate coordinate 
  //                       axis not parallel to  | fact that the top                                              // signs
  //                       z axis                | and bottom lines
  //                                             | of the trapezium
  //                                             | don't have their
  //                                             | center at the same 
  //                                             | x value (alpha !=0)                                            
  corners[0] = GlobalPoint(-dz*tan_theta_cos_phi - h1 * tan_alpha1       - bl1, -dz*tan_theta_sin_phi - h1 , -dz); // (-,-,-)
  corners[1] = GlobalPoint(-dz*tan_theta_cos_phi + h1 * tan_alpha1       - tl1, -dz*tan_theta_sin_phi + h1 , -dz); // (-,+,-)
  corners[2] = GlobalPoint(-dz*tan_theta_cos_phi + h1 * tan_alpha1       + tl1, -dz*tan_theta_sin_phi + h1 , -dz); // (+,+,-)
  corners[3] = GlobalPoint(-dz*tan_theta_cos_phi - h1 * tan_alpha1       + bl1, -dz*tan_theta_sin_phi - h1 , -dz); // (+,-,-)
                                                                             
  corners[4] = GlobalPoint(dz*tan_theta_cos_phi  - h2 * tan_alpha2       - bl2,  dz*tan_theta_sin_phi - h2  , dz); // (-,-,+)
  corners[5] = GlobalPoint(dz*tan_theta_cos_phi  + h2 * tan_alpha2       - tl2,  dz*tan_theta_sin_phi + h2  , dz); // (-,+,+)
  corners[6] = GlobalPoint(dz*tan_theta_cos_phi  + h2 * tan_alpha2       + tl2,  dz*tan_theta_sin_phi + h2  , dz); // (+,+,+)
  corners[7] = GlobalPoint(dz*tan_theta_cos_phi  - h2 * tan_alpha2       + bl2,  dz*tan_theta_sin_phi - h2  , dz); // (+,-,+)

  // corners[0],corners[3] and corners[1],corners[2] make the parallel lines of the -dz trapezium
  // corners[4],corners[7] and corners[5],corners[6] make the parallel lines of the +dz trapezium

  // determine which one is the front face

  
  //--------------------

  thetaAxis = theta;
  phiAxis = phi;

  // create the boundary planes, the first boundary plane will be the
  // front side

  if (frontSideIsPositiveZ)
    {
      boundaries.push_back(HepPlane3D(HepGeom::Point3D<float>(corners[4].x(),corners[4].y(),corners[4].z()), HepGeom::Point3D<float>(corners[5].x(),corners[5].y(),corners[5].z()), HepGeom::Point3D<float>(corners[6].x(),corners[6].y(),corners[6].z())));
      boundaries.push_back(HepPlane3D(HepGeom::Point3D<float>(corners[0].x(),corners[0].y(),corners[0].z()),HepGeom::Point3D<float>(corners[1].x(),corners[1].y(),corners[1].z()),HepGeom::Point3D<float>(corners[2].x(),corners[2].y(),corners[2].z())));
    }
  else
    {
      boundaries.push_back(HepPlane3D(HepGeom::Point3D<float>(corners[0].x(),corners[0].y(),corners[0].z()),HepGeom::Point3D<float>(corners[1].x(),corners[1].y(),corners[1].z()),HepGeom::Point3D<float>(corners[2].x(),corners[2].y(),corners[2].z())));
      boundaries.push_back(HepPlane3D(HepGeom::Point3D<float>(corners[4].x(),corners[4].y(),corners[4].z()), HepGeom::Point3D<float>(corners[5].x(),corners[5].y(),corners[5].z()), HepGeom::Point3D<float>(corners[6].x(),corners[6].y(),corners[6].z())));
    }
  
  // now generate the other boundaries
  boundaries.push_back(HepPlane3D(HepGeom::Point3D<float>(corners[0].x(),corners[0].y(),corners[0].z()),HepGeom::Point3D<float>(corners[1].x(),corners[1].y(),corners[1].z()),HepGeom::Point3D<float>(corners[4].x(),corners[4].y(),corners[4].z())));
  boundaries.push_back(HepPlane3D(HepGeom::Point3D<float>(corners[1].x(),corners[1].y(),corners[1].z()),HepGeom::Point3D<float>(corners[2].x(),corners[2].y(),corners[2].z()),HepGeom::Point3D<float>(corners[5].x(),corners[5].y(),corners[5].z())));
  boundaries.push_back(HepPlane3D(HepGeom::Point3D<float>(corners[0].x(),corners[0].y(),corners[0].z()),HepGeom::Point3D<float>(corners[3].x(),corners[3].y(),corners[3].z()),HepGeom::Point3D<float>(corners[4].x(),corners[4].y(),corners[4].z())));
  boundaries.push_back(HepPlane3D(HepGeom::Point3D<float>(corners[2].x(),corners[2].y(),corners[2].z()),HepGeom::Point3D<float>(corners[3].x(),corners[3].y(),corners[3].z()),HepGeom::Point3D<float>(corners[6].x(),corners[6].y(),corners[6].z())));
  
  // normalize the planes
  for (unsigned i=0; i < boundaries.size(); ++i)
    boundaries[i] = boundaries[i].normalize();
}

//----------------------------------------------------------------------

bool TruncatedPyramid::inside(const GlobalPoint & Point) const
{

  const HepGeom::Point3D<float> testPoint(Point.x(),Point.y(),Point.z());
  const HepGeom::Point3D<float> testPointb(Point.x(),Point.y(),Point.z());
  // with ordered planes this method becomes simpler
  if ((testPoint-boundaries[0].point(testPointb))*(testPoint-boundaries[1].point(testPointb))>0){
   return false ;
  } 
  if ((testPoint-boundaries[2].point(testPointb))*(testPoint-boundaries[5].point(testPointb))>0){   
   return false ;
  } 
  if ((testPoint-boundaries[3].point(testPointb))*(testPoint-boundaries[4].point(testPointb))>0){ 
   return false ;
  }
  return true;
}


const std::vector<GlobalPoint> & TruncatedPyramid::getCorners() const
{ return corners ; }

//----------------------------------------------------------------------
void 
TruncatedPyramid::hepTransform(const HepTransform3D &transformation)
{


  unsigned int i;

  //Updating corners
  for (i=0; i<corners.size(); ++i)
    {
      HepGeom::Point3D<float> newCorner(corners[i].x(),corners[i].y(),corners[i].z());
      newCorner.transform(transformation);
      corners[i]=GlobalPoint(newCorner.x(),newCorner.y(),newCorner.z());
    }

  //Updating reference position
  const GlobalPoint& position_=CaloCellGeometry::getPosition();
  HepGeom::Point3D<float> newPosition(position_.x(),position_.y(),position_.z());
  newPosition.transform(transformation);
  setPosition(GlobalPoint(newPosition.x(),newPosition.y(),newPosition.z()));

  //Updating Bondaries
  for (i=0; i<boundaries.size(); ++i)
    boundaries[i].transform(transformation);

  //Updating Theta and Phi
  HepVector3D axe(1.,1.,1.); 
  axe.setMag(1.); // must do this first
  axe.setTheta(thetaAxis);
  axe.setPhi(phiAxis);
  axe.transform(transformation);
  thetaAxis = axe.getTheta();
  phiAxis = axe.getPhi();

}

//----------------------------------------------------------------------
double 
TruncatedPyramid::trapeziumArea(double halfHeight, double
                                halfTopLength, double
                                halfBottomLength)
{
  return (halfTopLength + halfTopLength) * halfHeight * 2;
}


const GlobalPoint 
TruncatedPyramid::getPosition(float depth) const 
{
  GlobalPoint point = CaloCellGeometry::getPosition();

  if (depth <= 0)
    return point;

  // only add the vector if depth is positive
  Hep3Vector move(1.,1.,1.); 
  move.setMag(depth); // must do this first
  move.setTheta(thetaAxis);
  move.setPhi(phiAxis);

  // Bart Van de Vyver 10/5/2002 explicit GlobalPoint constructor 
  // to avoid compiler warning
  point = point + GlobalVector(move.x(),move.y(),move.z());

  return point;
}

//----------------------------------------------------------------------

const GlobalPoint TruncatedPyramid::getPosition(float depth, GlobalVector dir) const 
{
  GlobalPoint point = CaloCellGeometry::getPosition();
  
  if (depth <= 0)
    return point;

  // only add the vector if depth is non-zero      
  Hep3Vector move(1.,1.,1.); 
  depth *= cos(thetaAxis - dir.theta()); // project onto crystal axis
  move.setMag(depth); // must do this first
  move.setTheta(thetaAxis);
  move.setPhi(phiAxis);
  // Bart Van de Vyver 10/5/2002 explicit GlobalPoint constructor 
  // to avoid compiler warning
  point = point + GlobalVector(move.x(),move.y(),move.z());
  return point;
}

//----------------------------------------------------------------------



void TruncatedPyramid::dump(const char * prefix="") const {

  // a crystal must have eight corners (not only the front face...)
  assert(getCorners().size() == 8);
  std::cout << prefix << "Center: " <<  CaloCellGeometry::getPosition() << std::endl;
  float thetaaxis_= getThetaAxis();
  float phiaxis_= getPhiAxis();
  std::cout << prefix << "Axis: " <<  thetaaxis_ << " " << phiaxis_ << std::endl;
  //  vector<HepPoint3D> xtCorners=getCorners();
  for ( unsigned int  ci=0; ci !=corners.size(); ci++) {
    std::cout << prefix << "Corner: " << corners[ci] << std::endl;
  }
}
//----------------------------------------------------------------------

std::ostream& operator<<(std::ostream& s,const TruncatedPyramid& cell) {
  assert(cell.getCorners().size() == 8);
  s  << "Center: " <<  cell.getPosition(0.) << std::endl;
  float thetaaxis_= cell.getThetaAxis();
  float phiaxis_= cell.getPhiAxis();
  s  << "Axis: " <<  thetaaxis_ << " " << phiaxis_ << std::endl;
  const std::vector<GlobalPoint>& corners=cell.getCorners(); 
  //  vector<HepPoint3D> xtCorners=getCorners();
  for ( unsigned int  ci=0; ci !=corners.size(); ci++) {
    s  << "Corner: " << corners[ci] << std::endl;
  }
  return s;
}
  
