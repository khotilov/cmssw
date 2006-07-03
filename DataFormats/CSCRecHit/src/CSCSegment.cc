/** \file CSCSegment.cc
 *
 *  $Date: 2006/06/29 16:44:31 $
 *  \author Matteo Sani
 */

#include <DataFormats/CSCRecHit/interface/CSCSegment.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h>
#include <iostream>

CSCSegment::CSCSegment(std::vector<CSCRecHit2D> proto_segment, LocalPoint origin, 
	LocalVector direction, AlgebraicSymMatrix errors, double chi2) : theCSCRecHits(proto_segment),
	theOrigin(origin), theLocalDirection(direction), theCovMatrix(errors), theChi2(chi2) {

	// Get CSCDetId from one of the rechits, but then remove the layer part so it's a _chamber_ id
	CSCDetId id = theCSCRecHits.begin()->cscDetId();
	theDetId = CSCDetId(id.endcap(),id.station(),id.ring(),id.chamber(), 0);
}

CSCSegment::~CSCSegment() {}

std::vector<const TrackingRecHit*> CSCSegment::recHits() const{
  std::vector<const TrackingRecHit*> pointersOfRecHits;
  for (std::vector<CSCRecHit2D>::const_iterator irh = theCSCRecHits.begin(); irh!=theCSCRecHits.end(); ++irh) {
    pointersOfRecHits.push_back(&(*irh));
  }
  return pointersOfRecHits;
}

std::vector<TrackingRecHit*> CSCSegment::recHits() {
  
  std::vector<TrackingRecHit*> pointersOfRecHits;
  for (std::vector<CSCRecHit2D>::iterator irh = theCSCRecHits.begin(); irh!=theCSCRecHits.end(); ++irh) {
    pointersOfRecHits.push_back(&(*irh));
  }
  return pointersOfRecHits;
}

LocalError CSCSegment::localPositionError() const {
  return LocalError(theCovMatrix[2][2], theCovMatrix[2][3], theCovMatrix[3][3]);
}

LocalError CSCSegment::localDirectionError() const {
  return LocalError(theCovMatrix[0][0], theCovMatrix[0][1], theCovMatrix[1][1]); 
}


AlgebraicVector CSCSegment::parameters() const {
  // (x,y,dx/dz,dy/dz)
  AlgebraicVector result(4);
  result[2] = theOrigin.x();
  result[3] = theOrigin.y();
  result[0] = theLocalDirection.x()/theLocalDirection.z();
  result[1] = theLocalDirection.y()/theLocalDirection.z();    
  return result;
}


AlgebraicMatrix CSCSegment::projectionMatrix() const {
  static AlgebraicMatrix theProjectionMatrix( 4, 5, 0);
  static bool isInitialized = false;
  if (!isInitialized) {
    theProjectionMatrix[0][1] = 1;
    theProjectionMatrix[1][2] = 1;
    theProjectionMatrix[2][3] = 1;
    theProjectionMatrix[3][4] = 1;
    isInitialized=true;
  }    
  return theProjectionMatrix;
}


void CSCSegment::print() const {
  std::cout << *this << std::endl;
}

std::ostream& operator<<(std::ostream& os, const CSCSegment& seg) {
  os << "CSCSegment: local pos = " << seg.localPosition() << 
    " dir = " << seg.localDirection() <<
    " chi2 = " << seg.chi2() << " #rechits = " << seg.specificRecHits().size();
  return os;  
}

/*
const CSCChamber* CSCSegment::chamber() const { return theChamber; }
*/
