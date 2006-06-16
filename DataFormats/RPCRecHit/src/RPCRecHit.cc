/*
 *  See header file for a description of this class.
 *
 *  $Date: 2006/06/16 08:08:21 $
 *  $Revision: 1.4 $
 *  \author M. Maggi -- INFN Bari
 */


#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"


RPCRecHit::RPCRecHit(const RPCDetId& rpcId, int bx) : 
  theRPCId(rpcId), theBx(bx),theFirstStrip(99),theClusterSize(99), theLocalPosition(), theLocalError() 
{
}

RPCRecHit::RPCRecHit() : 
  theRPCId(), theBx(99),theFirstStrip(99),theClusterSize(99), theLocalPosition(), theLocalError() 
{
}


RPCRecHit::RPCRecHit(const RPCDetId& rpcId, int bx, const LocalPoint& pos) : 
  theRPCId(rpcId), theBx(bx), theFirstStrip(99),theClusterSize(99), theLocalPosition(pos) 
{
  float stripResolution = 3.0 ; //cm  this sould be taken from trimmed cluster size times strip size 
                                 //    taken out from geometry service i.e. topology
  theLocalError =
    LocalError(stripResolution*stripResolution, 0., 0.); //FIXME: is it really needed?
}



// Constructor from a local position and error, wireId and digi time.
RPCRecHit::RPCRecHit(const RPCDetId& rpcId,
		     int bx,
		     const LocalPoint& pos,
		     const LocalError& err) :
  theRPCId(rpcId), theBx(bx),theFirstStrip(99), theClusterSize(99), theLocalPosition(pos), theLocalError(err) 
{
}


// Constructor from a local position and error, wireId, bx and cluster size.
RPCRecHit::RPCRecHit(const RPCDetId& rpcId,
		     int bx,
		     int firstStrip,
		     int clustSize,
		     const LocalPoint& pos,
		     const LocalError& err) :
  theRPCId(rpcId), theBx(bx),theFirstStrip(firstStrip), theClusterSize(clustSize), theLocalPosition(pos), theLocalError(err) 
{
}




// Destructor
RPCRecHit::~RPCRecHit()
{
}



RPCRecHit * RPCRecHit::clone() const {
  return new RPCRecHit(*this);
}



// Return the detId of the Det (a RPCLayer).
DetId RPCRecHit::geographicalId() const {
  return theRPCId;
}



// Access to component RecHits.
// No components rechits: it returns a null vector
std::vector<const TrackingRecHit*> RPCRecHit::recHits() const {
  std::vector<const TrackingRecHit*> nullvector;
  return nullvector; 
}



// Non-const access to component RecHits.
// No components rechits: it returns a null vector
std::vector<TrackingRecHit*> RPCRecHit::recHits() {
  std::vector<TrackingRecHit*> nullvector;
  return nullvector; 
}


// Comparison operator, based on the wireId and the digi time
bool RPCRecHit::operator==(const RPCRecHit& hit) const {
  return this->geographicalId() == hit.geographicalId(); 
}


// The ostream operator
std::ostream& operator<<(std::ostream& os, const RPCRecHit& hit) {
  os << "pos: " << hit.localPosition().x() ; 
  os << " +/- " << sqrt(hit.localPositionError().xx()) ;
  return os;
}
