#include "FastSimulation/Event/interface/FBaseSimEvent.h"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "FastSimulation/Event/interface/FSimVertex.h"

  /// Default constructor
FSimVertex::FSimVertex() : SimVertex(), mom_(0), id_(-1) {;}
  
  /// constructor from the embedded vertex index in the FBaseSimEvent
FSimVertex::FSimVertex(const XYZTLorentzVector& v, int im, int id, FBaseSimEvent* mom) : 
  //    SimVertex(Hep3Vector(v.vect(),v.e(),im), mom_(mom), id_(id) 
  SimVertex(Hep3Vector(v.X(),v.Y(),v.Z()),v.T(),im), mom_(mom), id_(id),
  position_(v) {;}

std::ostream& operator <<(std::ostream& o , const FSimVertex& t) {
  return o << t;
}
