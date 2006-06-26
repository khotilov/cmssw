//FAMOS Headers
#include "FastSimulation/Event/interface/FSimEvent.h"

//C++ Headers

FSimEvent::FSimEvent(const edm::ParameterSet& vtx,
		     const edm::ParameterSet& kine) 
    : FBaseSimEvent(vtx,kine), id_(edm::EventID(0,0)), weight_(0)
{}
 
FSimEvent::~FSimEvent()
{}

void 
FSimEvent::fill(const HepMC::GenEvent& hev, edm::EventID& Id) { 
  FBaseSimEvent::fill(hev); 
  id_ = Id;
}
    
edm::EventID 
FSimEvent::id() const { 
  return id_; 
}
   
float FSimEvent::weight() const { 
  return weight_; 
}

unsigned int 
FSimEvent::nTracks() const {
  return FBaseSimEvent::nTracks();
}

unsigned int 
FSimEvent::nVertices() const { 
  return FBaseSimEvent::nVertices();
}

unsigned int 
FSimEvent::nGenParts() const {
  return FBaseSimEvent::nGenParts();
}

void 
FSimEvent::load(edm::SimTrackContainer & c) const
{
  for (unsigned int i=0; i<nTracks(); ++i) {
    //    SimTrack t = SimTrack(ip,p,iv,ig);
    c.push_back(embdTrack(i));
  }
}

void 
FSimEvent::load(edm::SimVertexContainer & c) const
{
  for (unsigned int i=0; i<nVertices(); ++i) {
    //    SimTrack t = SimTrack(ip,p,iv,ig);
    c.push_back(embdVertex(i));
  }
}









