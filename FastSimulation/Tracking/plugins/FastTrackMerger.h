#ifndef FastSimulation_Tracking_FastTrackMerger_h
#define FastSimulation_Tracking_FastTrackMerger_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include <vector>

namespace edm { 
  class ParameterSet;
  class Event;
  class EventSetup;
}

namespace reco { 
  class Track;
}

class FastTrackMerger : public edm::EDProducer
{
 public:
  
  explicit FastTrackMerger(const edm::ParameterSet& conf);
  
  virtual ~FastTrackMerger() {}
  
  virtual void beginJob (edm::EventSetup const & es) {}
  
  virtual void produce(edm::Event& e, const edm::EventSetup& es);
  
 private:

  int findId(const reco::Track& aTrack) const;

 private:

  std::vector<edm::InputTag> trackProducers;
  std::vector<edm::InputTag> removeTrackProducers;
  bool tracksOnly;
  double pTMin2;
  unsigned minHits;

};

#endif
