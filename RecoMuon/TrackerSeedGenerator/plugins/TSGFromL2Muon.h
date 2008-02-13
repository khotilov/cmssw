#ifndef RecoMuon_TrackerSeedGenerator_TSGFromL2Muon_H
#define RecoMuon_TrackerSeedGenerator_TSGFromL2Muon_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace edm { class ParameterSet; class Event; class EventSetup; }

class MuonServiceProxy;
class TrackerSeedGenerator;
class MuonTrackingRegionBuilder;

//
// generate seeds corresponding to L2 muons
//

class TSGFromL2Muon : public edm::EDProducer {
public:
  TSGFromL2Muon(const edm::ParameterSet& cfg);
  virtual ~TSGFromL2Muon();
  virtual void beginJob(const edm::EventSetup& es);
  virtual void produce(edm::Event& ev, const edm::EventSetup& es);
 
private:
  edm::ParameterSet theConfig;
  edm::InputTag theL2CollectionLabel;

  MuonServiceProxy* theService;
  double thePtCut;
  MuonTrackingRegionBuilder* theRegionBuilder;
  TrackerSeedGenerator* theTkSeedGenerator;

};
#endif
