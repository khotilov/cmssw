#ifndef RecoMuon_TrackerSeedGenerator_TSGFromL1Muon_H
#define RecoMuon_TrackerSeedGenerator_TSGFromL1Muon_H

/** \class TSGFromL1Muon
 * Description: 
 * EDPRoducer to generate L3MuonTracjectorySeed from L1MuonParticles
 * \author Marcin Konecki
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace edm { class Event; class EventSetup; }
class L1MuonRegionProducer;
class L1MuonPixelTrackFitter;
class OrderedHitsGenerator;
class PixelTrackFilter;
class L1MuonSeedsMerger;


class TSGFromL1Muon : public edm::EDProducer {
public:
  TSGFromL1Muon(const edm::ParameterSet& cfg);
  virtual ~TSGFromL1Muon();
  virtual void beginRun(edm::Run & run, const edm::EventSetup&es);
  virtual void produce(edm::Event& ev, const edm::EventSetup& es);
private:
 
private:
  edm::ParameterSet theConfig;
  edm::InputTag theSourceTag;


  L1MuonRegionProducer * theRegionProducer;
  OrderedHitsGenerator * theHitGenerator;
  L1MuonPixelTrackFitter * theFitter;
  PixelTrackFilter * theFilter;
  L1MuonSeedsMerger * theMerger;

};
#endif
