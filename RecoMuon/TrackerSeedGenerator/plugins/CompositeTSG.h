#ifndef RecoMuon_TrackerSeedGenerator_CompositeTSG_H
#define RecoMuon_TrackerSeedGenerator_CompositeTSG_H

#include "RecoMuon/TrackerSeedGenerator/interface/TrackerSeedGenerator.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include <TH1.h>

class TrackingRegion;
class MuonServiceProxy;

#include "FWCore/MessageLogger/interface/MessageLogger.h"

class CompositeTSG : public TrackerSeedGenerator {

public:
  typedef std::vector<TrajectorySeed> BTSeedCollection;  
  typedef std::pair<const Trajectory*, reco::TrackRef> TrackCand;

  CompositeTSG(const edm::ParameterSet &pset);
  virtual ~CompositeTSG();

  void init(const MuonServiceProxy *service);
  void setEvent(const edm::Event &event);

  virtual void trackerSeeds(const TrackCand&, const TrackingRegion&, BTSeedCollection &) =0;

 protected:
  bool useTFileService_;

  uint nTSGs() { return theTSGs.size();}
  std::vector<TrackerSeedGenerator*> theTSGs;
  std::vector<std::string> theNames;
  std::vector<TH1F*> theHistos;
  std::string theCategory;

  const MuonServiceProxy * theProxyService;
};


#endif 
