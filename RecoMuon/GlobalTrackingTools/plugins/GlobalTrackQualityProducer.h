#ifndef GlobalTrackingTools_GlobalTrackQualityProducer_h
#define GlobalTrackingTools_GlobalTrackQualityProducer_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonRefitter.h"

#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

class GlobalMuonRefitter;

class GlobalTrackQualityProducer : public edm::EDProducer {
 public:
  explicit GlobalTrackQualityProducer(const edm::ParameterSet& iConfig);

  virtual ~GlobalTrackQualityProducer(); // {}
  
 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual std::pair<double,double> kink(Trajectory& muon) const ;
  virtual std::pair<double,double> newChi2(Trajectory& muon) const;
  
  edm::InputTag inputCollection_;
  MuonServiceProxy* theService;
  GlobalMuonRefitter* theGlbRefitter;
  MeasurementEstimator *theEstimator;
  //muon::SelectionType selectionType_;
};
#endif
