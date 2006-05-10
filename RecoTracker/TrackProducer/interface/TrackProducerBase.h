#ifndef TrackProducerBase_h
#define TrackProducerBase_h

//
// Package:    RecoTracker/TrackProducer
// Class:      TrackProducerBase
// 
//
// Description: Base Class To Produce Tracks
//
//
// Original Author:  Giuseppe Cerati
//         Created:  Wed May  10 12:29:31 CET 2006
//

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoTracker/TrackProducer/interface/TrackProducerAlgorithm.h"

class Propagator;
class TrajectoryStateUpdator;
class MeasurementEstimator;
class TrackerGeometry;

class TrackProducerBase {
public:
  TrackProducerBase(){}

  virtual ~TrackProducerBase();
  
  virtual void getFromES(const edm::EventSetup&,
			 edm::ESHandle<TrackerGeometry>& ,
			 edm::ESHandle<MagneticField>& ,
			 edm::ESHandle<TrajectoryFitter>& ,
			 edm::ESHandle<Propagator>& );

  virtual void getFromEvt(edm::Event&, edm::Handle<TrackCandidateCollection>&);
  virtual void getFromEvt(edm::Event&, edm::Handle<reco::TrackCollection>&);

  virtual void putInEvt(edm::Event&,
			std::auto_ptr<TrackingRecHitCollection>&,
			std::auto_ptr<reco::TrackCollection>&,
			std::auto_ptr<reco::TrackExtraCollection>&,
			AlgoProductCollection&);

  virtual void produce(edm::Event&, const edm::EventSetup&) = 0;

  void setConf(edm::ParameterSet conf){conf_=conf;}
  //edm::ParameterSet conf(){return conf;}
  void setSrc(std::string src){src_=src;}
private:
  edm::ParameterSet conf_;
  std::string src_;

};

#endif
