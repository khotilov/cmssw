#ifndef TrackProducerBase_h
#define TrackProducerBase_h

/** \class TrackProducerBase
 *  Base Class To Produce Tracks
 *
 *  $Date: 2008/02/18 17:16:12 $
 *  $Revision: 1.14 $
 *  \author cerati
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackerRecHit2D/interface/ClusterRemovalInfo.h"

class Propagator;
class TrajectoryStateUpdator;
class MeasurementEstimator;
class TrackerGeometry;
class TrajectoryFitter;
class TransientTrackingRecHitBuilder;

template <class T>
class TrackProducerBase {
public:
  typedef std::vector<T> TrackCollection;
  typedef std::pair<Trajectory*, std::pair<T*,PropagationDirection> > AlgoProduct;
  typedef std::vector< AlgoProduct >  AlgoProductCollection;
public:
  /// Constructor
  TrackProducerBase(bool trajectoryInEvent = false):
     trajectoryInEvent_(trajectoryInEvent),
        rekeyClusterRefs_(false) {}

  /// Destructor
  virtual ~TrackProducerBase();
  
  /// Get needed services from the Event Setup
  virtual void getFromES(const edm::EventSetup&,
			 edm::ESHandle<TrackerGeometry>& ,
			 edm::ESHandle<MagneticField>& ,
			 edm::ESHandle<TrajectoryFitter>& ,
			 edm::ESHandle<Propagator>& ,
			 edm::ESHandle<TransientTrackingRecHitBuilder>& );

  /// Get TrackCandidateCollection from the Event (needed by TrackProducer)
  virtual void getFromEvt(edm::Event&, edm::Handle<TrackCandidateCollection>&, reco::BeamSpot&);
  /// Get TrackCollection from the Event (needed by TrackRefitter)
  virtual void getFromEvt(edm::Event&, edm::Handle<TrackCollection>&, reco::BeamSpot&);

  /// Method where the procduction take place. To be implemented in concrete classes
  virtual void produce(edm::Event&, const edm::EventSetup&) = 0;

  /// Set parameter set
  void setConf(edm::ParameterSet conf){conf_=conf;}

  /// set label of source collection
  void setSrc(edm::InputTag src, edm::InputTag bsSrc){src_=src;bsSrc_=bsSrc;}

  /// set the aliases of produced collections
  void setAlias(std::string alias){
    alias.erase(alias.size()-6,alias.size());
    alias_=alias;
  }

  /// Sets the information on cluster removal, and turns it on
  void setClusterRemovalInfo(const edm::InputTag &clusterRemovalInfo) {
    rekeyClusterRefs_ = true;
    clusterRemovalInfo_ = clusterRemovalInfo;
  }

  const edm::ParameterSet& getConf() const {return conf_;}
 private:
  edm::ParameterSet conf_;
  edm::InputTag src_;
 protected:
  std::string alias_;
  bool trajectoryInEvent_;
  edm::OrphanHandle<TrackCollection> rTracks_;
  edm::InputTag bsSrc_;

  bool rekeyClusterRefs_;
  edm::InputTag clusterRemovalInfo_;
};

#include "RecoTracker/TrackProducer/interface/TrackProducerBase.icc"

#endif
