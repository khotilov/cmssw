#ifndef TrackingTools_TrackRefitter_TracksToTrajectories_H
#define TrackingTools_TrackRefitter_TracksToTrajectories_H

/** \class TracksToTrajectories
 *  This class, which is a EDProducer, takes a reco::TrackCollection from the Event and refits the rechits 
 *  strored in the reco::Tracks. The final result is a std::vector of Trajectories (objs of the type "Trajectory"), 
 *  which is loaded into the Event in a transient way
 *
 *  $Date: 2006/11/23 11:34:38 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - INFN Torino <riccardo.bellan@cern.ch>
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

namespace edm {class ParameterSet; class Event; class EventSetup;}
class TrackTransformer;

class TracksToTrajectories: public edm::EDProducer{
public:

  /// Constructor
  TracksToTrajectories(const edm::ParameterSet&);

  /// Destructor
  virtual ~TracksToTrajectories();
  
  // Operations

  /// Convert a reco::TrackCollection into std::vector<Trajectory>
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
 protected:
  
 private:
  
  edm::InputTag theTracksLabel;
  TrackTransformer *theTrackTransformer;
};
#endif

