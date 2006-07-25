#ifndef RecoMuon_TrackingTools_MuonTrajectoryBuilder_H
#define RecoMuon_TrackingTools_MuonTrajectoryBuilder_H

/** \class MuonTrajectoryBuilder
 *  Base class for the Muon reco Trajectory Builder 
 *
 *  $Date: 2006/07/25 12:22:29 $
 *  $Revision: 1.9 $
 *  \author R. Bellan - INFN Torino
 */

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoMuon/TrackingTools/interface/MuonCandidate.h"
#include <vector>

namespace edm {class ParameterSet; class EventSetup; class Event;}

class TrajectorySeed;

class MuonTrajectoryBuilder {

  public:
  
    typedef MuonCandidate::TrajectoryContainer TrajectoryContainer;
    typedef MuonCandidate::CandidateContainer CandidateContainer;

    /// constructor
    MuonTrajectoryBuilder() {}
  
    /// constructor with Parameter set
    MuonTrajectoryBuilder(const edm::ParameterSet&) {}

    /// destructor
    virtual ~MuonTrajectoryBuilder() {}

    /// return a container of the reconstructed trajectories compatible with a given seed
    virtual TrajectoryContainer trajectories(const TrajectorySeed&) = 0;

    /// return a container reconstructed muons starting from a given track
    virtual CandidateContainer trajectories(const reco::TrackRef&) = 0;

    /// pass the Event Setup to the algo at each event
    virtual void setES(const edm::EventSetup& setup) = 0;

    /// pass the Event to the algo at each event
    virtual void setEvent(const edm::Event& event) = 0;
  
 private:

};
#endif
