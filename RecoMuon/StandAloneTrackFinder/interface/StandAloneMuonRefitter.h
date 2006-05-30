#ifndef RecoMuon_StandAloneTrackFinder_StandAloneMuonRefitter_H
#define RecoMuon_StandAloneTrackFinder_StandAloneMuonRefitter_H

/** \class StandAloneMuonRefitter
 *  The inward-outward fitter (starts from seed state).
 *
 *  $Date: 2006/05/29 17:56:45 $
 *  $Revision: 1.9 $
 *  \author R. Bellan - INFN Torino
 */

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

// FIXME
// change FreeTrajectoryState in TSOS

class Propagator;
class DetLayer;
class MuonBestMeasurementFinder;
class MuonTrajectoryUpdator;
class Trajectory;

namespace edm {class ParameterSet; class EventSetup; class Event;}

class StandAloneMuonRefitter {
public:
  /// Constructor
  StandAloneMuonRefitter(const edm::ParameterSet& par);

  /// Destructor
  virtual ~StandAloneMuonRefitter();

  // Operations
  
  /// Perform the inner-outward fitting
  void refit(TrajectoryStateOnSurface& initialState, const DetLayer*, Trajectory& trajectory);

  /// the last free trajectory state
  FreeTrajectoryState lastUpdatedFTS() const {return *theLastUpdatedTSOS.freeTrajectoryState();}

  /// the last but one free trajectory state
  FreeTrajectoryState lastButOneUpdatedFTS() const {return *theLastButOneUpdatedTSOS.freeTrajectoryState();}
  
  /// the Trajectory state on the last surface of the fitting
  TrajectoryStateOnSurface lastUpdatedTSOS() const {return theLastUpdatedTSOS;}

  /// the Trajectory state on the last surface of the fitting
  TrajectoryStateOnSurface lastButOneUpdatedTSOS() const {return theLastButOneUpdatedTSOS;}

  void reset();

  /// Pass the Event Setup to the algo at each event
  virtual void setES(const edm::EventSetup& setup);

  /// Pass the Event to the algo at each event
  virtual void setEvent(const edm::Event& event);

  int getTotalChamberUsed() const {return totalChambers;}
  int getDTChamberUsed() const {return dtChambers;}
  int getCSCChamberUsed() const {return cscChambers;}
  int getRPCChamberUsed() const {return rpcChambers;}

  /// Return the propagation direction
  PropagationDirection propagationDirection() const {return thePropagationDirection;}

protected:

private:

  /// Set the last TSOS
  void setLastUpdatedTSOS(TrajectoryStateOnSurface tsos) { theLastUpdatedTSOS = tsos;}
  
  /// Set the last but one TSOS
  void setLastButOneUpdatedTSOS(TrajectoryStateOnSurface tsos) { theLastButOneUpdatedTSOS = tsos;}
  
  /// Increment the DT,CSC,RPC counters
  void incrementChamberCounters(const DetLayer *layer);

  /// I have to use this method since I have to cope with two propagation direction
  void vectorLimits(std::vector<const DetLayer*> &vect,
		    std::vector<const DetLayer*>::const_iterator &vector_begin,
		    std::vector<const DetLayer*>::const_iterator &vector_end) const;
  /// I have to use this method since I have to cope with two propagation direction
  void incrementIterator(std::vector<const DetLayer*>::const_iterator &iter) const;

  /// the trajectory state on the last available surface
  TrajectoryStateOnSurface theLastUpdatedTSOS;
  /// the trajectory state on the last but one available surface
  TrajectoryStateOnSurface theLastButOneUpdatedTSOS;

  /// The Measurement extractor
  MuonDetLayerMeasurements theMeasurementExtractor;
  
  /// The propagator
  Propagator *thePropagator;
  
  /// access at the propagator
  Propagator *propagator() const {return thePropagator;}

  /// The Estimator
  MeasurementEstimator *theEstimator;

  /// access at the estimator
  MeasurementEstimator *estimator() const {return theEstimator;}

  /// the bestMeasurementFinder
  MuonBestMeasurementFinder *theBestMeasurementFinder;

  /// access at the bestMeasurementFinder
  MuonBestMeasurementFinder *bestMeasurementFinder() const {return theBestMeasurementFinder;}

  /// the muon updator (it doesn't inhert from an updator, but it has one!)
  MuonTrajectoryUpdator *theMuonUpdator;

  /// access at the muon updator
  MuonTrajectoryUpdator *updator() const {return theMuonUpdator;}

  // The max allowed chi2 to accept a rechit in the fit
  double theMaxChi2;
  // The errors of the trajectory state are multiplied by nSigma 
  // to define acceptance of BoundPlane and maximalLocalDisplacement
  double theNSigma;

  // the propagation direction
  PropagationDirection thePropagationDirection;

  int totalChambers;
  int dtChambers;
  int cscChambers;
  int rpcChambers;
};
#endif

