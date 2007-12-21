#ifndef Alignment_KalmanAlignmentAlgorithm_KalmanAlignmentTracklet_h
#define Alignment_KalmanAlignmentAlgorithm_KalmanAlignmentTracklet_h

#include "DataFormats/GeometrySurface/interface/ReferenceCounted.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmBase.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

/// Smart container for purely transient trajectory/track-pairs and, if available, an
/// external measurement (stored as TrajectoryStateOnSurface).


class KalmanAlignmentTracklet : public ReferenceCounted
{

public:

  typedef std::pair<Trajectory*, reco::Track*> TrajTrackPair;
  typedef std::vector< TrajTrackPair > TrajTrackPairCollection;
  typedef AlignmentAlgorithmBase::ConstTrajTrackPair ConstTrajTrackPair;
  typedef AlignmentAlgorithmBase::ConstTrajTrackPairCollection ConstTrajTrackPairCollection;

  typedef ReferenceCountingPointer< KalmanAlignmentTracklet > TrackletPtr;

  /// Contructor. NOTE: The container gains the ownership of the trajectory/track at construction time.
  KalmanAlignmentTracklet( TrajTrackPair& trajTrackPair,
			   const TrajectoryStateOnSurface& external );

  KalmanAlignmentTracklet( TrajTrackPair& trajTrackPair );

  /// Destructor.
  ~KalmanAlignmentTracklet( void );

  inline const Trajectory* trajectory( void ) const { return theTrajTrackPair.first; }
  inline const reco::Track* track( void ) const { return theTrajTrackPair.second; }
  inline const ConstTrajTrackPair trajTrackPair( void ) const { return theTrajTrackPair; }

  inline const TrajectoryStateOnSurface externalPrediction( void ) const { return theExternalPrediction; }

  inline bool externalPredictionAvailable( void ) const { return theExternalPredictionFlag; }

private:

  ConstTrajTrackPair theTrajTrackPair;
  TrajectoryStateOnSurface theExternalPrediction;
  bool theExternalPredictionFlag;
};

#endif
