#ifndef RecoMuon_MuonSeedGenerator_SETFilter_H
#define RecoMuon_MuonSeedGenerator_SETFilter_H

/** \class SETFilter
    I. Bloch, E. James, S. Stoynev
 */

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/DetLayers/interface/NavigationDirection.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"

#include "CLHEP/Matrix/DiagMatrix.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "FWCore/Framework/interface/ESHandle.h"


// used in the SET algorithm
struct seedSet{
  MuonTransientTrackingRecHit::MuonRecHitContainer theSet;
  Hep3Vector momentum;
  int charge;
  double weight;
  std::vector <TrajectoryMeasurement> trajectoryMeasurementsInTheSet;
};

class DetLayer;
class Trajectory;
//class MuonServiceProxy;
class TrajectoryFitter;

namespace edm {class ParameterSet; class EventSetup; class Event;}

class SETFilter {

 public:
    /// Constructor
  SETFilter(const edm::ParameterSet& par, const MuonServiceProxy* service);

  /// Destructor
  virtual ~SETFilter();

  // Operations
  
  /// Perform the inner-outward fitting
  void refit(const TrajectoryStateOnSurface& initialState, const DetLayer*, Trajectory& trajectory);

  /// the last free trajectory state
  FreeTrajectoryState lastUpdatedFTS() const {return *theLastUpdatedTSOS.freeTrajectoryState();}

  /// the Trajectory state on the last surface of the fitting
  TrajectoryStateOnSurface lastUpdatedTSOS() const {return theLastUpdatedTSOS;}

  /// Perform the SET inner-outward fitting
  bool fwfit_SET( std::vector < seedSet> & validSegmentsSet,
                  std::vector < TrajectoryMeasurement > & trajectoryMeasurementsInTheSet);


  ///  transforms "segment trajectory" to "rechit container"  
  bool transform(Trajectory::DataContainer &measurements_segments,
                 TransientTrackingRecHit::ConstRecHitContainer & hitContainer, 
		 TrajectoryStateOnSurface & firstTSOS);

  ///  transforms "segment trajectory" to "segment container" 
  bool transformLight(Trajectory::DataContainer &measurements_segments,
		      TransientTrackingRecHit::ConstRecHitContainer & hitContainer, 
		      TrajectoryStateOnSurface & firstTSOS);



  void reset();

  /// Pass the Event to the algo at each event
  virtual void setEvent(const edm::Event& event);

  int getTotalChamberUsed() const {return totalChambers;}
  int getDTChamberUsed() const {return dtChambers;}
  int getCSCChamberUsed() const {return cscChambers;}
  int getRPCChamberUsed() const {return rpcChambers;}

  inline bool goodState() const {return totalChambers >= 2 && 
				   ((dtChambers + cscChambers) >0 );}
  
  /// return the layer used for the refit
  std::vector<const DetLayer*> layers() const {return theDetLayers;}

  /// return the last det layer
  const DetLayer* lastDetLayer() const {return theDetLayers.back();}

  /// Return the propagation direction
  PropagationDirection propagationDirection() const;

protected:

private:

  /// Set the last TSOS
  void setLastUpdatedTSOS(TrajectoryStateOnSurface tsos) { theLastUpdatedTSOS = tsos;}
  
  /// Set the last but one TSOS
  void setLastButOneUpdatedTSOS(TrajectoryStateOnSurface tsos) { theLastButOneUpdatedTSOS = tsos;}
  
  /// Increment the DT,CSC,RPC counters
  void incrementChamberCounters(const DetLayer *layer);
 
  /// access at the propagator
  const Propagator *propagator() const;
  

  //---- SET
  // FTS <-> parameters
  void getFromFTS(const FreeTrajectoryState& fts,
                  Hep3Vector& p3, Hep3Vector& r3,
                  int& charge, AlgebraicSymMatrix66& cov);

  FreeTrajectoryState getFromCLHEP(const Hep3Vector& p3, const Hep3Vector& r3,
                                   int charge, const AlgebraicSymMatrix66& cov,
                                   const MagneticField* field);

  // chi2 functions (calculate chi2)
  double findChi2(double pX, double pY, double pZ,
                    const Hep3Vector& r3T,
                    seedSet & muonCandidate,
                    TrajectoryStateOnSurface  &lastUpdatedTSOS,
                    std::vector < TrajectoryMeasurement > & trajectoryMeasurementsInTheSet,
                    bool detailedOutput);

  double findMinChi2(unsigned int iSet, const Hep3Vector& r3T,
                 seedSet & muonCandidate,
                 std::vector < TrajectoryStateOnSurface > &lastUpdatedTSOS_Vect,
                 std::vector < TrajectoryMeasurement > & trajectoryMeasurementsInTheSet);

  double chi2AtSpecificStep(Hep3Vector &foot,
                            const Hep3Vector& r3T,
                            seedSet & muonCandidate,
                            TrajectoryStateOnSurface  &lastUpdatedTSOS,
                            std::vector < TrajectoryMeasurement > & trajectoryMeasurementsInTheSet,
                            bool detailedOutput);

  // find initial points for the SIMPLEX minimization
  std::vector <Hep3Vector> find3MoreStartingPoints(Hep3Vector &key_foot,
                                                   const Hep3Vector& r3T,
                                                   seedSet & muonCandidate);

  std::pair <double,double> findParabolaMinimum(std::vector <double> &quadratic_var,
                                                std::vector <double> &quadratic_chi2);

  // SIMPLEX minimization functions
  void pickElemets(std::vector <double> &chi2Feet,
                   unsigned int & high, unsigned int & second_high, unsigned int & low);

  Hep3Vector reflectFoot(std::vector <Hep3Vector> & feet,
                         unsigned int key_foot, double scale );

  void nDimContract(std::vector <Hep3Vector> & feet, unsigned int low);
  //---- SET
  
  /// the propagator name
  std::string thePropagatorName;

  /// the propagation direction
  NavigationDirection theFitDirection;

  /// the trajectory state on the last available surface
  TrajectoryStateOnSurface theLastUpdatedTSOS;
  /// the trajectory state on the last but one available surface
  TrajectoryStateOnSurface theLastButOneUpdatedTSOS;

  /// the det layer used in the reconstruction
  std::vector<const DetLayer*> theDetLayers;

  int totalChambers;
  int dtChambers;
  int cscChambers;
  int rpcChambers;

  bool useSegmentsInTrajectory;

  /// used in the SET BW fit
  edm::ESHandle<TrajectoryFitter> theBWLightFitter;
  std::string theBWLightFitterName;

  const MuonServiceProxy *theService;
  //bool theOverlappingChambersFlag;
};
#endif

