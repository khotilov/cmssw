#ifndef RecoEGAMMA_ConversionTrack_ConversionTrackFinder_h
#define RecoEGAMMA_ConversionTrack_ConversionTrackFinder_h

/** \class ConversionTrackFinder
 **  
 **
 **  $Id: ConversionTrackFinder.h,v 1.6 2007/05/29 10:11:45 elmer Exp $ 
 **  $Date: 2007/05/29 10:11:45 $ 
 **  $Revision: 1.6 $
 **  \author Nancy Marinelli, U. of Notre Dame, US
 **
 ***/
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
//
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
//
#include "RecoEgamma/EgammaPhotonAlgos/interface/ConversionSeedFinder.h"

#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "TrackingTools/PatternTools/interface/TrajectoryBuilder.h"

// C/C++ headers
#include <string>
#include <vector>

class TransientInitialStateEstimator;
class ConversionTrackFinder {

 public:
  
  ConversionTrackFinder( const edm::EventSetup& es,
			 const edm::ParameterSet& config );
                       
  
  virtual ~ConversionTrackFinder();
 
  
  virtual std::vector<Trajectory> tracks(const TrajectorySeedCollection seeds , TrackCandidateCollection &candidate) const =0;

  /// Initialize EventSetup objects at each event
  void setEventSetup( const edm::EventSetup& es ) ; 
  void setEvent(const  edm::Event& e ) ; 


 private:

   



 protected: 
  
  edm::ParameterSet conf_;
  const MagneticField* theMF_;
  
  const MeasurementTracker*     theMeasurementTracker_;
  const TrajectoryBuilder*  theCkfTrajectoryBuilder_;

  TransientInitialStateEstimator* theInitialState_;  
  const TrackerGeometry* theTrackerGeom_;
  KFUpdator*                          theUpdator_;


struct ExtractNumOfHits {
  typedef int result_type;
  result_type operator()(const Trajectory& t) const {return t.foundHits();}
  result_type operator()(const Trajectory* t) const {return t->foundHits();}
};


struct ExtractChi2 {
  typedef float result_type;
  result_type operator()(const Trajectory& t) const {return t.chiSquared();}
  result_type operator()(const Trajectory* t) const {return t->chiSquared();}
};


 

};

#endif
