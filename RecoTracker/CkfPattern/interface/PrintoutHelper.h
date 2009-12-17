#ifndef RecoTracker_CkfPattern_PrintoutHelper_h
#define RecoTracker_CkfPattern_PrintoutHelper_h

#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/PatternTools/interface/TrajectoryBuilder.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/PatternTools/interface/bqueue.h"

class TrackerGeometry;

class PrintoutHelper{
 public:
  template< class collection > static std::string dumpCandidates( collection & candidates);
  static std::string dumpMeasurements(const std::vector<TrajectoryMeasurement> & v) ;
  static std::string dumpMeasurements(const cmsutils::bqueue<TrajectoryMeasurement> & v);
  static std::string dumpMeasurement(const TrajectoryMeasurement & tm);
  static std::string regressionTest(const TrackerGeometry & tracker,std::vector<Trajectory> & unsmoothedResult);
};



template< class collection > 
std::string PrintoutHelper::dumpCandidates( collection & candidates) {
  std::stringstream buffer;
  uint ic=0;
  typename collection::const_iterator traj=candidates.begin();
  for (;traj!=candidates.end(); traj++) {  
    buffer<<ic++<<"] ";
    if (!traj->measurements().empty()){
      const TrajectoryMeasurement & last = traj->lastMeasurement();

      buffer<<"with: "<<traj->measurements().size()<<" measurements."<< traj->lostHits() << " lost, " << traj->foundHits()<<" found, chi2="<<traj->chiSquared()<<"\n";
      if (last.updatedState().isValid()) {
	const TrajectoryStateOnSurface & tsos = last.updatedState();
	buffer <<"Last [Updated] state\n x: "<<tsos.globalPosition()<<"\n p: "<<tsos.globalMomentum()<<"\n";
      } else if(last.forwardPredictedState().isValid()){
	const TrajectoryStateOnSurface & tsos = last.forwardPredictedState();
	buffer <<"Last [fwdPredicted] state\n x: "<<tsos.globalPosition()<<"\n p: "<<tsos.globalMomentum()<<"\n";
      } else if (last.predictedState().isValid()){
	const TrajectoryStateOnSurface & tsos = last.predictedState();
	buffer <<"Last [Predicted] state\n x: "<<tsos.globalPosition()<<"\n p: "<<tsos.globalMomentum()<<"\n";
      }
      buffer <<" hit is: "<<(last.recHit()->isValid()?"valid":"invalid")<<"\n";
      if (last.recHit()->isValid())
	buffer <<"on detId: "<<last.recHit()->geographicalId().rawId()<<"\n";
    }
    else{
      buffer<<" no measurement. \n";}
  }
  return buffer.str();
}


#endif
