#ifndef TrackingTools_TrackRefitter_TrajectoryReader_H
#define TrackingTools_TrackRefitter_TrajectoryReader_H

/** \class TrajectoryReader
 *  No description available.
 *
 *  $Date: 2006/11/22 18:37:21 $
 *  $Revision: 1.1 $
 *  \author R. Bellan - INFN Torino <riccardo.bellan@cern.ch>
 */
// Base Class Headers
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

namespace reco {class Track;}

class TFile;
class TH1F;
class TH2F;
class Trajectory;


#include <vector>

class TrajectoryReader: public edm::EDAnalyzer {

 public:
  typedef std::vector<Trajectory> Trajectories;

 public:
  /// Constructor
  TrajectoryReader(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~TrajectoryReader();

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);

  // Operations
  void beginJob(const edm::EventSetup&);
  void endJob();
  
protected:
  void printTrajectoryRecHits(const Trajectory &, edm::ESHandle<GlobalTrackingGeometry>) const;
  void printTrackRecHits(const reco::Track &, edm::ESHandle<GlobalTrackingGeometry>) const;
  
 private:
  
  edm::InputTag theInputLabel;
  TFile *theFile;
  std::string theRootFileName;
  
  TH1F *hDPtIn;
  TH1F *hDPtOut;
};
#endif
