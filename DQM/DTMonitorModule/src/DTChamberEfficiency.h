#ifndef DTCHAMBEREFFICIENCY_H
#define DTCHAMBEREFFICIENCY_H

/** \class DTChamberEfficiency
 *
 * Description:
 *  
 * This class provides the histograms for the calculation of the
 * efficiency of muons reconstruction in the DTs. It is applicable
 * both in presence or absence of a magnetic field.
 * Histos are 2D Sector vs Chamber plots for each wheel
 *
 * \author : Mario Pelliccioni - INFN Torino <pellicci@cern.ch>
 * $date   : 05/12/2008 16:51:04 CET $
 *
 * Modification:
 *
 */

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h" 
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h" 
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include <iosfwd>
#include <bitset>
#include <vector>

class MuonDetLayerMeasurements;
class TransientTrack;
class DQMStore;
class MonitorElement;

class DTChamberEfficiency : public edm::EDAnalyzer
{

 public:
  //Constructor 
  DTChamberEfficiency(const edm::ParameterSet& pset) ;

  //Destructor
  ~DTChamberEfficiency() ;

  //Operations
  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  void beginJob(const edm::EventSetup&);
  void beginRun(const edm::Run& , const edm::EventSetup&);
  void endJob();

 private:

  //functions
  std::vector<const DetLayer*> compatibleLayers(const DetLayer *initialLayer,
						const FreeTrajectoryState& fts, PropagationDirection propDir);


  void bookHistos();
  MeasurementContainer segQualityCut(const MeasurementContainer seg_list) const;
  bool chamberSelection(const DetId& idDetLay, reco::TransientTrack& trans_track) const;
  inline edm::ESHandle<Propagator> propagator() const;

  //data members
  bool debug;

  edm::InputTag theTracksLabel;

  edm::InputTag labelRPCRecHits;
  edm::InputTag thedt4DSegments;
  edm::InputTag thecscSegments;

  double theMaxChi2;
  double theNSigma;
  int theMinNrec;

  std::string theNavigationType;

  edm::ESHandle<DTGeometry> dtGeom;

  DQMStore* theDbe;

  MuonServiceProxy* theService;
  MuonDetLayerMeasurements* theMeasurementExtractor;
  Chi2MeasurementEstimator* theEstimator;

  edm::ESHandle<MagneticField> magfield;
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;

  std::map<DTChamberId, std::vector<MonitorElement*> > histosPerW;

 protected:

};

#endif // DTANALYZER_H
