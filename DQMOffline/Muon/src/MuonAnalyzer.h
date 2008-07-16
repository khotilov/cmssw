#ifndef MuonAnalyzer_H
#define MuonAnalyzer_H


/** \class MuonAnalyzer
 *
 *  DQM muon analysis monitoring
 *
 *  $Date: 2008/06/09 13:41:20 $
 *  $Revision: 1.11 $
 *  \author G. Mila - INFN Torino
 */


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

class MuonEnergyDepositAnalyzer;
class MuonSeedsAnalyzer;
class MuonRecoAnalyzer;
class SegmentTrackAnalyzer;
class DQMStore;
class MuonServiceProxy;

class MuonAnalyzer : public edm::EDAnalyzer {
 public:

  /// Constructor
  MuonAnalyzer(const edm::ParameterSet&);
  
  /// Destructor
  virtual ~MuonAnalyzer();
  
  /// Inizialize parameters for histo binning
  void beginJob(edm::EventSetup const& iSetup);

  /// Get the analysis
  void analyze(const edm::Event&, const edm::EventSetup&);

  /// Save the histos
  void endJob(void);

 private:
  // ----------member data ---------------------------
  
  DQMStore* theDbe;
  edm::ParameterSet parameters;
  MuonServiceProxy *theService;
  // Switch for verbosity
  std::string metname;

  // Muon Label
  edm::InputTag theMuonCollectionLabel;
  // Glb Muon Track label
  edm::InputTag theGlbMuTrackCollectionLabel;
  // Sta Muon Track label
  edm::InputTag theStaMuTrackCollectionLabel;
  // Seed Label
  edm::InputTag theSeedsCollectionLabel;
  
  bool theMuEnergyAnalyzerFlag;
  bool theSeedsAnalyzerFlag;
  bool theMuonRecoAnalyzerFlag;
  bool theMuonSegmentsAnalyzerFlag;

  // the muon energy analyzer
  MuonEnergyDepositAnalyzer * theMuEnergyAnalyzer;
  // the seeds analyzer
  MuonSeedsAnalyzer * theSeedsAnalyzer;
  // the muon reco analyzer
  MuonRecoAnalyzer * theMuonRecoAnalyzer;
  // the track segments analyzer for glb muons
  SegmentTrackAnalyzer * theGlbMuonSegmentsAnalyzer;
  // the track segments analyzer for sta muons
  SegmentTrackAnalyzer * theStaMuonSegmentsAnalyzer;

};
#endif  
