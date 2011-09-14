#ifndef METAnalyzer_H
#define METAnalyzer_H


/** \class METAnalyzer
 *
 *  DQM monitoring source for MET (Mu corrected/TcMET)
 *
 *  $Date: 2010/09/22 19:40:34 $
 *  $Revision: 1.21 $
 *  \author A.Apresyan - Caltech
 */


#include <memory>
#include <fstream>
#include "TMath.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMOffline/JetMET/interface/METAnalyzerBase.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
//
#include "CommonTools/TriggerUtils/interface/GenericTriggerEventFlag.h"
//
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "RecoMET/METAlgorithms/interface/HcalNoiseRBXArray.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include "RecoJets/JetProducers/interface/JetIDHelper.h"

#include "DataFormats/MuonReco/interface/MuonMETCorrectionData.h"
#include "DataFormats/Common/interface/ValueMap.h"  
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonMETCorrectionData.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DQMOffline/JetMET/interface/JetMETDQMDCSFilter.h"

class METAnalyzer : public METAnalyzerBase {
 public:

  /// Constructor
  METAnalyzer(const edm::ParameterSet&);
  
  /// Destructor
  virtual ~METAnalyzer();

  /// Inizialize parameters for histo binning
  void beginJob(DQMStore * dbe);

  /// Finish up a job
  void endJob();

  // Book MonitorElements
  void bookMESet(std::string);
  void bookMonitorElement(std::string, bool);

  /// Get the analysis
  void analyze(const edm::Event&, const edm::EventSetup&, 
               const edm::TriggerResults&);

  /// Initialize run-based parameters
  void beginRun(const edm::Run&,  const edm::EventSetup&);

  /// Finish up a run
  void endRun(const edm::Run& iRun, const edm::EventSetup& iSetup, DQMStore *dbe);

  // Fill MonitorElements
  void fillMESet(const edm::Event&, std::string, const reco::MET&);
  void fillMonitorElement(const edm::Event&, std::string, std::string, const reco::MET&, bool);
  void makeRatePlot(std::string, double);

  bool selectHighPtJetEvent(const edm::Event&);
  bool selectLowPtJetEvent(const edm::Event&);
  bool selectWElectronEvent(const edm::Event&);
  bool selectWMuonEvent(const edm::Event&);

  void setSource(std::string source) {
    _source = source;
  }

  int evtCounter;

 private:
  // ----------member data ---------------------------
  
  edm::ParameterSet parameters;
  // Switch for verbosity
  int _verbose;

  std::string metname;
  std::string _source;

  std::string _FolderName;

  edm::InputTag theMETCollectionLabel;
  edm::InputTag HcalNoiseRBXCollectionTag;
  edm::InputTag HcalNoiseSummaryTag;
  edm::InputTag theJetCollectionLabel;
  edm::InputTag thePfJetCollectionLabel;
  edm::InputTag TcCandidatesTag;
  edm::InputTag BeamHaloSummaryTag;
  edm::InputTag vertexTag;
  edm::InputTag gtTag;

  edm::InputTag inputTrackLabel;
  edm::InputTag inputMuonLabel;
  edm::InputTag inputElectronLabel;
  edm::InputTag inputBeamSpotLabel;


  // list of Jet or MB HLT triggers
  std::vector<std::string > HLTPathsJetMBByName_;

  GenericTriggerEventFlag * _HighPtJetEventFlag;
  GenericTriggerEventFlag * _LowPtJetEventFlag;
  GenericTriggerEventFlag * _MinBiasEventFlag;
  GenericTriggerEventFlag * _HighMETEventFlag;
  GenericTriggerEventFlag * _LowMETEventFlag;
  GenericTriggerEventFlag * _EleEventFlag;
  GenericTriggerEventFlag * _MuonEventFlag;

  std::string _hlt_HighPtJet;
  std::string _hlt_LowPtJet;
  std::string _hlt_MinBias;
  std::string _hlt_HighMET;
  std::string _hlt_LowMET;
  std::string _hlt_Ele;
  std::string _hlt_Muon;

  edm::ParameterSet theCleaningParameters;
  std::string _hlt_PhysDec;

  std::vector<unsigned > _techTrigsAND;
  std::vector<unsigned > _techTrigsOR;
  std::vector<unsigned > _techTrigsNOT;

  bool _doPVCheck;
  bool _doHLTPhysicsOn;

  bool     _tightBHFiltering;
  int      _tightJetIDFiltering;
  bool     _tightHcalFiltering;

  int _nvtx_min;
  int _nvtxtrks_min;
  int _vtxndof_min;
  double _vtxchi2_max;
  double _vtxz_max;
  
  int _trig_JetMB;
  int _trig_HighPtJet;
  int _trig_LowPtJet;
  int _trig_MinBias;
  int _trig_HighMET;
  int _trig_LowMET;
  int _trig_Ele;
  int _trig_Muon;
  int _trig_PhysDec;


  double _highPtJetThreshold;
  double _lowPtJetThreshold;
  double _highMETThreshold;
  double _lowMETThreshold;

  // Et threshold for MET plots
  double _etThreshold;

  // HF calibration factor (in 31X applied by TcProducer)
  double hfCalibFactor_;  //

  // JetID helper
  reco::helper::JetIDHelper *jetID;


  // DCS filter
  JetMETDQMDCSFilter *DCSFilter;

  //
  bool _allhist;
  bool _allSelection;
  bool _cleanupSelection;

  //
  std::vector<std::string> _FolderNames;

  //
  math::XYZPoint bspot;

  edm::Handle< edm::ValueMap<reco::MuonMETCorrectionData> > tcMet_ValueMap_Handle;
  edm::Handle< reco::MuonCollection > muon_h;
  edm::Handle< edm::View<reco::Track> > track_h;
  edm::Handle< edm::View<reco::GsfElectron > > electron_h;
  edm::Handle< reco::BeamSpot > beamSpot_h;

  //
  DQMStore *_dbe;

  //trigger histos
  MonitorElement* hTriggerName_HighPtJet;
  MonitorElement* hTriggerName_LowPtJet;
  MonitorElement* hTriggerName_MinBias;
  MonitorElement* hTriggerName_HighMET;
  MonitorElement* hTriggerName_LowMET;
  MonitorElement* hTriggerName_Ele;
  MonitorElement* hTriggerName_Muon;

  //the histos
  MonitorElement* hMETRate;

  MonitorElement* hmetME;
  //removed for optimizations//MonitorElement* hNevents;
  MonitorElement* hMEx;
  MonitorElement* hMEy;
  //removed for optimizations//MonitorElement* hEz;
  MonitorElement* hMETSig;
  MonitorElement* hMET;
  MonitorElement* hMETPhi;
  MonitorElement* hSumET;

  MonitorElement* hMET_logx;
  MonitorElement* hSumET_logx;

  //removed for optimizations//MonitorElement* hMETIonFeedbck;
  //removed for optimizations//MonitorElement* hMETHPDNoise;
  //removed for optimizations//MonitorElement* hMETRBXNoise;

  MonitorElement* hMExLS;
  MonitorElement* hMEyLS;

  MonitorElement* htrkPt;
  MonitorElement* htrkEta;
  MonitorElement* htrkNhits;
  MonitorElement* htrkChi2;
  MonitorElement* htrkD0;
  MonitorElement* helePt;
  MonitorElement* heleEta;
  MonitorElement* heleHoE;
  MonitorElement* hmuPt;
  MonitorElement* hmuEta;
  MonitorElement* hmuNhits;
  MonitorElement* hmuChi2;
  MonitorElement* hmuD0;
  
  MonitorElement* hMExCorrection;
  MonitorElement* hMEyCorrection;
  MonitorElement* hMuonCorrectionFlag;


};
#endif
