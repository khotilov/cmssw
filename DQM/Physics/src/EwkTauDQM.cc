#include "DQM/Physics/src/EwkTauDQM.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"

const std::string dqmSeparator = "/";

std::string dqmDirectoryName(const std::string& dqmRootDirectory, const std::string& dqmSubDirectory)
{
//--- concatenate names of dqmRootDirectory and dqmSubDirectory;
//    add "/" separator inbetween if necessary
  std::string dirName = dqmRootDirectory;
  if ( dirName != "" && dirName.find_last_of(dqmSeparator) != (dirName.length() - 1) )  dirName.append(dqmSeparator);
  dirName.append(dqmSubDirectory);
  return dirName;
}

EwkTauDQM::EwkTauDQM(const edm::ParameterSet& cfg)
  : dqmDirectory_(cfg.getParameter<std::string>("dqmDirectory")),
    dqmError_(0)
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("EwkTauDQM") << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!";
    dqmError_ = 1;
    return;
  }

  DQMStore* dqmStore = &(*edm::Service<DQMStore>());

  edm::ParameterSet cfgChannels = cfg.getParameter<edm::ParameterSet>("channels");

  edm::ParameterSet cfgElecTauChannel = cfgChannels.getParameter<edm::ParameterSet>("elecTauChannel");
  std::string dqmSubDirectoryElecTauChannel = cfgElecTauChannel.getParameter<std::string>("dqmSubDirectory");
  cfgElecTauChannel.addParameter<std::string>("dqmDirectory", dqmDirectoryName(dqmDirectory_, dqmSubDirectoryElecTauChannel));
  elecTauHistManager_ = new EwkElecTauHistManager(cfgElecTauChannel, dqmStore);

  edm::ParameterSet cfgMuTauChannel = cfgChannels.getParameter<edm::ParameterSet>("muTauChannel");
  std::string dqmSubDirectoryMuTauChannel = cfgMuTauChannel.getParameter<std::string>("dqmSubDirectory");
  cfgMuTauChannel.addParameter<std::string>("dqmDirectory", dqmDirectoryName(dqmDirectory_, dqmSubDirectoryMuTauChannel));
  muTauHistManager_ = new EwkMuTauHistManager(cfgMuTauChannel, dqmStore);
}

EwkTauDQM::~EwkTauDQM()
{
  delete elecTauHistManager_;
  delete muTauHistManager_;
}

void EwkTauDQM::beginJob()
{
  if ( dqmError_ ) return;
 
  elecTauHistManager_->bookHistograms();
  muTauHistManager_->bookHistograms();
}

void EwkTauDQM::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  if ( dqmError_ ) return;

  elecTauHistManager_->fillHistograms(evt, es);
  muTauHistManager_->fillHistograms(evt, es);
}

void EwkTauDQM::endJob()
{
  if ( dqmError_ ) return;

  elecTauHistManager_->finalizeHistograms();
  muTauHistManager_->finalizeHistograms();
}

//-------------------------------------------------------------------------------
// code specific to Z --> e + tau-jet channel
//-------------------------------------------------------------------------------

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TMath.h"

#include <iostream>
#include <iomanip>

EwkElecTauHistManager::EwkElecTauHistManager(const edm::ParameterSet& cfg, DQMStore* dqmStore) 
  : dqmStore_(dqmStore),
    dqmDirectory_(cfg.getParameter<std::string>("dqmDirectory")),
    numEventsAnalyzed_(0),
    numEventsSelected_(0),
    cfgError_(0)
{
  triggerResultsSource_ = cfg.getParameter<edm::InputTag>("triggerResultsSource");
  vertexSource_ = cfg.getParameter<edm::InputTag>("vertexSource");
  beamSpotSource_ = cfg.getParameter<edm::InputTag>("beamSpotSource");
  electronSource_ = cfg.getParameter<edm::InputTag>("electronSource");
  tauJetSource_ = cfg.getParameter<edm::InputTag>("tauJetSource");
  caloMEtSource_ = cfg.getParameter<edm::InputTag>("caloMEtSource");
  pfMEtSource_ = cfg.getParameter<edm::InputTag>("pfMEtSource");

  tauDiscrByLeadTrackFinding_ = cfg.getParameter<edm::InputTag>("tauDiscrByLeadTrackFinding");
  tauDiscrByLeadTrackPtCut_ = cfg.getParameter<edm::InputTag>("tauDiscrByLeadTrackPtCut");
  tauDiscrByTrackIso_ = cfg.getParameter<edm::InputTag>("tauDiscrByTrackIso");
  tauDiscrByEcalIso_ = cfg.getParameter<edm::InputTag>("tauDiscrByEcalIso");
  tauDiscrAgainstElectrons_ = cfg.getParameter<edm::InputTag>("tauDiscrAgainstElectrons");
  tauDiscrAgainstMuons_ = cfg.getParameter<edm::InputTag>("tauDiscrAgainstMuons");

  hltPaths_ = cfg.getParameter<vstring>("hltPaths");

  electronEtaCut_ = cfg.getParameter<double>("electronEtaCut");
  electronPtCut_ = cfg.getParameter<double>("electronPtCut");
  electronTrackIsoCut_ = cfg.getParameter<double>("electronTrackIsoCut");
  electronEcalIsoCut_ = cfg.getParameter<double>("electronEcalIsoCut");
  std::string electronIsoMode_string = cfg.getParameter<std::string>("electronIsoMode");
  electronIsoMode_ = getIsoMode(electronIsoMode_string, cfgError_);

  tauJetEtaCut_ = cfg.getParameter<double>("tauJetEtaCut");
  tauJetPtCut_ = cfg.getParameter<double>("tauJetPtCut");

  visMassCut_ = cfg.getParameter<double>("visMassCut");
}

void EwkElecTauHistManager::bookHistograms()
{
  dqmStore_->setCurrentFolder(dqmDirectory_);

  //hNumIdElectrons_ = dqmStore_->book1D("NumIdElectronsMuons" , "Num. id. Muons", 5, -0.5, 4.5);
  hElectronPt_ = dqmStore_->book1D("ElectronPt" , "P_{T}^{e}", 20, 0., 100.);
  hElectronEta_ = dqmStore_->book1D("ElectronEta" , "#eta_{e}", 20, -4.0, +4.0);
  hElectronPhi_ = dqmStore_->book1D("ElectronPhi" , "#phi_{e}", 20, -TMath::Pi(), +TMath::Pi());
  hElectronTrackIsoPt_ = dqmStore_->book1D("ElectronTrackIsoPt" , "Electron Track Iso.", 20, -0.01, 10.);
  hElectronEcalIsoPt_ = dqmStore_->book1D("ElectronEcalIsoPt" , "Electron Ecal Iso.", 20, -0.01, 10.);

  //hTauJetPt_ = dqmStore_->book1D("TauJetPt" , "P_{T}^{#tau-Jet}", 20, 0., 100.);
  //hTauJetEta_ = dqmStore_->book1D("TauJetEta" , "#eta_{#tau-Jet}", 20, -4.0, +4.0);
  //hTauJetPhi_ = dqmStore_->book1D("TauJetPhi" , "#phi_{#tau-Jet}", 20, -TMath::Pi(), +TMath::Pi());
  //hTauLeadTrackPt_ = dqmStore_->book1D("TauLeadTrackPt" , "P_{T}^{#tau-Jet}", 20, 0., 50.);
  //hTauTrackIsoPt_ = dqmStore_->book1D("TauTrackIsoPt" , "Tau Track Iso.", 20, -0.01, 40.);
  //hTauEcalIsoPt_ = dqmStore_->book1D("TauEcalIsoPt" , "Tau Ecal Iso.", 10, -0.01, 10.);
  //hTauDiscrAgainstElectrons_ = dqmStore_->book1D("TauDiscrAgainstElectrons" , "Tau Discr. against Electrons", 2, -0.5, +1.5);
  //hTauDiscrAgainstMuons_ = dqmStore_->book1D("TauDiscrAgainstMuons" , "Tau Discr. against Muons", 2, -0.5, +1.5);
  //hTauJetCharge_ = dqmStore_->book1D("TauJetCharge" , "Q_{#tau-Jet}", 11, -5.5, +5.5);
  //hTauJetNumSignalTracks_ = dqmStore_->book1D("TauJetNumSignalTracks" , "Num. Tau signal Cone Tracks", 20, -0.5, +19.5);
  //hTauJetNumIsoTracks_ = dqmStore_->book1D("TauJetNumIsoTracks" , "Num. Tau isolation Cone Tracks", 20, -0.5, +19.5);
  
  hVisMass_ = dqmStore_->book1D("VisMass", "e + #tau-Jet visible Mass", 20, 20., 120.);
  //hMtElecCaloMEt_ = dqmStore_->book1D("MtElecCaloMEt", "e + E_{T}^{miss} (Calo) transverse Mass", 20, 20., 120.);
  hMtElecPFMEt_ = dqmStore_->book1D("MtElecPFMEt", "e + E_{T}^{miss} (PF) transverse Mass", 20, 20., 120.);
  //hPzetaCaloMEt_ = dqmStore_->book1D("PzetaCaloMEt", "P_{#zeta} - 1.5*P_{#zeta}^{vis} (Calo)", 20, -40., 40.);
  //hPzetaPFMEt_ = dqmStore_->book1D("PzetaPFMEt", "P_{#zeta} - 1.5*P_{#zeta}^{vis} (PF)", 20, -40., 40.);
  hElecTauAcoplanarity_ = dqmStore_->book1D("ElecTauAcoplanarity", "#Delta #phi_{e #tau-Jet}", 20, -TMath::Pi(), +TMath::Pi());
  //hElecTauCharge_ = dqmStore_->book1D("ElecTauCharge" , "Q_{e + #tau-Jet}", 11, -5.5, +5.5);

  //hVertexChi2_ = dqmStore_->book1D("VertexChi2", "Event Vertex #chi^{2} / n.d.o.f.", 20, 0., 2.0);
  hVertexZ_ = dqmStore_->book1D("VertexZ", "Event Vertex z-Position", 20, -25., +25.);
  //hVertexD0_ = dqmStore_->book1D("VertexD0", "Event Vertex d_{0}", 20, -0.0001, 0.05);

  hCaloMEtPt_ = dqmStore_->book1D("CaloMEtPt", "E_{T}^{miss} (Calo)", 20, 0., 100.);
  //hCaloMEtPhi_ = dqmStore_->book1D("CaloMEtPhi", "#phi^{miss} (Calo)", 20, -TMath::Pi(), +TMath::Pi());

  hPFMEtPt_ = dqmStore_->book1D("PFMEtPt", "E_{T}^{miss} (PF)", 20, 0., 100.);
  //hPFMEtPhi_ = dqmStore_->book1D("PFMEtPhi", "#phi^{miss} (PF)", 20, -TMath::Pi(), +TMath::Pi());

  hCutFlowSummary_ = dqmStore_->book1D("CutFlowSummary", "Cut-flow Summary", 11, 0.5, 11.5);
  hCutFlowSummary_->setBinLabel(kPassedPreselection, "Preselection");
  hCutFlowSummary_->setBinLabel(kPassedTrigger, "HLT");
  hCutFlowSummary_->setBinLabel(kPassedElectronId, "e ID");
  hCutFlowSummary_->setBinLabel(kPassedElectronTrackIso, "e Trk Iso.");
  hCutFlowSummary_->setBinLabel(kPassedElectronEcalIso, "e Ecal Iso.");
  hCutFlowSummary_->setBinLabel(kPassedTauLeadTrack, "#tau lead. Track");
  hCutFlowSummary_->setBinLabel(kPassedTauLeadTrackPt, "#tau lead. Track P_{T}");
  hCutFlowSummary_->setBinLabel(kPassedTauTrackIso, "#tau Track Iso.");
  hCutFlowSummary_->setBinLabel(kPassedTauEcalIso, "#tau Ecal Iso.");
  hCutFlowSummary_->setBinLabel(kPassedTauDiscrAgainstElectrons, "#tau anti-e Discr.");
  hCutFlowSummary_->setBinLabel(kPassedTauDiscrAgainstMuons, "#tau anti-#mu Discr.");
}

void EwkElecTauHistManager::fillHistograms(const edm::Event& evt, const edm::EventSetup& es)
{
  if ( cfgError_ ) return;

  //-----------------------------------------------------------------------------
  // access event-level information
  //-----------------------------------------------------------------------------

//--- get decision of high-level trigger for the event
  edm::Handle<edm::TriggerResults> hltDecision;
  if ( !evt.getByLabel(triggerResultsSource_, hltDecision) ) {
    edm::LogWarning ("EwkElecTauHistManager") << "Failed to access Trigger results !!";
    return;
  }
  
  edm::TriggerNames triggerNames;
  triggerNames.init(*hltDecision);
   
  bool isTriggered = false;
  for ( vstring::const_iterator hltPath = hltPaths_.begin();
	hltPath != hltPaths_.end(); ++hltPath ) {
    unsigned int index = triggerNames.triggerIndex(*hltPath);
    if ( index < triggerNames.size() ) {
      if ( hltDecision->accept(index) ) isTriggered = true;
    } else {
      edm::LogWarning ("EwkElecTauHistManager") << " Undefined HLT path = " << (*hltPath) << " !!";
      continue;
    }
  }
  
//--- get reconstructed primary event vertex of the event
//   (take as "the" primary event vertex the first entry in the collection
//    of vertex objects, corresponding to the vertex associated to the highest Pt sum of tracks)
  edm::Handle<reco::VertexCollection> vertexCollection;
  if ( !evt.getByLabel(vertexSource_, vertexCollection) ) {
    edm::LogWarning ("EwkElecTauHistManager") << "Failed to access Vertex collection !!";
    return;
  }
  const reco::Vertex* theEventVertex = ( vertexCollection->size() > 0 ) ? &(vertexCollection->at(0)) : 0;

//--- get beam-spot (expected vertex position) for the event
  edm::Handle<reco::BeamSpot> beamSpot;
  if ( !evt.getByLabel(beamSpotSource_, beamSpot) ) {
    edm::LogWarning ("EwkElecTauHistManager") << "Failed to access Beam-spot !!";
    return;
  }
  
//--- get collections of reconstructed electrons from the event
  edm::Handle<reco::GsfElectronCollection> electrons;
  if ( !evt.getByLabel(electronSource_, electrons) ) {
    edm::LogWarning ("EwkElecTauHistManager") << "Failed to access Electron collection !!";
    return;
  }

  const reco::GsfElectron* theElectron = getTheElectron(*electrons, electronEtaCut_, electronPtCut_);

  double theElectronTrackIsoPt = 1.e+3;
  double theElectronEcalIsoPt = 1.e+3;
  if ( theElectron ) {
    theElectronTrackIsoPt = theElectron->dr03TkSumPt();
    theElectronEcalIsoPt = theElectron->dr03EcalRecHitSumEt();

    if ( electronIsoMode_ == kRelativeIso && theElectron->pt() > 0. ) {
      theElectronTrackIsoPt /= theElectron->pt();
      theElectronEcalIsoPt /= theElectron->pt();
    }
  }

//--- get collections of reconstructed tau-jets from the event
  edm::Handle<reco::PFTauCollection> tauJets;
  if ( !evt.getByLabel(tauJetSource_, tauJets) ) {
    edm::LogWarning ("EwkElecTauHistManager") << "Failed to access Tau-jet collection !!";
    return;
  }

//--- get collections of tau-jet discriminators for those tau-jets
  edm::Handle<reco::PFTauDiscriminator> tauDiscrByLeadTrackFinding;
  evt.getByLabel(tauDiscrByLeadTrackFinding_, tauDiscrByLeadTrackFinding);
  
  edm::Handle<reco::PFTauDiscriminator> tauDiscrByLeadTrackPtCut;
  evt.getByLabel(tauDiscrByLeadTrackPtCut_, tauDiscrByLeadTrackPtCut);

  edm::Handle<reco::PFTauDiscriminator> tauDiscrByTrackIso;
  evt.getByLabel(tauDiscrByTrackIso_, tauDiscrByTrackIso);

  edm::Handle<reco::PFTauDiscriminator> tauDiscrByEcalIso;
  evt.getByLabel(tauDiscrByEcalIso_, tauDiscrByEcalIso);

  edm::Handle<reco::PFTauDiscriminator> tauDiscrAgainstElectrons;
  evt.getByLabel(tauDiscrAgainstElectrons_, tauDiscrAgainstElectrons);

  edm::Handle<reco::PFTauDiscriminator> tauDiscrAgainstMuons;
  evt.getByLabel(tauDiscrAgainstMuons_, tauDiscrAgainstMuons);

  int theTauJetIndex = -1;
  const reco::PFTau* theTauJet = getTheTauJet(*tauJets, tauJetEtaCut_, tauJetPtCut_, theTauJetIndex);

  double theTauDiscrByLeadTrackFinding = -1.;
  double theTauDiscrByLeadTrackPtCut = -1.;
  double theTauDiscrByTrackIso = -1.;
  double theTauDiscrByEcalIso = -1.;
  double theTauDiscrAgainstElectrons = -1.;
  double theTauDiscrAgainstMuons = -1.;
  if ( theTauJetIndex != -1 ) {
    reco::PFTauRef theTauJetRef(tauJets, theTauJetIndex);
    theTauDiscrByLeadTrackFinding = (*tauDiscrByLeadTrackFinding)[theTauJetRef];
    theTauDiscrByLeadTrackPtCut = (*tauDiscrByLeadTrackPtCut)[theTauJetRef];
    theTauDiscrByTrackIso = (*tauDiscrByTrackIso)[theTauJetRef];
    theTauDiscrByEcalIso = (*tauDiscrByEcalIso)[theTauJetRef];
    theTauDiscrAgainstElectrons = (*tauDiscrAgainstElectrons)[theTauJetRef];
    theTauDiscrAgainstMuons = (*tauDiscrAgainstMuons)[theTauJetRef];
  }

//--- get missing transverse momentum
//    measured by calorimeters/reconstructed by particle-flow algorithm
  edm::Handle<reco::CaloMETCollection> caloMEtCollection;
  if ( !evt.getByLabel(caloMEtSource_, caloMEtCollection) ) {
    edm::LogWarning ("EwkElecTauHistManager") << "Failed to access calo. MET collection !!";
    return;
  }
  const reco::CaloMET& caloMEt = caloMEtCollection->at(0);
  
  edm::Handle<reco::PFMETCollection> pfMEtCollection;
  if ( !evt.getByLabel(pfMEtSource_, pfMEtCollection) ) {
    edm::LogWarning ("EwkElecTauHistManager") << "Failed to access pf. MET collection !!";
    return;
  }
  const reco::PFMET& pfMEt = pfMEtCollection->at(0);

  if ( !(theElectron && theTauJet && theTauJetIndex != -1) ) return;

  //-----------------------------------------------------------------------------
  // compute EWK tau analysis specific quantities
  //-----------------------------------------------------------------------------

  double dPhiElecTau = calcDeltaPhi(theElectron->phi(), theTauJet->phi());

  double mElecTau = (theElectron->p4() + theTauJet->p4()).M();

  //double mtElecCaloMEt = calcMt(theElectron->px(), theElectron->px(), caloMEt.px(), caloMEt.py());
  double mtElecPFMEt = calcMt(theElectron->px(), theElectron->px(), pfMEt.px(), pfMEt.py());

  //double pZetaCaloMEt = calcPzeta(theElectron->p4(), theTauJet->p4(), caloMEt.px(), caloMEt.py());
  //double pZetaPFMEt = calcPzeta(theElectron->p4(), theTauJet->p4(), pfMEt.px(), pfMEt.py());

  //-----------------------------------------------------------------------------
  // apply selection criteria; fill histograms
  //-----------------------------------------------------------------------------

//--- fill electron multiplicity histogram
  unsigned numIdElectrons = 0;
  for ( reco::GsfElectronCollection::const_iterator electron = electrons->begin();
	electron != electrons->end(); ++electron ) {
    if ( passesElectronId(*electron) ) {
      ++numIdElectrons;
    }
  }

  //hNumIdElectrons_->Fill(numIdElectrons);

  ++numEventsAnalyzed_;

  bool isSelected = false;
  int cutFlowStatus = -1;

  if ( mElecTau > visMassCut_ ) {
    cutFlowStatus = kPassedPreselection;
  }
  if ( cutFlowStatus == kPassedPreselection && (isTriggered || hltPaths_.size() == 0) ) {
    cutFlowStatus = kPassedTrigger;
  }
  if ( cutFlowStatus == kPassedTrigger && passesElectronId(*theElectron) ) {
    cutFlowStatus = kPassedElectronId;
    hElectronTrackIsoPt_->Fill(theElectronTrackIsoPt);
  }
  if ( cutFlowStatus == kPassedElectronId && theElectronTrackIsoPt < electronTrackIsoCut_ ) {
    cutFlowStatus = kPassedElectronTrackIso;
    hElectronEcalIsoPt_->Fill(theElectronEcalIsoPt);
  }
  if ( cutFlowStatus == kPassedElectronTrackIso && theElectronEcalIsoPt < electronEcalIsoCut_ ) {
    cutFlowStatus = kPassedElectronEcalIso;
  }
  if ( cutFlowStatus == kPassedElectronEcalIso && theTauDiscrByLeadTrackFinding > 0.5 ) {
    cutFlowStatus = kPassedTauLeadTrack;
    if ( theTauJet->leadTrack().isAvailable() ) {
      //hTauLeadTrackPt_->Fill(theTauJet->leadTrack()->pt());
    } else {
      edm::LogWarning ("EwkElecTauHistManager") << "Failed to access lead. Track of tau-jet !!";
    }
  }
  if ( cutFlowStatus == kPassedTauLeadTrack && theTauDiscrByLeadTrackPtCut > 0.5 ) {
    cutFlowStatus = kPassedTauLeadTrackPt;
    //hTauTrackIsoPt_->Fill(theTauJet->isolationPFChargedHadrCandsPtSum());
  }
  if ( cutFlowStatus == kPassedTauLeadTrackPt && theTauDiscrByTrackIso > 0.5 ) {
    cutFlowStatus = kPassedTauTrackIso;
    //hTauEcalIsoPt_->Fill(theTauJet->isolationPFGammaCandsEtSum());
  }
  if ( cutFlowStatus == kPassedTauTrackIso && theTauDiscrByEcalIso > 0.5 ) {
    cutFlowStatus = kPassedTauEcalIso;
    //hTauDiscrAgainstElectrons_->Fill(theTauDiscrAgainstElectrons);
  }
  if ( cutFlowStatus == kPassedTauEcalIso && theTauDiscrAgainstElectrons > 0.5 ) {
    cutFlowStatus = kPassedTauDiscrAgainstElectrons;
    //hTauDiscrAgainstMuons_->Fill(theTauDiscrAgainstMuons);
  }
  if ( cutFlowStatus == kPassedTauDiscrAgainstElectrons && theTauDiscrAgainstMuons > 0.5 ) {
    cutFlowStatus = kPassedTauDiscrAgainstMuons;
    isSelected = true;
  }

  for ( int iCut = 1; iCut <= cutFlowStatus; ++iCut ) {
    hCutFlowSummary_->Fill(iCut);
  }

  if ( isSelected ) {
    hElectronPt_->Fill(theElectron->pt());
    hElectronEta_->Fill(theElectron->eta());
    hElectronPhi_->Fill(theElectron->phi());

    //hTauJetPt_->Fill(theTauJet->pt());
    //hTauJetEta_->Fill(theTauJet->eta());
    //hTauJetPhi_->Fill(theTauJet->phi());

    //hTauJetCharge_->Fill(theTauJet->charge());
    if ( theTauJet->signalTracks().isAvailable() ) {
      //hTauJetNumSignalTracks_->Fill(theTauJet->signalTracks().size());
    } else {
      edm::LogWarning ("EwkElecTauHistManager") << "Failed to access signal Tracks associated to tau-jet !!";
    }
    if ( theTauJet->isolationTracks().isAvailable() ) {
      //hTauJetNumIsoTracks_->Fill(theTauJet->isolationTracks().size());
    } else {
      edm::LogWarning ("EwkElecTauHistManager") << "Failed to access isolation Tracks associated to tau-jet !!";
    }
  
    hVisMass_->Fill(mElecTau);
    //hMtElecCaloMEt_->Fill(mtElecCaloMEt);
    hMtElecPFMEt_->Fill(mtElecPFMEt);
    //hPzetaCaloMEt_->Fill(pZetaCaloMEt);
    //hPzetaPFMEt_->Fill(pZetaPFMEt);
    hElecTauAcoplanarity_->Fill(dPhiElecTau);
    //hElecTauCharge_->Fill(theElectron->charge() + theTauJet->charge());

    if ( theEventVertex ) {
      //hVertexChi2_->Fill(theEventVertex->normalizedChi2());
      hVertexZ_->Fill(theEventVertex->z());
      //hVertexD0_->Fill(getVertexD0(*theEventVertex, *beamSpot));
    }
    
    hCaloMEtPt_->Fill(caloMEt.pt());
    //hCaloMEtPhi_->Fill(caloMEt.phi());

    hPFMEtPt_->Fill(pfMEt.pt());
    //hPFMEtPhi_->Fill(pfMEt.phi());
  }

  if ( isSelected ) {
    ++numEventsSelected_;
    LogTrace ("EwkElecTauHistManager") << "Event accepted.";
  } else {
    LogTrace ("EwkElecTauHistManager") << "Event rejected.";
  }
}

void EwkElecTauHistManager::finalizeHistograms()
{
  edm::LogInfo ("EwkElecTauHistManager") 
    << "Filter-Statistics Summary:" << std::endl
    << " Events analyzed = " << numEventsAnalyzed_ << std::endl
    << " Events selected = " << numEventsSelected_;
  if ( numEventsAnalyzed_ > 0 ) {
    double eff = numEventsSelected_/(double)numEventsAnalyzed_;
    edm::LogInfo ("") 
      << "Overall efficiency = " << std::setprecision(4) << eff*100. 
      << " +/- " << std::setprecision(4) << TMath::Sqrt(eff*(1 - eff)/numEventsAnalyzed_)*100. << ")%";
  }
}

//-------------------------------------------------------------------------------
// code specific to Z --> mu + tau-jet channel
//-------------------------------------------------------------------------------

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TMath.h"

#include <iostream>
#include <iomanip>

EwkMuTauHistManager::EwkMuTauHistManager(const edm::ParameterSet& cfg, DQMStore* dqmStore) 
  : dqmStore_(dqmStore),
    dqmDirectory_(cfg.getParameter<std::string>("dqmDirectory")),
    numEventsAnalyzed_(0),
    numEventsSelected_(0),
    cfgError_(0)
{
  triggerResultsSource_ = cfg.getParameter<edm::InputTag>("triggerResultsSource");
  vertexSource_ = cfg.getParameter<edm::InputTag>("vertexSource");
  beamSpotSource_ = cfg.getParameter<edm::InputTag>("beamSpotSource");
  muonSource_ = cfg.getParameter<edm::InputTag>("muonSource");
  tauJetSource_ = cfg.getParameter<edm::InputTag>("tauJetSource");
  caloMEtSource_ = cfg.getParameter<edm::InputTag>("caloMEtSource");
  pfMEtSource_ = cfg.getParameter<edm::InputTag>("pfMEtSource");

  tauDiscrByLeadTrackFinding_ = cfg.getParameter<edm::InputTag>("tauDiscrByLeadTrackFinding");
  tauDiscrByLeadTrackPtCut_ = cfg.getParameter<edm::InputTag>("tauDiscrByLeadTrackPtCut");
  tauDiscrByTrackIso_ = cfg.getParameter<edm::InputTag>("tauDiscrByTrackIso");
  tauDiscrByEcalIso_ = cfg.getParameter<edm::InputTag>("tauDiscrByEcalIso");
  tauDiscrAgainstMuons_ = cfg.getParameter<edm::InputTag>("tauDiscrAgainstMuons");

  hltPaths_ = cfg.getParameter<vstring>("hltPaths");

  muonEtaCut_ = cfg.getParameter<double>("muonEtaCut");
  muonPtCut_ = cfg.getParameter<double>("muonPtCut");
  muonTrackIsoCut_ = cfg.getParameter<double>("muonTrackIsoCut");
  muonEcalIsoCut_ = cfg.getParameter<double>("muonEcalIsoCut");
  std::string muonIsoMode_string = cfg.getParameter<std::string>("muonIsoMode");
  muonIsoMode_ = getIsoMode(muonIsoMode_string, cfgError_);

  tauJetEtaCut_ = cfg.getParameter<double>("tauJetEtaCut");
  tauJetPtCut_ = cfg.getParameter<double>("tauJetPtCut");

  visMassCut_ = cfg.getParameter<double>("visMassCut");
}

void EwkMuTauHistManager::bookHistograms()
{
  dqmStore_->setCurrentFolder(dqmDirectory_);

  //hNumGlobalMuons_ = dqmStore_->book1D("NumGlobalMuons" , "Num. global Muons", 5, -0.5, 4.5);
  hMuonPt_ = dqmStore_->book1D("MuonPt" , "P_{T}^{#mu}", 20, 0., 100.);
  hMuonEta_ = dqmStore_->book1D("MuonEta" , "#eta_{#mu}", 20, -4.0, +4.0);
  hMuonPhi_ = dqmStore_->book1D("MuonPhi" , "#phi_{#mu}", 20, -TMath::Pi(), +TMath::Pi());
  hMuonTrackIsoPt_ = dqmStore_->book1D("MuonTrackIsoPt" , "Muon Track Iso.", 20, -0.01, 10.);
  hMuonEcalIsoPt_ = dqmStore_->book1D("MuonEcalIsoPt" , "Muon Ecal Iso.", 20, -0.01, 10.);

  hTauJetPt_ = dqmStore_->book1D("TauJetPt" , "P_{T}^{#tau-Jet}", 20, 0., 100.);
  hTauJetEta_ = dqmStore_->book1D("TauJetEta" , "#eta_{#tau-Jet}", 20, -4.0, +4.0);
  hTauJetPhi_ = dqmStore_->book1D("TauJetPhi" , "#phi_{#tau-Jet}", 20, -TMath::Pi(), +TMath::Pi());
  hTauLeadTrackPt_ = dqmStore_->book1D("TauLeadTrackPt" , "P_{T}^{#tau-Jet}", 20, 0., 50.);
  hTauTrackIsoPt_ = dqmStore_->book1D("TauTrackIsoPt" , "Tau Track Iso.", 20, -0.01, 40.);
  hTauEcalIsoPt_ = dqmStore_->book1D("TauEcalIsoPt" , "Tau Ecal Iso.", 10, -0.01, 10.);
  hTauDiscrAgainstMuons_ = dqmStore_->book1D("TauDiscrAgainstMuons" , "Tau Discr. against Muons", 2, -0.5, +1.5);
  //hTauJetCharge_ = dqmStore_->book1D("TauJetCharge" , "Q_{#tau-Jet}", 11, -5.5, +5.5);
  hTauJetNumSignalTracks_ = dqmStore_->book1D("TauJetNumSignalTracks" , "Num. Tau signal Cone Tracks", 20, -0.5, +19.5);
  hTauJetNumIsoTracks_ = dqmStore_->book1D("TauJetNumIsoTracks" , "Num. Tau isolation Cone Tracks", 20, -0.5, +19.5);
  
  hVisMass_ = dqmStore_->book1D("VisMass", "#mu + #tau-Jet visible Mass", 20, 20., 120.);
  //hMtMuCaloMEt_ = dqmStore_->book1D("MtMuCaloMEt", "#mu + E_{T}^{miss} (Calo) transverse Mass", 20, 20., 120.);
  hMtMuPFMEt_ = dqmStore_->book1D("MtMuPFMEt", "#mu + E_{T}^{miss} (PF) transverse Mass", 20, 20., 120.);
  //hPzetaCaloMEt_ = dqmStore_->book1D("PzetaCaloMEt", "P_{#zeta} - 1.5*P_{#zeta}^{vis} (Calo)", 20, -40., 40.);
  //hPzetaPFMEt_ = dqmStore_->book1D("PzetaPFMEt", "P_{#zeta} - 1.5*P_{#zeta}^{vis} (PF)", 20, -40., 40.);
  hMuTauAcoplanarity_ = dqmStore_->book1D("MuTauAcoplanarity", "#Delta #phi_{#mu #tau-Jet}", 20, -TMath::Pi(), +TMath::Pi());
  //hMuTauCharge_ = dqmStore_->book1D("MuTauCharge" , "Q_{#mu + #tau-Jet}", 11, -5.5, +5.5);

  //hVertexChi2_ = dqmStore_->book1D("VertexChi2", "Event Vertex #chi^{2} / n.d.o.f.", 20, 0., 2.0);
  hVertexZ_ = dqmStore_->book1D("VertexZ", "Event Vertex z-Position", 20, -25., +25.);
  //hVertexD0_ = dqmStore_->book1D("VertexD0", "Event Vertex d_{0}", 20, -0.0001, 0.05);

  hCaloMEtPt_ = dqmStore_->book1D("CaloMEtPt", "E_{T}^{miss} (Calo)", 20, 0., 100.);
  //hCaloMEtPhi_ = dqmStore_->book1D("CaloMEtPhi", "#phi^{miss} (Calo)", 20, -TMath::Pi(), +TMath::Pi());

  hPFMEtPt_ = dqmStore_->book1D("PFMEtPt", "E_{T}^{miss} (PF)", 20, 0., 100.);
  //hPFMEtPhi_ = dqmStore_->book1D("PFMEtPhi", "#phi^{miss} (PF)", 20, -TMath::Pi(), +TMath::Pi());

  hCutFlowSummary_ = dqmStore_->book1D("CutFlowSummary", "Cut-flow Summary", 10, 0.5, 10.5);
  hCutFlowSummary_->setBinLabel(kPassedPreselection, "Preselection");
  hCutFlowSummary_->setBinLabel(kPassedTrigger, "HLT");
  hCutFlowSummary_->setBinLabel(kPassedMuonId, "#mu ID");
  hCutFlowSummary_->setBinLabel(kPassedMuonTrackIso, "#mu Trk Iso.");
  hCutFlowSummary_->setBinLabel(kPassedMuonEcalIso, "#mu Ecal Iso.");
  hCutFlowSummary_->setBinLabel(kPassedTauLeadTrack, "#tau lead. Track");
  hCutFlowSummary_->setBinLabel(kPassedTauLeadTrackPt, "#tau lead. Track P_{T}");
  hCutFlowSummary_->setBinLabel(kPassedTauTrackIso, "#tau Track Iso.");
  hCutFlowSummary_->setBinLabel(kPassedTauEcalIso, "#tau Ecal Iso.");
  hCutFlowSummary_->setBinLabel(kPassedTauDiscrAgainstMuons, "#tau anti-#mu Discr.");
}

void EwkMuTauHistManager::fillHistograms(const edm::Event& evt, const edm::EventSetup& es)
{
  if ( cfgError_ ) return;

  //-----------------------------------------------------------------------------
  // access event-level information
  //-----------------------------------------------------------------------------

//--- get decision of high-level trigger for the event
  edm::Handle<edm::TriggerResults> hltDecision;
  if ( !evt.getByLabel(triggerResultsSource_, hltDecision) ) {
    edm::LogWarning ("EwkMuTauHistManager") << "Failed to access Trigger results !!";
    return;
  }
  
  edm::TriggerNames triggerNames;
  triggerNames.init(*hltDecision);
   
  bool isTriggered = false;
  for ( vstring::const_iterator hltPath = hltPaths_.begin();
	hltPath != hltPaths_.end(); ++hltPath ) {
    unsigned int index = triggerNames.triggerIndex(*hltPath);
    if ( index < triggerNames.size() ) {
      if ( hltDecision->accept(index) ) isTriggered = true;
    } else {
      edm::LogWarning ("EwkMuTauHistManager") << " Undefined HLT path = " << (*hltPath) << " !!";
      continue;
    }
  }
  
//--- get reconstructed primary event vertex of the event
//   (take as "the" primary event vertex the first entry in the collection
//    of vertex objects, corresponding to the vertex associated to the highest Pt sum of tracks)
  edm::Handle<reco::VertexCollection> vertexCollection;
  if ( !evt.getByLabel(vertexSource_, vertexCollection) ) {
    edm::LogWarning ("EwkMuTauHistManager") << "Failed to access Vertex collection !!";
    return;
  }
  const reco::Vertex* theEventVertex = ( vertexCollection->size() > 0 ) ? &(vertexCollection->at(0)) : 0;

//--- get beam-spot (expected vertex position) for the event
  edm::Handle<reco::BeamSpot> beamSpot;
  if ( !evt.getByLabel(beamSpotSource_, beamSpot) ) {
    edm::LogWarning ("EwkMuTauHistManager") << "Failed to access Beam-spot !!";
    return;
  }
  
//--- get collections of reconstructed muons from the event
  edm::Handle<reco::MuonCollection> muons;
  if ( !evt.getByLabel(muonSource_, muons) ) {
    edm::LogWarning ("EwkMuTauHistManager") << "Failed to access Muon collection !!";
    return;
  }

  const reco::Muon* theMuon = getTheMuon(*muons, muonEtaCut_, muonPtCut_);

  double theMuonTrackIsoPt = 1.e+3;
  double theMuonEcalIsoPt = 1.e+3;
  if ( theMuon ) {
    theMuonTrackIsoPt = theMuon->isolationR05().sumPt;
    theMuonEcalIsoPt = theMuon->isolationR05().emEt;

    if ( muonIsoMode_ == kRelativeIso && theMuon->pt() > 0. ) {
      theMuonTrackIsoPt /= theMuon->pt();
      theMuonEcalIsoPt /= theMuon->pt();
    }
  }

//--- get collections of reconstructed tau-jets from the event
  edm::Handle<reco::PFTauCollection> tauJets;
  if ( !evt.getByLabel(tauJetSource_, tauJets) ) {
    edm::LogWarning ("EwkMuTauHistManager") << "Failed to access Tau-jet collection !!";
    return;
  }

//--- get collections of tau-jet discriminators for those tau-jets
  edm::Handle<reco::PFTauDiscriminator> tauDiscrByLeadTrackFinding;
  evt.getByLabel(tauDiscrByLeadTrackFinding_, tauDiscrByLeadTrackFinding);
  
  edm::Handle<reco::PFTauDiscriminator> tauDiscrByLeadTrackPtCut;
  evt.getByLabel(tauDiscrByLeadTrackPtCut_, tauDiscrByLeadTrackPtCut);

  edm::Handle<reco::PFTauDiscriminator> tauDiscrByTrackIso;
  evt.getByLabel(tauDiscrByTrackIso_, tauDiscrByTrackIso);

  edm::Handle<reco::PFTauDiscriminator> tauDiscrByEcalIso;
  evt.getByLabel(tauDiscrByEcalIso_, tauDiscrByEcalIso);

  edm::Handle<reco::PFTauDiscriminator> tauDiscrAgainstMuons;
  evt.getByLabel(tauDiscrAgainstMuons_, tauDiscrAgainstMuons);

  int theTauJetIndex = -1;
  const reco::PFTau* theTauJet = getTheTauJet(*tauJets, tauJetEtaCut_, tauJetPtCut_, theTauJetIndex);

  double theTauDiscrByLeadTrackFinding = -1.;
  double theTauDiscrByLeadTrackPtCut = -1.;
  double theTauDiscrByTrackIso = -1.;
  double theTauDiscrByEcalIso = -1.;
  double theTauDiscrAgainstMuons = -1.;
  if ( theTauJetIndex != -1 ) {
    reco::PFTauRef theTauJetRef(tauJets, theTauJetIndex);
    theTauDiscrByLeadTrackFinding = (*tauDiscrByLeadTrackFinding)[theTauJetRef];
    theTauDiscrByLeadTrackPtCut = (*tauDiscrByLeadTrackPtCut)[theTauJetRef];
    theTauDiscrByTrackIso = (*tauDiscrByTrackIso)[theTauJetRef];
    theTauDiscrByEcalIso = (*tauDiscrByEcalIso)[theTauJetRef];
    theTauDiscrAgainstMuons = (*tauDiscrAgainstMuons)[theTauJetRef];
  }

//--- get missing transverse momentum
//    measured by calorimeters/reconstructed by particle-flow algorithm
  edm::Handle<reco::CaloMETCollection> caloMEtCollection;
  if ( !evt.getByLabel(caloMEtSource_, caloMEtCollection) ) {
    edm::LogWarning ("EwkMuTauHistManager") << "Failed to access calo. MET collection !!";
    return;
  }
  const reco::CaloMET& caloMEt = caloMEtCollection->at(0);
  
  edm::Handle<reco::PFMETCollection> pfMEtCollection;
  if ( !evt.getByLabel(pfMEtSource_, pfMEtCollection) ) {
    edm::LogWarning ("EwkMuTauHistManager") << "Failed to access pf. MET collection !!";
    return;
  }
  const reco::PFMET& pfMEt = pfMEtCollection->at(0);

  if ( !(theMuon && theTauJet && theTauJetIndex != -1) ) return;

  //-----------------------------------------------------------------------------
  // compute EWK tau analysis specific quantities
  //-----------------------------------------------------------------------------

  double dPhiMuTau = calcDeltaPhi(theMuon->phi(), theTauJet->phi());

  double mMuTau = (theMuon->p4() + theTauJet->p4()).M();

  //double mtMuCaloMEt = calcMt(theMuon->px(), theMuon->px(), caloMEt.px(), caloMEt.py());
  double mtMuPFMEt = calcMt(theMuon->px(), theMuon->px(), pfMEt.px(), pfMEt.py());

  //double pZetaCaloMEt = calcPzeta(theMuon->p4(), theTauJet->p4(), caloMEt.px(), caloMEt.py());
  //double pZetaPFMEt = calcPzeta(theMuon->p4(), theTauJet->p4(), pfMEt.px(), pfMEt.py());

  //-----------------------------------------------------------------------------
  // apply selection criteria; fill histograms
  //-----------------------------------------------------------------------------

//--- fill muon multiplicity histogram
  unsigned numGlobalMuons = 0;
  for ( reco::MuonCollection::const_iterator muon = muons->begin();
	muon != muons->end(); ++muon ) {
    if ( muon->isGlobalMuon() ) {
      ++numGlobalMuons;
    }
  }

  //hNumGlobalMuons_->Fill(numGlobalMuons);

  ++numEventsAnalyzed_;

  bool isSelected = false;
  int cutFlowStatus = -1;

  if ( mMuTau > visMassCut_ ) {
    cutFlowStatus = kPassedPreselection;
  }
  if ( cutFlowStatus == kPassedPreselection && (isTriggered || hltPaths_.size() == 0) ) {
    cutFlowStatus = kPassedTrigger;
  }
  if ( cutFlowStatus == kPassedTrigger && theMuon->isGlobalMuon() ) {
    cutFlowStatus = kPassedMuonId;
    hMuonTrackIsoPt_->Fill(theMuonTrackIsoPt);
  }
  if ( cutFlowStatus == kPassedMuonId && theMuonTrackIsoPt < muonTrackIsoCut_ ) {
    cutFlowStatus = kPassedMuonTrackIso;
    hMuonEcalIsoPt_->Fill(theMuonEcalIsoPt);
  }
  if ( cutFlowStatus == kPassedMuonTrackIso && theMuonEcalIsoPt < muonEcalIsoCut_ ) {
    cutFlowStatus = kPassedMuonEcalIso;
  }
  if ( cutFlowStatus == kPassedMuonEcalIso && theTauDiscrByLeadTrackFinding > 0.5 ) {
    cutFlowStatus = kPassedTauLeadTrack;
    if ( theTauJet->leadTrack().isAvailable() ) {
      hTauLeadTrackPt_->Fill(theTauJet->leadTrack()->pt());
    } else {
      edm::LogWarning ("EwkMuTauHistManager") << "Failed to access lead. Track of tau-jet !!";
    }
  }
  if ( cutFlowStatus == kPassedTauLeadTrack && theTauDiscrByLeadTrackPtCut > 0.5 ) {
    cutFlowStatus = kPassedTauLeadTrackPt;
    hTauTrackIsoPt_->Fill(theTauJet->isolationPFChargedHadrCandsPtSum());
  }
  if ( cutFlowStatus == kPassedTauLeadTrackPt && theTauDiscrByTrackIso > 0.5 ) {
    cutFlowStatus = kPassedTauTrackIso;
    hTauEcalIsoPt_->Fill(theTauJet->isolationPFGammaCandsEtSum());
  }
  if ( cutFlowStatus == kPassedTauTrackIso && theTauDiscrByEcalIso > 0.5 ) {
    cutFlowStatus = kPassedTauEcalIso;
    hTauDiscrAgainstMuons_->Fill(theTauDiscrAgainstMuons);
  }
  if ( cutFlowStatus == kPassedTauEcalIso && theTauDiscrAgainstMuons > 0.5 ) {
    cutFlowStatus = kPassedTauDiscrAgainstMuons;
    isSelected = true;
  }

  for ( int iCut = 1; iCut <= cutFlowStatus; ++iCut ) {
    hCutFlowSummary_->Fill(iCut);
  }

  if ( isSelected ) {
    hMuonPt_->Fill(theMuon->pt());
    hMuonEta_->Fill(theMuon->eta());
    hMuonPhi_->Fill(theMuon->phi());

    hTauJetPt_->Fill(theTauJet->pt());
    hTauJetEta_->Fill(theTauJet->eta());
    hTauJetPhi_->Fill(theTauJet->phi());

    //hTauJetCharge_->Fill(theTauJet->charge());
    if ( theTauJet->signalTracks().isAvailable() ) {
      hTauJetNumSignalTracks_->Fill(theTauJet->signalTracks().size());
    } else {
      edm::LogWarning ("EwkMuTauHistManager") << "Failed to access signal Tracks associated to tau-jet !!";
    }
    if ( theTauJet->isolationTracks().isAvailable() ) {
      hTauJetNumIsoTracks_->Fill(theTauJet->isolationTracks().size());
    } else {
      edm::LogWarning ("EwkMuTauHistManager") << "Failed to access isolation Tracks associated to tau-jet !!";
    }
  
    hVisMass_->Fill(mMuTau);
    //hMtMuCaloMEt_->Fill(mtMuCaloMEt);
    hMtMuPFMEt_->Fill(mtMuPFMEt);
    //hPzetaCaloMEt_->Fill(pZetaCaloMEt);
    //hPzetaPFMEt_->Fill(pZetaPFMEt);
    hMuTauAcoplanarity_->Fill(dPhiMuTau);
    //hMuTauCharge_->Fill(theMuon->charge() + theTauJet->charge());

    if ( theEventVertex ) {
      //hVertexChi2_->Fill(theEventVertex->normalizedChi2());
      hVertexZ_->Fill(theEventVertex->z());
      //hVertexD0_->Fill(getVertexD0(*theEventVertex, *beamSpot));
    }
    
    hCaloMEtPt_->Fill(caloMEt.pt());
    //hCaloMEtPhi_->Fill(caloMEt.phi());

    hPFMEtPt_->Fill(pfMEt.pt());
    //hPFMEtPhi_->Fill(pfMEt.phi());
  }

  if ( isSelected ) {
    ++numEventsSelected_;
    LogTrace ("EwkMuTauHistManager") << "Event accepted.";
  } else {
    LogTrace ("EwkMuTauHistManager") << "Event rejected.";
  }
}

void EwkMuTauHistManager::finalizeHistograms()
{
  edm::LogInfo ("EwkMuTauHistManager") 
    << "Filter-Statistics Summary:" << std::endl
    << " Events analyzed = " << numEventsAnalyzed_ << std::endl
    << " Events selected = " << numEventsSelected_;
  if ( numEventsAnalyzed_ > 0 ) {
    double eff = numEventsSelected_/(double)numEventsAnalyzed_;
    edm::LogInfo ("") 
      << "Overall efficiency = " << std::setprecision(4) << eff*100. 
      << " +/- " << std::setprecision(4) << TMath::Sqrt(eff*(1 - eff)/numEventsAnalyzed_)*100. << ")%";
  }
}

//-------------------------------------------------------------------------------
// common auxiliary functions used by different channels
//-------------------------------------------------------------------------------

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <TMath.h>

int getIsoMode(const std::string& isoMode_string, int& error)
{
  int isoMode_int;
  if ( isoMode_string == "absoluteIso" ) {
    isoMode_int = kAbsoluteIso;
  } else if ( isoMode_string == "relativeIso" ) {
    isoMode_int = kRelativeIso;
  } else { 
    edm::LogError ("getIsoMode") << " Failed to decode isoMode string = " << isoMode_string << " !!";
    isoMode_int = kUndefinedIso;
    error = 1;
  }
  return isoMode_int;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

double calcDeltaPhi(double phi1, double phi2) 
{
  double deltaPhi = phi1 - phi2;

  if ( deltaPhi < 0. ) deltaPhi = -deltaPhi;

  if ( deltaPhi > TMath::Pi() ) deltaPhi = 2*TMath::Pi() - deltaPhi;

  return deltaPhi;
}

double calcMt(double px1, double py1, double px2, double py2) 
{
  double pt1 = TMath::Sqrt(px1*px1 + py1*py1);
  double pt2 = TMath::Sqrt(px2*px2 + py2*py2);

  double p1Dotp2 = px1*px2 + py1*py2;
  double cosAlpha = p1Dotp2/(pt1*pt2);

  return TMath::Sqrt(2*pt1*pt2*(1 - cosAlpha));
}

double calcPzeta(const reco::Candidate::LorentzVector& p1,const reco::Candidate::LorentzVector& p2, double pxMEt, double pyMEt)
{
  double cosPhi1 = cos(p1.phi());
  double sinPhi1 = sin(p1.phi());
  double cosPhi2 = cos(p2.phi());
  double sinPhi2 = sin(p2.phi());
  double zetaX = cosPhi1 + cosPhi2;
  double zetaY = sinPhi1 + sinPhi2;
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) {
    zetaX /= zetaR;
    zetaY /= zetaR;
  }

  double pxVis = p1.px() + p2.px();
  double pyVis = p1.py() + p2.py();
  double pZetaVis = pxVis*zetaX + pyVis*zetaY;

  double px = pxVis + pxMEt;
  double py = pyVis + pyMEt;
  double pZeta = px*zetaX + py*zetaY;

  return pZeta - 1.5*pZetaVis;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

bool passesElectronPreId(const reco::GsfElectron& electron)
{
  if ( (TMath::Abs(electron.eta()) < 1.479 || TMath::Abs(electron.eta()) > 1.566) && // cut ECAL barrel/endcap crack
       electron.deltaPhiSuperClusterTrackAtVtx() < 0.58 &&
       electron.deltaEtaSuperClusterTrackAtVtx() < 0.01 &&
       electron.sigmaIetaIeta() < 0.027 ) {
    return true;
  } else {
    return false;
  }
}

bool passesElectronId(const reco::GsfElectron& electron)
{
  if ( passesElectronPreId(electron) &&
       ((TMath::Abs(electron.eta()) > 1.566 && // electron reconstructed in ECAL barrel
	 electron.sigmaEtaEta() < 0.028 && electron.hcalOverEcal() < 0.1 && 
	 TMath::Abs(electron.deltaEtaSuperClusterTrackAtVtx()) < 0.0066 && 
	 TMath::Abs(electron.deltaPhiSuperClusterTrackAtVtx()) < 0.020 ) ||
	(TMath::Abs(electron.eta()) < 1.479 && // electron reconstructed in ECAL endcap
	 electron.sigmaEtaEta() < 0.0099 && electron.hcalOverEcal() < 0.1 && 
	 TMath::Abs(electron.deltaEtaSuperClusterTrackAtVtx()) < 0.004 && 
	 TMath::Abs(electron.deltaPhiSuperClusterTrackAtVtx()) < 0.025)) ) {
    return true;
  } else {
    return false;
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

const reco::GsfElectron* getTheElectron(const reco::GsfElectronCollection& electrons, double electronEtaCut, double electronPtCut)
{
  const reco::GsfElectron* theElectron = 0;
  
  for ( reco::GsfElectronCollection::const_iterator electron = electrons.begin();
	electron != electrons.end(); ++electron ) {
    if ( TMath::Abs(electron->eta()) < electronEtaCut && electron->pt() > electronPtCut && passesElectronPreId(*electron) ) {
      if ( theElectron == 0 || electron->pt() > theElectron->pt() ) theElectron = &(*electron);
    }
  }
  
  return theElectron;
}

const reco::Muon* getTheMuon(const reco::MuonCollection& muons, double muonEtaCut, double muonPtCut)
{
  const reco::Muon* theMuon = 0;

  for ( reco::MuonCollection::const_iterator muon = muons.begin();
	muon != muons.end(); ++muon ) {
    if ( TMath::Abs(muon->eta()) < muonEtaCut && muon->pt() > muonPtCut ) {
      if ( theMuon == 0 || muon->pt() > theMuon->pt() ) theMuon = &(*muon);
    }
  }

  return theMuon;
}

const reco::PFTau* getTheTauJet(const reco::PFTauCollection& tauJets, double tauJetEtaCut, double tauJetPtCut, int& theTauJetIndex)
{
  const reco::PFTau* theTauJet = 0;
  theTauJetIndex = -1;

  int numTauJets = tauJets.size();
  for ( int iTauJet = 0; iTauJet < numTauJets; ++iTauJet ) {
    const reco::PFTau& tauJet = tauJets.at(iTauJet);

    if ( TMath::Abs(tauJet.eta()) < tauJetEtaCut && tauJet.pt() > tauJetPtCut ) {
      if ( theTauJet == 0 || tauJet.pt() > theTauJet->pt() ) {
	theTauJet = &tauJet;
	theTauJetIndex = iTauJet;
      }
    }
  }

  return theTauJet;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

double getVertexD0(const reco::Vertex& vertex, const reco::BeamSpot& beamSpot)
{
  double dX = vertex.x() - beamSpot.x0();
  double dY = vertex.y() - beamSpot.y0();
  return TMath::Sqrt(dX*dX + dY*dY);
}
