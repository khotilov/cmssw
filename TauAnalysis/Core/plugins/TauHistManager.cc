#include "TauAnalysis/Core/plugins/TauHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/IsolationAlgos/interface/IsoDepositVetoFactory.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"
#include "TauAnalysis/Core/interface/eventAuxFunctions.h"
#include "TauAnalysis/GenSimTools/interface/genParticleAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>

#include <stdlib.h>

bool matchesGenTau(const pat::Tau& patTau)
{
  bool isGenTauMatched = false;
  std::vector<reco::GenParticleRef> associatedGenParticles = patTau.genParticleRefs();
  for ( std::vector<reco::GenParticleRef>::const_iterator it = associatedGenParticles.begin();
	it != associatedGenParticles.end(); ++it ) {
    if ( isValidRef(*it) ) {
      const reco::GenParticleRef& genParticle = (*it);
      if ( genParticle->pdgId() == -15 || genParticle->pdgId() == +15 ) isGenTauMatched = true;
    }
  }
  return isGenTauMatched;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TauHistManager::TauHistManager(const edm::ParameterSet& cfg)
  : HistManagerBase(cfg)
{
  //std::cout << "<TauHistManager::TauHistManager>:" << std::endl;

  tauSrc_ = cfg.getParameter<edm::InputTag>("tauSource");
  //std::cout << " tauSrc = " << tauSrc_.label() << std::endl;

  vertexSrc_ = ( cfg.exists("vertexSource") ) ? cfg.getParameter<edm::InputTag>("vertexSource") : edm::InputTag();
  if ( vertexSrc_.label() == "" ) {
    edm::LogWarning("TauHistManager")
      << " Configuration parameter 'vertexSource' not specified"
      << " --> Impact Parameter histograms will NOT be plotted !!";
  }
  //std::cout << " vertexSrc = " << vertexSrc_.label() << std::endl;

  jetSrc_ = cfg.getParameter<edm::InputTag>("jetSource");
  //std::cout << " jetSrc = " << jetSrc_ << std::endl;

  genParticleSrc_ = ( cfg.exists("genParticleSource") ) ? cfg.getParameter<edm::InputTag>("genParticleSource") : edm::InputTag();
  //std::cout << " genParticleSrc = " << genParticleSrc_ << std::endl;

  tauJetWeightExtractors_ = getTauJetWeightExtractors<pat::Tau>(cfg, "tauJetWeightSource");

  if ( cfg.exists("tauIndicesToPlot") ) {
    std::string tauIndicesToPlot_string = cfg.getParameter<std::string>("tauIndicesToPlot");
    if ( tauIndicesToPlot_string != "all" ) {
      size_t posTauIndexBegin = 0;
      size_t posTauIndexEnd = posTauIndexBegin;
      do {
	size_t posNextSeparator = tauIndicesToPlot_string.find(",", posTauIndexBegin);
	posTauIndexEnd = ( posNextSeparator < std::string::npos ) ? posNextSeparator : tauIndicesToPlot_string.length();

	std::string tauIndex_string(tauIndicesToPlot_string, posTauIndexBegin, posTauIndexEnd - posTauIndexBegin);
	int tauIndex_int = atoi(tauIndex_string.data());
	//std::cout << "--> adding tauIndex_int = " << tauIndex_int << std::endl;
	tauIndicesToPlot_.push_back(tauIndex_int);

	posTauIndexBegin = posTauIndexEnd + 1;
      } while ( posTauIndexEnd < tauIndicesToPlot_string.length() );
    }
  }

  requireGenTauMatch_ = cfg.getParameter<bool>("requireGenTauMatch");
  //std::cout << " requireGenTauMatch = " << requireGenTauMatch_ << std::endl;

  skipPdgIdsGenParticleMatch_ = cfg.getParameter<vint>("skipPdgIdsGenParticleMatch");

  useHPSpTaNCalgorithm_ = cfg.getParameter<bool>("useHPSpTaNCalgorithm");
  useHPSclassicAlgorithm_ = cfg.getParameter<bool>("useHPSclassicAlgorithm");
  if(useHPSpTaNCalgorithm_ && useHPSclassicAlgorithm_ )
      std::cout << "TauHistManager::TauHistManager --> must NOT set both 'useHPSpTaNCalgorithm' and 'useHPSclassicAlgorithm' true" << std::endl
          << " Using 'useHPSpTaNCalgorithm'! " << std::endl;
    

  std::string normalization_string = cfg.getParameter<std::string>("normalization");
  normMethod_ = getNormMethod(normalization_string, "taus");

  makeIsoPtCtrlHistograms_ = ( cfg.exists("makeIsoPtCtrlHistograms") ) ?
    cfg.getParameter<bool>("makeIsoPtCtrlHistograms") : false;

  makeIsoPtConeSizeDepHistograms_ = ( cfg.exists("makeIsoPtConeSizeDepHistograms") ) ?
    cfg.getParameter<bool>("makeIsoPtConeSizeDepHistograms") : false;

  checkWeightConsistency_ = ( cfg.exists("checkWeightConsistency") ) ?
    cfg.getParameter<bool>("checkWeightConsistency") : false;
  //std::cout << " checkWeightConsistency = " << checkWeightConsistency_ << std::endl;

  numTauIsoConeSizes_ = 15;
  tauIsoConeSizeIncr_ = 0.1;
  numTauIsoPtThresholds_ = 4;
  tauIsoPtThresholdIncr_ = 0.5;

//--- create "veto" objects for computation of IsoDeposit sums
  tauParticleFlowIsoParam_.push_back(IsoDepositVetoFactory::make("0.0"));
  tauParticleFlowIsoParam_.push_back(IsoDepositVetoFactory::make("Threshold(0.5)"));

  makeTauIdEfficiencyHistograms_ = ( cfg.exists("makeTauIdEfficiencyHistograms") ) ?
    cfg.getParameter<vstring>("makeTauIdEfficiencyHistograms") : vstring();
  makeTauFakeRateHistograms_ = ( cfg.exists("makeTauFakeRateHistograms") ) ?
    cfg.getParameter<vstring>("makeTauFakeRateHistograms") : vstring();
}

TauHistManager::~TauHistManager()
{
//--- delete "veto" objects for computation of IsoDeposit sums
  clearIsoParam(tauParticleFlowIsoParam_);

  for ( std::vector<FakeRateJetWeightExtractor<pat::Tau>*>::iterator it = tauJetWeightExtractors_.begin();
	it != tauJetWeightExtractors_.end(); ++it ) {
    delete (*it);
  }
}

void TauHistManager::bookHistogramsImp()
{
  //std::cout << "<TauHistManager::bookHistogramsImp>:" << std::endl;

//--- book histogram for number of tau-jets in each event
  hNumTaus_ = book1D("NumTaus", "Tau Jets in Event", 10, -0.5, 9.5);

//--- book histograms for Pt, eta and phi distributions
//    of tau-jets passing all id. and isolation selections
  bookTauHistograms(hTauPt_, hTauEta_, hTauPhi_, "Tau");
  hTauPtVsEta_ = book2D("TauPtVsEta", "TauPtVsEta", 24, -3., +3., 30, 0., 150.);
  hTauCharge_ = book1D("TauCharge", "Tau Charge (#Sigma Tracks in Signal Cone)", 11, -5.5, +5.5);
  hTauJetRadius_ = book1D("TauJetRadius", "Tau jet-Radius", 51, -0.005, +0.505);
  hTauJetRadiusPtProfile_ = bookProfile1D("TauJetRadiusPtProfile", "Tau jet-Radius vs. P_{T}", 30, 0., 150.);
  hTauJetRadiusEnProfile_ = bookProfile1D("TauJetRadiusEnProfile", "Tau jet-Radius vs. Energy", 50, 0., 250.);

  hTauVisMass_ = book1D("TauVisMass", "TauVisMass", 100, 0., 2.5);
  hTauVisMassRes_ = book1D("TauVisMassRes", "TauVisMassRes", 100, -1.25, 1.25);
  hTauVisMassOneProngOnePi0_ = book1D("TauVisMassOneProngOnePi0", "TauVisMassOneProngOnePi0", 100, 0., 2.5);
  hTauVisMassResOneProngOnePi0_ = book1D("TauVisMassResOneProngOnePi0", "TauVisMassResOneProngOnePi0", 100, -1.25, 1.25);
  hTauVisMassOneProngTwoPi0s_ = book1D("TauVisMassOneProngTwoPi0s", "TauVisMassOneProngTwoPi0s", 100, 0., 2.5);
  hTauVisMassResOneProngTwoPi0s_ = book1D("TauVisMassResOneProngTwoPi0s", "TauVisMassResOneProngTwoPi0s", 100, -1.25, 1.25);
  hTauVisMassThreeProngNoPi0s_ = book1D("TauVisMassThreeProngNoPi0s", "TauVisMassThreeProngNoPi0s", 100, 0., 2.5);
  hTauVisMassResThreeProngNoPi0s_ = book1D("TauVisMassResThreeProngNoPi0s", "TauVisMassResThreeProngNoPi0s", 100, -1.25, 1.25);
  hTauVisMassThreeProngOnePi0_ = book1D("TauVisMassThreeProngOnePi0", "TauVisMassThreeProngOnePi0", 100, 0., 2.5);
  hTauVisMassResThreeProngOnePi0_ = book1D("TauVisMassResThreeProngOnePi0", "TauVisMassResThreeProngOnePi0", 100, -1.25, 1.25);

  hDistPionEnResOneProngOnePi0_ = book1D("DistPionEnResOneProngOnePi0", "DistPionEnResOneProngOnePi0", 100, -1.25, 1.25);
  hDistPionEnResOneProngTwoPi0s_ = book1D("DistPionEnResOneProngTwoPi0s", "DistPionEnResOneProngTwoPi0s", 100, -1.25, 1.25);
  hDistPionEnResThreeProngNoPi0s_ = book1D("DistPionEnResThreeProngNoPi0s", "DistPionEnResThreeProngNoPi0s", 100, -1.25, 1.25);

  bookWeightHistograms(*dqmStore_, "TauJetWeight", "Tau Weight",
		       hTauJetWeightPosLog_, hTauJetWeightNegLog_, hTauJetWeightZero_,
		       hTauJetWeightLinear_);

  bookTauHistograms(hGenTauPt_, hGenTauEta_, hGenTauPhi_, "GenTau");
  hGenTauX_ = book1D("GenTauX", "GenTauX", 50, 0, 1);
  bookTauHistograms(hGenVisTauPt_, hGenVisTauEta_, hGenVisTauPhi_, "GenVisTau");

  hTauEnCompToGen_ = book1D("TauEnCompToGen", "RECO-GEN #Delta E", 100, -2.50, +2.50);
  hTauThetaCompToGen_ = book1D("TauThetaCompToGen", "RECO-GEN #Delta#theta", 200, -0.250, +0.250);
  hTauPhiCompToGen_ = book1D("TauPhiCompToGen", "RECO-GEN #Delta#phi", 100, -0.050, +0.050);

  hTauMatchingGenParticlePdgId_ = book1D("TauMatchingGenParticlePdgId", "matching gen. Particle PdgId", 26, -1.5, 24.5);
  hTauMatchingFinalStateGenParticlePdgId_ = book1D("TauMatchingFinalStateGenParticlePdgId", "matching final state gen. Particle PdgId", 26, -1.5, 24.5);
  hTauMatchingGenTauDecayMode_ = book1D("TauMatchingGenTauDecayMode", "matching gen. Tau decay mode", 20, -0.5, 19.5);
  setAxisLabelsGenTauDecayMode(hTauMatchingGenTauDecayMode_->getTH1()->GetXaxis());

  hTauLeadPFChargedHadCandRefValidity_ = book1D("LeadPfChargedHadValidity","Ref for Leading PF charged hadron is valid",2,-0.5,1.5);
  hTauLeadKfTrkRefValidity_ = book1D("LeadKfTrkRefIsValid","KF Track Ref of lead trk is Valid",2,-0.5,1.5);
  hTauLeadGsfTrkRefValidity_ = book1D("LeadGsfTrkRefIsValid","GSF Track Ref of lead trk is Valid",2,-0.5,1.5);

  hTauNumTracksSignalCone_ = book1D("TauNumTracksSignalCone", "Tracks in Signal Cone", 10, -0.5, 9.5);
  hTauNumTracksIsoCone_ = book1D("TauNumTracksIsoCone", "Tracks in Isolation Cone", 20, -0.5, 19.5);

  bookTauHistograms(hTauLeadTrkPt_, hTauLeadTrkEta_, hTauLeadTrkPhi_, "TauLeadTrk");
  hTauLeadTrkMatchDist_ = book1D("TauLeadTrkMatchDist", "TauLeadTrkMatchDist", 100, -0.500, 0.500);
  hTauLeadTrkIPxy_ = book1D("TauLeadTrkIPxy", "Lead Track Impact Parameter (xy)", 100, -0.100, 0.100);
  hTauLeadTrkIPz_ = book1D("TauLeadTrkIPz", "Lead Track Impact Parameter (z)", 100, -1.0, 1.0);
  hTauLeadTrkNumHits_ = book1D("TauLeadTrkNumHits", "Lead Track Number of Pixel + Strip Hits", 25, -0.5, 24.5);
  hTauLeadTrkNumPixelHits_ = book1D("TauLeadTrkNumPixelHits", "Lead Track Number of Pixel Hits", 5, -0.5, 4.5);
  hTauLeadTrkNumStripHits_ = book1D("TauLeadTrkNumStripHits", "Lead Track Number of Strip Hits", 20, -0.5, 19.5);

  hTauDiscriminatorByIsolation_ = book1D("TauDiscriminatorByIsolation",
					 "Discriminator by Isolation (Track and ECAL)", 2, -0.5, 1.5);
  hTauDiscriminatorByTrackIsolation_ = book1D("TauDiscriminatorByTrackIsolation",
					      "Discriminator by Track Isolation", 2, -0.5, 1.5);
  hTauDiscriminatorByEcalIsolation_ = book1D("TauDiscriminatorByEcalIsolation",
					     "Discriminator by ECAL Isolation", 2, -0.5, 1.5);

  hTauDiscriminatorAgainstElectronsLoose_ = 
    book1D("TauDiscriminatorAgainstElectronsLoose", "Discriminator against Electrons (loose)", 2, -0.5, 1.5);
  hTauDiscriminatorAgainstElectronsMedium_ = 
    book1D("TauDiscriminatorAgainstElectronsMedium", "Discriminator against Electrons (medium)", 2, -0.5, 1.5);
  hTauDiscriminatorAgainstElectronsTight_ = 
    book1D("TauDiscriminatorAgainstElectronsTight", "Discriminator against Electrons (tight)", 2, -0.5, 1.5);
  hTauPFElectronMVA_ = book1D("TauPFElectronMVA", "TauPFElectronMVA", 40, -1.01, +1.01);
  hTauEmFraction_ = book1D("TauEmFraction", "TauEmFraction", 101, -0.01, 2.01);
  hTauHcalTotOverPLead_ = book1D("TauHcalTotOverPLead", "TauHcalTotOverPLead", 101, -0.01, 2.01);
  hTauHcalMaxOverPLead_ = book1D("TauHcalMaxOverPLead", "TauHcalMaxOverPLead", 101, -0.01, 2.01);
  hTauHcal3x3OverPLead_ = book1D("TauHcal3x3OverPLead", "TauHcal3x3OverPLead", 101, -0.01, 2.01);
  hTauEcalStripSumEOverPLead_ = book1D("TauEcalStripSumEOverPLead", "TauEcalStripSumEOverPLead", 101, -0.01, 2.01);
  hTauBremsRecoveryEOverPLead_ = book1D("TauBremsRecoveryEOverPLead", "TauBremsRecoveryEOverPLead", 101, -0.01, 2.01);
  hTauCaloEOverPLead_ = book1D("TauCaloEOverPLead", "TauCaloEOverPLead", 101, -0.01, 2.01);

  hTauDiscriminatorAgainstMuonsLoose_ = 
    book1D("TauDiscriminatorAgainstMuonsLoose", "Discriminator against Muons (loose)", 2, -0.5, 1.5);
  hTauDiscriminatorAgainstMuonsTight_ = 
    book1D("TauDiscriminatorAgainstMuonsTight", "Discriminator against Muons (tight)", 2, -0.5, 1.5);

  hTauRecDecayMode_ = book1D("TauRecDecayMode", "rec. Tau decay mode", 20, -0.5, 19.5);
  setAxisLabelsRecTauDecayMode(hTauRecDecayMode_->getTH1()->GetXaxis());
  hTauRecVsGenDecayMode_ = book2D("TauRecVsGenDecayMode", "rec. vs. gen. Tau decay mode", 20, -0.5, 19.5, 20, -0.5, 19.5);
  setAxisLabelsGenTauDecayMode(hTauRecVsGenDecayMode_->getTH1()->GetXaxis());
  setAxisLabelsRecTauDecayMode(hTauRecVsGenDecayMode_->getTH1()->GetYaxis());

  hTauTaNCoutputOneProngNoPi0s_ = book1D("TauTaNCoutputOneProngNoPi0s",
					 "TauTaNCoutputOneProngNoPi0s", 102, -0.01, 1.01);
  hTauTaNCoutputOneProngOnePi0_ = book1D("TauTaNCoutputOneProngOnePi0",
					 "TauTaNCoutputOneProngOnePi0", 102, -0.01, 1.01);
  hTauTaNCoutputOneProngTwoPi0s_ = book1D("TauTaNCoutputOneProngTwoPi0s",
					  "TauTaNCoutputOneProngTwoPi0s", 102, -0.01, 1.01);
  hTauTaNCoutputThreeProngNoPi0s_ = book1D("TauTaNCoutputThreeProngNoPi0s",
					   "TauTaNCoutputThreeProngNoPi0s", 102, -0.01, 1.01);
  hTauTaNCoutputThreeProngOnePi0_ = book1D("TauTaNCoutputThreeProngOnePi0",
					   "TauTaNCoutputThreeProngOnePi0", 102, -0.01, 1.01);
  hTauTaNCoutputTransform_ = book1D("TauTaNCoutputTransform",
				    "TauTaNCoutputTransform", 102, -0.01, 1.01);

  hTauDiscriminatorTaNCfrOnePercent_ = book1D("TauDiscriminatorTaNCfrOnePercent",
					      "TauDiscriminatorTaNCfrOnePercent", 2, -0.5, 1.5);
  hTauDiscriminatorTaNCfrHalfPercent_ = book1D("TauDiscriminatorTaNCfrHalfPercent",
					       "TauDiscriminatorTaNCfrHalfPercent", 2, -0.5, 1.5);
  hTauDiscriminatorTaNCfrQuarterPercent_ = book1D("TauDiscriminatorTaNCfrQuarterPercent",
						  "TauDiscriminatorTaNCfrQuarterPercent", 2, -0.5, 1.5);
  hTauDiscriminatorTaNCfrTenthPercent_ = book1D("TauDiscriminatorTaNCfrTenthPercent",
						"TauDiscriminatorTaNCfrTenthPercent", 2, -0.5, 1.5);

  hTauDiscriminatorCombinedTaNCvloose_ = book1D("TauDiscriminatorHPSpTaNCwTaNCvloose",
					"TauDiscriminatorHPS+TaNC: TaNCvloose", 2, -0.5, 1.5);
  hTauDiscriminatorCombinedTaNCloose_ = book1D("TauDiscriminatorHPSpTaNCwTaNCloose",
				       "TauDiscriminatorHPS+TaNC: TaNCloose", 2, -0.5, 1.5);
  hTauDiscriminatorCombinedTaNCmedium_ = book1D("TauDiscriminatorHPSpTaNCwTaNCmedium",
					"TauDiscriminatorHPS+TaNC: TaNCmedium", 2, -0.5, 1.5);
  hTauDiscriminatorCombinedTaNCtight_ = book1D("TauDiscriminatorHPSpTaNCwTaNCtight",
				       "TauDiscriminatorHPS+TaNC: TaNCtight", 2, -0.5, 1.5);

  hTauDiscriminatorCombinedHPSvloose_ = book1D("TauDiscriminatorHPSpTaNCwHPSvloose",
				      "TauDiscriminatorHPS+TaNC:HPSvloose", 2, -0.5, 1.5);
  hTauDiscriminatorCombinedHPSloose_ = book1D("TauDiscriminatorHPSpTaNCwHPSloose",
				      "TauDiscriminatorHPS+TaNC:HPSloose", 2, -0.5, 1.5);
  hTauDiscriminatorCombinedHPSmedium_ = book1D("TauDiscriminatorHPSpTaNCwHPSmedium",
				       "TauDiscriminatorHPS+TaNC:HPSmedium", 2, -0.5, 1.5);
  hTauDiscriminatorCombinedHPStight_ = book1D("TauDiscriminatorHPSpTaNCwHPStight",
				      "TauDiscriminatorHPS+TaNC:HPStight", 2, -0.5, 1.5);

  hTauDiscriminatorHPSvloose_ = book1D("TauDiscriminatorHPSvloose",
				      "TauDiscriminatorHPS: vloose", 2, -0.5, 1.5);
  hTauDiscriminatorHPSloose_ = book1D("TauDiscriminatorHPSloose",
				      "TauDiscriminatorHPS: loose", 2, -0.5, 1.5);
  hTauDiscriminatorHPSmedium_ = book1D("TauDiscriminatorHPSmedium",
				       "TauDiscriminatorHPS: medium", 2, -0.5, 1.5);
  hTauDiscriminatorHPStight_ = book1D("TauDiscriminatorHPStight",
				      "TauDiscriminatorHPS: tight", 2, -0.5, 1.5);

  hTauDiscriminatorHPSvlooseDeltaB_ = book1D("TauDiscriminatorHPSvlooseDeltaBeta",
				      "TauDiscriminatorHPSdeltaBeta corr.: vloose", 2, -0.5, 1.5);
  hTauDiscriminatorHPSlooseDeltaB_ = book1D("TauDiscriminatorHPSlooseDeltaBeta",
				      "TauDiscriminatorHPSdeltaBeta corr.: loose", 2, -0.5, 1.5);
  hTauDiscriminatorHPSmediumDeltaB_ = book1D("TauDiscriminatorHPSmediumDeltaBeta",
				       "TauDiscriminatorHPSdeltaBeta corr.: medium", 2, -0.5, 1.5);
  hTauDiscriminatorHPStightDeltaB_ = book1D("TauDiscriminatorHPStightDeltaBeta",
				      "TauDiscriminatorHPSdeltaBeta corr.: tight", 2, -0.5, 1.5);
  
  hTauDiscriminatorHPSvlooseCombDeltaB_ = book1D("TauDiscriminatorHPSvlooseCombDeltaBeta",
				      "TauDiscriminatorHPScombDeltaBeta corr.: vloose", 2, -0.5, 1.5);
  hTauDiscriminatorHPSlooseCombDeltaB_ = book1D("TauDiscriminatorHPSlooseCombDeltaBeta",
				      "TauDiscriminatorHPScombDeltaBeta corr.: loose", 2, -0.5, 1.5);
  hTauDiscriminatorHPSmediumCombDeltaB_ = book1D("TauDiscriminatorHPSmediumCombDeltaBeta",
				       "TauDiscriminatorHPScombDeltaBeta corr.: medium", 2, -0.5, 1.5);
  hTauDiscriminatorHPStightCombDeltaB_ = book1D("TauDiscriminatorHPStightCombDeltaBeta",
				      "TauDiscriminatorHPScombDeltaBeta corr.: tight", 2, -0.5, 1.5);
  
  hTauTrkIsoPt_ = book1D("TauTrkIsoPt", "Track Isolation P_{T}", 100, 0., 10.);
  hTauEcalIsoPt_ = book1D("TauEcalIsoPt", "ECAL Isolation P_{T}", 100, 0., 10.);
  hTauHcalIsoPt_ = book1D("TauHcalIsoPt", "HCAL Isolation P_{T}", 100, 0., 10.);
  hTauIsoSumPt_ = book1D("TauIsoSumPt", "Isolation Sum(P_{T})", 100, 0., 10.);

  hTauDeltaRnearestJet_ = book1D("TauDeltaRnearestJet", "#DeltaR(nearest Jet)", 102, -0.1, 10.1);

  hTauParticleFlowIsoPt_ = book1D("TauParticleFlowIsoPt", "Particle Flow Isolation P_{T}", 100, 0., 10.);
  hTauPFChargedHadronIsoPt_ = book1D("TauPFChargedHadronIsoPt", "Particle Flow (Charged Hadron) Isolation P_{T}", 100, 0., 10.);
  hTauPFNeutralHadronIsoPt_ = book1D("TauPFNeutralHadronIsoPt", "Particle Flow (Neutral Hadron) Isolation P_{T}", 100, 0., 10.);
  hTauPFGammaIsoPt_ = book1D("TauPFGammaIsoPt", "Particle Flow (Photon) Isolation P_{T}", 100, 0., 10.);

  hTauNumSignalPFChargedHadrons_ = book1D("TauNumSignalPFChargedHadrons", "PF charged hadrons in signal cone", 10, -0.5, 9.5);
  hTauNumIsoPFChargedHadrons_ = book1D("TauNumIsoPFChargedHadrons", "PF charged hadrons in isolation cone", 10, -0.5, 9.5);
  hTauNumPFChargedHadrons_ = book1D("TauNumPFChargedHadrons", "PF charged hadrons in Tau jet", 20, -0.5, 19.5);

  hTauNumSignalPFGammas_ = book1D("TauNumSignalPFGammas", "PF photons in signal cone", 10, -0.5, 9.5);
  hTauNumIsoPFGammas_ = book1D("TauNumIsoPFGammas", "PF photons in isolation cone", 10, -0.5, 9.5);
  hTauNumPFGammas_= book1D("TauNumPFGammas", "PF photons in Tau jet", 20, -0.5, 19.5);

 //--- book "control" histograms to check agreement between tau isolation variables
//    computed by PAT-level IsoDeposits with the values computed by reco::PFTau producer
  if ( makeIsoPtCtrlHistograms_ ) {
    hTauPFChargedHadronIsoPtCtrl_ = book2D("TauPFChargedHadronIsoPtCtrl", "Particle Flow (Charged Hadron) Isolation P_{T} (reco::PFTau vs. IsoDeposit)", 40, 0., 20., 40, 0., 20.);
    hTauPFGammaIsoPtCtrl_ = book2D("TauPFGammaIsoPtCtrl", "Particle Flow (Photon) Isolation P_{T} (reco::PFTau vs. IsoDeposit)", 40, 0., 20., 40, 0., 20.);
  } else {
    hTauPFChargedHadronIsoPtCtrl_ = 0;
    hTauPFGammaIsoPtCtrl_ = 0;
  }

  hTauTrkIsoEnProfile_ = book1D("TauTrkIsoEnProfile", "All Isolation Tracks #Delta P", 100, 0., 10.);
  hTauTrkIsoPtProfile_ = book1D("TauTrkIsoPtProfile", "All Isolation Tracks #P_{T}", 100, 0., 10.);
  hTauTrkIsoEtaDistProfile_ = book1D("TauTrkIsoEtaDistProfile", "All Isolation Tracks |#Delta#eta|", 15, 0., 1.5);
  hTauTrkIsoPhiDistProfile_ = book1D("TauTrkIsoPhiDistProfile", "All Isolation Tracks |#Delta#phi|", 15, 0., 1.5);

  if ( makeIsoPtConeSizeDepHistograms_ ) bookTauIsoConeSizeDepHistograms();

//--- book "control" histograms to check parametrization of
//    tau id. efficiency and fake-rate values stored in the pat::Tau
  bookTauIdEfficiencyHistograms(hTauIdEfficiencies_, "IdEfficiency", makeTauIdEfficiencyHistograms_);
  bookTauIdEfficiencyHistograms(hTauFakeRates_, "FakeRate", makeTauFakeRateHistograms_);
}

void TauHistManager::fillTauDiscriminatorHistogram(MonitorElement* h, const pat::Tau& patTau, const char* discrName,
						   std::map<std::string, bool>& discrAvailability_hasBeenChecked, double weight)
{
//--- tau id. discriminators not available for all kinds of taus
//    (in particular those based on TaNC are available only for shrinking signal cone PFTaus so far),
//    so need to check whether a given discriminator is available before filling histogram,
//    in order to avoid pat::Tau::tauID method from triggering an exception;
//    availability is checked only once and a warning is printed
//    in case the discriminator given as function argument is unavailable
  if ( !discrAvailability_hasBeenChecked[discrName] ) {
    if ( !patTau.isTauIDAvailable(discrName) ) {
      edm::LogWarning("TauHistManager")
	<< " Discriminator = " << discrName << " unavailable for pat::Tau collection = " << tauSrc_.label()
	<< " --> skipping filling of histogram = " << h->getName() << " !!";
    }

    discrAvailability_hasBeenChecked[discrName] = true;
  }

  if ( patTau.isTauIDAvailable(discrName) ) h->Fill(patTau.tauID(discrName), weight);
}

bool isIndexed(int patTauIndex, const std::vector<int>& tauIndicesToPlot)
{
  bool isIndexed = true;

  if ( tauIndicesToPlot.size() > 0 ) {
    isIndexed = false;
    for ( std::vector<int>::const_iterator tauIndexToPlot = tauIndicesToPlot.begin();
	  tauIndexToPlot != tauIndicesToPlot.end(); ++tauIndexToPlot ) {
      if ( (*tauIndexToPlot) == patTauIndex ) isIndexed = true;
    }
  }

  return isIndexed;
}

void fillTauVisMassHistogram(MonitorElement* hVisMass, MonitorElement* hVisMassRes, const std::string& selectedtauDecayMode,
			     double recVisMass, double genVisMass, const std::string& genTauDecayMode, double weight)
{
  if ( genTauDecayMode == selectedtauDecayMode || selectedtauDecayMode == "all" ) {
    hVisMass->Fill(recVisMass, weight);
    hVisMassRes->Fill(recVisMass - genVisMass, weight);
  }
}

void TauHistManager::fillHistogramsImp(const edm::Event& evt, const edm::EventSetup& es, double evtWeight)
{
  //std::cout << "<TauHistManager::fillHistogramsImp>:" << std::endl;

  edm::Handle<pat::TauCollection> patTaus;
  getCollection(evt, tauSrc_, patTaus);

  edm::Handle<pat::JetCollection> patJets;
  getCollection(evt, jetSrc_, patJets);

  edm::Handle<reco::GenParticleCollection> genParticles;
  if ( genParticleSrc_.label() != "" ) evt.getByLabel(genParticleSrc_, genParticles);

  double tauJetWeightSum = 0.;
  int patTauIndex_run1 = 0;
  for ( std::vector<pat::Tau>::const_iterator patTau = patTaus->begin();
	patTau != patTaus->end(); ++patTau, ++patTauIndex_run1 ) {

    if ( tauIndicesToPlot_.size() > 0 && (!isIndexed(patTauIndex_run1, tauIndicesToPlot_)) ) continue;

    if ( requireGenTauMatch_ && !matchesGenTau(*patTau) ) continue;

    tauJetWeightSum += getTauJetWeight<pat::Tau>(*patTau, tauJetWeightExtractors_);
  }

  //std::cout << " tauJetWeightSum = " << tauJetWeightSum << std::endl;
  //std::cout << " evtWeight = " << evtWeight << std::endl;

  if ( checkWeightConsistency_ && evtWeight != 1. && TMath::Abs(evtWeight - tauJetWeightSum) > 1.e-4 ) {
    edm::LogWarning("TauHistManager")
      << " Mismatch between tauJetWeightSum = " << tauJetWeightSum << " and evtWeight = " << evtWeight << " !!";
  }

  //std::cout << " patTaus.size = " << patTaus->size() << std::endl;
  hNumTaus_->Fill(patTaus->size(), evtWeight);

  int patTauIndex_run2 = 0;
  for ( std::vector<pat::Tau>::const_iterator patTau = patTaus->begin();
	patTau != patTaus->end(); ++patTau, ++patTauIndex_run2 ) {

    if ( tauIndicesToPlot_.size() > 0 && (!isIndexed(patTauIndex_run2, tauIndicesToPlot_)) ) continue;

    //bool isGenTauMatched = matchesGenTau(*patTau);
    //std::cout << " Pt = " << patTau->pt() << ", eta = " << patTau->eta() << ", phi = " << patTau->phi() << std::endl;
    //std::cout << " isGenTauMatched = " << isGenTauMatched << std::endl;

    if ( requireGenTauMatch_ && !matchesGenTau(*patTau) ) continue;

    double tauJetWeight = getTauJetWeight<pat::Tau>(*patTau, tauJetWeightExtractors_);
    double weight = getWeight(evtWeight, tauJetWeight, tauJetWeightSum);

    fillTauHistograms(*patTau, hTauPt_, hTauEta_, hTauPhi_, weight);
    hTauPtVsEta_->Fill(patTau->eta(), patTau->pt(), weight);
    hTauCharge_->Fill(patTau->charge(), weight);
    double jetRadius = TMath::Sqrt(patTau->etaetaMoment() + patTau->phiphiMoment());
    if ( !TMath::IsNaN(jetRadius) ) {
      hTauJetRadius_->Fill(jetRadius, weight);
/*

  CV: temporary work-around until MonitorElement::Fill(double, double, double) is fixed for TProfiles

      hTauJetRadiusPtProfile_->Fill(patTau->pt(), jetRadius, weight);
      hTauJetRadiusEnProfile_->Fill(patTau->energy(), jetRadius, weight);
 */
      hTauJetRadiusPtProfile_->getTProfile()->Fill(patTau->pt(), jetRadius, weight);
      hTauJetRadiusEnProfile_->getTProfile()->Fill(patTau->energy(), jetRadius, weight);
    }

    if ( patTau->genJet() != 0 ) {
      double recVisMass = patTau->mass();
      double genVisMass = patTau->genJet()->mass();

      const reco::GenParticle* genTau = 0;
      if ( genParticles.isValid() ) genTau = findGenParticle(patTau->p4(), *genParticles);
      std::string genTauDecayMode;
      if ( genTau ) genTauDecayMode = getGenTauDecayMode(genTau);
      else {
	edm::LogWarning("TauHistManager")
	  << " No genTau associated to pat::Tau --> calling deprecated JetMCTagUtils::genTauDecayMode function !!";
	genTauDecayMode = JetMCTagUtils::genTauDecayMode(*patTau->genJet());
      }

      fillTauVisMassHistogram(hTauVisMass_, hTauVisMassRes_, "all",
			      recVisMass, genVisMass, genTauDecayMode, weight);
      fillTauVisMassHistogram(hTauVisMassOneProngOnePi0_, hTauVisMassResOneProngOnePi0_, "oneProng1Pi0",
			      recVisMass, genVisMass, genTauDecayMode, weight);
      fillTauVisMassHistogram(hTauVisMassOneProngTwoPi0s_, hTauVisMassResOneProngTwoPi0s_, "oneProng2Pi0",
			      recVisMass, genVisMass, genTauDecayMode, weight);
      fillTauVisMassHistogram(hTauVisMassThreeProngNoPi0s_, hTauVisMassResThreeProngNoPi0s_, "threeProng0Pi0",
			      recVisMass, genVisMass, genTauDecayMode, weight);
      fillTauVisMassHistogram(hTauVisMassThreeProngOnePi0_, hTauVisMassResThreeProngOnePi0_, "threeProng1Pi0",
			      recVisMass, genVisMass, genTauDecayMode, weight);

      if ( (genTauDecayMode == "oneProng1Pi0"   ||
	    genTauDecayMode == "oneProng2Pi0"   ||
	    genTauDecayMode == "threeProng0Pi0") &&
	   (patTau->signalPFChargedHadrCands().size() == 1 ||
	    patTau->signalPFChargedHadrCands().size() == 3) ) {
	const reco::Candidate* recDistPion = getDistPion(*patTau);
	const reco::Candidate* genDistPion = getDistPion(*patTau->genJet());
	if ( recDistPion && genDistPion ) {
	  double distPionPtRes = (recDistPion->pt() - genDistPion->pt())/patTau->pt();
	  if      ( genTauDecayMode == "oneProng1Pi0"   ) hDistPionEnResOneProngOnePi0_->Fill(distPionPtRes, weight);
	  else if ( genTauDecayMode == "oneProng2Pi0"   ) hDistPionEnResOneProngTwoPi0s_->Fill(distPionPtRes, weight);
	  else if ( genTauDecayMode == "threeProng0Pi0" ) hDistPionEnResThreeProngNoPi0s_->Fill(distPionPtRes, weight);
	} else {
	  // CV: disable warning for now...
	  //edm::LogWarning("TauHistManager::fillHistogramsImp")
	  //  << " Failed to identify 'distinguishable' pion !!";
	}
      }
    }

    fillWeightHistograms(hTauJetWeightPosLog_, hTauJetWeightNegLog_, hTauJetWeightZero_,
			 hTauJetWeightLinear_, tauJetWeight);

//--- compare reconstructed tau-jet
//    to visible decay products on generator level;
//    normalize difference between reconstructed and generated energy
//    to expected energy dependence of resolution
    if ( patTau->genJet() ) {
      hGenVisTauPt_->Fill(patTau->genJet()->pt(), weight);
      hGenVisTauEta_->Fill(patTau->genJet()->eta(), weight);
      hGenVisTauPhi_->Fill(patTau->genJet()->phi(), weight);

      // Find mother tau
      const reco::GenParticle* firstDaughter =
        dynamic_cast<const reco::GenParticle*>(patTau->genJet()->daughter(0));
      if (firstDaughter) {
        const reco::GenParticle* genTau = findMotherWithPdgId(firstDaughter, 15);
        if (genTau) {
          hGenTauPt_->Fill(genTau->pt(), weight);
          hGenTauEta_->Fill(genTau->eta(), weight);
          hGenTauPhi_->Fill(genTau->phi(), weight);
          hGenTauX_->Fill(patTau->genJet()->energy()/genTau->energy(), weight);
        }
      }
      hTauEnCompToGen_->Fill((patTau->energy() - patTau->genJet()->energy())/patTau->genJet()->energy(), weight);
      hTauThetaCompToGen_->Fill(patTau->theta() - patTau->genJet()->theta(), weight);
      hTauPhiCompToGen_->Fill(patTau->phi() - patTau->genJet()->phi(), weight);
    }

    // get PDG IDs of matching generator particles
    int matchingGenParticlePdgId = 0;
    int matchingFinalStateGenParticlePdgId = 0;
    if ( genParticles.isValid() ) {
      matchingGenParticlePdgId = getMatchingGenParticlePdgId(patTau->p4(), *genParticles, &skipPdgIdsGenParticleMatch_, true);
      matchingFinalStateGenParticlePdgId = getMatchingGenParticlePdgId(patTau->p4(), *genParticles, &skipPdgIdsGenParticleMatch_, false);
    }

    if ( matchingGenParticlePdgId == 0 ) {
      hTauMatchingGenParticlePdgId_->Fill(-1, weight);
    } else if ( abs(matchingGenParticlePdgId) > 22 ) {
      hTauMatchingGenParticlePdgId_->Fill(24, weight);
    } else {
      hTauMatchingGenParticlePdgId_->Fill(TMath::Abs(matchingGenParticlePdgId), weight);
    }

    if ( matchingFinalStateGenParticlePdgId == 0 ) {
      hTauMatchingFinalStateGenParticlePdgId_->Fill(-1, weight);
    } else if ( abs(matchingFinalStateGenParticlePdgId) > 22 ) {
      hTauMatchingFinalStateGenParticlePdgId_->Fill(24, weight);
    } else {
      hTauMatchingFinalStateGenParticlePdgId_->Fill(TMath::Abs(matchingFinalStateGenParticlePdgId), weight);
    }

    std::string genTauDecayMode = "";
    if ( patTau->genJet() != 0 ) {
      genTauDecayMode = JetMCTagUtils::genTauDecayMode(*patTau->genJet());
    } else if ( TMath::Abs(matchingGenParticlePdgId) == 15 ) { // special handling of tau --> electron/muon decays
      if      ( TMath::Abs(matchingFinalStateGenParticlePdgId) == 11 ) genTauDecayMode = "electron";
      else if ( TMath::Abs(matchingFinalStateGenParticlePdgId) == 13 ) genTauDecayMode = "muon";
    }

    if ( genTauDecayMode != "" ) hTauMatchingGenTauDecayMode_->getTH1()->Fill(genTauDecayMode.data(), weight);

    hTauNumTracksSignalCone_->Fill(patTau->signalPFChargedHadrCands().size(), weight);
    hTauNumTracksIsoCone_->Fill(patTau->isolationTracks().size(), weight);

    hTauLeadPFChargedHadCandRefValidity_ ->Fill( isValidRef( patTau->leadPFChargedHadrCand() ) );
    if( isValidRef( patTau->leadPFChargedHadrCand() ) ) {
        hTauLeadKfTrkRefValidity_->Fill( isValidRef( patTau->leadPFChargedHadrCand()->trackRef() ) );
        hTauLeadGsfTrkRefValidity_->Fill( isValidRef( patTau->leadPFChargedHadrCand()->gsfTrackRef() ) );
    }

    if ( isValidRef(patTau->leadPFChargedHadrCand()) && isValidRef(patTau->leadPFChargedHadrCand()->trackRef()) ) {
      hTauLeadTrkPt_->Fill(patTau->leadPFChargedHadrCand()->trackRef()->pt(), weight);
      hTauLeadTrkEta_->Fill(patTau->leadPFChargedHadrCand()->trackRef()->eta(), weight);
      hTauLeadTrkPhi_->Fill(patTau->leadPFChargedHadrCand()->trackRef()->phi(), weight);

      hTauLeadTrkMatchDist_->Fill(reco::deltaR(patTau->leadPFChargedHadrCand()->trackRef()->momentum(), patTau->p4()), weight);

      if ( vertexSrc_.label() != "" ) {
	edm::Handle<std::vector<reco::Vertex> > recoVertices;
	evt.getByLabel(vertexSrc_, recoVertices);
	if ( recoVertices->size() >= 1 ) {
	  const reco::Vertex& thePrimaryEventVertex = (*recoVertices->begin());
	  hTauLeadTrkIPxy_->Fill(patTau->leadPFChargedHadrCand()->trackRef()->dxy(thePrimaryEventVertex.position()), weight);
	  hTauLeadTrkIPz_->Fill(patTau->leadPFChargedHadrCand()->trackRef()->dz(thePrimaryEventVertex.position()), weight);
	}
      }

      const reco::HitPattern& hitPattern = patTau->leadPFChargedHadrCand()->trackRef()->hitPattern();
      hTauLeadTrkNumHits_->Fill(hitPattern.numberOfValidTrackerHits(), weight);
      hTauLeadTrkNumPixelHits_->Fill(hitPattern.numberOfValidPixelHits(), weight);
      hTauLeadTrkNumStripHits_->Fill(hitPattern.numberOfValidStripHits(), weight);
    }

    static std::map<std::string, bool> discrAvailability_hasBeenChecked;

    if ( !useHPSpTaNCalgorithm_ && !useHPSclassicAlgorithm_ ) {
      fillTauDiscriminatorHistogram(hTauDiscriminatorByIsolation_, *patTau, "byIsolation",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorByTrackIsolation_, *patTau, "trackIsolation",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorByEcalIsolation_, *patTau, "ecalIsolation",
				    discrAvailability_hasBeenChecked, weight);
    }

    fillTauDiscriminatorHistogram(hTauDiscriminatorAgainstElectronsLoose_, *patTau, "againstElectronLoose",
				  discrAvailability_hasBeenChecked, weight);
    fillTauDiscriminatorHistogram(hTauDiscriminatorAgainstElectronsMedium_, *patTau, "againstElectronMedium",
				  discrAvailability_hasBeenChecked, weight);
    fillTauDiscriminatorHistogram(hTauDiscriminatorAgainstElectronsTight_, *patTau, "againstElectronTight",
				  discrAvailability_hasBeenChecked, weight);
    if ( patTau->leadPFCand().isNonnull() ) {
	double pfElectronMVA = patTau->leadPFCand()->mva_e_pi();
	if ( pfElectronMVA > +1.0 ) pfElectronMVA = +1.0;
	if ( pfElectronMVA < -1.0 ) pfElectronMVA = -1.0;
	hTauPFElectronMVA_->Fill(pfElectronMVA, weight);
    }
    hTauEmFraction_->Fill(patTau->emFraction(), weight);
    hTauHcalTotOverPLead_->Fill(patTau->hcalTotOverPLead(), weight);
    hTauHcalMaxOverPLead_->Fill(patTau->hcalMaxOverPLead(), weight);
    hTauHcal3x3OverPLead_->Fill(patTau->hcal3x3OverPLead(), weight);
    hTauEcalStripSumEOverPLead_->Fill(patTau->ecalStripSumEOverPLead(), weight);
    hTauBremsRecoveryEOverPLead_->Fill(patTau->bremsRecoveryEOverPLead(), weight);
    hTauCaloEOverPLead_->Fill(patTau->ecalStripSumEOverPLead() + patTau->hcalTotOverPLead(), weight);

    fillTauDiscriminatorHistogram(hTauDiscriminatorAgainstMuonsLoose_, *patTau, "againstMuonLoose",
				  discrAvailability_hasBeenChecked, weight);
    fillTauDiscriminatorHistogram(hTauDiscriminatorAgainstMuonsTight_, *patTau, "againstMuonTight",
				  discrAvailability_hasBeenChecked, weight);

    int recTauDecayMode = patTau->decayMode();
    hTauRecDecayMode_->Fill(recTauDecayMode, weight);
    if ( genTauDecayMode != "" ) {
      TH2* tauRecVsGenDecayMode_th2 = dynamic_cast<TH2*>(hTauRecVsGenDecayMode_->getTH1());
      tauRecVsGenDecayMode_th2->Fill(genTauDecayMode.data(), recTauDecayMode, weight);
    }

    if ( useHPSpTaNCalgorithm_ ) {
      fillTauDiscriminatorHistogram(hTauTaNCoutputTransform_, *patTau, "byTaNCtransform",
				    discrAvailability_hasBeenChecked, weight);

      fillTauDiscriminatorHistogram(hTauDiscriminatorCombinedTaNCvloose_, *patTau, "byTaNCvloose",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorCombinedTaNCloose_, *patTau, "byTaNCloose",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorCombinedTaNCmedium_, *patTau, "byTaNCmedium",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorCombinedTaNCtight_, *patTau, "byTaNCtight",
				    discrAvailability_hasBeenChecked, weight);

      fillTauDiscriminatorHistogram(hTauDiscriminatorCombinedHPSvloose_, *patTau, "byHPSvloose",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorCombinedHPSloose_, *patTau, "byHPSloose",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorCombinedHPSmedium_, *patTau, "byHPSmedium",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorCombinedHPStight_, *patTau, "byHPStight",
				    discrAvailability_hasBeenChecked, weight);
    
    } else if (useHPSclassicAlgorithm_) {
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPSvloose_, *patTau, "byVLooseIsolation",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPSloose_, *patTau, "byLooseIsolation",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPSmedium_, *patTau, "byMediumIsolation",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPStight_, *patTau, "byTightIsolation",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPSvlooseDeltaB_, *patTau, "byVLooseIsolationDeltaBetaCorr",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPSlooseDeltaB_, *patTau, "byLooseIsolationDeltaBetaCorr",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPSmediumDeltaB_, *patTau, "byMediumIsolationDeltaBetaCorr",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPStightDeltaB_, *patTau, "byTightIsolationDeltaBetaCorr",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPSvlooseCombDeltaB_, *patTau, "byVLooseCombinedIsolationDeltaBetaCorr",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPSlooseCombDeltaB_, *patTau, "byLooseCombinedIsolationDeltaBetaCorr",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPSmediumCombDeltaB_, *patTau, "byMediumCombinedIsolationDeltaBetaCorr",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorHPStightCombDeltaB_, *patTau, "byTightCombinedIsolationDeltaBetaCorr",
				    discrAvailability_hasBeenChecked, weight);
    
    } else {
      fillTauDiscriminatorHistogram(hTauDiscriminatorTaNCfrOnePercent_, *patTau, "byTaNCfrOnePercent",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorTaNCfrHalfPercent_, *patTau, "byTaNCfrHalfPercent",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorTaNCfrQuarterPercent_, *patTau, "byTaNCfrQuarterPercent",
				    discrAvailability_hasBeenChecked, weight);
      fillTauDiscriminatorHistogram(hTauDiscriminatorTaNCfrTenthPercent_, *patTau, "byTaNCfrTenthPercent",
				    discrAvailability_hasBeenChecked, weight);
    }

    MonitorElement* hTauTaNCoutput = 0;
    if ( recTauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion0PiZero ) {
      hTauTaNCoutput = hTauTaNCoutputOneProngNoPi0s_;
    } else if ( recTauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero ) {
      hTauTaNCoutput = hTauTaNCoutputOneProngOnePi0_;
    } else if ( recTauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion2PiZero ) {
      hTauTaNCoutput = hTauTaNCoutputOneProngTwoPi0s_;
    } else if ( recTauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion0PiZero ) {
      hTauTaNCoutput = hTauTaNCoutputThreeProngNoPi0s_;
    } else if ( recTauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion1PiZero ) {
      hTauTaNCoutput = hTauTaNCoutputThreeProngOnePi0_;
    }
    if ( hTauTaNCoutput ) {
      fillTauDiscriminatorHistogram(hTauTaNCoutput, *patTau, "byTaNC",
				    discrAvailability_hasBeenChecked, weight);
    }

    fillTauIsoHistograms(*patTau, weight);
    hTauDeltaRnearestJet_->Fill(getDeltaRnearestJet(patTau->p4(), patJets), weight);
    if ( makeIsoPtConeSizeDepHistograms_ ) fillTauIsoConeSizeDepHistograms(*patTau, weight);

    fillTauIdEfficiencyHistograms(*patTau, weight, hTauIdEfficiencies_, makeTauIdEfficiencyHistograms_);
    fillTauIdEfficiencyHistograms(*patTau, weight, hTauFakeRates_, makeTauFakeRateHistograms_);
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void TauHistManager::bookTauHistograms(MonitorElement*& hTauPt, MonitorElement*& hTauEta, MonitorElement*& hTauPhi, const char* histoSetName)
{
  std::string hTauPtName = std::string(histoSetName).append("Pt");
  hTauPt = book1D(hTauPtName, hTauPtName, 75, 0., 150.);

  std::string hTauEtaName = std::string(histoSetName).append("Eta");
  hTauEta = book1D(hTauEtaName, hTauEtaName, 60, -3., +3.);

  std::string hTauPhiName = std::string(histoSetName).append("Phi");
  hTauPhi = book1D(hTauPhiName, hTauPhiName, 36, -TMath::Pi(), +TMath::Pi());
}

void TauHistManager::bookTauIsoConeSizeDepHistograms()
{
  for ( unsigned iConeSize = 1; iConeSize <= numTauIsoConeSizes_; ++iConeSize ) {
    std::ostringstream iConeSizeString;
    iConeSizeString << std::setfill('0') << std::setw(2) << iConeSize;

    std::string hTauParticleFlowIsoPtConeSizeDepName_i
      = std::string("TauParticleFlowIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hTauParticleFlowIsoPtConeSizeDep_.push_back(book1D(hTauParticleFlowIsoPtConeSizeDepName_i,
						       hTauParticleFlowIsoPtConeSizeDepName_i, 40, 0., 10.));
    std::string hTauPFChargedHadronIsoPtConeSizeDepName_i
      = std::string("TauChargedHadronIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hTauPFChargedHadronIsoPtConeSizeDep_.push_back(book1D(hTauPFChargedHadronIsoPtConeSizeDepName_i,
							  hTauPFChargedHadronIsoPtConeSizeDepName_i, 40, 0., 10.));
    std::string hTauPFNeutralHadronIsoPtConeSizeDepName_i
      = std::string("TauPFNeutralHadronIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hTauPFNeutralHadronIsoPtConeSizeDep_.push_back(book1D(hTauPFNeutralHadronIsoPtConeSizeDepName_i,
							  hTauPFNeutralHadronIsoPtConeSizeDepName_i, 40, 0., 10.));
    std::string hTauPFGammaIsoPtConeSizeDepName_i
      = std::string("TauPFGammaIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hTauPFGammaIsoPtConeSizeDep_.push_back(book1D(hTauPFGammaIsoPtConeSizeDepName_i,
						  hTauPFGammaIsoPtConeSizeDepName_i, 40, 0., 10.));
  }
}

void TauHistManager::bookTauIdEfficiencyHistograms(std::vector<MonitorElement*>& histograms,
						   const char* type, const vstring& labels)
{
  for ( vstring::const_iterator label = labels.begin();
	label != labels.end(); ++label ) {
    std::string histogramName = std::string("Tau").append(type).append(" (").append(*label).append(")");
    histograms.push_back(book1D(histogramName, histogramName, 102, -0.01, 1.01));
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void TauHistManager::fillTauHistograms(const pat::Tau& patTau, MonitorElement* hTauPt, MonitorElement* hTauEta, MonitorElement* hTauPhi,
				       double weight)
{
  //std::cout << "<TauHistManager::fillTauHistograms>:" << std::endl;

  hTauPt->Fill(patTau.pt(), weight);
  hTauEta->Fill(patTau.eta(), weight);
  hTauPhi->Fill(patTau.phi(), weight);
}

void TauHistManager::fillTauIsoHistograms(const pat::Tau& patTau, double weight)
{
  //std::cout << "<TauHistManager::fillTauIsoHistograms>:" << std::endl;

  hTauTrkIsoPt_->Fill(patTau.trackIso(), weight);
  hTauEcalIsoPt_->Fill(patTau.ecalIso(), weight);
  hTauHcalIsoPt_->Fill(patTau.hcalIso(), weight);
  //hTauIsoSumPt_->Fill(patTau.trackIso() + patTau.ecalIso() + patTau.hcalIso(), weight);
  hTauIsoSumPt_->Fill(patTau.trackIso() + patTau.ecalIso(), weight);

  //std::cout << " particleIso = " << patTau.particleIso() << std::endl;
  //std::cout << " chargedHadronIso = " << patTau.chargedHadronIso() << std::endl;
  //std::cout << " neutralHadronIso = " << patTau.neutralHadronIso() << std::endl;
  //std::cout << " photonIso = " << patTau.photonIso() << std::endl;

  hTauParticleFlowIsoPt_->Fill(patTau.particleIso(), weight);
  hTauPFChargedHadronIsoPt_->Fill(patTau.chargedHadronIso(), weight);
  hTauPFNeutralHadronIsoPt_->Fill(patTau.neutralHadronIso(), weight);
  hTauPFGammaIsoPt_->Fill(patTau.photonIso(), weight);

  hTauNumSignalPFChargedHadrons_->Fill(patTau.signalPFChargedHadrCands().size(), weight);
  hTauNumIsoPFChargedHadrons_->Fill(patTau.isolationPFChargedHadrCands().size(), weight);
  hTauNumPFChargedHadrons_->Fill(patTau.signalPFChargedHadrCands().size() + patTau.isolationPFChargedHadrCands().size(), weight);

  hTauNumSignalPFGammas_->Fill(patTau.signalPFGammaCands().size(), weight);
  hTauNumIsoPFGammas_->Fill(patTau.isolationPFGammaCands().size(), weight);
  hTauNumPFGammas_->Fill(patTau.signalPFGammaCands().size() + patTau.isolationPFGammaCands().size(), weight);

  if ( makeIsoPtCtrlHistograms_ ) {
    double sumPtIsolationConePFChargedHadrons = 0.;
    for ( reco::PFCandidateRefVector::const_iterator pfChargedHadron = patTau.isolationPFChargedHadrCands().begin();
	  pfChargedHadron != patTau.isolationPFChargedHadrCands().end(); ++pfChargedHadron ) {
      if ( (*pfChargedHadron)->pt() > 1.0 ) sumPtIsolationConePFChargedHadrons += (*pfChargedHadron)->pt();
    }
    hTauPFChargedHadronIsoPtCtrl_->Fill(sumPtIsolationConePFChargedHadrons, patTau.chargedHadronIso(), weight);

    double sumPtIsolationConePFGammas = 0.;
    for ( reco::PFCandidateRefVector::const_iterator pfGamma = patTau.isolationPFGammaCands().begin();
	  pfGamma != patTau.isolationPFGammaCands().end(); ++pfGamma ) {
      if ( (*pfGamma)->pt() > 1.5 ) sumPtIsolationConePFGammas += (*pfGamma)->pt();
    }
    hTauPFGammaIsoPtCtrl_->Fill(sumPtIsolationConePFGammas, patTau.photonIso(), weight);
  }

  for ( reco::TrackRefVector::const_iterator isolationTrack = patTau.isolationTracks().begin();
	isolationTrack != patTau.isolationTracks().end(); ++isolationTrack ) {
    hTauTrkIsoEnProfile_->Fill((*isolationTrack)->p(), weight);
    hTauTrkIsoPtProfile_->Fill((*isolationTrack)->pt(), weight);
    hTauTrkIsoEtaDistProfile_->Fill(TMath::Abs(patTau.eta() - (*isolationTrack)->eta()), weight);
    hTauTrkIsoPhiDistProfile_->Fill(TMath::Abs(patTau.phi() - (*isolationTrack)->phi()), weight);
  }
}

void TauHistManager::fillTauIsoConeSizeDepHistograms(const pat::Tau& patTau, double weight)
{
  //std::cout << "<TauHistManager::fillTauIsoConeSizeDepHistograms>:" << std::endl;

  for ( unsigned iConeSize = 1; iConeSize <= numTauIsoConeSizes_; ++iConeSize ) {
    double isoConeSize_i = iConeSize*tauIsoConeSizeIncr_;

    if ( patTau.isoDeposit(pat::PfAllParticleIso) ) {
      double tauParticleFlowIsoDeposit_i
	= patTau.isoDeposit(pat::PfAllParticleIso)->countWithin(isoConeSize_i, tauParticleFlowIsoParam_, false);
      hTauParticleFlowIsoPtConeSizeDep_[iConeSize - 1]->Fill(tauParticleFlowIsoDeposit_i, weight);
    }

    if ( patTau.isoDeposit(pat::PfChargedHadronIso) ) {
      double tauPFChargedHadronIsoDeposit_i
	= patTau.isoDeposit(pat::PfChargedHadronIso)->countWithin(isoConeSize_i, tauParticleFlowIsoParam_, false);
      hTauPFChargedHadronIsoPtConeSizeDep_[iConeSize - 1]->Fill(tauPFChargedHadronIsoDeposit_i, weight);
    }

    if ( patTau.isoDeposit(pat::PfNeutralHadronIso) ) {
      double tauPFNeutralHadronIsoDeposit_i
	= patTau.isoDeposit(pat::PfNeutralHadronIso)->countWithin(isoConeSize_i, tauParticleFlowIsoParam_, false);
      hTauPFNeutralHadronIsoPtConeSizeDep_[iConeSize - 1]->Fill(tauPFNeutralHadronIsoDeposit_i, weight);
    }

    if ( patTau.isoDeposit(pat::PfGammaIso) ) {
      double tauPFGammaIsoDeposit_i
	= patTau.isoDeposit(pat::PfGammaIso)->countWithin(isoConeSize_i, tauParticleFlowIsoParam_, false);
      hTauPFGammaIsoPtConeSizeDep_[iConeSize - 1]->Fill(tauPFGammaIsoDeposit_i, weight);
    }
  }
}

void TauHistManager::fillTauIdEfficiencyHistograms(const pat::Tau& patTau, double weight,
						   std::vector<MonitorElement*>& histograms, const vstring& labels)
{
  assert(histograms.size() == labels.size());

  unsigned numHistograms = histograms.size();
  for ( unsigned iHistogram = 0; iHistogram < numHistograms; ++iHistogram ) {
    histograms[iHistogram]->Fill(patTau.efficiency(labels[iHistogram].data()).value(), weight);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, TauHistManager, "TauHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, TauHistManager, "TauHistManager");

#include "TauAnalysis/Core/interface/HistManagerAdapter.h"

typedef HistManagerAdapter<TauHistManager> TauAnalyzer;

DEFINE_FWK_MODULE(TauAnalyzer);
