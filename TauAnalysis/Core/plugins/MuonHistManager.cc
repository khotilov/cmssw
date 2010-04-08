#include "TauAnalysis/Core/plugins/MuonHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/IsolationAlgos/interface/IsoDepositVetoFactory.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"
#include "TauAnalysis/Core/interface/eventAuxFunctions.h"
#include "TauAnalysis/GenSimTools/interface/genParticleAuxFunctions.h"

#include <TMath.h>

#include <stdlib.h>

bool matchesGenMuon(const pat::Muon& patMuon)
{
  //std::cout << "<matchesGenMuon>:" << std::endl;

  bool isGenMuonMatched = false;
  std::vector<reco::GenParticleRef> associatedGenParticles = patMuon.genParticleRefs();
  for ( std::vector<reco::GenParticleRef>::const_iterator it = associatedGenParticles.begin(); 
	it != associatedGenParticles.end(); ++it ) {
    if ( isValidRef(*it) ) {
      const reco::GenParticleRef& genParticle = (*it);
      if ( genParticle->pdgId() == -13 || genParticle->pdgId() == +13 ) isGenMuonMatched = true;
    } else {
      edm::LogWarning("matchesGenMuon") << " edm::Ref of genParticle associated to pat::Muon is invalid !!";
    }
  }
  return isGenMuonMatched;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

MuonHistManager::MuonHistManager(const edm::ParameterSet& cfg)
  : HistManagerBase(cfg)
{
  //std::cout << "<MuonHistManager::MuonHistManager>:" << std::endl;

  muonSrc_ = cfg.getParameter<edm::InputTag>("muonSource");
  //std::cout << " muonSrc = " << muonSrc_ << std::endl;

  vertexSrc_ = ( cfg.exists("vertexSource") ) ? cfg.getParameter<edm::InputTag>("vertexSource") : edm::InputTag();
  if ( vertexSrc_.label() == "" ) {
    edm::LogWarning("MuonHistManager") << " Configuration parameter 'vertexSource' not specified" 
				       << " --> Impact Parameter histograms will NOT be plotted !!";
  }
  //std::cout << " vertexSrc = " << vertexSrc_ << std::endl;

  jetSrc_ = cfg.getParameter<edm::InputTag>("jetSource");
  //std::cout << " jetSrc = " << jetSrc_ << std::endl;

  genParticleSrc_ = ( cfg.exists("genParticleSource") ) ? cfg.getParameter<edm::InputTag>("genParticleSource") : edm::InputTag();
  if ( genParticleSrc_.label() == "" ) {
    edm::LogWarning("MuonHistManager") 
      << " Configuration parameter 'genParticleSource' not specified" 
      << " --> matching gen. Particle PdgId histogram will NOT be plotted !!";
  }
  //std::cout << " genParticleSrc = " << genParticleSrc_ << std::endl;

  requireGenMuonMatch_ = cfg.getParameter<bool>("requireGenMuonMatch");
  //std::cout << " requireGenMuonMatch = " << requireGenMuonMatch_ << std::endl;

  skipPdgIdsGenParticleMatch_ = cfg.getParameter<vint>("skipPdgIdsGenParticleMatch");

  std::string normalization_string = cfg.getParameter<std::string>("normalization");
  normMethod_ = getNormMethod(normalization_string, "muons");

  makeIsoPtCtrlHistograms_ = ( cfg.exists("makeIsoPtCtrlHistograms") ) ? 
    cfg.getParameter<bool>("makeIsoPtCtrlHistograms") : false;

  makeIsoPtConeSizeDepHistograms_ = ( cfg.exists("makeIsoPtConeSizeDepHistograms") ) ? 
    cfg.getParameter<bool>("makeIsoPtConeSizeDepHistograms") : false;

  numMuonIsoConeSizes_ = 15;
  muonIsoConeSizeIncr_ = 0.1;
  numMuonIsoPtThresholds_ = 4;
  muonIsoPtThresholdIncr_ = 0.5;

//--- create "veto" objects for computation of IsoDeposit sums
  muonTrkIsoParam_.push_back(IsoDepositVetoFactory::make("0.02"));
  muonTrkIsoParam_.push_back(IsoDepositVetoFactory::make("Threshold(1.0)"));
  
  muonEcalIsoParam_.push_back(IsoDepositVetoFactory::make("0.0"));
  muonEcalIsoParam_.push_back(IsoDepositVetoFactory::make("Threshold(0.0)"));

  muonHcalIsoParam_.push_back(IsoDepositVetoFactory::make("0.0"));
  muonHcalIsoParam_.push_back(IsoDepositVetoFactory::make("Threshold(0.0)"));

  muonParticleFlowIsoParam_.push_back(IsoDepositVetoFactory::make("0.0"));
  muonParticleFlowIsoParam_.push_back(IsoDepositVetoFactory::make("Threshold(0.5)"));
}

MuonHistManager::~MuonHistManager()
{
//--- delete "veto" objects for computation of IsoDeposit sums
  clearIsoParam(muonTrkIsoParam_);
  clearIsoParam(muonEcalIsoParam_);
  clearIsoParam(muonHcalIsoParam_);
  clearIsoParam(muonParticleFlowIsoParam_);
}

void MuonHistManager::bookHistogramsImp()
{
  //std::cout << "<MuonHistManager::bookHistogramsImp>:" << std::endl;

//--- book histogram for number of muons in each event
  hNumMuons_ = book1D("NumMuons", "NumMuons", 10, -0.5, 9.5);

//--- book histograms for Pt, eta and phi distributions
//    of muons passing all id. and isolation selections
  bookMuonHistograms(hMuonPt_, hMuonEta_, hMuonPhi_, "Muon");
  hMuonPtVsEta_ = book2D("MuonPtVsEta", "MuonPtVsEta", 24, -3., +3., 30, 0., 150.);
  float muonAbsEtaBinning2d[] = { 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1 };
  float muonPtBinning2d[] = { 15., 20., 25., 30., 40., 60., 120. };
  hMuonPtVsAbsEta_ = book2D("MuonPtVsAbsEta", "MuonPtVsAbsEta", 7, muonAbsEtaBinning2d, 6, muonPtBinning2d);
  hMuonCharge_ = book1D("MuonCharge", "Muon Charge", 3, -1.5, +1.5);
  
  bookWeightHistograms(*dqmStore_, "MuonWeight", "Muon Weight", 
		       hMuonWeightPosLog_, hMuonWeightNegLog_, hMuonWeightZero_, 
		       hMuonWeightLinear_);
  
  hMuonPtCompToGen_ = book1D("MuonPtCompToGen", "MuonPtCompToGen", 200, -0.10, +0.10);
  hMuonThetaCompToGen_ = book1D("MuonThetaCompToGen", "MuonThetaCompToGen", 200, -0.010, +0.010);
  hMuonPhiCompToGen_ = book1D("MuonPhiCompToGen", "MuonPhiCompToGen", 200, -0.010, +0.010);
  
  hMuonMatchingGenParticlePdgId_ = book1D("MuonMatchingGenParticlePdgId", "matching gen. Particle PdgId", 26, -1.5, 24.5);

  hMuonTrackIPxy_ = book1D("MuonTrackIPxy", "MuonTrackIPxy", 100, -0.100, 0.100);
  hMuonTrackIPz_ = book1D("MuonTrackIPz", "MuonTrackIPz", 100, -1.0, 1.0);
  
  hMuonEcalDeposits_ = book1D("MuonEcalDeposits", "MuonEcalDeposits", 100, 0., 10.);
  hMuonHcalDeposits_ = book1D("MuonHcalDeposits", "MuonHcalDeposits", 100, 0., 10.);
  hMuonCaloDeposits_ = book1D("MuonCaloDeposits", "MuonCaloDeposits", 100, 0., 10.);
  hMuonCaloCompatibility_ = book1D("MuonCaloCompatibility", "MuonCaloCompatibility", 102, -0.01, 1.01);
  
  hMuonNumberOfChambers_ = book1D("MuonNumberOfChambers", "MuonNumberOfChambers", 25, -0.5, 24.5);
  hMuonSegmentCompatibility_ = book1D("MuonSegmentCompatibility", "MuonSegmentCompatibility", 102, -0.01, 1.01);
  
  hMuonTrkIsoPt_ = book1D("MuonTrkIsoPt", "MuonTrkIsoPt", 100, 0., 10.);    
  hMuonTrkIsoPtVsMuonPt_ = book2D("MuonTrkIsoPtVsMuonPt", "MuonTrkIsoPtVsMuonPt", 20, 0., 100., 20, 0., 10.); 
  hMuonEcalIsoPt_ = book1D("MuonEcalIsoPt", "MuonEcalIsoPt", 100, 0., 10.);
  hMuonEcalIsoPtVsMuonPt_ = book2D("MuonEcalIsoPtVsMuonPt", "MuonEcalIsoPtVsMuonPt", 20, 0., 100., 20, 0., 10.); 
  hMuonHcalIsoPt_ = book1D("MuonHcalIsoPt", "MuonHcalIsoPt", 100, 0., 10.);
  hMuonIsoSumPt_ = book1D("MuonIsoSumPt", "MuonIsoSumPt", 100, 0., 10.);
  hMuonIsoSumPtVsMuonPt_ = book2D("MuonIsoSumPtVsMuonPt", "MuonIsoSumPtVsMuonPt", 20, 0., 100., 20, 0., 10.); 
  hMuonTrkIsoPtRel_ = book1D("MuonTrkIsoPtRel", "MuonTrkIsoPtRel", 100, 0., 1.);    
  hMuonEcalIsoPtRel_ = book1D("MuonEcalIsoPtRel", "MuonEcalIsoPtRel", 100, 0., 1.);
  hMuonHcalIsoPtRel_ = book1D("MuonHcalIsoPtRel", "MuonHcalIsoPtRel", 100, 0., 1.);
  hMuonIsoSumPtRel_ = book1D("MuonIsoSumPtRel", "MuonIsoSumPtRel", 100, 0., 1.);

  hMuonDeltaRnearestJet_ = book1D("MuonDeltaRnearestJet", "#DeltaR(nearest Jet)", 102, -0.1, 10.1);

  hMuonParticleFlowIsoPt_ = book1D("MuonParticleFlowIsoPt", "MuonParticleFlowIsoPt", 100, 0., 10.);    
  hMuonPFChargedHadronIsoPt_ = book1D("MuonPFChargedHadronIsoPt", "MuonPFChargedHadronIsoPt", 100, 0., 10.);   
  hMuonPFNeutralHadronIsoPt_ = book1D("MuonPFNeutralHadronIsoPt", "MuonPFNeutralHadronIsoPt", 100, 0., 10.);   
  hMuonPFGammaIsoPt_ = book1D("MuonPFGammaIsoPt", "MuonPFGammaIsoPt", 100, 0., 10.);  
  
//--- book "control" histograms to check agreement between muon isolation variables
//    computed by PAT-level IsoDeposits based on particle flow 
//    with values computed on AOD level, based on ECAL recHits/CaloTowers and reco::Tracks
  if ( makeIsoPtCtrlHistograms_ ) {
    hMuonPFChargedHadronIsoPtCtrl_ = book2D("MuonPFChargedHadronIsoPtCtrl", "PFChargedHadron vs. reco::Track Isolation P_{T}", 40, 0., 20., 40, 0., 20.);
    hMuonPFGammaIsoPtCtrl_ = book2D("MuonPFGammaIsoPtCtrl", "PFGamma vs. ECAL recHit/CaloTower Isolation P_{T}", 40, 0., 20., 40, 0., 20.);
  } else {
    hMuonPFChargedHadronIsoPtCtrl_ = 0;
    hMuonPFGammaIsoPtCtrl_ = 0;
  } 

  hMuonTrkIsoValProfile_ = book1D("MuonTrkIsoValProfile", "MuonTrkIsoValProfile", 100, 0., 10.);
  hMuonTrkIsoEtaDistProfile_ = book1D("MuonTrkIsoEtaDistProfile", "MuonTrkIsoEtaDistProfile", 15, 0., 1.5);
  hMuonTrkIsoPhiDistProfile_ = book1D("MuonTrkIsoPhiDistProfile", "MuonTrkIsoPhiDistProfile", 15, 0., 1.5);
  
  hMuonEcalIsoValProfile_ = book1D("MuonEcalIsoValProfile", "MuonEcalIsoValProfile", 100, 0., 10.);
  hMuonEcalIsoEtaDistProfile_ = book1D("MuonEcalIsoEtaDistProfile", "MuonEcalIsoEtaDistProfile", 15, 0., 1.5);
  hMuonEcalIsoPhiDistProfile_ = book1D("MuonEcalIsoPhiDistProfile", "MuonEcalIsoPhiDistProfile", 15, 0., 1.5);
  
  hMuonHcalIsoValProfile_ = book1D("MuonHcalIsoValProfile", "MuonHcalIsoValProfile", 100, 0., 10.);
  hMuonHcalIsoEtaDistProfile_ = book1D("MuonHcalIsoEtaDistProfile", "MuonHcalIsoEtaDistProfile", 15, 0., 1.5);
  hMuonHcalIsoPhiDistProfile_  = book1D("MuonHcalIsoPhiDistProfile", "MuonHcalIsoPhiDistProfile", 15, 0., 1.5);
  
  if ( makeIsoPtConeSizeDepHistograms_ ) bookMuonIsoConeSizeDepHistograms();
}

double MuonHistManager::getMuonWeight(const pat::Muon& patMuon)
{
  return 1.;
}

void MuonHistManager::fillHistogramsImp(const edm::Event& evt, const edm::EventSetup& es, double evtWeight)
{  
  //std::cout << "<MuonHistManager::fillHistogramsImp>:" << std::endl; 

  edm::Handle<pat::MuonCollection> patMuons;
  getCollection(evt, muonSrc_, patMuons);

  edm::Handle<pat::JetCollection> patJets;
  getCollection(evt, jetSrc_, patJets);

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(genParticleSrc_, genParticles);

  //std::cout << " patMuons.size = " << patMuons->size() << std::endl;
  hNumMuons_->Fill(patMuons->size(), evtWeight);

  double muonWeightSum = 0.;
  for ( std::vector<pat::Muon>::const_iterator patMuon = patMuons->begin(); 
	patMuon != patMuons->end(); ++patMuon ) {
    if ( requireGenMuonMatch_ && !matchesGenMuon(*patMuon) ) continue;

    muonWeightSum += getMuonWeight(*patMuon);
  }

  for ( std::vector<pat::Muon>::const_iterator patMuon = patMuons->begin(); 
	patMuon != patMuons->end(); ++patMuon ) {
  
    //bool isGenMuonMatched = matchesGenMuon(*patMuon);
    //std::cout << " Pt = " << patMuon->pt() << ", eta = " << patMuon->eta() << ", phi = " << patMuon->phi() << std::endl;
    //std::cout << " isGenMuonMatched = " << isGenMuonMatched << std::endl;
    
    if ( requireGenMuonMatch_ && !matchesGenMuon(*patMuon) ) continue;
    
    double muonWeight = getMuonWeight(*patMuon);
    double weight = getWeight(evtWeight, muonWeight, muonWeightSum);

    fillMuonHistograms(*patMuon, hMuonPt_, hMuonEta_, hMuonPhi_, weight);
    hMuonPtVsEta_->Fill(patMuon->eta(), patMuon->pt(), weight);
    hMuonPtVsAbsEta_->Fill(TMath::Abs(patMuon->eta()), patMuon->pt(), weight); 
    hMuonCharge_->Fill(patMuon->charge(), weight);

    fillWeightHistograms(hMuonWeightPosLog_, hMuonWeightNegLog_, hMuonWeightZero_, 
			 hMuonWeightLinear_, muonWeight);

//--- compare reconstructed muon to generator level one;
//    normalize difference between reconstructed and generated Pt
//    to expected Pt dependence of resolution
    if ( patMuon->genLepton() ) {
      hMuonPtCompToGen_->Fill((patMuon->pt() - patMuon->genLepton()->pt())/patMuon->genLepton()->pt(), weight);
      hMuonThetaCompToGen_->Fill(patMuon->theta() - patMuon->genLepton()->theta(), weight);
      hMuonPhiCompToGen_->Fill(patMuon->phi() - patMuon->genLepton()->phi(), weight);
    }

    int matchingGenParticlePdgId = getMatchingGenParticlePdgId(patMuon->p4(), genParticles, &skipPdgIdsGenParticleMatch_);
    if ( matchingGenParticlePdgId == -1 ) {
      hMuonMatchingGenParticlePdgId_->Fill(-1, weight);
    } else if ( abs(matchingGenParticlePdgId) > 22 ) {
      hMuonMatchingGenParticlePdgId_->Fill(24, weight);
    } else {
      hMuonMatchingGenParticlePdgId_->Fill(abs(matchingGenParticlePdgId), weight);
    }

    if ( vertexSrc_.label() != "" && isValidRef(patMuon->track()) ) {
      edm::Handle<std::vector<reco::Vertex> > recoVertices;
      evt.getByLabel(vertexSrc_, recoVertices);
      if ( recoVertices->size() >= 1 ) {
	const reco::Vertex& thePrimaryEventVertex = (*recoVertices->begin());
	hMuonTrackIPxy_->Fill(patMuon->track()->dxy(thePrimaryEventVertex.position()), weight);
	hMuonTrackIPz_->Fill(patMuon->track()->dz(thePrimaryEventVertex.position()), weight);
      }
    }

    if ( patMuon->isEnergyValid() ) {
      hMuonEcalDeposits_->Fill(patMuon->calEnergy().em, weight);
      hMuonHcalDeposits_->Fill(patMuon->calEnergy().had + patMuon->calEnergy().ho, weight);
      hMuonCaloDeposits_->Fill(patMuon->calEnergy().em + patMuon->calEnergy().had + patMuon->calEnergy().ho, weight);
    }
    if ( patMuon->isCaloCompatibilityValid() ) hMuonCaloCompatibility_->Fill(patMuon->caloCompatibility(), weight);
      
    hMuonNumberOfChambers_->Fill(patMuon->numberOfChambers(), weight);
    double segmentCompatibility = muon::segmentCompatibility(*patMuon);
    hMuonSegmentCompatibility_->Fill(segmentCompatibility, weight);

    fillMuonIsoHistograms(*patMuon, weight);
    hMuonDeltaRnearestJet_->Fill(getDeltaRnearestJet(patMuon->p4(), patJets), weight);
    if ( makeIsoPtConeSizeDepHistograms_ ) fillMuonIsoConeSizeDepHistograms(*patMuon, weight);
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void MuonHistManager::bookMuonHistograms(MonitorElement*& hMuonPt, MonitorElement*& hMuonEta, MonitorElement*& hMuonPhi, const char* histoSetName)
{
  std::string hMuonPtName = std::string(histoSetName).append("Pt");
  hMuonPt = book1D(hMuonPtName, hMuonPtName, 75, 0., 150.);
  
  std::string hMuonEtaName = std::string(histoSetName).append("Eta");
  hMuonEta = book1D(hMuonEtaName, hMuonEtaName, 60, -3., +3.);
  
  std::string hMuonPhiName = std::string(histoSetName).append("Phi");
  hMuonPhi = book1D(hMuonPhiName, hMuonPhiName, 36, -TMath::Pi(), +TMath::Pi());
}

void MuonHistManager::bookMuonIsoConeSizeDepHistograms()
{
  for ( unsigned iConeSize = 1; iConeSize <= numMuonIsoConeSizes_; ++iConeSize ) {
    std::ostringstream iConeSizeString;
    iConeSizeString << std::setfill('0') << std::setw(2) << iConeSize;
    
    std::string hMuonTrkIsoPtConeSizeDepName_i 
      = std::string("MuonTrkIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonTrkIsoPtConeSizeDep_.push_back(book1D(hMuonTrkIsoPtConeSizeDepName_i, 
					       hMuonTrkIsoPtConeSizeDepName_i, 40, 0., 10.));
    std::string hMuonEcalIsoPtConeSizeDepName_i 
      = std::string("MuonEcalIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonEcalIsoPtConeSizeDep_.push_back(book1D(hMuonEcalIsoPtConeSizeDepName_i, 
						hMuonEcalIsoPtConeSizeDepName_i, 40, 0., 10.));
    
    std::string hMuonHcalIsoPtConeSizeDepName_i 
      = std::string("MuonHcalIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonHcalIsoPtConeSizeDep_.push_back(book1D(hMuonHcalIsoPtConeSizeDepName_i, 
						hMuonHcalIsoPtConeSizeDepName_i, 40, 0., 10.));
    
    std::string hMuonParticleFlowIsoPtConeSizeDepName_i 
      = std::string("MuonParticleFlowIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonParticleFlowIsoPtConeSizeDep_.push_back(book1D(hMuonParticleFlowIsoPtConeSizeDepName_i, 
							hMuonParticleFlowIsoPtConeSizeDepName_i, 40, 0., 10.));
    
    std::string hMuonPFChargedHadronIsoPtConeSizeDepName_i 
      = std::string("MuonChargedHadronIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonPFChargedHadronIsoPtConeSizeDep_.push_back(book1D(hMuonPFChargedHadronIsoPtConeSizeDepName_i, 
							   hMuonPFChargedHadronIsoPtConeSizeDepName_i, 40, 0., 10.));
    
    std::string hMuonPFNeutralHadronIsoPtConeSizeDepName_i 
      = std::string("MuonPFNeutralHadronIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonPFNeutralHadronIsoPtConeSizeDep_.push_back(book1D(hMuonPFNeutralHadronIsoPtConeSizeDepName_i, 
							   hMuonPFNeutralHadronIsoPtConeSizeDepName_i, 40, 0., 10.));
    
    std::string hMuonPFGammaIsoPtConeSizeDepName_i 
      = std::string("MuonPFGammaIsoPtConeSizeDep").append("_").append(iConeSizeString.str());
    hMuonPFGammaIsoPtConeSizeDep_.push_back(book1D(hMuonPFGammaIsoPtConeSizeDepName_i, 
						   hMuonPFGammaIsoPtConeSizeDepName_i, 40, 0., 10.));
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void MuonHistManager::fillMuonHistograms(const pat::Muon& patMuon, 
					 MonitorElement* hMuonPt, MonitorElement* hMuonEta, MonitorElement* hMuonPhi, 
					 double weight)
{
  //std::cout << "<MuonHistManager::fillMuonHistograms>:" << std::endl;

  hMuonPt->Fill(patMuon.pt(), weight);
  hMuonEta->Fill(patMuon.eta(), weight);
  hMuonPhi->Fill(patMuon.phi(), weight);
}

void MuonHistManager::fillMuonIsoHistograms(const pat::Muon& patMuon, double weight)
{
  //std::cout << "<MuonHistManager::fillMuonIsoHistograms>:" << std::endl;

  double muonTrackIso = patMuon.userIsolation(pat::TrackIso);
  double muonEcalIso = patMuon.userIsolation(pat::EcalIso);
  double muonHcalIso = patMuon.userIsolation(pat::HcalIso);

  hMuonTrkIsoPt_->Fill(muonTrackIso, weight);
  hMuonTrkIsoPtVsMuonPt_->Fill(patMuon.pt(), muonTrackIso, weight);
  hMuonEcalIsoPt_->Fill(muonEcalIso, weight);
  hMuonEcalIsoPtVsMuonPt_->Fill(patMuon.pt(), muonEcalIso, weight);
  hMuonHcalIsoPt_->Fill(muonHcalIso, weight);
  //hMuonIsoSumPt_->Fill(muonTrackIso + muonEcalIso + muonHcalIso, weight);
  hMuonIsoSumPt_->Fill(muonTrackIso + muonEcalIso, weight);
  hMuonIsoSumPtVsMuonPt_->Fill(patMuon.pt(), muonTrackIso + muonEcalIso, weight);
  hMuonTrkIsoPtRel_->Fill(muonTrackIso/patMuon.pt(), weight);
  hMuonEcalIsoPtRel_->Fill(muonEcalIso/patMuon.pt(), weight);
  hMuonHcalIsoPtRel_->Fill(muonHcalIso/patMuon.pt(), weight);
  //hMuonIsoSumPtRel_->Fill((muonTrackIso + muonEcalIso + muonHcalIso)/patMuon.pt(), weight);
  hMuonIsoSumPtRel_->Fill((muonTrackIso + muonEcalIso)/patMuon.pt(), weight);

  //std::cout << " particleIso = " << patMuon.particleIso() << std::endl;
  //std::cout << " chargedHadronIso = " << patMuon.chargedHadronIso() << std::endl;
  //std::cout << " neutralHadronIso = " << patMuon.neutralHadronIso() << std::endl;
  //std::cout << " photonIso = " << patMuon.photonIso() << std::endl;
  
  hMuonParticleFlowIsoPt_->Fill(patMuon.particleIso(), weight);
  hMuonPFChargedHadronIsoPt_->Fill(patMuon.chargedHadronIso(), weight);
  hMuonPFNeutralHadronIsoPt_->Fill(patMuon.neutralHadronIso(), weight);
  hMuonPFGammaIsoPt_->Fill(patMuon.photonIso(), weight);
  
  if ( makeIsoPtCtrlHistograms_ ) {
    hMuonPFChargedHadronIsoPtCtrl_->Fill(muonTrackIso, patMuon.chargedHadronIso(), weight);
    hMuonPFGammaIsoPtCtrl_->Fill(muonEcalIso, patMuon.photonIso(), weight);
  }

  fillLeptonIsoDepositHistograms(patMuon.trackIsoDeposit(), 
				 hMuonTrkIsoValProfile_, hMuonTrkIsoEtaDistProfile_, hMuonTrkIsoPhiDistProfile_, weight);
  fillLeptonIsoDepositHistograms(patMuon.ecalIsoDeposit(), 
				 hMuonEcalIsoValProfile_, hMuonEcalIsoEtaDistProfile_, hMuonEcalIsoPhiDistProfile_, weight);
  fillLeptonIsoDepositHistograms(patMuon.hcalIsoDeposit(), 
				 hMuonHcalIsoValProfile_, hMuonHcalIsoEtaDistProfile_, hMuonHcalIsoPhiDistProfile_, weight);
}

void MuonHistManager::fillMuonIsoConeSizeDepHistograms(const pat::Muon& patMuon, double weight)
{
  //std::cout << "<MuonHistManager::fillMuonIsoConeSizeDepHistograms>:" << std::endl;

  for ( unsigned iConeSize = 1; iConeSize <= numMuonIsoConeSizes_; ++iConeSize ) {
    double isoConeSize_i = iConeSize*muonIsoConeSizeIncr_;
    
    double muonTrkIsoDeposit_i = patMuon.trackIsoDeposit()->countWithin(isoConeSize_i, muonTrkIsoParam_, false);
    hMuonTrkIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonTrkIsoDeposit_i, weight);
    
    double muonEcalIsoDeposit_i = patMuon.ecalIsoDeposit()->countWithin(isoConeSize_i, muonEcalIsoParam_, false);
    hMuonEcalIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonEcalIsoDeposit_i, weight);
    
    double muonHcalIsoDeposit_i = patMuon.hcalIsoDeposit()->countWithin(isoConeSize_i, muonHcalIsoParam_, false);
    hMuonHcalIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonHcalIsoDeposit_i, weight);

    if ( patMuon.isoDeposit(pat::PfAllParticleIso) ) {
      double muonParticleFlowIsoDeposit_i 
	= patMuon.isoDeposit(pat::PfAllParticleIso)->countWithin(isoConeSize_i, muonParticleFlowIsoParam_, false);
      hMuonParticleFlowIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonParticleFlowIsoDeposit_i, weight);
    }
    
    if ( patMuon.isoDeposit(pat::PfChargedHadronIso) ) {
      double muonPFChargedHadronIsoDeposit_i 
	= patMuon.isoDeposit(pat::PfChargedHadronIso)->countWithin(isoConeSize_i, muonParticleFlowIsoParam_, false);
      hMuonPFChargedHadronIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonPFChargedHadronIsoDeposit_i, weight);
    }
    
    if ( patMuon.isoDeposit(pat::PfNeutralHadronIso) ) {
      double muonPFNeutralHadronIsoDeposit_i 
	= patMuon.isoDeposit(pat::PfNeutralHadronIso)->countWithin(isoConeSize_i, muonParticleFlowIsoParam_, false);
      hMuonPFNeutralHadronIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonPFNeutralHadronIsoDeposit_i, weight);
    }

    if ( patMuon.isoDeposit(pat::PfGammaIso) ) {
      double muonPFGammaIsoDeposit_i 
	= patMuon.isoDeposit(pat::PfGammaIso)->countWithin(isoConeSize_i, muonParticleFlowIsoParam_, false);
      hMuonPFGammaIsoPtConeSizeDep_[iConeSize - 1]->Fill(muonPFGammaIsoDeposit_i, weight);
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, MuonHistManager, "MuonHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, MuonHistManager, "MuonHistManager");

#include "TauAnalysis/Core/interface/HistManagerAdapter.h"

typedef HistManagerAdapter<MuonHistManager> MuonAnalyzer;

DEFINE_ANOTHER_FWK_MODULE(MuonAnalyzer);
