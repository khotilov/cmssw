#ifndef VgKit_test_VgAnalyzerKit_h
#define VgKit_test_VgAnalyzerKit_h

#include <memory>
#include <fstream>

#include "ElectroWeakAnalysis/MultiBosons/interface/PatKitHelper.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"

#include "TTree.h"
#include "TH1F.h"

using namespace edm;

const int maxP = 500;

class VgAnalyzerKit : public EDProducer {

public:

  explicit VgAnalyzerKit(const ParameterSet&);
  virtual ~VgAnalyzerKit();

protected:

  virtual void beginJob();
  virtual void produce( Event &, const EventSetup & );
  virtual void endJob();
  double eT(double pt1, double pt2) const;
  double massT(double pt1, double pt2, double wpx, double wpy) const;
  double acop(double phi1, double phi2) const;

  bool saveHistograms_;
  bool saveHLTInfo_;
  bool doGenParticles_;

  InputTag vtxlabel_;
  InputTag gtdigilabel_;
  InputTag tracklabel_;
  InputTag trgResults_;
  InputTag caloTowerlabel_;
  InputTag tcMETlabel_;
  InputTag pfMETlabel_;
  InputTag ebReducedRecHitCollection_;
  InputTag eeReducedRecHitCollection_;
  InputTag beamSpotCollection_;

  int verbosity_;
  int leadingElePtCut_;
  int leadingMuPtCut_;
  int leadingPhoPtCut_;

  pat::PatKitHelper helper_;
  TTree *tree_;
  TH1F  *hEvents_; 

  int    run_;
  int    event_;
  int    orbit_;
  int    bx_;
  int    lumis_;
  bool   isData_;
  int    ttbit0_;
  int    ttbit34_;
  int    ttbit40_;
  int    ttbit41_;
  int    ttbitBH_;
  float  pdf_[7];
  float  processID_;
  int    nHLT_;
  int    HLT_[maxP];
  int    nHFTowersP_;
  int    nHFTowersN_;
  int    nVtx_;
  float  vtx_[50][3];
  int    vtxNTrk_[50];
  int    vtxNDF_[50];
  float  vtxD0_[50];
  int    IsVtxGood_;
  int    nTrk_;
  int    nGoodTrk_;
  int    IsTracksGood_;
  // genParticle
  int    nMC_;
  int    mcPID[maxP];
  float  mcPt[maxP];
  float  mcMass[maxP];
  float  mcEta[maxP];
  float  mcPhi[maxP];
  int    mcGMomPID[maxP];
  int    mcMomPID[maxP];
  float  mcMomPt[maxP];
  float  mcMomMass[maxP];
  float  mcMomEta[maxP];
  float  mcMomPhi[maxP];
  int    mcIndex[maxP];
  int    mcDecayType[maxP];
  // Gen & Reco MET
  float  genMET_;
  float  genMETx_;
  float  genMETy_;
  float  genMETPhi_;
  float  MET_;
  float  METx_;
  float  METy_;
  float  METPhi_;
  float  METsumEt_;
  float  uncorrMET_[3];
  float  uncorrMETPhi_[3];
  float  uncorrMETSumEt_[3];
  // tcMET
  float  tcMET_;
  float  tcMETx_;
  float  tcMETy_;
  float  tcMETPhi_;
  float  tcMETsumEt_;
  float  tcMETmEtSig_;
  float  tcMETSig_;
  // pfMET
  float  pfMET_;
  float  pfMETx_;
  float  pfMETy_;
  float  pfMETPhi_;
  float  pfMETsumEt_;
  float  pfMETmEtSig_;
  float  pfMETSig_;
  // Electron
  int    nEle_;
  int    eleID_[maxP][12];
  int    eleClass_[maxP];
  int    eleCharge_[maxP];
  float  eleEn_[maxP];
  float  eleSCEn_[maxP];
  float  eleESEn_[maxP];
  float  eleSCEtaWidth_[maxP];
  float  eleSCPhiWidth_[maxP];
  float  eleVtx_[maxP][3];
  float  eleCaloPos_[maxP][3];
  float  eleSCPos_[maxP][3];
  float  elePt_[maxP];
  float  elePz_[maxP];
  float  eleEta_[maxP];
  float  elePhi_[maxP];
  float  eleSCEta_[maxP];
  float  eleSCPhi_[maxP];
  float  eleSCRawEn_[maxP];
  float  eleHoverE_[maxP];
  float  eleHoverE1_[maxP];
  float  eleHoverE2_[maxP];
  float  eleEoverP_[maxP];
  float  elePin_[maxP];
  float  elePout_[maxP];
  float  eleBrem_[maxP];
  float  eledEtaAtVtx_[maxP];
  float  eledPhiAtVtx_[maxP];
  float  eleSigmaEtaEta_[maxP];
  float  eleSigmaIEtaIEta_[maxP];
  float  eleEMax_[maxP];
  float  eleE2nd_[maxP];
  float  eleE2x2_[maxP];
  float  eleE3x2_[maxP];
  float  eleE3x3_[maxP];
  float  eleE4x4_[maxP];
  float  eleE1x4_[maxP];
  float  eleE1x5_[maxP];
  float  eleE5x5_[maxP];
  float  eleE2x5Right_[maxP];
  float  eleE2x5Left_[maxP];
  float  eleE2x5Top_[maxP];
  float  eleE2x5Bottom_[maxP];
  float  eleE2x5Max_[maxP];
  float  eleERight_[maxP];
  float  eleELeft_[maxP];
  float  eleETop_[maxP];
  float  eleEBottom_[maxP];
  int    eleGenIndex_[maxP];
  int    eleGenGMomPID_[maxP];
  int    eleGenMomPID_[maxP];
  float  eleGenMomPt_[maxP];
  float  eleIsoTrkDR03_[maxP];
  float  eleIsoEcalDR03_[maxP];
  float  eleIsoHcalDR03_[maxP];
  float  eleIsoTrkDR04_[maxP];
  float  eleIsoEcalDR04_[maxP];
  float  eleIsoHcalDR04_[maxP];
  float  eleChi2NDF_[maxP];
  float  eleD0_[maxP];
  int    eleNumberOfValidHits_[maxP];
  int    eleValidHitInFirstPXB_[maxP];
  int    eleTrkExpectHitsInner_[maxP];
  float  eleDist_[maxP];
  float  eleDcot_[maxP];
  float  eleConvRadius_[maxP];
  float  eleConvPoint_[maxP][3];

  // Photon
  int    nPho_;
  bool   phoIsPhoton_[maxP];
  float  phoE_[maxP];
  float  phoEt_[maxP];
  float  phoPz_[maxP];
  float  phoEta_[maxP];
  float  phoPhi_[maxP];
  float  phoR9_[maxP];
  float  phoTrkIsoSolidDR03_[maxP];
  float  phoTrkIsoHollowDR03_[maxP];
  int    phoNTrkSolidDR03_[maxP];
  int    phoNTrkHollowDR03_[maxP];
  float  phoEcalIsoDR03_[maxP];
  float  phoHcalIsoDR03_[maxP];
  float  phoTrkIsoSolidDR04_[maxP];
  float  phoTrkIsoHollowDR04_[maxP];
  int    phoNTrkSolidDR04_[maxP];
  int    phoNTrkHollowDR04_[maxP];
  float  phoEcalIsoDR04_[maxP];
  float  phoHcalIsoDR04_[maxP];
  float  phoHoverE_[maxP];
  float  phoHoverE1_[maxP];
  float  phoHoverE2_[maxP];
  float  phoSeedTime_[maxP];
  float  phoSigmaEtaEta_[maxP];
  float  phoSigmaIEtaIEta_[maxP];
  int    phoPos_[maxP];
  float  phoEMax_[maxP];
  float  phoE2nd_[maxP];
  float  phoE2x2_[maxP];
  float  phoE3x2_[maxP];
  float  phoE3x3_[maxP];
  float  phoE4x4_[maxP];
  float  phoE1x5_[maxP];
  float  phoE5x5_[maxP];
  float  phoE2x5Right_[maxP];
  float  phoE2x5Left_[maxP];
  float  phoE2x5Top_[maxP];
  float  phoE2x5Bottom_[maxP];
  float  phoE2x5Max_[maxP];
  float  phoERight_[maxP];
  float  phoELeft_[maxP];
  float  phoETop_[maxP];
  float  phoEBottom_[maxP];
  float  phoRoundness_[maxP];
  float  phoAngle_[maxP];
  int    phoGenIndex_[maxP];
  int    phoGenGMomPID[maxP];
  int    phoGenMomPID[maxP];
  float  phoGenMomPt[maxP];
  float  phoSCE_[maxP];
  float  phoSCEt_[maxP];
  float  phoSCEta_[maxP];
  float  phoSCPhi_[maxP];
  float  phoSCEtaWidth_[maxP];
  float  phoSCPhiWidth_[maxP];
  int    phoOverlap_[maxP];
  int    phohasPixelSeed_[maxP];
  // Muon
  int    nMu_;
  float  muEta_[maxP];
  float  muPhi_[maxP];
  int    muCharge_[maxP];
  float  muPt_[maxP];
  float  muPz_[maxP];
  int    muGenIndex_[maxP];
  float  muIsoTrk_[maxP];
  float  muIsoCalo_[maxP];
  float  muIsoEcal_[maxP];
  float  muIsoHcal_[maxP];
  float  muChi2NDF_[maxP];
  float  muEmVeto_[maxP];
  float  muHadVeto_[maxP];
  int    muType_[maxP];
  bool   muID_[maxP][6];
  float  muD0_[maxP];
  int    muNumberOfValidTrkHits_[maxP];
  int    muNumberOfValidPixelHits_[maxP];
  int    muNumberOfValidMuonHits_[maxP];
  // Zee candidate
  int    nZee_;
  float  ZeeMass_[maxP];
  float  ZeePt_[maxP];
  float  ZeeEta_[maxP];
  float  ZeePhi_[maxP];
  int    ZeeLeg1Index_[maxP];
  int    ZeeLeg2Index_[maxP];
  // Zmumu candidate
  int    nZmumu_;
  float  ZmumuMass_[maxP];
  float  ZmumuPt_[maxP];
  float  ZmumuEta_[maxP];
  float  ZmumuPhi_[maxP];
  int    ZmumuLeg1Index_[maxP];
  int    ZmumuLeg2Index_[maxP];

  // Wenu candidate
  int    nWenu_;
  float  WenuMassTCaloMET_[maxP];
  float  WenuEtCaloMET_[maxP];
  float  WenuACopCaloMET_[maxP];
  float  WenuMassTTcMET_[maxP];
  float  WenuEtTcMET_[maxP];
  float  WenuACopTcMET_[maxP];
  float  WenuMassTPfMET_[maxP];
  float  WenuEtPfMET_[maxP];
  float  WenuACopPfMET_[maxP];
  int    WenuEleIndex_[maxP];

  // Wmunu candidate
  int    nWmunu_;
  float  WmunuMassTCaloMET_[maxP];
  float  WmunuEtCaloMET_[maxP];
  float  WmunuACopCaloMET_[maxP];
  float  WmunuMassTTcMET_[maxP];
  float  WmunuEtTcMET_[maxP];
  float  WmunuACopTcMET_[maxP];
  float  WmunuMassTPfMET_[maxP];
  float  WmunuEtPfMET_[maxP];
  float  WmunuACopPfMET_[maxP];
  int    WmunuMuIndex_[maxP];

  // Physics objects handles
  Handle<View<pat::Muon> >                  muonHandle_;
  Handle<View<pat::Electron> >              electronHandle_;
  Handle<View<pat::Tau> >                   tauHandle_;
  Handle<View<pat::Jet> >                   jetHandle_;
  Handle<View<pat::MET> >                   METHandle_;
  Handle<View<pat::Photon> >                photonHandle_;
  Handle<std::vector<reco::RecoChargedCandidate> > trackHandle_;
  Handle<std::vector<reco::GenParticle> >          genParticlesHandle_;
  Handle<reco::CandidateView>               zmumuHandle_;

};

#endif
