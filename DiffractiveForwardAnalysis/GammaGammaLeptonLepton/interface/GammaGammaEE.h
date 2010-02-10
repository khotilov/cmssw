#ifndef DiffractiveForwardAnalysis_GammaGammaEE
#define DiffractiveForwardAnalysis_GammaGammaEE

// user include files
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/EDFilter.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/Event.h>
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Common/interface/TriggerResults.h" 
#include "FWCore/Framework/interface/TriggerNames.h" 

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

class GammaGammaEE : public edm::EDAnalyzer {
 public:
  explicit GammaGammaEE(const edm::ParameterSet&);
  ~GammaGammaEE();
  
  
 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  
  edm::InputTag recTrackLabel;
  edm::InputTag theGLBMuonLabel;
  edm::InputTag thePixelGsfELabel;
  edm::InputTag theJetLabel;
  edm::InputTag theMetLabel;
  edm::InputTag thePhotonLabel;
  edm::InputTag theCaloTowLabel;
  edm::InputTag recCastorTowerLabel;  

  double eldetmax;
  double eldphimin;
  double drisocalo;

  std::string rootfilename;

  TFile *thefile;
  TTree *thetree;

  int nEleCand;
  int ELEMAX;// used to set maximum of arrays
  double EleCand_px[4];
  double EleCand_py[4];
  double EleCand_pz[4];
  double EleCand_p[4];
  double EleCand_e[4];
  double EleCand_et[4];
  double EleCand_phi[4];
  double EleCand_eta[4];
  double EleCandTrack_p[100];
  int EleCand_charge[4];
  int EleCand_looseid[4];
  double EleCand_likelihoodid[4];
  int EleCand_robustid[4];

  int nPU; 

  int nHLTEleCand; 
  double HLTEleCand_pt[10]; 
  double HLTEleCand_eta[10]; 
  double HLTEleCand_phi[10]; 
  int HLTEleCand_charge[10]; 

  double ElEl_mass;
  double ElEl_dphi;
  double ElEl_vtxx; 
  double ElEl_vtxy; 
  double ElEl_vtxz; 
  double ElEl_vtxchi2dof; 
  int ElEl_vtxisvalid; 
  int ElEl_extratracks1mm; 
  int ElEl_extratracks2mm; 
  int ElEl_extratracks3mm; 
  int ElEl_extratracks4mm; 
  int ElEl_extratracks5mm; 
  int ElEl_extratracks1cm; 
  int ElEl_extratracks2cm; 
  int ElEl_extratracks5cm; 
  int ElEl_extratracks10cm; 

  int nJetCand;
  int JETMAX;// used to set maximum of arrays
  double JetCand_px[30];
  double JetCand_py[30];
  double JetCand_pz[30];
  double JetCand_e[30];
  double JetCand_eta[30];
  double JetCand_phi[30];
  double HighestJet_e;
  double HighestJet_eta;
  double HighestJet_phi;
  double SumJet_e;

  int nPFlowCand;
  int PFlowCandIds[500];

  int HitInZDC; 
  int HitInCastor; 

  double Etmiss;

  int nCaloCand;
  int nExtraCaloTowersE1, nExtraCaloTowersE2, nExtraCaloTowersE3, nExtraCaloTowersE4, nExtraCaloTowersE5, nExtraCaloTowersE6, nExtraCaloTowersE7, nExtraCaloTowersE8, nExtraCaloTowersE9;  
  int nExtraCaloTowersEt0pt1, nExtraCaloTowersEt0pt2, nExtraCaloTowersEt0pt5, nExtraCaloTowersEt1, nExtraCaloTowersEt2, nExtraCaloTowersEt3, nExtraCaloTowersEt4;  
  int nExtraCaloTowersE0hf, nExtraCaloTowersE1hf, nExtraCaloTowersE2hf;
  int nExtraCaloTowersE1he, nExtraCaloTowersE2he, nExtraCaloTowersE3he;
  int nExtraCaloTowersE2hb, nExtraCaloTowersE3hb, nExtraCaloTowersE4hb;
  double CaloTower_e[1000];
  double CaloTower_et[1000];
  double CaloTower_eta[1000];
  double CaloTower_phi[1000];
  double CaloTower_dr[1000];
  double HighestCaloTower_e;
  double HighestCaloTower_eta;
  double HighestCaloTower_phi;
  double HighestCaloTower_dr;
  double HighestEtCaloTower_et;
  double HighestEtCaloTower_eta;
  double HighestEtCaloTower_phi;
  double HighestEtCaloTower_dr;
  double SumCalo_e;

  int nCastorTowerCand;  
  double CastorTower_e[1000];  
  double CastorTower_eta[1000];   
  double CastorTower_phi[1000];  
  double CastorTower_emratio[1000];  
  double HighestCastorTowerFwd_e;  
  double HighestCastorTowerBwd_e;  
  double SumCastorFwd_e;
  double SumCastorBwd_e;

  int nTrackCand;
  int TRACKMAX;
  int nExtraTrackCand; 
  double TrackCand_px[100];
  double TrackCand_py[100];
  double TrackCand_pz[100];
  double TrackCand_p[100];
  double TrackCand_eta[100];
  double TrackCand_pt[100];
  double TrackCand_phi[100];
  int TrackCand_charge[100];
  double TrackCand_vtxdxyz[100];
  double ClosestExtraTrack_vtxdxyz; 

  double evweight;

  int HLT2Electron5_L1R_NI;
  int HLT2ElectronExclusive;

  edm::TriggerNames trigNames ;

};
#endif
