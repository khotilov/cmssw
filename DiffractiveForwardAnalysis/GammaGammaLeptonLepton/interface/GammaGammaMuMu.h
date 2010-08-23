#ifndef DiffractiveForwardAnalysis_GammaGammaMuMu
#define DiffractiveForwardAnalysis_GammaGammaMuMu

// user include files
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/EDFilter.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h> 
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>
#include <FWCore/ParameterSet/interface/ParameterDescriptionNode.h> 
#include <FWCore/Framework/interface/Event.h>
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Common/interface/TriggerResults.h" 
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h" 
#include "FWCore/Common/interface/TriggerNames.h"

#include "MuonAnalysis/TagAndProbe/interface/MuonPerformanceReadback.h" 
#include "MuonAnalysis/TagAndProbe/interface/MuonPerformance.h"  

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AcceptanceTableHelper.h"

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TLorentzVector.h>

class GammaGammaMuMu : public edm::EDAnalyzer {
 public:
  explicit GammaGammaMuMu(const edm::ParameterSet&);
  ~GammaGammaMuMu();
  
  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
  
 private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  MuonPerformanceReadback *effreader; 
  std::vector<std::string> algonames; 
  
  edm::InputTag recTrackLabel;
  edm::InputTag recVertexLabel;
  edm::InputTag theGLBMuonLabel;
  edm::InputTag thePixelGsfELabel;
  edm::InputTag theJetLabel;
  edm::InputTag theMetLabel;
  edm::InputTag thePhotonLabel;
  edm::InputTag theCaloTowLabel;
  edm::InputTag recCastorTowerLabel;   
  edm::InputTag recZDCRecHitsLabel;
  edm::InputTag recCastorRecHitsLabel;
  std::string hltMenuLabel;

  double mudptmax;
  double mudphimin;
  double drisocalo; 
  bool keepsamesign;

  std::string rootfilename;

  TFile *thefile;
  TTree *thetree;

  int BX;
  int Run;
  int LumiSection;
  int EventNum;

  int L1TechnicalTriggers[128]; 

  int nEvt;
  int nMuonCand;
  int MUONMAX;// used to set maximum of arrays
  double MuonCand_px[10];
  double MuonCand_py[10];
  double MuonCand_pz[10];
  double MuonCand_vtxx[10]; 
  double MuonCand_vtxy[10]; 
  double MuonCand_vtxz[10]; 
  double MuonCand_p[10];
  double MuonCand_eta[10];
  double MuonCand_pt[10];
  double MuonCand_phi[10];
  double MuonCand_e[10];
  double MuonCandTrack_p[100];
  double MuonCand_efficiency[10];
  int MuonCand_charge[10];
  int MuonCand_tmlsloosemuonid[10];
  int MuonCand_tmlsOptLowPtloosemuonid[10];
  int MuonCand_tm2dloosemuid[10];
  int MuonCand_arbmuid[10];
  int MuonCand_tmlsAngloosemuonid[10];
  int MuonCand_tmlsAngtightmuonid[10]; 
  int MuonCand_tmosAngloosemuonid[10]; 
  int MuonCand_tmosAngtightmuonid[10]; 
  int MuonCand_isglobal[10];
  int MuonCand_istracker[10];
  int MuonCand_isstandalone[10];
  double MuonCand_ecalisor3[10]; 
  double MuonCand_hcalisor3[10]; 
  double MuonCand_trkisor3[10];
  double MuonCand_hoisor3[10]; 
  double MuonCand_ecalisor5[10];  
  double MuonCand_hcalisor5[10];  
  double MuonCand_trkisor5[10]; 
  double MuonCand_hoisor5[10];
  double MuonCand_timein[10]; 	 
  double MuonCand_timeout[10]; 	 
  double MuonCand_timeouterr[10]; 	 
  double MuonCand_timeinerr[10]; 	 
  int MuonCand_validtrackhits[10];
  int MuonCand_validhits[10];
  double MuonCand_normchi2[10];
  double MuonCand_normtrackchi2[10];

  int MuonPairCand[2];

  int nHLTMu3MuonCand;
  double HLT_Mu3_MuonCand_pt[10];
  double HLT_Mu3_MuonCand_eta[10];
  double HLT_Mu3_MuonCand_phi[10];
  int HLT_Mu3_MuonCand_charge[10];

  int nHLTDiMu3MuonCand; 
  double HLT_DoubleMu3_MuonCand_pt[10]; 
  double HLT_DoubleMu3_MuonCand_eta[10]; 
  double HLT_DoubleMu3_MuonCand_phi[10]; 
  int HLT_DoubleMu3_MuonCand_charge[10]; 

  int nHLTDiMu0MuonCand;  
  double HLT_DoubleMu0_MuonCand_pt[10];  
  double HLT_DoubleMu0_MuonCand_eta[10];  
  double HLT_DoubleMu0_MuonCand_phi[10];  
  int HLT_DoubleMu0_MuonCand_charge[10];  

  double MuMu_mass;
  double MuMu_dphi;
  double MuMu_dpt;
  double MuMu_vtxx;
  double MuMu_vtxy;
  double MuMu_vtxz;
  double MuMu_vtxT;
  double MuMu_vtxchi2dof;
  int MuMu_vtxisvalid;
  int MuMu_extratracks1mm;
  int MuMu_extratracks3mm;
  int MuMu_extratracks5mm;
  int MuMu_extratracks1cm;
  int MuMu_extratracks3cm;
  int MuMu_extratracks5cm;
  int MuMu_extratracks10cm;
  double MuMuGamma_mass[50];

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

  int HitInZDC; 
  int HitInCastor; 

  int nGenPhotCand;
  int GENPHOTONMAX;
  double GenPhotCand_pt[5];
  double GenPhotCand_eta[5];
  double GenPhotCand_phi[5];

  int nGenMuonCand;
  int GENMUONMAX;
  double GenMuonCand_px[10]; 
  double GenMuonCand_py[10]; 
  double GenMuonCand_pz[10]; 
  double GenMuMu_eta; 
  double GenMuMu_pt; 


  double Etmiss;

  int nCaloCand;
  int nExtraCaloTowersE0pt6eb, nExtraCaloTowersE2pt45ee, nExtraCaloTowersE1pt25hb, nExtraCaloTowersE1pt9he, nExtraCaloTowersE4pt5hfp, nExtraCaloTowersE4pt0hfm;
  int nExtraCaloTowersE1, nExtraCaloTowersE2, nExtraCaloTowersE3, nExtraCaloTowersE4, nExtraCaloTowersE5, nExtraCaloTowersE6, nExtraCaloTowersE7, nExtraCaloTowersE8, nExtraCaloTowersE9; 
  int nExtraCaloTowersEt0pt1, nExtraCaloTowersEt0pt2, nExtraCaloTowersEt0pt5, nExtraCaloTowersEt1, nExtraCaloTowersEt2, nExtraCaloTowersEt3, nExtraCaloTowersEt4; 

  double CaloTower_e[1000];
  double CaloTower_et[1000];
  double CaloTower_eta[1000];
  double CaloTower_phi[1000];
  double CaloTower_dr[1000];
  double CaloTower_emE[1000];
  double CaloTower_hadE[1000];
  double CaloTower_outE[1000];
  int CaloTower_ID[1000];
  double CaloTower_x[1000];
  double CaloTower_y[1000];
  double CaloTower_z[1000];
  double CaloTower_t[1000];
  int CaloTower_badhcalcells[1000];
  int CaloTower_problemhcalcells[1000];
  int CaloTower_badecalcells[1000];
  int CaloTower_problemecalcells[1000];

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
  int nCastorTowerCandE3;
  double CastorTower_e[1000];   
  double CastorTower_eta[1000];    
  double CastorTower_phi[1000];   
  double CastorTower_emratio[1000];   
  double HighestCastorTowerFwd_e;   
  double HighestCastorTowerBwd_e;   
  double SumCastorFwd_e; 
  double SumCastorBwd_e; 

  int nZDChitCand;
  int ZDChit_section[5000];
  double ZDChit_energy[5000];
  double ZDChit_time[5000];
  int ZDChit_side[5000];
  double ZDCsumEMplus;
  double ZDCsumHADplus;
  double ZDCsumEMminus;
  double ZDCsumHADminus;

  double CASTORsumRecHitsE;

  int nPrimVertexCand;
  double PrimVertexCand_x[10];
  double PrimVertexCand_y[10];
  double PrimVertexCand_z[10];
  int PrimVertexCand_tracks[10];
  double PrimVertexCand_chi2[10];
  double PrimVertexCand_ndof[10];

  int nTrackCand;
  int nQualityTrackCand;
  int TRACKMAX;
  double TrackCand_purity[500];
  int TrackCand_nhits[500];
  double TrackCand_px[500];
  double TrackCand_py[500];
  double TrackCand_pz[500];
  double TrackCand_p[500];
  double TrackCand_eta[500];
  double TrackCand_pt[500];
  double TrackCand_phi[500];
  double TrackCand_vtxdxyz[500];
  double TrackCand_vtxT[500];
  double TrackCand_vtxZ[500];
  double TrackCand_X[500];
  double TrackCand_Y[500];
  double TrackCand_Z[500];
  int TrackCand_charge[500];
  double ClosestExtraTrack_vtxdxyz;

  int nPFPhotonCand;  
  int PHOTONMAX;
  double PFPhotonCand_pt[50];
  double PFPhotonCand_eta[50];
  double PFPhotonCand_phi[50];
  double PFPhotonCand_drtrue[50];

  double evweight;
  
  int HLT_DoubleMu3;
  int HLT_Mu3;
  int HLT_DoubleMu0;
  int HLT_L2Mu0;
  int HLT_L1DoubleMuOpen;
  int HLT_DoubleMu3_Prescl; 
  int HLT_Mu3_Prescl; 
  int HLT_DoubleMu0_Prescl; 
  int HLT_L2Mu0_Prescl; 
  int HLT_L1DoubleMuOpen_Prescl; 

  double LowPt_pt[10];
  double LowPt_eta[10];

  int nPU;
  HLTConfigProvider hltConfig_;  

  AcceptanceTableHelper helper420beam1;   
  AcceptanceTableHelper helper420beam2;   
  AcceptanceTableHelper helper220beam1;   
  AcceptanceTableHelper helper220beam2;   
  AcceptanceTableHelper helper420a220beam1;   
  AcceptanceTableHelper helper420a220beam2;   

  edm::TriggerNames trigNames ;
};
#endif
