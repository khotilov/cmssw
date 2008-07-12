// -*- C++ -*-
//
// Package:    HcalProm
// Class:      HcalProm
// 
/**\class HcalProm HcalProm.cc RecoHcal/HcalProm/src/HcalProm.cc

 SIMPLE EXAMPLE CODE demonstrating the access to officially reconstructed global run data.
 The config file is under the test directory.
 The config file should read the officially reconstructed root file in castor or 
 The output from Configuration/GlobalRuns/data/recoT0.cfg (official reco cfg). 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Efe Yazgan
// Updated        :  Taylan Yetkin (2008/05/08)
//         Created:  Wed Apr 16 10:03:18 CEST 2008
// $Id: HcalProm.h,v 1.22 2008/07/12 08:35:59 tyetkin Exp $
//
//


// system include files
#include <memory>
#include <fstream>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
//TFile Service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

//geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
//point
#include "DataFormats/GeometryVector/interface/GlobalPoint.h" 
//vector
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

//HBHERecHits
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

//trigger
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtPsbWord.h"

//propagation
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FastSimulation/CalorimeterProperties/interface/Calorimeter.h"
#include "FastSimulation/CaloGeometryTools/interface/CaloGeometryHelper.h"

//Ntuple creation
#include "TTree.h"

class TH1F;
class TH2F;
class TStyle;
class TTree;

#include "TH1I.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
//
// class decleration
//
class HcalProm : public edm::EDAnalyzer {
   public:
      explicit HcalProm(const edm::ParameterSet&);
      ~HcalProm();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      //adding HTML output functions (very similar to DQM) 
      void htmlOutput(void);
      void CaloTowerHTMLOutput(std::string startTime, std::string htmlDir, std::string htmlName);
      void DigiHTMLOutput(std::string startTime, std::string htmlDir, std::string htmlName);
      void RecHitHTMLOutput(std::string startTime, std::string htmlDir, std::string htmlName);
      void JetMetHTMLOutput(std::string startTime, std::string htmlDir, std::string htmlName);
      void MuonHTMLOutput(std::string startTime, std::string htmlDir, std::string htmlName);
      void HPDNoiseHTMLOutput(std::string startTime, std::string htmlDir, std::string htmlName);
      void histoHTML(TH1F * hist, const char *xlab, const char *ylab, int width, std::ofstream & htmlFile, std::string htmlDir);
      void histoHTML2(TH2F * hist, const char *xlab, const char *ylab, int width, std::ofstream & htmlFile, std::string htmlDir, bool color = false);
      std::string getIMG(TH1F* hist, int size, std::string htmlDir, const char* xlab, const char* ylab);
      std::string getIMG2(TH2F* hist, int size, std::string htmlDir, const char* xlab, const char* ylab,bool color);
      void parseString(std::string& title);
      void cleanString(std::string& title);
      float deltaR(float eta0, float phi0, float eta, float phi);
      float getNxNClusterEnergy(int gridSize,int ieta, int iphi, const HBHERecHitCollection*  HBHERecHits);
      float get1xNClusterEnergy(int gridSize,int ieta, int iphi, const HBHERecHitCollection*  HBHERecHits);
      float getEtaCluster(int gridSize,int ieta, int iphi, const HBHERecHitCollection*  HBHERecHits);
      float getPhiCluster(int gridSize,int ieta, int iphi, const HBHERecHitCollection*  HBHERecHits);
      float getEnergyEtaPhiDepth(int ieta, int iphi, int depth, const HFRecHitCollection * HFRecHits);
      //muon extrapolation
      void getHBmax(
          // inputs
          int tphi, const HBHERecHitCollection*  HBHERecHits,
          // outputs
          int *ietaMax, int *iphiMax, double *Emax, double *E2max
          );

      bool Extrapolate(
          // inputs
	  double ox, double oy, double oz, double px, double py, double pz, 
	  double radius, double zwidth, double ImpPar,
	  // outputs
	  double *xp_out, double *yp_out, double *zp_out,
	  double *xm_out, double *ym_out, double *zm_out,
	  double *thetap_out, double *phip_out,
	  double *thetam_out, double *phim_out
	  );
      reco::Track bestTrack(const reco::TrackCollection&) const;
      HcalDetId getClosestCell(float dR,float eta, float phi, const HBHERecHitCollection*  HBHERecHits, const CaloGeometry* geo, float &time, bool &found);
  //------- histogram bookers ----------
  void bookHistograms ();
  TH1F* book1DHistogram (TFileDirectory& fDir, const std::string& fName, const std::string& fTitle, 
			 int fNbins, double fXmin, double fXmax) const;
  TH2F* book2DHistogram (TFileDirectory& fDir, const std::string& fName, const std::string& fTitle, 
			 int fNbinsX, double fXmin, double fXmax,
			 int fNbinsY, double fYmin, double fYmax) const;
  TProfile* bookProfileHistogram (TFileDirectory & fDir, const std::string & fName, const std::string & fTitle, int fNbins, double fXmin, double fXmax, double fYmin, double fYmax) const;

      // ----------member data ---------------------------
      edm::Service<TFileService> fs;

      TH1F* h_global_trigger_bit;

      TH1F* h_hbhe_rechit_energy;
      TH1F* h_maxhbherec;
      TH1F* h_maxhbMinusrec;
      TH1F* h_maxhbPlusrec;
      TH2F* h_hbhe_eta_phi;	
      TH1F* h_hf_rechit_energy;
      TH2F* h_hf_eta_phi;
      TH1F* h_ho_rechit_energy;
      TH2F* h_ho_eta_phi;

      TH1F* h_eb_rechit_energy;
      TH1F* h_maxebeerec;
      TH1F* h_maxhorec;
//->kropiv
      TH1F* h_Emuon1hbMinus;
      TH1F* h_Emuon1hbPlus;
      TH1F* h_ETAmuon1hbMinus;
      TH1F* h_ETAmuon1hbPlus;
      TH1F* h_PHImuon1hbMinus;
      TH1F* h_PHImuon1hbPlus;
      TH1F* h_NumbermuonDT;
      TH1F* h_NumberHBTowersmuon; 
      TH1F* h_RphiinPosMuon;
      TH1F* h_RphioutPosMuon;
      TH1F* h_DeltaRphiPosMuon;
      TH2F* h_EtaPhiCorrelMuonPhiPlane;
      TH1F* h_Trigger;
// ntuple creation:
//     TFile*      hOutputFile ;
     TTree * myTree;
     int run, event, TriggerBit[4], NumMuonHBphiPlane, NumHBTowersMuon[50];
     int IdTowerPhiMuonIn[50], IdTowerPhiMuonOut[50], IdTowerEtaMuonIn[50], IdTowerEtaMuonOut[50];
     float TimeAvMuonHB[50];
     float PHIoutMuonHB[50], PHIinMuonHB[50],  ETAoutMuonHB[50], ETAinMuonHB[50];
     float ImpXYmuonHB[50], ZImpXYmuonHB[50];
//     float ImpXYmuon[50], ZImpXYmuon[50],
     float PHIinTowerHB[50];
     float EmuonHB[50];
//
     float XinPosMuonDT[50], YinPosMuonDT[50], ZinPosMuonDT[50];
     float XoutPosMuonDT[50], YoutPosMuonDT[50], ZoutPosMuonDT[50];
     float LengthMuonHB[50], LengthMuonDT[50]; 
     int NumHitsMuonDT[50], NumHitsMuonDTall[50];
//
     //int  NumHBTowersMuon2[50];
     int IdTowerPhiMuonIn2[50], IdTowerPhiMuonOut2[50], IdTowerEtaMuonIn2[50], IdTowerEtaMuonOut2[50];
     float TimeAvMuonHB2[50],  TimeAvMuonHBcut2[50], TimeAvMuonHBcutwide2[50], LengthMuonHB2[50];
     float EmuonHB2[50], EmuonHBcut2[50], EmuonHBcutwide2[50];
// make initialisation of Muons that not belong to the same Phi Plane in HB
     int NumMuonHBnoPhiPlane; 
     int IdTowerMuonNoPlanePhiIn[50], IdTowerMuonNoPlanePhiOut[50];
     int IdTowerMuonNoPlaneEtaIn[50], IdTowerMuonNoPlaneEtaOut[50];
     float XinPosMuonNoPlaneDT[50], YinPosMuonNoPlaneDT[50], ZinPosMuonNoPlaneDT[50];
     float XoutPosMuonNoPlaneDT[50], YoutPosMuonNoPlaneDT[50], ZoutPosMuonNoPlaneDT[50];
     float ImpXYmuonNoPlaneDT[50], ZImpXYmuonNoPlaneDT[50], LengthMuonNoPlaneDT[50];
     int NumHitsMuonNoPlaneDTall[50];
// end ntuple creation
//<-kropiv

//----

//->kropiv-efe
// 2nd (Partner) Muon is good: in one Phi Plane in HB
   //matching between 2 muons
  TH1F *hAngleMuonHB2DT;
  TH1F *hDeltaImpXYHB2DT;
  TH1F *hDeltaZImpXYHB2DT;
  TH1F *hPhiDeltaTowerHB2DT;
//Energy
  TH1F *hEmuonHB2DTTopPlus;
  TH1F *hEmuonHB2DTTopMinus;
  TH1F *hEmuonHB2DTBotPlus;
  TH1F *hEmuonHB2DTBotMinus;
   //Time
  TH1F *hTimeMuonHB2DTTopPlus;
  TH1F *hTimeMuonHB2DTTopMinus;
  TH1F *hTimeMuonHB2DTBotPlus;
  TH1F *hTimeMuonHB2DTBotMinus;
   //Number Towers in Eta Plane passing by good muon
  TH1F *hNumTowerMuonHB2DTTopPlus;
  TH1F *hNumTowerMuonHB2DTTopMinus;
  TH1F *hNumTowerMuonHB2DTBotPlus;
  TH1F *hNumTowerMuonHB2DTBotMinus;
  //idPhiTower for eta+ and eta-
  TH1F *hIdPhiTowerHB2DTMinus;
  TH1F *hIdPhiTowerHB2DTPlus;
   // Muon Track quality in DT
  TH1F *hNumHitsHB2DT;
  TH1F *hNumHitsHB2DT2;
  TH1F *hLengthMuonDTHB2DT;
  TH1F *hLengthMuonDTHB2DT2;
   // energy for time- and time+
  TH1F *hEmuonHB2DTTopPlusTimePlus;
  TH1F *hEmuonHB2DTTopMinusTimePlus;
  TH1F *hEmuonHB2DTBotPlusTimePlus;
  TH1F *hEmuonHB2DTBotMinusTimePlus;
  TH1F *hEmuonHB2DTTopPlusTimeMinus;
  TH1F *hEmuonHB2DTTopMinusTimeMinus;
  TH1F *hEmuonHB2DTBotPlusTimeMinus;
  TH1F *hEmuonHB2DTBotMinusTimeMinus;
  // in different time:
  TH1F *hEmuonHB2DTBotPlusTimePlus1;
  TH1F *hEmuonHB2DTBotPlusTimePlus2;
  TH1F *hEmuonHB2DTBotPlusTimePlus3;
  TH1F *hEmuonHB2DTBotPlusTimePlus4;
  TH1F *hEmuonHB2DTBotPlusTimeMinus1;
  TH1F *hEmuonHB2DTBotPlusTimeMinus2;
  TH1F *hEmuonHB2DTBotPlusTimeMinus3;
  TH1F *hEmuonHB2DTBotPlusTimeMinus4;
   // IdPhi vs Emuon
  TH2F *hIdPhiPlusVsE;
  TH2F *hIdPhiMinusVsE;
   // IdEta vs Emuon
  TH2F *hIdEtaTopVsE;
  TH2F *hIdEtaBotVsE;
 
   // control plots
  TH1F *hImpXYHB2DT;
  TH1F *hZImpXYHB2DT;
  TH2F *hLmuonDTImpXY;
  TH2F *hImpXYvsZ;

   //energy of selected muon in different part of phi Towers
  TH1F *hEmuonPhiDetaTower1;
  TH1F *hEmuonPhiDetaTower2;
  TH1F *hEmuonPhiDetaTower3;
  TH1F *hEmuonPhiDetaTower4;
  TH1F *hEmuonPhiDetaTower5;
  TH1F *hEmuonPhiDetaTower6;

   // fill TProfile for Time as iphi
  TProfile *hProfTimeAsIdPhiMinus;
  TProfile *hProfTimeAsIdPhiPlus;
  TProfile *hProfTimeAsIdEtaTop;
  TProfile *hProfTimeAsIdEtaBot;

//->kropiv-efe
      TH1F* h_maxebee_plus_maxhbhe;

      TH1F* h_hbplustiming;
      TH1F* h_hbminustiming;
      TH1F* h_heplustiming;
      TH1F* h_heminustiming;
      TH1F* h_hfplustiming;
      TH1F* h_hfminustiming;
      TH1F* h_ecal_cluster_energy;
      TH1F* h_ecal_cluster_eta;
      TH1F* h_ecal_cluster_phi;
      TH2F* h_ecal_vs_hcal_X;
      TH2F* h_ecal_vs_hcal_Y; 
      
      TH1F* h_HFplus_energy_CSC_trig;
      TH1F* h_HFminus_energy_CSC_trig;
      TH2F* h_HFplus_energy_vs_RecHitTime_CSC_trig;
      TH2F* h_HFminus_energy_vs_RecHitTime_CSC_trig;
      TH2F* h_HF_CSC_Trig_Eta_Phi;	
      TH1F* h_S_over_L_plus_S_normal;
      TH1F* h_S_over_L_plus_S_abnormal;
      TH2F* h_Long_vs_Short;

      TH1F* h_calo_tower_energy;

      TH1F* h_Qiesum;
      TH2F* h_NoiseChan[9];
      TH1F* h_NoisePhi[9];
	      

      TH1F *h_jet_multiplicity;
      TH1F *h_jet_Pt;
      TH1F *h_jet_Eta;
      TH1F *h_jet_Phi;
      
      TH1F *h_leadJet_Pt;
      TH1F *h_leadJet_Eta;
      TH1F *h_leadJet_Phi;
      
      TH1F *h_caloMet_Met;
      TH1F *h_caloMet_Phi;
      TH1F *h_caloMet_SumEt;
      TH1F *h_HT;
      TH1F *h_MHT;

      TH1F* h_muon_vertex_x;
      TH1F* h_muon_px;
      TH1F* h_muon_p;
      
      TH2F* h_ecalx_vs_muonx;
      TH2F* h_ecaly_vs_muony;
      TH2F* h_jetphi_vs_muonphi;
      TH1F* h_impact_diff;
      
      TH2F*  DT_HCAL_eta_correlation;
      TH2F*  DT_HCAL_eta_correlation_all;
      TH2F*  DT_HCAL_phi_correlation;
      TH2F*  DT_HCAL_phi_correlation_all;
      
      TH1F*  HCAL_energy_correlation;
      TH1F*  HCAL_energy_correlation_all;

      TH1F* h_AllTracks;
      TH1F* h_ct1;
      TH1F* h_ct2;
      TH1F* h_ct3;
      TH1F* h_ct4;
      TH1F* h_ct5;
      TH1F* h_ct6;
      TH1F* h_ct7;
      TH1F* h_ct8;
      TH1F* h_ct9;
      TH1F* h_NTracks1;
      TH1F* h_NTracks2;
      TH1F* h_NPTracks;
      TH1F* h_NMTracks;
      TH1F* h_NMuons;
      TH1F* h_Pt;
      TH1F* h_Px;
      TH1F* h_Py;
      TH1F* h_d0;
      TH1F* h_Chi2;
      TH1F* h_Vtx_X;
      TH1F* h_Vtx_Y;
      TH1F* h_Vtx_Z;
      TH1F* h_Vtx_R1;
      TH1F* h_Vtx_R2;
      TH1F* h_diff;
      TH1F* h_InnerPos_X;
      TH1F* h_InnerPos_Y;
      TH1F* h_InnerPos_Z;
      TH1F* h_OuterPos_X;
      TH1F* h_OuterPos_Y;
      TH1F* h_OuterPos_Z;
      TH1F* h_Angle1;
      TH1F* h_Angle2;
      TH1F* h_Pt_Top;
      TH1F* h_Pt_Bottom;
      TH1F* h_Eta_Top;
      TH1F* h_Eta_Bottom;
      TH1F* h_Phi_Top;
      TH1F* h_Phi_Bottom;
      TH1F* h_dR_Top;
      TH1F* h_dR_Bottom;
      
      TH1F* h_HBP_e1x1;
      TH1F* h_HBP_e3x3;
      TH1F* h_HBP_e5x5;
      TH1F* h_HBM_e1x1;
      TH1F* h_HBM_e3x3;
      TH1F* h_HBM_e5x5;

      TH1F* h_HBP_e1[72];
      TH1F* h_HBP_e3[72];
      TH1F* h_HBP_e5[72];
      TH1F* h_HBM_e1[72];
      TH1F* h_HBM_e3[72];
      TH1F* h_HBM_e5[72];
      
      TH1F* h_HBTopMin_e1x1;
      TH1F* h_HBTopMin_e3x3cr;
      TH1F* h_HBTopMin_e1x3;
      TH1F* h_HBTopMin_e1x3cr;
      TH1F* h_HBTopMin_e5x5;

      TH1F* h_HBTopPlu_e1x1;
      TH1F* h_HBTopPlu_e3x3cr;
      TH1F* h_HBTopPlu_e1x3;
      TH1F* h_HBTopPlu_e1x3cr;
      TH1F* h_HBTopPlu_e5x5;
      
      TH1F* h_HBBottomMin_e1x1;
      TH1F* h_HBBottomMin_e3x3cr;
      TH1F* h_HBBottomMin_e1x3;
      TH1F* h_HBBottomMin_e1x3cr;
      TH1F* h_HBBottomMin_e5x5;
      
      TH1F* h_HBBottomPlu_e1x1;
      TH1F* h_HBBottomPlu_e3x3cr;
      TH1F* h_HBBottomPlu_e1x3;
      TH1F* h_HBBottomPlu_e1x3cr;
      TH1F* h_HBBottomPlu_e5x5;
      TH1F* h_corr;
      //timing for HB rechits
      TH1F *h_HBBottomPlu_timing[4];
      TH1F *h_HBBottomMin_timing[4];
      TH1F *h_HBTopPlu_timing[4];
      TH1F *h_HBTopMin_timing[4];
      //timing for muon signal
      TH1F *h_HBBottomPlu_Muontiming[4];
      TH1F *h_HBBottomMin_Muontiming[4];
      TH1F *h_HBTopPlu_Muontiming[4];
      TH1F *h_HBTopMin_Muontiming[4];


      TH1F* PROJ_eta;
      TH1F* PROJ_phi;
      TH1F*  FIN_zdist;
      TH2F*  h_EtaPhi;
      
      TH1F* h_HB_maxE;
      TH1F* h_HB_maxE2;
      TH1F* h_HB_maxEta;
      TH1F* h_HB_maxPhi;
      TH1F* h_HB_maxIEta;
      TH1F* h_HB_maxIPhi;
      TH2F* h_HB_maxIEtaIPhi;
      TH1F* h_HB_IEta;
      TH1F* h_HB_IPhi;
      TH1F* h_HB_DEta;
      TH1F* h_HB_DPhi;
      TH1F* h_HBEnt_IPhi;
      TH1F* h_HBEnt_IEta;
      TH1F* h_HBExt_IPhi;
      TH1F* h_HBExt_IEta;
      TH2F* h_HBCorrEta;
      TH2F* h_HBCorrPhi;
      TH2F* h_HBCorrEta1;
      TH2F* h_HBCorrPhi1;
      TH2F* h_HBCorrEta2;
      TH2F* h_HBCorrPhi2;
      TH1F* h_Trk1_x;
      TH1F* h_Trk1_y;
      TH1F* h_Trk1_z;
      TH1F* h_Trk2_x;
      TH1F* h_Trk2_y;
      TH1F* h_Trk2_z;

      TH1F* h_InnerPos_x;
      TH1F* h_InnerPos_y;
      TH1F* h_InnerPos_z;
      TH1F* h_InnerPos_eta;
      TH1F* h_InnerPos_phi;
      
      TH1F* h_OuterPos_x;
      TH1F* h_OuterPos_y;
      TH1F* h_OuterPos_z;
      TH1F* h_OuterPos_eta;
      TH1F* h_OuterPos_phi;

      TH1F* h_distance2D;
      TH1F* h_distance3D;
      TH2F* h_distance2D_vs_HB;
      TH2F* h_distance3D_vs_HB;
      TH2F* h_distance2D_vs_nHits;
      TH2F* h_distance2D_vs_chi2;
      
      TH1F* h_muonNHits;


  TH1F* Chi2;
  TH1F* Chi2overDoF;
  TH1F* DistanceOfClosestApproach;
  TH2F* DistanceOfClosestApproachVsTheta;
  TH2F* DistanceOfClosestApproachVsPhi;
  TH2F* DistanceOfClosestApproachVsEta;
  TH1F* xPointOfClosestApproach;
  TH1F* yPointOfClosestApproach;
  TH1F* zPointOfClosestApproach;

      int trigDT;
      std::string baseHtmlDir_;
      long runBegin,lumibegin,lumiend,evtNo;
      std::string startTime;
      bool prompt_htmlPrint;
      bool doDigiHTML;
      bool doCaloTowerHTML;
      bool doRecHitHTML;
      bool doJetMetHTML;
      bool doMuonHTML;
      bool doHPDNoiseHTML;

      CaloGeometryHelper *myCalorimeter_;
      int NDTT;  
      int NRPCT;  
      int NCSCT;  
      int NHcalT;  
      float fIPDistance;
      int fAlgoType; 
      double fIP_z;
      double fIP_r;
      double fHBThreshold;
      int NTotal;
      int NAccepted;
      int NAccepted2;
      int top_ihpi_low;
      int top_ihpi_high;
      int bottom_ihpi_low;
      int bottom_ihpi_high;
      int t1;
      int t2;
      int t3;
      int t4;
      int t5;
      int t6;
      int t7;

};
