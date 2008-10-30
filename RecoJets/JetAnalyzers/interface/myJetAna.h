#ifndef RecoExamples_myJetAna_h
#define RecoExamples_myJetAna_h
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>

/* \class myJetAna
 *
 * \author Frank Chlebana
 *
 * \version 1
 *
 */
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

// class TFile;

/****
class RBX {
  RBX();

 private:
  int nTowers;
  int ieta;
  int iphi;
  float energy;
  float time;
};

class RBXCollection {

  RBXCollection();
  void addRBX(RBX r)  {rbx_.push_back(r)};

 private:
  std::vector<RBX> rbx_;

};
*****/


class myJetAna : public edm::EDAnalyzer {

public:
  myJetAna( const edm::ParameterSet & );

private:
  void beginJob( const edm::EventSetup & );
  void analyze ( const edm::Event& , const edm::EventSetup& );
  void endJob();

  std::string CaloJetAlgorithm;
  std::string GenJetAlgorithm;
  edm::InputTag theTriggerResultsLabel;
  std::string JetCorrectionService;


  // --- Passed selection cuts
  TH1F *h_pt;
  TH1F *h_ptTower;
  TH1F *h_ptRBX;
  TH1F *h_ptHPD;
  TH1F *h_et;
  TH1F *h_eta;
  TH1F *h_phi;
  // ---
  
  // --- RecHits
  TH1F *HBEne;
  TH1F *HBTime;
  TH1F *HEEne;
  TH1F *HETime;
  TH1F *HFEne;
  TH1F *HFTime;
  TH1F *HOEne;
  TH1F *HOTime;
  TH1F *EBEne;
  TH1F *EBTime;
  TH1F *EEEne;
  TH1F *EETime;

  TH1F *RBX_et;
  TH1F *RBX_hadEnergy;
  TH1F *RBX_hcalTime;
  TH1F *RBX_nTowers;
  TH1F *RBX_N;

  TH1F *HPD_et;
  TH1F *HPD_hadEnergy;
  TH1F *HPD_hcalTime;
  TH1F *HPD_nTowers;
  TH1F *HPD_N;

  // --- from reco calomet
  TH1F *SumEt;
  TH1F *MET;
  TH1F *MET_Tower;
  TH1F *MET_RBX;
  TH1F *MET_HPD;
  TH1F *METSig;
  TH1F *MEx;
  TH1F *MEy;
  TH1F *METPhi;
  // ---

  // --- from reco vertexs
  TH1F *h_Vx;
  TH1F *h_Vy;
  TH1F *h_Vz;
  TH1F *h_VNTrks;
  // ---

  // --- from reco tracks
  TH1F *h_Trk_pt;
  TH1F *h_Trk_NTrk;
  // ---
 
  TH1F *hf_sumTowerAllEx; 
  TH1F *hf_sumTowerAllEy;
  TH1F *hf_TowerJetEt;

  TH1F *ETime; 
  TH1F *HTime; 

  TH1F *nTowers1; 
  TH1F *nTowers2; 
  TH1F *nTowers3; 
  TH1F *nTowers4;
  TH1F *nTowersLeadJetPt1; 
  TH1F *nTowersLeadJetPt2; 
  TH1F *nTowersLeadJetPt3; 
  TH1F *nTowersLeadJetPt4;

  TH1F *totEneLeadJetEta1;
  TH1F *totEneLeadJetEta2; 
  TH1F *totEneLeadJetEta3;
  TH1F *hadEneLeadJetEta1; 
  TH1F *hadEneLeadJetEta2; 
  TH1F *hadEneLeadJetEta3;
  TH1F *emEneLeadJetEta1;  
  TH1F *emEneLeadJetEta2;  
  TH1F *emEneLeadJetEta3;

  TH1F *hadFracEta1; 
  TH1F *hadFracEta2; 
  TH1F *hadFracEta3;

  TH1F *tMassGen;

  TH1F *dijetMass;

  TH1F *h_nCalJets;
  TH1F *h_nGenJets;

  TH1F *h_ptCal;
  TH1F *h_etaCal;
  TH1F *h_phiCal;

  TH1F *h_ptGen; 
  TH1F *h_etaGen; 
  TH1F *h_phiGen;

  TH1F *h_ptGenL;
  TH1F *h_etaGenL;
  TH1F *h_phiGenL;

  TH1F *h_jetEt;

  TH1F *h_UnclusteredEt;
  TH1F *h_UnclusteredEts;
  TH1F *h_TotalUnclusteredEt;

  TH1F *h_UnclusteredE;
  TH1F *h_TotalUnclusteredE;

  TH1F *h_ClusteredE;
  TH1F *h_TotalClusteredE;

  TH1F *h_jet1Pt;
  TH1F *h_jet2Pt;
  TH1F *h_jet1PtHLT;

  TH1F *EMFraction;
  TH1F *NTowers;

  TH2F *h_EmEnergy;
  TH2F *h_HadEnergy;

};

#endif
