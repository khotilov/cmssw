#ifndef JetAnalyzer_H
#define JetAnalyzer_H


/** \class JetAnalyzer
 *
 *  DQM monitoring source for Calo Jets
 *
 *  $Date: 2009/06/30 13:48:04 $
 *  $Revision: 1.1 $
 *  \author F. Chlebana - Fermilab
 */


#include <memory>
#include <fstream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMOffline/JetMET/interface/JetAnalyzerBase.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
//
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "RecoJets/JetAlgorithms/interface/JetIDHelper.h"

#include <string>
using namespace std;
using namespace edm;

class JetAnalyzer : public JetAnalyzerBase {
 public:

  /// Constructor
  //  JetAnalyzer(const edm::ParameterSet&, JetServiceProxy *theService);
  JetAnalyzer(const edm::ParameterSet&);
  
  /// Destructor
  virtual ~JetAnalyzer();

  /// Inizialize parameters for histo binning
  void beginJob(edm::EventSetup const& iSetup, DQMStore *dbe);

  /// Finish up a job
  void endJob();

  /// Get the analysis
    //  void analyze(const edm::Event&, const edm::EventSetup&, const edm::TriggerResults&,
    //	       const reco::CaloJet& caloJet);
  void analyze(const edm::Event&, const edm::EventSetup&, 
	       const reco::CaloJet& caloJet);

  void setSource(std::string source) {
    _source = source;
  }

  void setLeadJetFlag(int flag) {
    _leadJetFlag = flag;
  }
  int getLeadJetFlag() {
    return  _leadJetFlag;
  }
  void setNJets(int njets) {
    _NJets = njets;
  }
  int getNJets() {
    return  _NJets;
  }
  void setJetLoPass(int pass) {
    _JetLoPass = pass;
  }

  void setJetHiPass(int pass) {
    _JetHiPass = pass;
  }
  void setDPhi(double dphi) {
    _DPhi = dphi;
  }
  double getDPhi() {
    return  _DPhi;
  }


 private:
  // ----------member data ---------------------------
  
  edm::ParameterSet parameters;
  // Switch for verbosity
  std::string jetname;
  std::string _source;
  // Calo Jet Label
  edm::InputTag theCaloJetCollectionLabel;

  int   _JetLoPass;
  int   _JetHiPass;
  int   _leadJetFlag;
  int   _NJets;
  double _ptThreshold;
  double _fHPDMax;
  double _resEMFMin;
  int _n90HitsMin;
  double _DPhi;

  //histo binning parameters
  int    etaBin;
  double etaMin;
  double etaMax;

  int    phiBin;
  double phiMin;
  double phiMax;

  int    ptBin;
  double ptMin;
  double ptMax;

  int    eBin;
  double eMin;
  double eMax;

  int    pBin;
  double pMin;
  double pMax;

  //the histos
  MonitorElement* jetME;

  // JetID helper
  reco::helper::JetIDHelper *jetID;

  // Calo Jets

  //  std::vector<MonitorElement*> etaCaloJet;
  //  std::vector<MonitorElement*> phiCaloJet;
  //  std::vector<MonitorElement*> ptCaloJet;
  //  std::vector<MonitorElement*> qGlbTrack;

  //  MonitorElement* etaCaloJet;
  //  MonitorElement* phiCaloJet;
  //  MonitorElement* ptCaloJet;

  // Generic Jet Parameters

  // --- Used for Data Certification
  MonitorElement* mPt;
  MonitorElement* mPt_1;
  MonitorElement* mPt_2;
  MonitorElement* mPt_3;
  MonitorElement* mEta;
  MonitorElement* mPhi;
  MonitorElement* mConstituents;
  MonitorElement* mHFrac;
  MonitorElement* mEFrac;
  MonitorElement* mPhiVSEta;

  MonitorElement* mPt_Barrel;
  MonitorElement* mEta_Barrel;
  MonitorElement* mPhi_Barrel;

  MonitorElement* mPt_EndCap;
  MonitorElement* mEta_EndCap;
  MonitorElement* mPhi_EndCap;

  MonitorElement* mPt_Forward;
  MonitorElement* mEta_Forward;
  MonitorElement* mPhi_Forward;

  MonitorElement* mPt_Barrel_Lo;
  MonitorElement* mEta_Barrel_Lo;
  MonitorElement* mPhi_Barrel_Lo;
  MonitorElement* mConstituents_Barrel_Lo;
  MonitorElement* mHFrac_Barrel_Lo;
  MonitorElement* mPt_EndCap_Lo;
  MonitorElement* mEta_EndCap_Lo;
  MonitorElement* mPhi_EndCap_Lo;
  MonitorElement* mConstituents_EndCap_Lo;
  MonitorElement* mHFrac_EndCap_Lo;
  MonitorElement* mPt_Forward_Lo;
  MonitorElement* mEta_Forward_Lo;
  MonitorElement* mPhi_Forward_Lo;
  MonitorElement* mConstituents_Forward_Lo;
  MonitorElement* mHFrac_Forward_Lo;

  MonitorElement* mPt_Barrel_Hi;
  MonitorElement* mEta_Barrel_Hi;
  MonitorElement* mPhi_Barrel_Hi;
  MonitorElement* mConstituents_Barrel_Hi;
  MonitorElement* mHFrac_Barrel_Hi;
  MonitorElement* mPt_EndCap_Hi;
  MonitorElement* mEta_EndCap_Hi;
  MonitorElement* mPhi_EndCap_Hi;
  MonitorElement* mConstituents_EndCap_Hi;
  MonitorElement* mHFrac_EndCap_Hi;
  MonitorElement* mPt_Forward_Hi;
  MonitorElement* mEta_Forward_Hi;
  MonitorElement* mPhi_Forward_Hi;
  MonitorElement* mConstituents_Forward_Hi;
  MonitorElement* mHFrac_Forward_Hi;
  // ---


  MonitorElement* mE_Barrel;
  MonitorElement* mE_EndCap;
  MonitorElement* mE_Forward;

  MonitorElement* mE;
  MonitorElement* mP;
  MonitorElement* mMass;
  MonitorElement* mNJets;
  MonitorElement* mDPhi;

  // Leading Jet Parameters
  MonitorElement* mEtaFirst;
  MonitorElement* mPhiFirst;
  MonitorElement* mEFirst;
  MonitorElement* mPtFirst;


  // CaloJet specific
  MonitorElement* mMaxEInEmTowers;
  MonitorElement* mMaxEInHadTowers;
  MonitorElement* mHadEnergyInHO;
  MonitorElement* mHadEnergyInHB;
  MonitorElement* mHadEnergyInHF;
  MonitorElement* mHadEnergyInHE;
  MonitorElement* mEmEnergyInEB;
  MonitorElement* mEmEnergyInEE;
  MonitorElement* mEmEnergyInHF;
  //  MonitorElement* mEnergyFractionHadronic;
  //  MonitorElement* mEnergyFractionEm;
  MonitorElement* mN90Hits;
  MonitorElement* mfHPD;
  MonitorElement* mresEMF;


  // Events passing the jet triggers
  MonitorElement* mEta_Lo;
  MonitorElement* mPhi_Lo;
  MonitorElement* mPt_Lo;

  MonitorElement* mEta_Hi;
  MonitorElement* mPhi_Hi;
  MonitorElement* mPt_Hi;

};
#endif
