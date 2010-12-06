#ifndef EwkMuLumiMonitorDQM_H
#define EwkMuLumiMonitorDQM_H

/** \class EwkMuLumiMonitorDQM
 *
 *  DQM offline for EWK MuLumiMonitor: intended for luminosity purposes using Z/W  
 *  \authors:  Michele de Gruttola, INFN Naples - Maria Cepeda, CIEMAT 
 *  on behalf EWK-Muon group 
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"

class DQMStore;
class MonitorElement;
class EwkMuLumiMonitorDQM : public edm::EDAnalyzer {
public:
  EwkMuLumiMonitorDQM (const edm::ParameterSet &);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&);

  void init_histograms();
  double muIso(const reco::Muon &);
  double tkIso(reco::Track, edm::Handle<reco::TrackCollection>, edm::Handle<CaloTowerCollection>);
  bool IsMuMatchedToHLTMu ( const reco::Muon & , std::vector<reco::Particle> ,double ,double );

private:

  edm::InputTag trigTag_; 
  edm::InputTag trigEv_;
  edm::InputTag muonTag_;
  edm::InputTag trackTag_;
  edm::InputTag caloTowerTag_;
  edm::InputTag metTag_;
  bool  metIncludesMuons_;

  // const std::string hltPath_;  
  //  const std::string L3FilterName_;    
 
  double ptMuCut_;
  double etaMuCut_;

  bool isRelativeIso_;
  bool isCombinedIso_;
  double isoCut03_;
 
  double deltaRTrk_; 
  double ptThreshold_; 
  double deltaRVetoTrk_; 
  double maxDPtRel_; 
  double maxDeltaR_;
  double mtMin_;
  double mtMax_;
  double acopCut_; 
  double  dxyCut_;


  

  DQMStore* theDbe;

  
  MonitorElement* mass2HLT_;
  MonitorElement* highMass2HLT_;
  // MonitorElement* highest_mupt2HLT_;
  //MonitorElement* lowest_mupt2HLT_;

  MonitorElement* mass1HLT_;
  MonitorElement* highMass1HLT_;
  //  MonitorElement* highest_mupt1HLT_;
  //MonitorElement* lowest_mupt1HLT_;

  MonitorElement* massNotIso_;
  MonitorElement* highMassNotIso_;
  //MonitorElement* highest_muptNotIso_;
  // MonitorElement* lowest_muptNotIso_;

  MonitorElement* massGlbSta_;
  MonitorElement* highMassGlbSta_;
  // MonitorElement* highest_muptGlbSta_;
  //MonitorElement* lowest_muptGlbSta_;

  MonitorElement* massGlbTrk_;
  MonitorElement* highMassGlbTrk_;
  //  MonitorElement* highest_muptGlbTrk_;
  //MonitorElement* lowest_muptGlbTrk_;

  MonitorElement* TMass_;

  MonitorElement* massIsBothGlbTrkThanW_;
  MonitorElement* highMassIsBothGlbTrkThanW_;
  

  unsigned int nall;
  unsigned int nEvWithHighPtMu; 
  unsigned int nInKinRange;
  unsigned int nsel;
  unsigned int niso;
  unsigned int nhlt;
  unsigned int n1hlt;
  unsigned int n2hlt;
  unsigned int nNotIso;
  unsigned int nGlbSta;
  unsigned int nGlbTrk;
  unsigned int nTMass;
  unsigned int nW;

  bool isZGolden1HLT_;  
  bool isZGolden2HLT_; 
  bool isZGoldenNoIso_;
  bool isZGlbSta_;
  bool isZGlbTrk_;
  bool isW_;

  bool isValidHltConfig_;
  HLTConfigProvider  hltConfigProvider_;


};

#endif
