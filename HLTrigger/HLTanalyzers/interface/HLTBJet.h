#ifndef HLTrigger_HLTanalyzers_HLTBJet_h
#define HLTrigger_HLTanalyzers_HLTBJet_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

class TTree;

class HLTBJet {
public:
  HLTBJet();
  ~HLTBJet();
  
  void setup(const edm::ParameterSet & config, TTree * tree);
  void clear(void);
  void analyze(const edm::Event & event, const edm::EventSetup & setup, TTree * tree);

private:
  void analyseJets(
      const edm::View<reco::Jet>   & jets);
  
  void analyseCorrectedJets(
      const edm::View<reco::Jet>   & jets);
  
  void analyseLifetime(
      const edm::View<reco::Jet>   & jets, 
      const reco::JetTagCollection & tagsL25, 
      const reco::JetTagCollection & tagsL3);

  void analyseLifetimeLoose(
      const edm::View<reco::Jet>   & jets, 
      const reco::JetTagCollection & tagsL25, 
      const reco::JetTagCollection & tagsL3);

  void analyseSoftmuon(
      const edm::View<reco::Jet>   & jets, 
      const reco::JetTagCollection & tagsL25, 
      const reco::JetTagCollection & tagsL3);

  void analysePerformance(
      const edm::View<reco::Jet>   & jets, 
      const reco::JetTagCollection & tagsL25, 
      const reco::JetTagCollection & tagsL3);

  // labels for input HLT collections
  edm::InputTag m_jets;
  edm::InputTag m_correctedJets;
  edm::InputTag m_lifetimeBJetsL25;
  edm::InputTag m_lifetimeBJetsL3;
  edm::InputTag m_lifetimeBJetsL25Relaxed;
  edm::InputTag m_lifetimeBJetsL3Relaxed;
  edm::InputTag m_softmuonBJetsL25;
  edm::InputTag m_softmuonBJetsL3;
  edm::InputTag m_performanceBJetsL25;
  edm::InputTag m_performanceBJetsL3;

  // set of variables for uncorrected L2 jets
  int NohBJetL2;
  float * ohBJetL2Energy;
  float * ohBJetL2Et;
  float * ohBJetL2Pt;
  float * ohBJetL2Eta;
  float * ohBJetL2Phi;
              
  // set of variables for corrected L2 jets
  int NohBJetL2Corrected;
  float * ohBJetL2CorrectedEnergy;
  float * ohBJetL2CorrectedEt;
  float * ohBJetL2CorrectedPt;
  float * ohBJetL2CorrectedEta;
  float * ohBJetL2CorrectedPhi;
  
  // set of variables for lifetime-based b-tag
  float * ohBJetIPL25Tag;
  float * ohBJetIPL3Tag;
  
  // set of variables for lifetime-based relaxed b-tag
  float * ohBJetIPLooseL25Tag;
  float * ohBJetIPLooseL3Tag;
  
  // set of variables for soft-muon-based b-tag
  int   * ohBJetMuL25Tag;           // do not optimize
  float * ohBJetMuL3Tag;
  
  // set of variables for b-tagging performance measurements
  int   * ohBJetPerfL25Tag;         // do not optimize 
  int   * ohBJetPerfL3Tag;          // do not optimize
};

#endif // HLTrigger_HLTanalyzers_HLTBJet_h
