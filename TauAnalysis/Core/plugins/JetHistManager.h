#ifndef TauAnalysis_Core_JetHistManager_h  
#define TauAnalysis_Core_JetHistManager_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TauAnalysis/Core/interface/HistManagerBase.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TauAnalysis/RecoTools/interface/PATJetAlphaExtractor.h"

#include <vector>
#include <string>

class JetHistManager : public HistManagerBase 
{
 public:  
  explicit JetHistManager(const edm::ParameterSet&);
  ~JetHistManager();
  
 private: 
//--- histogram booking and filling functions 
//    inherited from HistManagerBase class
  void bookHistogramsImp();
  void fillHistogramsImp(const edm::Event&, const edm::EventSetup&, double);

//--- auxiliary functions
  void bookJetHistograms(MonitorElement*&, MonitorElement*&, MonitorElement*&, const char*);
  
  double getJetWeight(const pat::Jet&);

  void fillJetHistograms(const pat::Jet&, MonitorElement*, MonitorElement*, MonitorElement*, double);
  void fillNumCentralJetsToBeVetoesHistograms(const std::vector<pat::Jet>&, MonitorElement*, double, double, double, double);

//--- configuration parameters
  edm::InputTag jetSrc_;
  edm::InputTag genParticleSrc_;

  bool requireGenJetMatch_;

  typedef std::vector<int> vint;
  vint skipPdgIdsGenParticleMatch_;

  typedef std::vector<double> vdouble;
  vdouble centralJetsToBeVetoedEtMin_;
  vdouble centralJetsToBeVetoedEtaMax_;
  vdouble centralJetsToBeVetoedAlphaMin_;

  std::vector<std::string>  bTaggingDiscriminators_;
  std::vector<double> bTaggingDiscriminatorThresholds_;

//--- auxiliary class to compute quantity alpha,
// defined as ratio of sum of charged particle transverse momenta 
// to sum of charged plus neutral particle transverse momenta,
// for a jet
  PATJetAlphaExtractor jetAlphaExtractor_;

//--- histograms
  MonitorElement* hNumJets_;
  MonitorElement* hSumEtJets_;

  MonitorElement* hJetPt_;
  MonitorElement* hJetEta_;
  MonitorElement* hJetPtVsEta_;
  MonitorElement* hJetPhi_;

  MonitorElement* hJetWeightPosLog_;
  MonitorElement* hJetWeightNegLog_;
  MonitorElement* hJetWeightZero_;
  MonitorElement* hJetWeightLinear_;

  MonitorElement* hJetMatchingGenParticlePdgId_;
  
  MonitorElement* hJetAlpha_;
  MonitorElement* hJetNumTracks_;
  MonitorElement* hJetTrkPt_;
  MonitorElement* hJetLeadTrkPt_;

  std::vector<MonitorElement*> hNumBtags_;
  std::vector<MonitorElement*> hPtBtags_;
  std::vector<MonitorElement*> hBtagDisc_;

  std::vector<MonitorElement*> hNumCentralJetsToBeVetoed_;
};

#endif  


