#ifndef HLTriggerOffline_Egamma_EmDQM_H
#define HLTriggerOffline_Egamma_EmDQM_H


// Base Class Headers
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <vector>
#include "TDirectory.h"
#include "HepMC/GenParticle.h"
#include "PhysicsTools/Utilities/interface/PtComparator.h"

class EmDQM : public edm::EDAnalyzer{
public:
  /// Constructor
  explicit EmDQM(const edm::ParameterSet& pset);

  /// Destructor
  ~EmDQM();

  // Operations

  void analyze(const edm::Event & event, const edm::EventSetup&);
  void beginJob(const edm::EventSetup&);
  void endJob();

private:
  // Input from cfg file
  std::vector<edm::InputTag> theHLTCollectionLabels;  
  unsigned int numOfHLTCollectionLabels;  // Will be size of above vector
  edm::InputTag theL1Seed;
  std::vector<int> theHLTOutputTypes;
  std::vector<bool> plotiso;
  std::vector<std::vector<edm::InputTag> > isoNames; // there has to be a better solution
  std::vector<std::pair<double,double> > plotBounds; 
  std::string theHltName;

  ////////////////////////////////////////////////////////////
  //          Read from configuration file                  //
  ////////////////////////////////////////////////////////////
  // paramters for generator study
  unsigned int reqNum;
  int   pdgGen;
  double genEtaAcc;
  double genEtAcc;
  // plotting paramters
  double plotEtaMax;
  double plotPtMin ;
  double plotPtMax ;
  unsigned int plotBins ;
  // preselction cuts
  edm::InputTag gencutCollection_;
  unsigned int gencut_;


  ////////////////////////////////////////////////////////////
  //          Create Histograms                             //
  ////////////////////////////////////////////////////////////
  // Et & eta distributions
  std::vector<MonitorElement*> etahist;
  std::vector<MonitorElement*> ethist;
  std::vector<MonitorElement*> etahistmatch;
  std::vector<MonitorElement*> ethistmatch;
  std::vector<MonitorElement*> histEtOfHltObjMatchToGen;
  std::vector<MonitorElement*> histEtaOfHltObjMatchToGen;
  // Isolation distributions
  std::vector<MonitorElement*> etahistiso;
  std::vector<MonitorElement*> ethistiso;
  std::vector<MonitorElement*> etahistisomatch;
  std::vector<MonitorElement*> ethistisomatch;
  std::vector<MonitorElement*> histEtIsoOfHltObjMatchToGen;
  std::vector<MonitorElement*> histEtaIsoOfHltObjMatchToGen;
  // Plots of efficiency per step
  MonitorElement* total;
  MonitorElement* totalmatch;
  //generator histograms
  MonitorElement* etgen;
  MonitorElement* etagen;

  // interface to DQM framework
  DQMStore * dbe;
  std::string dirname_;

  template <class T> void fillHistos(edm::Handle<trigger::TriggerEventWithRefs>& ,const edm::Event& ,unsigned int, std::vector<reco::Particle>& );
  GreaterByPt<reco::Particle> pTComparator_;
  
};
#endif
