#ifndef McSelector_h
#define McSelector_h

/*  \class McSelector
*
*  Class to apply analysis cuts in the TriggerValidation Code
*
*  Author: Massimiliano Chiorboli      Date: August 2007
//         Maurizio Pierini
//         Maria Spiropulu
*
*/
#include <memory>
#include <string>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"


class McSelector {

 public:
  
  //Constructor
  McSelector(edm::ParameterSet userCut_params);
  //Destructor
  virtual ~McSelector(){};

  //Methods
  void handleObjects(const edm::Event&);
  bool isSelected(const edm::Event&);

 private:
  
  // Define the parameters
  std::string m_genSrc;
  std::string m_genJetSrc;
  std::string m_genMetSrc;
  double mc_ptElecMin;
  double mc_ptMuonMin;
  double mc_ptPhotMin;
  double mc_ptJetMin;
  double mc_metMin;
  int    mc_nElec;
  int    mc_nMuon;
  int    mc_nPhot;
  int    mc_nJet;



  const reco::GenParticleCollection* theGenParticleCollection;
  const reco::GenJetCollection*    theGenJetCollection;
  const reco::GenMETCollection*    theGenMETCollection;

};

#endif
