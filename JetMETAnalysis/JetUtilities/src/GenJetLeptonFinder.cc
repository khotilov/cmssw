////////////////////////////////////////////////////////////////////////////////
//
// GenJetLeptonFinder_cc
// ---------------------
//
// -> given a GenJet, this class digs for lepton (mu/el)
//    and associated neutrino constituents
//    
//
//                   16.11.2009 Hauke Held <held@cern.ch>
//
////////////////////////////////////////////////////////////////////////////////


#include "JetMETAnalysis/JetUtilities/interface/GenJetLeptonFinder.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include <iostream>
#include <memory>
#include <vector>
#include <string>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// -------------------- constructor && destructor ------------------------------
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
GenJetLeptonFinder::GenJetLeptonFinder(const reco::Candidate& genJet) :
  genJet_(genJet),
  lepton_(0),
  neutrino_(0)
{

}

//______________________________________________________________________________
GenJetLeptonFinder::~GenJetLeptonFinder()
{

}


////////////////////////////////////////////////////////////////////////////////
// -------------------- member function implementation -------------------------
////////////////////////////////////////////////////////////////////////////////

bool GenJetLeptonFinder::run()
{
  // get the GenJet constituents (but sticking to the Candidate interface)
  unsigned njetdaus = genJet_.numberOfDaughters();
  vector<const reco::Candidate*> jetdaus;
  for (unsigned i=0;i<njetdaus;i++) jetdaus.push_back(genJet_.daughter(i));
  
  // loop over constituents until lepton is found, then search neutrino
  bool foundLep = false;
  for (unsigned itd=0;itd<njetdaus&&!foundLep;itd++){
    
    // retrieve constituent
    const reco::Candidate* dau = jetdaus[itd];
    int dauid = dau->pdgId();
    
    // decide wether lepton or not
    if (abs(dauid)!=11 && abs(dauid)!=13)continue;
    
    // if dau is lepton save as lepton_
    lepton_  = dau;
    foundLep = true;
    
    // check for mother
    const reco::Candidate* mo  = lepton_->mother();
    assert (0!=mo);
    unsigned nmodaus = mo->numberOfDaughters();
    bool foundNu = false;
    for (unsigned itcand=0;itcand<nmodaus&&!foundNu;itcand++){
      const reco::Candidate* nucand = mo->daughter(itcand);
      if (lepton_==nucand)continue;
      if ( (lepton_->pdgId())*(nucand->pdgId())==-132 ||
	   (lepton_->pdgId())*(nucand->pdgId())==-182 ){
	neutrino_ = nucand;
	foundNu   = true;
      }
    }
  }
  return true;
}
