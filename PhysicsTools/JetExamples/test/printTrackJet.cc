// system include files
#include <memory>
#include <string>
#include <iostream>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

using namespace std;
using namespace reco;
using namespace edm;

class printTrackJet : public edm::EDAnalyzer {
  public:
    explicit printTrackJet(const edm::ParameterSet & );
    ~printTrackJet() {};
    void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
     
  private:

    edm::InputTag source_;
    edm::Handle<reco::CandidateCollection> trackJets;
};

printTrackJet::printTrackJet(const edm::ParameterSet& iConfig)
{
  source_  = iConfig.getParameter<InputTag> ("src");
}

void printTrackJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  cout << "[printTrackJet] analysing event " << iEvent.id() << endl;
  
  try {
    iEvent.getByLabel (source_ ,trackJets);
  } catch(std::exception& ce) {
    cerr << "[printTrackJet] caught std::exception " << ce.what() << endl;
    return;
  }

  cout << "************************" << endl;
  cout << "* TrackJetCollection  *" << endl;
  cout << "************************" << endl;
  for( CandidateCollection::const_iterator f  = trackJets->begin();
                                           f != trackJets->end();
                                           f++) {

     printf("[printTrackJet] (pt,eta,phi) = %7.3f %6.3f %6.3f |\n",
              f->et(),
              f->eta(),
              f->phi()  );

     for( Candidate::const_iterator c  = f->begin();   
                                    c != f->end();   
                                    c ++) {  
       printf("        [Constituents] (pt,eta,phi) = %6.2f %5.2f %5.2f|\n",
               c->et(),                                                                       
               c->eta(),                                                                      
               c->phi() );
     }                                                                                          
  }
}

DEFINE_FWK_MODULE( printTrackJet );
