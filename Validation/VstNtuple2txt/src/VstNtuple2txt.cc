// -*- C++ -*-
//
// Package:    Ntuple2txt
// Class:      Ntuple2txt
// 
/**\class Ntuple2txt Ntuple2txt.cc Validation/VstNtuple2txt/src/Ntuple2txt.cc

 Description: Converts HepMC info to (text file)
   Vista file

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Stephen Mrenna
//         Created:  Wed Oct 10 12:52:50 CDT 2007
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/DataViewImpl.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimDataFormats/HepMCProduct/interface/GenInfoProduct.h"

//#include "HepMC/GenParticle.h"
//#include "HepMC/GenEvent.h"
//#include "HepMC/GenVertex.h"

#include "Validation/VstTurboSim/interface/QuaeroPartonObject.hh"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
//#include "CLHEP/Vector/ThreeVector.h"
#include <iostream>
#include <fstream>
//#include "Validation/VstTurboSim/interface/hepevt2parton.hh"
//#include "Validation/VstQuaeroUtils/interface/HepevtBlock.hh"
#include <sstream>


//
// class declaration
//

//void hepevt2partonEvent(ofstream& , const HepevtEvent& );
//void stdhep2partonWrapSimple(std::string, const HepevtEvent& );



#include <vector>
#include <map>
#include <set>
#include <algorithm>
//namespace edm { class ParameterSet; }
//namespace HepMC { class GenParticle; class GenEvent; }
//namespace reco { class GenParticleCandidate; }

///using namespace edm;
//using namespace reco;

const static double ptMinShower=1.0;
const static double etaMaxShower=5.0;


class Ntuple2txt : public edm::EDAnalyzer {
public:
  explicit Ntuple2txt(const edm::ParameterSet&);
  ~Ntuple2txt();


private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&);
  std::string fileName_;
  std::ofstream outputFile_;
  int numberOfEvents_;
  edm::RunNumber_t thisRun_;
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Ntuple2txt::Ntuple2txt(const edm::ParameterSet& iConfig) :
  fileName_(iConfig.getParameter<std::string>("fileName_"))
,
  outputFile_(fileName_.c_str()) 
{
  if(!outputFile_.good()) {
    std::cerr << " Cannot open file " << std::endl;
  }
}


Ntuple2txt::~Ntuple2txt()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
Ntuple2txt::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  numberOfEvents_++;
  Handle<reco::CandidateCollection> genParticles;
  iEvent.getByLabel("genParticleCandidates",genParticles);
  Handle<reco::CandidateCollection> particles;
  iEvent.getByLabel ("mergerVst",particles);

  Handle<reco::CandidateCollection> cands;
  Handle<reco::CandidateCollection> jetParticles;
  Handle<reco::GenJetCollection> generatorJet;
  Handle<reco::GenJetCollection> tauJet;
  iEvent.getByLabel ("VistaJet",generatorJet);
  iEvent.getByLabel ("VistaTauJet",tauJet);
  iEvent.getByLabel ("VistaJetClone",jetParticles);
  iEvent.getByLabel ("genVista","otherStableVst",cands);
  Handle<reco::CandMatchMap> particleJetMap;
  iEvent.getByLabel("particleJetMatch", particleJetMap);

  // "fix up" the jets
  // used matched particles to get jet mass
  // used parton ids in jets to get jet type
  size_t psize = jetParticles->size();
  std::vector<HepLorentzVector> newJets;
  std::vector<int> jetId;
  for( reco::GenJetCollection::const_iterator it = generatorJet->begin(); it!= generatorJet->end(); ++it) {
    //create a 0 vector for all jets
    newJets.push_back(HepLorentzVector(0,0,0,0));
    //get mc particles in jet
    std::vector<const reco::GenParticleCandidate*> mcparts = it->getConstituents();
    std::vector<const reco::GenParticleCandidate*>::const_iterator mcbegin = mcparts.begin(); 
    std::vector<const reco::GenParticleCandidate*>::const_iterator mcend = mcparts.end(); 
    // pick highest-pt parton to tag the jet type
    int tempId = (*mcbegin)->pdgId();
    if( mcparts.size() > 1 ) {
      double ptold = (*mcbegin)->pt();
      for( std::vector<const reco::GenParticleCandidate*>::const_iterator jt = mcbegin; jt!=mcend; ++jt) {
	if ( (*jt)->pt() > ptold ) {
	  ptold = (*jt)->pt();
	  tempId = (*jt)->pdgId();
	}
      }
    }
    jetId.push_back(tempId);
  }

  //accumulate the particle level 4 vectors in the jets
  for( size_t c = 0; c != cands->size(); ++ c ) {
    reco::CandidateRef cref( cands, c );
    reco::CandMatchMap::const_iterator f = particleJetMap->find( cref );
    if( f != particleJetMap->end() ) {
      /*      int id = ( f->val->pdgId() );
      std::cout << c << ") found match with id: " << id << " " << f->val.key() << " " 
		<< "; pt, eta, phi = " << cref->pt() << ", " << cref->eta() << ", " << cref->phi()
		<< " mass = " << cref->mass();
		std::cout << " [*]"; */
      newJets[f->val.key()] += HepLorentzVector(cref->px(),cref->py(),cref->pz(),cref->energy());
    }
    //    std::cout << std::endl;
  }

  //string to accumulate Vista information
  std::stringstream hepEvtStream;

  std::string printFormat = "(m)-pt-eta-phi(deg)";

  HepLorentzVector rootS = HepLorentzVector(0,0,0,0);
  double genSumPt = 0.0;
  for( reco::CandidateCollection::const_iterator imc = genParticles->begin(); imc != genParticles->end(); ++imc) {
    if(imc->status()==1) {
      rootS += HepLorentzVector(imc->px(),imc->py(),imc->pz(),imc->energy());
      genSumPt += imc->pt();
    }
  }
  double energyCMS = rootS.e();
  if( fabs(energyCMS-14000) < .01 ) energyCMS = 14000;

  std::stringstream ss;
  std::string collisionType = "pp", generatedSumPt = "";
  ss << genSumPt;
  ss >> generatedSumPt;

  double wt = 1, zvtx = 0;

  hepEvtStream << "sig " << iEvent.id().event() << " " << wt << " " 
	    << collisionType << " " << energyCMS << " " << generatedSumPt << " ";

  hepEvtStream << zvtx << " ";

  for( reco::CandidateCollection::const_iterator p  = particles->begin();
       p != particles->end(); ++p) {

    std::string objectType="j";
    // Particle Name
    int idabs = abs(p->pdgId());
    // Particle Index
    if( idabs==11 ) {
      objectType="e";
      objectType += ( p->pdgId() >0 ? "-" : "+" );		  
    } else if ( idabs==13 ) {
      objectType="mu";
      objectType += ( p->pdgId() >0 ? "-" : "+" );		  
    } else if ( idabs==22 ) {
      objectType="ph";
    }

    hepEvtStream << QuaeroPartonObject(objectType,HepLorentzVector(p->px(),p->py(),p->pz(),p->energy())).print(printFormat) << " ";
  }

  // add this in later --> maybe just treat as a jet!

  //edit out treatment of taus
  /*  for( reco::GenJetCollection::const_iterator p  = tauJet->begin();
       p != tauJet->end(); ++p) {
    std::string objectType="tau";
    hepEvtStream << QuaeroPartonObject(objectType,HepLorentzVector(p->px(),p->py(),p->pz(),p->energy())).print(printFormat) << " ";
    }*/


  for( size_t c = 0; c != psize; ++ c ) {
    std::string objectType="j";
    if( abs(jetId[c])== 5 ) objectType="b";
    if( newJets[c].perp()>5.0) hepEvtStream << QuaeroPartonObject(objectType,newJets[c]).print(printFormat) << " ";
  }

  hepEvtStream << ";" << std::endl;
  outputFile_ << hepEvtStream.str();
}


// ------------ method called once each job just before starting event loop  ------------
void 
Ntuple2txt::beginJob(const edm::EventSetup& )
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Ntuple2txt::endJob() {

}

void Ntuple2txt::beginRun(const edm::Run & r, const edm::EventSetup & e) {

  using namespace edm;

  numberOfEvents_=0;

  thisRun_ = r.run();
  std::cout << "Run number is " << thisRun_ << std::endl;

  Handle< GenInfoProduct > gi;
  if( r.getByLabel("source", gi) ) {

    double auto_cross_section = gi->cross_section(); 
// calculated at the end of RUN --  units in nb
    double external_cross_section = gi->external_cross_section(); 
// is the precalculated one written in the cfg file -- units is pb

    std::cout << auto_cross_section << " " << external_cross_section << std::endl;
 }

}



void Ntuple2txt::endRun(const edm::Run & r, const edm::EventSetup & e) {
 
  std::cout << "Number of events in run " << numberOfEvents_ << std::endl;

}




//define this as a plug-in
DEFINE_FWK_MODULE(Ntuple2txt);

