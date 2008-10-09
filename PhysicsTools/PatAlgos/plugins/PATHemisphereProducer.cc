// -*- C++ -*-
//
// Package:    PatShapeAna
// Class:      PatShapeAna
// 
/**\class PatShapeAna PatShapeAna.cc PhysicsTools/PatShapeAna/src/PatShapeAna.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Tanja Rommerskirchen
//         Created:  Sat Mar 22 12:58:04 CET 2008
// $Id: PATHemisphereProducer.cc,v 1.6 2008/09/29 09:53:00 gpetrucc Exp $
//
//


//system
#include <vector>
#include <memory>
//PAT
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Hemisphere.h"
//DataFormats
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/Particle.h"
//User
#include  "PhysicsTools/PatAlgos/plugins/PATHemisphereProducer.h"


using namespace pat;


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PATHemisphereProducer::PATHemisphereProducer(const edm::ParameterSet& iConfig) :
  _patJets       ( iConfig.getParameter<edm::InputTag>( "patJets" ) ),
  _patMuons      ( iConfig.getParameter<edm::InputTag>( "patMuons" ) ),
  _patElectrons  ( iConfig.getParameter<edm::InputTag>( "patElectrons" ) ),
  _patPhotons    ( iConfig.getParameter<edm::InputTag>( "patPhotons" ) ),
  _patTaus       ( iConfig.getParameter<edm::InputTag>( "patTaus" ) ),

  _minJetEt       ( iConfig.getParameter<double>("minJetEt") ),
  _minMuonEt       ( iConfig.getParameter<double>("minMuonEt") ),
  _minElectronEt       ( iConfig.getParameter<double>("minElectronEt") ),
  _minTauEt       ( iConfig.getParameter<double>("minTauEt") ), 
  _minPhotonEt       ( iConfig.getParameter<double>("minPhotonEt") ),

  _maxJetEta       ( iConfig.getParameter<double>("maxJetEta") ),
  _maxMuonEta       ( iConfig.getParameter<double>("maxMuonEta") ),
  _maxElectronEta       ( iConfig.getParameter<double>("maxElectronEta") ),
  _maxTauEta       ( iConfig.getParameter<double>("maxTauEta") ), 
  _maxPhotonEta       ( iConfig.getParameter<double>("maxPhotonEta") ),

  _seedMethod    ( iConfig.getParameter<int>("seedMethod") ),
  _combinationMethod ( iConfig.getParameter<int>("combinationMethod") )
  //  _EJselectionCfg(iConfig.getParameter<edm::ParameterSet>("ElectronJetCrossCleaning")),    
  // _ElectronJetCC(reco::modules::make<ElectronJetCrossCleaner>(_EJselectionCfg))
{
  //produces<std::vector<std::double
  ///produces cross-cleaned collections of above objects
  //Alternative: produce cross-cleaning decision & MET correction per object
//    produces<HemiAxis>("hemi1"); //hemisphere 1 axis
//    produces<HemiAxis>("hemi2"); //hemisphere 1 axis
 
//   produces<std::vector<pat::Jet> >();
//   produces<std::vector<pat::MET> >();
//   produces<std::vector<pat::Muon> >();
//   produces<std::vector<pat::Electron> >();
//   produces<std::vector<pat::Photon> >();
//   produces<std::vector<pat::Tau> >();

  produces< std::vector<pat::Hemisphere> >();
}


PATHemisphereProducer::~PATHemisphereProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PATHemisphereProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   //Jets   
   Handle<reco::CandidateView> pJets;
   iEvent.getByLabel(_patJets,pJets);

   //Muons   
   Handle<reco::CandidateView> pMuons;
   iEvent.getByLabel(_patMuons,pMuons);

   //Electrons   
   Handle<reco::CandidateView> pElectrons;
   iEvent.getByLabel(_patElectrons,pElectrons);

   //Photons   
   Handle<reco::CandidateView> pPhotons;
   iEvent.getByLabel(_patPhotons,pPhotons);

   //Taus   
   Handle<reco::CandidateView> pTaus;
   iEvent.getByLabel(_patTaus,pTaus);


   //fill e,p vector with information from all objects (hopefully cleaned before)
   for(int i = 0; i < (int) (*pJets).size() ; i++){
     if((*pJets)[i].pt() <  _minJetEt || fabs((*pJets)[i].eta()) >  _maxJetEta) continue;
     vPx.push_back((*pJets)[i].px());
     vPy.push_back((*pJets)[i].py());
     vPz.push_back((*pJets)[i].pz());
     vE.push_back((*pJets)[i].energy());
     componentPtrs_.push_back(pJets->ptrAt(i));
   }

   for(int i = 0; i < (int) (*pMuons).size() ; i++){
     if((*pMuons)[i].pt() <  _minMuonEt || fabs((*pMuons)[i].eta()) >  _maxMuonEta) continue; 
     vPx.push_back((*pMuons)[i].px());
     vPy.push_back((*pMuons)[i].py());
     vPz.push_back((*pMuons)[i].pz());
     vE.push_back((*pMuons)[i].energy());
     componentPtrs_.push_back(pMuons->ptrAt(i));
   }
  
   for(int i = 0; i < (int) (*pElectrons).size() ; i++){
     if((*pElectrons)[i].pt() <  _minElectronEt || fabs((*pElectrons)[i].eta()) >  _maxElectronEta) continue;  
     vPx.push_back((*pElectrons)[i].px());
     vPy.push_back((*pElectrons)[i].py());
     vPz.push_back((*pElectrons)[i].pz());
     vE.push_back((*pElectrons)[i].energy());
     componentPtrs_.push_back(pElectrons->ptrAt(i));
   } 

   for(int i = 0; i < (int) (*pPhotons).size() ; i++){
     if((*pPhotons)[i].pt() <  _minPhotonEt || fabs((*pPhotons)[i].eta()) >  _maxPhotonEta) continue;   
     vPx.push_back((*pPhotons)[i].px());
     vPy.push_back((*pPhotons)[i].py());
     vPz.push_back((*pPhotons)[i].pz());
     vE.push_back((*pPhotons)[i].energy());
     componentPtrs_.push_back(pPhotons->ptrAt(i));
   } 

   //aren't taus included in jets?
   for(int i = 0; i < (int) (*pTaus).size() ; i++){
     if((*pTaus)[i].pt() <  _minTauEt || fabs((*pTaus)[i].eta()) >  _maxTauEta) continue;   
     vPx.push_back((*pTaus)[i].px());
     vPy.push_back((*pTaus)[i].py());
     vPz.push_back((*pTaus)[i].pz());
     vE.push_back((*pTaus)[i].energy());
     componentPtrs_.push_back(pTaus->ptrAt(i));
   }  

   // create product
   std::auto_ptr< std::vector<Hemisphere> > hemispheres(new std::vector<Hemisphere>);;
   hemispheres->reserve(2);

  //calls HemiAlgorithm for seed method 3 (transv. inv. Mass) and association method 3 (Lund algo)
  HemisphereAlgo myHemi(vPx,vPy,vPz,vE,_seedMethod,_combinationMethod);

  //get Hemisphere Axis 
  vA1 = myHemi.getAxis1();
  vA2 = myHemi.getAxis2();

  reco::Particle::LorentzVector p1(vA1[0]*vA1[3],vA1[1]*vA1[3],vA1[2]*vA1[3],vA1[4]);
  hemispheres->push_back(Hemisphere(p1));

  reco::Particle::LorentzVector p2(vA2[0]*vA2[3],vA2[1]*vA2[3],vA2[2]*vA2[3],vA2[4]);
  hemispheres->push_back(Hemisphere(p2));
 
  //get information to which Hemisphere each object belongs
  vgroups = myHemi.getGrouping(); 

  for ( unsigned int i=0; i<vgroups.size(); ++i ) {
    if ( vgroups[i]==1 ) {
      (*hemispheres)[0].addDaughter(componentPtrs_[i]);
    }
    else {
      (*hemispheres)[1].addDaughter(componentPtrs_[i]);
    }
  }

//   std::auto_ptr<HemiAxis > hemiAxis1(new HemiAxis(vA1));
//   std::auto_ptr<HemiAxis > hemiAxis2(new HemiAxis(vA2));

 

  //  hemi1->push_back(vA1);
  // hemi2->push_back(vA2);

//    iEvent.put(hemiAxis1,"hemi1");
//    iEvent.put(hemiAxis2,"hemi2");
  iEvent.put(hemispheres);

  //clean up
//     delete myHemi;
    vPx.clear();
    vPy.clear();
    vPz.clear();
    vE.clear();
    vgroups.clear();
    componentPtrs_.clear();
}



// ------------ method called once each job just before starting event loop  ------------
void 
PATHemisphereProducer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATHemisphereProducer::endJob() {
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATHemisphereProducer);
