// -*- C++ -*-
//
// Package:    CentralityProducer
// Class:      CentralityProducer
// 
/**\class CentralityProducer CentralityProducer.cc RecoHI/CentralityProducer/src/CentralityProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Yetkin Yilmaz, Young Soo Park
//         Created:  Wed Jun 11 15:31:41 CEST 2008
// $Id: CentralityProducer.cc,v 1.8 2009/08/19 11:45:05 yilmaz Exp $
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Ref.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

using namespace std;

//
// class declaration
//

class CentralityProducer : public edm::EDProducer {
   public:
      explicit CentralityProducer(const edm::ParameterSet&);
      ~CentralityProducer();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------

  bool recoLevel_;
  bool genLevel_;
  edm::InputTag   genLabel_;	

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
CentralityProducer::CentralityProducer(const edm::ParameterSet& iConfig)
{
   //register your products
  
  recoLevel_ = iConfig.getUntrackedParameter<bool>("recoLevel",true);
  genLevel_ = iConfig.getUntrackedParameter<bool>("genLevel",true);
  genLabel_ = iConfig.getParameter<edm::InputTag>("genLabel");

  if(recoLevel_){
    produces<reco::CentralityCollection>("recoBased");
  }
  if(genLevel_){
    produces<reco::CentralityCollection>("genBased");
  }

}


CentralityProducer::~CentralityProducer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CentralityProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
   using namespace HepMC;
  
   const GenEvent* evt;

   double eHF     =  0;        // Variable for computing total HF energy 
   double eCASTOR =  0;
   double eZDC    =  0;
   int cnt = 0;

   if(recoLevel_){
     Handle<HFRecHitCollection> hits;
     iEvent.getByLabel("hfreco",hits);
     for( size_t ihit = 0; ihit<hits->size(); ++ ihit){
       const HFRecHit & rechit = (*hits)[ ihit ];
       eHF = eHF + rechit.energy();
     }
     std::auto_ptr<CentralityCollection> centOutput(new CentralityCollection);
     Centrality creco(eHF,"HFTotalEnergy");
     centOutput->push_back(creco);
     iEvent.put(centOutput, "recoBased");

     /*
     To Do : 
     - Add other detectors into reconstruction.
     */

   }

   if(genLevel_){
     Handle<HepMCProduct> mc;
     iEvent.getByLabel(genLabel_,mc);   
     evt=mc->GetEvent();
     
     HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
     HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
     for(HepMC::GenEvent::particle_const_iterator it = begin;it!=end;++it)
       {
	 HepMC::GenParticle* p = *it;
	 //float status = p->status();
	 float pdg_id = p->pdg_id();
	 float eta = p->momentum().eta();
	 //float phi = p->momentum().phi();
	 //float pt = p->momentum().perp();
	 float e = p->momentum().e();
	 
	 
	 if(fabs(eta)>=3&&fabs(eta)<=5){
	   
	   eHF+=e;
	 }
	 if(fabs(eta)>=5.3&&fabs(eta)<=6.7){
	   
	   eCASTOR+=e;
	 }
	 if(fabs(eta)>=5.9){
	   if(pdg_id==2112){
	     
	     eZDC+=e;
	     cnt+=1;
	   }
	 }    
       }
     

     Centrality c1(eHF,"HFGenParticleEnergy");

     std::auto_ptr<CentralityCollection> centOutput(new CentralityCollection);
     centOutput->push_back(c1);

     iEvent.put(centOutput, "genBased");
     
   }
}

// ------------ method called once each job just before starting event loop  ------------
void 
CentralityProducer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CentralityProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CentralityProducer);
