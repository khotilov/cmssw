// -*- C++ -*-
//
// Package:    PFMuonCleaner
// Class:      PFMuonCleaner
// 
/**\class PFMuonCleaner PFMuonCleaner.cc MuonAnalysis/PFMuonCleaner/src/PFMuonCleaner.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michail Bachtis,32 3-B16,+41227678176,
//         Created:  Thu Mar 28 17:35:33 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"

//
// class declaration
//

class PFMuonUpdater : public edm::EDProducer {
   public:
      explicit PFMuonUpdater(const edm::ParameterSet&);
      ~PFMuonUpdater();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
  edm::InputTag pf_;
  edm::InputTag muons_;


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
PFMuonUpdater::PFMuonUpdater(const edm::ParameterSet& iConfig):
  pf_(iConfig.getParameter<edm::InputTag>("pfMuons")),
  muons_(iConfig.getParameter<edm::InputTag>("muons"))
{
   //register your products
  produces<reco::MuonCollection>();
}


PFMuonUpdater::~PFMuonUpdater()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PFMuonUpdater::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::MuonCollection> muonH;
   iEvent.getByLabel(muons_,muonH);
   
   Handle<reco::PFCandidateCollection> pfH;
   iEvent.getByLabel(pf_,pfH);
   
   std::auto_ptr<reco::MuonCollection> out(new reco::MuonCollection);
   

   //loop on the muons
   for( unsigned int i=0;i<muonH->size();++i) {
     //check if it is a PFMuon.If it is update the reco muon with the correct momentum
     reco::Muon mu = muonH->at(i);
     reco::MuonRef muonRef( muonH, i );

     for (unsigned int j=0;j<pfH->size();++j)
       if(pfH->at(j).muonRef() == muonRef) {
	 mu.setP4(pfH->at(j).p4());
	 break;
       }
     out->push_back(mu);
   }

   iEvent.put(out);
}

// ------------ method called once each job just before starting event loop  ------------
void 
PFMuonUpdater::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PFMuonUpdater::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
PFMuonUpdater::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PFMuonUpdater::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PFMuonUpdater::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PFMuonUpdater::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PFMuonUpdater::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFMuonUpdater);
