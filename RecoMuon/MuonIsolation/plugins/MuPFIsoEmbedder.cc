// -*- C++ -*-
//
// Package:    MuPFIsoEmbedder
// Class:      MuPFIsoEmbedder
// 
/**\class MuPFIsoEmbedder MuPFIsoEmbedder.cc RecoMuon/MuPFIsoEmbedder/src/MuPFIsoEmbedder.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michail Bachtis,32 3-B16,+41227675567,
//         Created:  Thu Jun  9 01:36:17 CEST 2011
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
#include "RecoMuon/MuonIsolation/interface/MuPFIsoHelper.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

//
// class declaration
//

class MuPFIsoEmbedder : public edm::EDProducer {
   public:
      explicit MuPFIsoEmbedder(const edm::ParameterSet&);
      ~MuPFIsoEmbedder();

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
  edm::InputTag muons_;

  MuPFIsoHelper *helper_;

};




//
MuPFIsoEmbedder::MuPFIsoEmbedder(const edm::ParameterSet& iConfig):
  muons_(iConfig.getParameter<edm::InputTag>("src"))
{

  helper_ = new MuPFIsoHelper(iConfig);
  produces<reco::MuonCollection>();
}


MuPFIsoEmbedder::~MuPFIsoEmbedder()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuPFIsoEmbedder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;


   helper_->beginEvent(iEvent);
  
   edm::Handle<reco::MuonCollection > muons;
   iEvent.getByLabel(muons_,muons);


   std::auto_ptr<MuonCollection> out(new MuonCollection);

   for(unsigned int i=0;i<muons->size();++i) {
     MuonRef muonRef(muons,i);
     Muon muon = muons->at(i);
     helper_->embedPFIsolation(muon,muonRef);
     out->push_back(muon);
   }

   iEvent.put(out);
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuPFIsoEmbedder::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuPFIsoEmbedder::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
MuPFIsoEmbedder::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MuPFIsoEmbedder::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MuPFIsoEmbedder::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MuPFIsoEmbedder::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuPFIsoEmbedder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuPFIsoEmbedder);
