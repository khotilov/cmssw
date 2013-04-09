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

class PFMuonCleaner : public edm::EDProducer {
   public:
      explicit PFMuonCleaner(const edm::ParameterSet&);
      ~PFMuonCleaner();

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
  PFMuonAlgo *muAlgo_;

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
PFMuonCleaner::PFMuonCleaner(const edm::ParameterSet& iConfig)
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  muAlgo_ = new PFMuonAlgo();
  muAlgo_->setParameters(iConfig);
  produces<reco::PFCandidateCollection>();

  
}


PFMuonCleaner::~PFMuonCleaner()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PFMuonCleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::MuonCollection> muonH;
   iEvent.getByLabel("muons",muonH);
   

   Handle<reco::PFCandidateCollection> pfH;
   iEvent.getByLabel("particleFlow",pfH);

   Handle<reco::VertexCollection> vertexH;
   iEvent.getByLabel("offlinePrimaryVertices",vertexH);
   
   std::auto_ptr<reco::PFCandidateCollection> out(new reco::PFCandidateCollection);

   muAlgo_->setInputsForCleaning(vertexH.product());
   
   //Loop on PF candidates
   for (unsigned int i =0 ;i<pfH->size();++i) 
     //if not PF muon just put in in the output collection
     if(abs(pfH->at(i).pdgId())!=13)
       out->push_back(pfH->at(i));
     else {
       //ok a muon we need to reset the momentum based on PFmuonAlgo+new TuneP
       //and then add it back
       reco::PFCandidate cand = pfH->at(i);
       if(muAlgo_->reconstructMuon(cand, cand.muonRef(),true)) {
	 out->push_back(cand);
	 
       }else {
	 printf("Rejected Muon with  pt=%f\n",cand.pt());

       }
     }

   //OK now perform the PFCleaning
   muAlgo_->postClean(out.get());
   
   //And finally add missing candidates
   muAlgo_->addMissingMuons(muonH,out.get());
   iEvent.put(out);
}

// ------------ method called once each job just before starting event loop  ------------
void 
PFMuonCleaner::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PFMuonCleaner::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
PFMuonCleaner::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PFMuonCleaner::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PFMuonCleaner::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PFMuonCleaner::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PFMuonCleaner::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFMuonCleaner);
