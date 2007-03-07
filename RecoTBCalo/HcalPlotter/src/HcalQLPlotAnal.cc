// -*- C++ -*-
//
// Package:    HcalQLPlotAnal
// Class:      HcalQLPlotAnal
// 
/**\class HcalQLPlotAnal HcalQLPlotAnal.cc RecoTBCalo/HcalQLPlotAnal/src/HcalQLPlotAnal.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Phillip R. Dudero
//         Created:  Tue Jan 16 21:11:37 CST 2007
// $Id: HcalQLPlotAnal.cc,v 1.2 2007/02/22 15:44:12 dudero Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBTriggerData.h"
#include "RecoTBCalo/HcalPlotter/src/HcalQLPlotAnalAlgos.h"
#include <string>
//
// class declaration
//

class HcalQLPlotAnal : public edm::EDAnalyzer {
   public:
      explicit HcalQLPlotAnal(const edm::ParameterSet&);
      ~HcalQLPlotAnal();


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  edm::InputTag hbheRHLabel_,hoRHLabel_,hfRHLabel_;
  edm::InputTag hcalDigiLabel_, hcalTrigLabel_;
  HcalQLPlotAnalAlgos * algo_;

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
HcalQLPlotAnal::HcalQLPlotAnal(const edm::ParameterSet& iConfig) :
  hbheRHLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hbheRHtag")),
  hoRHLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hoRHtag")),
  hfRHLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hfRHtag")),
  hcalDigiLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hcalDigiTag")),
  hcalTrigLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hcalTrigTag"))
{
  algo_ = new
    HcalQLPlotAnalAlgos(iConfig.getUntrackedParameter<std::string>("outputFilename").c_str(),
			iConfig.getParameter<edm::ParameterSet>("HistoParameters"));
}


HcalQLPlotAnal::~HcalQLPlotAnal()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HcalQLPlotAnal::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Step A/C: Get Inputs and process (repeatedly)
  try {
    edm::Handle<HcalTBTriggerData> trig;
    iEvent.getByLabel(hcalTrigLabel_,trig);
    algo_->SetEventType(*trig);
  } catch (std::exception& e) { // can't find it!
    edm::LogError("HcalQLPlotAnal::analyze") << "No Trigger Data found, skip event";
    return;
  }

  try {
    edm::Handle<HBHEDigiCollection> hbhedg;    iEvent.getByLabel(hcalDigiLabel_,hbhedg);
    edm::Handle<HBHERecHitCollection> hbherh;  iEvent.getByLabel(hbheRHLabel_,hbherh);
    algo_->processDigi(*hbhedg);
    algo_->processRH(*hbherh,*hbhedg);
  } catch (std::exception& e) { // can't find it!
    edm::LogWarning("HcalQLPlotAnal::analyze") << "One of HBHE Digis/RecHits not found";
  }

  try {
    edm::Handle<HODigiCollection> hodg;    iEvent.getByLabel(hcalDigiLabel_,hodg);
    edm::Handle<HORecHitCollection> horh;  iEvent.getByLabel(hoRHLabel_,horh);
    algo_->processDigi(*hodg);
    algo_->processRH(*horh,*hodg);
  } catch (std::exception& e) { // can't find it!
    edm::LogWarning("HcalQLPlotAnal::analyze") << "One of HO Digis/RecHits not found";
  }

  try {
    edm::Handle<HFDigiCollection> hfdg;    iEvent.getByLabel(hcalDigiLabel_,hfdg);
    edm::Handle<HFRecHitCollection> hfrh;  iEvent.getByLabel(hfRHLabel_,hfrh);
    algo_->processDigi(*hfdg);
    algo_->processRH(*hfrh,*hfdg);
  } catch (std::exception& e) { // can't find it!
    edm::LogWarning("HcalQLPlotAnal::analyze") << "One of HF Digis/RecHits not found";
  }

}


// ------------ method called once each job just after ending the event loop  ------------
void 
HcalQLPlotAnal::endJob()
{
  algo_->end();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalQLPlotAnal);
