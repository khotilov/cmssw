// -*- C++ -*-
//
// Package:    SiStripTools
// Class:      EventTimeDistribution
// 
/**\class EventTimeDistribution EventTimeDistribution.cc DPGAnalysis/SiStripTools/plugins/EventTimeDistribution.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Andrea Venturi
//         Created:  Tue Jul 19 11:56:00 CEST 2009
//
//


// system include files
#include <memory>

// user include files
#include <string>

#include "TH1F.h"
#include "TH2F.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DPGAnalysis/SiStripTools/interface/EventWithHistory.h"
#include "DPGAnalysis/SiStripTools/interface/APVCyclePhaseCollection.h"
//
// class decleration
//

class EventTimeDistribution : public edm::EDAnalyzer {
 public:
    explicit EventTimeDistribution(const edm::ParameterSet&);
    ~EventTimeDistribution();


   private:
      virtual void beginJob() ;
      virtual void beginRun(const edm::Run&, const edm::EventSetup&) ;
      virtual void endRun(const edm::Run&, const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

  const edm::InputTag _historyProduct;
  const edm::InputTag _apvphasecoll;
  const std::string _phasepart;
  unsigned int _nevents;
  const double _binsize;

  TH1F* _dbx;
  TH1F* _bxincycle;
  TH2F* _dbxvsbx;
  TH2F* _orbitvsbx;

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
EventTimeDistribution::EventTimeDistribution(const edm::ParameterSet& iConfig):
  _historyProduct(iConfig.getParameter<edm::InputTag>("historyProduct")),
  _apvphasecoll(iConfig.getParameter<edm::InputTag>("apvPhaseCollection")),
  _phasepart(iConfig.getUntrackedParameter<std::string>("phasePartition","None")),
  _nevents(0),
  _binsize(iConfig.getUntrackedParameter<double>("minBinSizeInSec",1.))
{
   //now do what ever initialization is needed

}


EventTimeDistribution::~EventTimeDistribution()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
EventTimeDistribution::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   _nevents++;

   edm::Handle<EventWithHistory> he;
   iEvent.getByLabel(_historyProduct,he);

   edm::Handle<APVCyclePhaseCollection> apvphase;
   iEvent.getByLabel(_apvphasecoll,apvphase);

   long long tbx = he->absoluteBX();
   if(apvphase.isValid() && !apvphase.failedToGet()) {
     const int thephase = apvphase->getPhase(_phasepart); 
     if(thephase!=APVCyclePhaseCollection::invalid &&
	thephase!=APVCyclePhaseCollection::multiphase &&
	thephase!=APVCyclePhaseCollection::nopartition)
       tbx -= thephase;
   }



   // improve the matchin between default and actual partitions
   
   _dbx->Fill(he->deltaBX());
   _bxincycle->Fill(tbx%70);
   _dbxvsbx->Fill(tbx%70,he->deltaBX());
   _orbitvsbx->Fill(tbx%70,iEvent.orbitNumber());


}

void 
EventTimeDistribution::beginRun(const edm::Run& iRun, const edm::EventSetup&)
{
  
  edm::Service<TFileService> tfserv;

  char dirname[300];
  sprintf(dirname,"run_%d",iRun.run());
  TFileDirectory subrun = tfserv->mkdir(dirname);

  _dbx = subrun.make<TH1F>("dbx","dbx",1000,-0.5,999.5);
  _bxincycle = subrun.make<TH1F>("bxcycle","bxcycle",70,-0.5,69.5);
  _dbxvsbx = subrun.make<TH2F>("dbxvsbx","dbxvsbx",70,-0.5,69.5,1000,-0.5,999.5);
  _orbitvsbx = subrun.make<TH2F>("orbitvsbx","orbitvsbx",70,-0.5,69.5,3600,0,11223*_binsize*3600);
  _orbitvsbx->SetBit(TH1::kCanRebin);

}

void 
EventTimeDistribution::endRun(const edm::Run& iRun, const edm::EventSetup&)
{
}


// ------------ method called once each job just before starting event loop  ------------
void 
EventTimeDistribution::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventTimeDistribution::endJob() {

  edm::LogInfo("EndOfJob") << _nevents << " analyzed events";

}


//define this as a plug-in
DEFINE_FWK_MODULE(EventTimeDistribution);
