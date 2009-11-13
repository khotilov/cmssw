// File: DataMixingModule.cc
// Description:  see DataMixingModule.h
// Author:  Mike Hildreth, University of Notre Dame
//
//--------------------------------------------

#include <map>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Framework/interface/ConstProductRegistry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/Provenance/interface/BranchDescription.h"
//
//
#include "DataMixingModule.h"


using namespace std;

namespace edm
{

  // Constructor 
  DataMixingModule::DataMixingModule(const edm::ParameterSet& ps) : BMixingModule(ps),
							    label_(ps.getParameter<std::string>("Label"))

  {                                                         // what's "label_"?

    // get the subdetector names
    this->getSubdetectorNames();  //something like this may be useful to check what we are supposed to do...

    // create input selector
    if (label_.size()>0){
      sel_=new Selector( ModuleLabelSelector(label_));
    }
    else {
      sel_=new Selector( MatchAllSelector());
    }

    // For now, list all of them here.  Later, make this selectable with input parameters
    // 

    // Check to see if we are working in Full or Fast Simulation

    DoFastSim_ = (ps.getParameter<std::string>("IsThisFastSim")=="YES");

    // Put Fast Sim Sequences here for Simplification: Fewer options!

    if(DoFastSim_) {

    // declare the products to produce

      //Ecal:

      EBRecHitCollectionDM_        = ps.getParameter<std::string>("EBRecHitCollectionDM");
      EERecHitCollectionDM_        = ps.getParameter<std::string>("EERecHitCollectionDM");
      ESRecHitCollectionDM_        = ps.getParameter<std::string>("ESRecHitCollectionDM");

      produces< EBRecHitCollection >(EBRecHitCollectionDM_);
      produces< EERecHitCollection >(EERecHitCollectionDM_);
      produces< ESRecHitCollection >(ESRecHitCollectionDM_);

      EMWorker_ = new DataMixingEMWorker(ps);

      //Hcal:

      HBHERecHitCollectionDM_ = ps.getParameter<std::string>("HBHERecHitCollectionDM");
      HORecHitCollectionDM_   = ps.getParameter<std::string>("HORecHitCollectionDM");
      HFRecHitCollectionDM_   = ps.getParameter<std::string>("HFRecHitCollectionDM");
      ZDCRecHitCollectionDM_  = ps.getParameter<std::string>("ZDCRecHitCollectionDM");

      produces< HBHERecHitCollection >(HBHERecHitCollectionDM_);
      produces< HORecHitCollection >(HORecHitCollectionDM_);
      produces< HFRecHitCollection >(HFRecHitCollectionDM_);
      produces< ZDCRecHitCollection >(ZDCRecHitCollectionDM_);

      HcalWorker_ = new DataMixingHcalWorker(ps);

      //Muons:

      DTDigiCollectionDM_  = ps.getParameter<std::string>("DTDigiCollectionDM");
      RPCDigiCollectionDM_ = ps.getParameter<std::string>("RPCDigiCollectionDM");
      CSCStripDigiCollectionDM_ = ps.getParameter<std::string>("CSCStripDigiCollectionDM");
      CSCWireDigiCollectionDM_  = ps.getParameter<std::string>("CSCWireDigiCollectionDM");
      CSCComparatorDigiCollectionDM_  = ps.getParameter<std::string>("CSCComparatorDigiCollectionDM");

      produces< DTDigiCollection >();
      produces< RPCDigiCollection >();
      produces< CSCStripDigiCollection >(CSCStripDigiCollectionDM_);
      produces< CSCWireDigiCollection >(CSCWireDigiCollectionDM_);
      produces< CSCComparatorDigiCollection >(CSCComparatorDigiCollectionDM_);

      MuonWorker_ = new DataMixingMuonWorker(ps);

      //Tracks:

      produces< reco::TrackCollection >();
      GeneralTrackWorker_ = new DataMixingGeneralTrackWorker(ps);

    }
    else{  // Full Simulation options

    // declare the products to produce
    // Start with EM

    MergeEMDigis_ = (ps.getParameter<std::string>("EcalMergeType")=="Digis");

    if(MergeEMDigis_) {
      EBDigiCollectionDM_        = ps.getParameter<std::string>("EBDigiCollectionDM");
      EEDigiCollectionDM_        = ps.getParameter<std::string>("EEDigiCollectionDM");
      ESDigiCollectionDM_        = ps.getParameter<std::string>("ESDigiCollectionDM");
      //   nMaxPrintout_            = ps.getUntrackedParameter<int>("nMaxPrintout",10);

      produces< EBDigiCollection >(EBDigiCollectionDM_);
      produces< EEDigiCollection >(EEDigiCollectionDM_);
      produces< ESDigiCollection >(ESDigiCollectionDM_);

      EMDigiWorker_ = new DataMixingEMDigiWorker(ps);
    }
    else { // merge RecHits 
      EBRecHitCollectionDM_        = ps.getParameter<std::string>("EBRecHitCollectionDM");
      EERecHitCollectionDM_        = ps.getParameter<std::string>("EERecHitCollectionDM");
      ESRecHitCollectionDM_        = ps.getParameter<std::string>("ESRecHitCollectionDM");
      //   nMaxPrintout_            = ps.getUntrackedParameter<int>("nMaxPrintout",10);

      produces< EBRecHitCollection >(EBRecHitCollectionDM_);
      produces< EERecHitCollection >(EERecHitCollectionDM_);
      produces< ESRecHitCollection >(ESRecHitCollectionDM_);

      EMWorker_ = new DataMixingEMWorker(ps);
    }
    // Hcal next

    MergeHcalDigis_ = (ps.getParameter<std::string>("HcalMergeType")=="Digis");

    if(MergeHcalDigis_){
      HBHEDigiCollectionDM_ = ps.getParameter<std::string>("HBHEDigiCollectionDM");
      HODigiCollectionDM_   = ps.getParameter<std::string>("HODigiCollectionDM");
      HFDigiCollectionDM_   = ps.getParameter<std::string>("HFDigiCollectionDM");
      ZDCDigiCollectionDM_  = ps.getParameter<std::string>("ZDCDigiCollectionDM");

      produces< HBHEDigiCollection >();
      produces< HODigiCollection >();
      produces< HFDigiCollection >();
      produces< ZDCDigiCollection >();

      MergeHcalDigisProd_ = (ps.getParameter<std::string>("HcalDigiMerge")=="FullProd");

      if(MergeHcalDigisProd_) {
	HcalDigiWorkerProd_ = new DataMixingHcalDigiWorkerProd(ps);
      }
      else {HcalDigiWorker_ = new DataMixingHcalDigiWorker(ps);
      }


    }
    else{
      HBHERecHitCollectionDM_ = ps.getParameter<std::string>("HBHERecHitCollectionDM");
      HORecHitCollectionDM_   = ps.getParameter<std::string>("HORecHitCollectionDM");
      HFRecHitCollectionDM_   = ps.getParameter<std::string>("HFRecHitCollectionDM");
      ZDCRecHitCollectionDM_  = ps.getParameter<std::string>("ZDCRecHitCollectionDM");

      produces< HBHERecHitCollection >(HBHERecHitCollectionDM_);
      produces< HORecHitCollection >(HORecHitCollectionDM_);
      produces< HFRecHitCollection >(HFRecHitCollectionDM_);
      produces< ZDCRecHitCollection >(ZDCRecHitCollectionDM_);

      HcalWorker_ = new DataMixingHcalWorker(ps);
    }

    // Muons

    DTDigiCollectionDM_  = ps.getParameter<std::string>("DTDigiCollectionDM");
    RPCDigiCollectionDM_ = ps.getParameter<std::string>("RPCDigiCollectionDM");
    CSCStripDigiCollectionDM_ = ps.getParameter<std::string>("CSCStripDigiCollectionDM");
    CSCWireDigiCollectionDM_  = ps.getParameter<std::string>("CSCWireDigiCollectionDM");
    CSCComparatorDigiCollectionDM_  = ps.getParameter<std::string>("CSCComparatorDigiCollectionDM");


    produces< DTDigiCollection >();
    produces< RPCDigiCollection >();
    produces< CSCStripDigiCollection >(CSCStripDigiCollectionDM_);
    produces< CSCWireDigiCollection >(CSCWireDigiCollectionDM_);
    produces< CSCComparatorDigiCollection >(CSCComparatorDigiCollectionDM_);

    MuonWorker_ = new DataMixingMuonWorker(ps);

    // Si-Strips

    SiStripDigiCollectionDM_  = ps.getParameter<std::string>("SiStripDigiCollectionDM");


    produces< edm::DetSetVector<SiStripDigi> > (SiStripDigiCollectionDM_);
    
    SiStripWorker_ = new DataMixingSiStripWorker(ps);

    // Pixels

    PixelDigiCollectionDM_  = ps.getParameter<std::string>("PixelDigiCollectionDM");

    produces< edm::DetSetVector<PixelDigi> > (PixelDigiCollectionDM_);

    SiPixelWorker_ = new DataMixingSiPixelWorker(ps);

    }

  }

  void DataMixingModule::getSubdetectorNames() {
    // get subdetector names
    // edm::Service<edm::ConstProductRegistry> reg;
    // Loop over provenance of products in registry.
    //for (edm::ProductRegistry::ProductList::const_iterator it = reg->productList().begin(); it != reg->productList().end(); ++it) {

      //  **** Check this out.... ****

      // See FWCore/Framework/interface/BranchDescription.h
      // BranchDescription contains all the information for the product.

      // This section not very backwards-compatible in terms of digi-merging.  Need to be able to specify here which data format
      // to look at...

      //      edm::BranchDescription desc = it->second;
      //if (!desc.friendlyClassName_.compare(0,9,"EBRecHitC")) {
      //	Subdetectors_.push_back(desc.productInstanceName_);
      //LogInfo("DataMixingModule") <<"Adding container "<<desc.productInstanceName_ <<" for pileup treatment";
      //}
      //else if (!desc.friendlyClassName_.compare(0,9,"EERecHitC")) {
	//      else if (!desc.friendlyClassName_.compare(0,9,"EErechitC") && desc.productInstanceName_.compare(0,11,"TrackerHits")) {
      //	Subdetectors_.push_back(desc.productInstanceName_);
      //LogInfo("DataMixingModule") <<"Adding container "<<desc.productInstanceName_ <<" for pileup treatment";
      //}
      //else if (!desc.friendlyClassName_.compare(0,9,"HBRecHitC")) {
      //	Subdetectors_.push_back(desc.productInstanceName_);
      //LogInfo("DataMixingModule") <<"Adding container "<<desc.productInstanceName_ <<" for pileup treatment";
      //}
      //else if (!desc.friendlyClassName_.compare(0,9,"HERecHitC")) {
      //	Subdetectors_.push_back(desc.productInstanceName_);
      //LogInfo("DataMixingModule") <<"Adding container "<<desc.productInstanceName_ <<" for pileup treatment";
      // }
	// and so on with other detector types...
    // }
  }       
	       

  void DataMixingModule::createnewEDProduct() {
  }
 

  // Virtual destructor needed.
  DataMixingModule::~DataMixingModule() { 
    delete sel_;
    if(MergeEMDigis_){ delete EMDigiWorker_;}
    else {delete EMWorker_;}
    if(MergeHcalDigis_) { 
      if(MergeHcalDigisProd_) { delete HcalDigiWorkerProd_;}
      else { delete HcalDigiWorker_; }}
    else {delete HcalWorker_;}
    delete MuonWorker_;
    delete SiStripWorker_;
    delete SiPixelWorker_;
  }  



  void DataMixingModule::addSignals(const edm::Event &e, const edm::EventSetup& ES) { 
    // fill in maps of hits

    LogDebug("DataMixingModule")<<"===============> adding MC signals for "<<e.id();

    // Ecal
    if(MergeEMDigis_) { EMDigiWorker_->addEMSignals(e, ES); }
    else{ EMWorker_->addEMSignals(e);}

    // Hcal
    if(MergeHcalDigis_) { 
      if(MergeHcalDigisProd_){
	HcalDigiWorkerProd_->addHcalSignals(e, ES);
      }
      else{
	HcalDigiWorker_->addHcalSignals(e, ES);
      }
    }
    else {HcalWorker_->addHcalSignals(e);}
    
    // Muon
    MuonWorker_->addMuonSignals(e);

    // SiStrips
    SiStripWorker_->addSiStripSignals(e);

    // SiPixels
    SiPixelWorker_->addSiPixelSignals(e);
    
  } // end of addSignals

  
  void DataMixingModule::checkSignal(const edm::Event &e){}

  void DataMixingModule::addPileups(const int bcr, EventPrincipal *ep, unsigned int eventNr, unsigned int worker, const edm::EventSetup& ES) {  


    LogDebug("DataMixingModule") <<"\n===============> adding pileups from event  "<<ep->id()<<" for bunchcrossing "<<bcr;

    // fill in maps of hits; same code as addSignals, except now applied to the pileup events

    // Ecal
    if(MergeEMDigis_) {    EMDigiWorker_->addEMPileups(bcr, ep, eventNr, ES);}
    else {EMWorker_->addEMPileups(bcr, ep, eventNr); }

    // Hcal
    if(MergeHcalDigis_) {    
      if(MergeHcalDigisProd_) {    
	HcalDigiWorkerProd_->addHcalPileups(bcr, ep, eventNr, ES);
      }
      else{
	HcalDigiWorker_->addHcalPileups(bcr, ep, eventNr, ES);}
    }
    else {HcalWorker_->addHcalPileups(bcr, ep, eventNr);}

    // Muon
    MuonWorker_->addMuonPileups(bcr, ep, eventNr);

    // SiStrips
    SiStripWorker_->addSiStripPileups(bcr, ep, eventNr);

    // SiPixels
    SiPixelWorker_->addSiPixelPileups(bcr, ep, eventNr);

  }


  void DataMixingModule::doPileUp(edm::Event &e, const edm::EventSetup& ES)
  {//                                                                                       
                                       
    for (int bunchCrossing=minBunch_;bunchCrossing<=maxBunch_;++bunchCrossing) {
      setBcrOffset();
      for (unsigned int isource=0;isource<maxNbSources_;++isource) {
        setSourceOffset(isource);
        if (doit_[isource]) {
          merge(bunchCrossing, (pileup_[isource])[bunchCrossing-minBunch_],1, ES);
        }
      }
    }
  }


  void DataMixingModule::put(edm::Event &e,const edm::EventSetup& ES) {

    // individual workers...

    // Ecal
    if(MergeEMDigis_) {EMDigiWorker_->putEM(e,ES);}
    else {EMWorker_->putEM(e);}

    // Hcal
    if(MergeHcalDigis_) {
      if(MergeHcalDigisProd_) {
	HcalDigiWorkerProd_->putHcal(e,ES);
      }
      else{
	HcalDigiWorker_->putHcal(e,ES);
      }
    }
    else {HcalWorker_->putHcal(e);}

    // Muon
    MuonWorker_->putMuon(e);

    // SiStrips
    SiStripWorker_->putSiStrip(e);

    // SiPixels
    SiPixelWorker_->putSiPixel(e);

  }

  void DataMixingModule::setBcrOffset() {
  }

  void DataMixingModule::setSourceOffset(const unsigned int is) {
  }

} //edm
