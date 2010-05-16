#include "DQMOffline/JetMET/interface/JetMETDQMDCSFilter.h" 
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

#include <iostream>
 
//
// -- Constructor
//
JetMETDQMDCSFilter::JetMETDQMDCSFilter( const edm::ParameterSet & pset ) {
   verbose_       = pset.getUntrackedParameter<bool>( "DebugOn", false );
   detectorTypes_ = pset.getUntrackedParameter<std::string>( "DetectorTypes", "ecal:hcal");
   filter_        = pset.getUntrackedParameter<bool>( "Filter", true );
   detectorOn_    = false;
   if (verbose_) std::cout << "JetMETDQMDCSFilter constructor: " << detectorTypes_ << std::endl;
}
//
// -- Destructor
//
JetMETDQMDCSFilter::~JetMETDQMDCSFilter() {
   if (verbose_) std::cout << "JetMETDQMDCSFilter destructor: " << std::endl;
}
 
bool JetMETDQMDCSFilter::filter(const edm::Event & evt, const edm::EventSetup & es) {
  
  detectorOn_ = true;

  if (!evt.isRealData()) return detectorOn_;
  if (!filter_) return detectorOn_;

  edm::Handle<DcsStatusCollection> dcsStatus;
  evt.getByLabel("scalersRawToDigi", dcsStatus);

  if (dcsStatus.isValid() && dcsStatus->size() != 0) {

  if (detectorTypes_.find("pixel") !=std::string::npos) {
      if ((*dcsStatus)[0].ready(DcsStatus::BPIX) && 
	  (*dcsStatus)[0].ready(DcsStatus::FPIX)) {
	if (verbose_) std::cout << "pixel on" << std::endl;
      } else detectorOn_ = false;
    }

    if (detectorTypes_.find("sistrip") !=std::string::npos){  
      if ((*dcsStatus)[0].ready(DcsStatus::TIBTID) &&
	  (*dcsStatus)[0].ready(DcsStatus::TOB) &&   
	  (*dcsStatus)[0].ready(DcsStatus::TECp) &&  
	  (*dcsStatus)[0].ready(DcsStatus::TECm)) {
	if (verbose_) std::cout << "sistrip on" << std::endl;
      } else detectorOn_ = false;
    }

    if (detectorTypes_.find("ecal") !=std::string::npos){  
      if ((*dcsStatus)[0].ready(DcsStatus::EBp) &&
	  (*dcsStatus)[0].ready(DcsStatus::EBm) &&   
	  (*dcsStatus)[0].ready(DcsStatus::EEp) &&  
	  (*dcsStatus)[0].ready(DcsStatus::EEm)) {
	if (verbose_) std::cout << "ecal on" << std::endl;
      } else detectorOn_ = false;
    }

    if (detectorTypes_.find("hbhe") !=std::string::npos){  
      if ((*dcsStatus)[0].ready(DcsStatus::HBHEa) &&
	  (*dcsStatus)[0].ready(DcsStatus::HBHEb) &&   
	  (*dcsStatus)[0].ready(DcsStatus::HBHEc)){
	if (verbose_) std::cout << "hbhe on" << std::endl;
      } else detectorOn_ = false;
    }

    if (detectorTypes_.find("hf") !=std::string::npos){  
      if ((*dcsStatus)[0].ready(DcsStatus::HF)){
	if (verbose_) std::cout << "hf on" << std::endl;
      } else detectorOn_ = false;
    }

    if (detectorTypes_.find("ho") !=std::string::npos){  
      if ((*dcsStatus)[0].ready(DcsStatus::HO)){
	if (verbose_) std::cout << "ho on" << std::endl;
      } else detectorOn_ = false;
    }

    if (detectorTypes_.find("es") !=std::string::npos){  
      if ((*dcsStatus)[0].ready(DcsStatus::ESp) &&
	  (*dcsStatus)[0].ready(DcsStatus::ESm)) {
	if (verbose_) std::cout << "es on" << std::endl;
      } else detectorOn_ = false;
    }

    if (detectorTypes_.find("muon") !=std::string::npos){  
      if ((*dcsStatus)[0].ready(DcsStatus::RPC)  &&
	  (*dcsStatus)[0].ready(DcsStatus::DT0)  &&   
	  (*dcsStatus)[0].ready(DcsStatus::DTp)  &&  
	  (*dcsStatus)[0].ready(DcsStatus::DTm)  &&
	  (*dcsStatus)[0].ready(DcsStatus::CSCp) &&  
	  (*dcsStatus)[0].ready(DcsStatus::CSCm)) {
	if (verbose_) std::cout << "muon on" << std::endl;
      } else detectorOn_ = false;
    }    

  }

  return detectorOn_;

}

//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(JetMETDQMDCSFilter);
