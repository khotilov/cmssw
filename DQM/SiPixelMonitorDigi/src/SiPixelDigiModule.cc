#include "DQM/SiPixelMonitorDigi/interface/SiPixelDigiModule.h"
#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
/// Framework
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
// STL
#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <stdlib.h>
//
// Constructors
//
SiPixelDigiModule::SiPixelDigiModule() : id_(0) { }
///
SiPixelDigiModule::SiPixelDigiModule(const uint32_t& id) : 
  id_(id),
  ncols_(416),
  nrows_(160)
{ 
}
///
SiPixelDigiModule::SiPixelDigiModule(const uint32_t& id, const int& ncols, const int& nrows) : 
  id_(id),
  ncols_(ncols),
  nrows_(nrows)
{ 
}
//
// Destructor
//
SiPixelDigiModule::~SiPixelDigiModule() {}
//
// Book histograms
//
void SiPixelDigiModule::book() {

  DaqMonitorBEInterface* theDMBE = edm::Service<DaqMonitorBEInterface>().operator->();

  char hkey[80];  
  // Number of digis
  sprintf(hkey, "ndigis_module_%i",id_);
  meNDigis_ = theDMBE->book1D(hkey,"Number of Digis",50,0.,50.);
  meNDigis_->setAxisTitle("Number of digis",1);
  // Charge in ADC counts
  sprintf(hkey, "adc_module_%i",id_);
  meADC_ = theDMBE->book1D(hkey,"Digi charge",500,0.,500.);
  meADC_->setAxisTitle("ADC counts",1);

  // 2D hit map
  sprintf(hkey, "pixdigis_module_%i",id_);
  mePixDigis_ = theDMBE->book2D(hkey,"Digis per four pixels",208,0.,float(ncols_),80,0.,float(nrows_));
  mePixDigis_->setAxisTitle("Columns",1);
  mePixDigis_->setAxisTitle("Rows",2);
}
//
// Fill histograms
//
void SiPixelDigiModule::fill(const edm::DetSetVector<PixelDigi>& input) {
  
  edm::DetSetVector<PixelDigi>::const_iterator isearch = input.find(id_); // search  digis of detid
  
  if( isearch != input.end() ) {  // Not at empty iterator
    
    unsigned int numberOfDigis = 0;
    
    // Look at digis now
    edm::DetSet<PixelDigi>::const_iterator  di;
    for(di = isearch->data.begin(); di != isearch->data.end(); di++) {
      numberOfDigis++;
      int adc = di->adc();    // charge
      int col = di->column(); // column 
      int row = di->row();    // row
      (mePixDigis_)->Fill((float)col,(float)row);
      (meADC_)->Fill((float)adc);
    }
    (meNDigis_)->Fill((float)numberOfDigis);
    //std::cout<<"number of digis="<<numberOfDigis<<std::endl;
      
  }
  
  
  //std::cout<<"number of detector units="<<numberOfDetUnits<<std::endl;
  
}
