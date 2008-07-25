// Package:    SiPixelMonitorTrack
// Class:      SiPixelTrackResidualModule
// 
// class SiPixelTrackResidualModule SiPixelTrackResidualModule.cc 
//       DQM/SiPixelMonitorTrack/src/SiPixelTrackResidualModule.cc
//
// Description: SiPixel hit-to-track residual data quality monitoring modules
// Implementation: prototype -> improved -> never final - end of the 1st step 
//
// Original Author: Shan-Huei Chuang
//         Created: Fri Mar 23 18:41:42 CET 2007
// $Id: SiPixelTrackResidualModule.cc,v 1.6 2008/03/01 20:19:52 lat Exp $


#include <string>
#include <iostream>

#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQM/SiPixelCommon/interface/SiPixelHistogramId.h"
#include "DQM/SiPixelMonitorTrack/interface/SiPixelTrackResidualModule.h"


using namespace std; 


SiPixelTrackResidualModule::SiPixelTrackResidualModule() : id_(0) {
}


SiPixelTrackResidualModule::SiPixelTrackResidualModule(uint32_t id) : id_(id) { 
}


SiPixelTrackResidualModule::~SiPixelTrackResidualModule() { 
}


void SiPixelTrackResidualModule::book(const edm::ParameterSet& iConfig) {
  DQMStore* dbe = edm::Service<DQMStore>().operator->();

  edm::InputTag src = iConfig.getParameter<edm::InputTag>("src");
  SiPixelHistogramId* theHistogramId = new SiPixelHistogramId(src.label());
  std::string hisID;

  hisID = theHistogramId->setHistoId("residualX",id_);
  meResidualX_ = dbe->book1D(hisID,"Hit-to-Track Residual in X",500,-5.,5.);
  meResidualX_->setAxisTitle("hit-to-track residual in x (cm)",1);

  hisID = theHistogramId->setHistoId("residualY",id_);
  meResidualY_ = dbe->book1D(hisID,"Hit-to-Track Residual in Y",500,-5.,5.);
  meResidualY_->setAxisTitle("hit-to-track residual in y (cm)",1);

  delete theHistogramId;
}


void SiPixelTrackResidualModule::fill(const Measurement2DVector& residual) {
  (meResidualX_)->Fill(residual.x()); 
  (meResidualY_)->Fill(residual.y()); 
}
