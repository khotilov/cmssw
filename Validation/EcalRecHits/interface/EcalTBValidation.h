#ifndef ECALTBVALIDATION_H
#define ECALTBVALIDATION_H

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

#include "TFile.h" 
#include "TH1.h"
#include "TH2.h"

class EcalTBValidation : public edm::EDAnalyzer {
 public:
  explicit EcalTBValidation( const edm::ParameterSet& );
  ~EcalTBValidation();
  
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  virtual void beginJob(edm::EventSetup const&);
  virtual void endJob();
  
 private:
  int data_;
  int xtalInBeam_;
  std::string rootfile_;
  std::string digiCollection_;
  std::string digiProducer_;
  std::string hitCollection_;
  std::string hitProducer_;
  std::string hodoRecInfoCollection_;
  std::string hodoRecInfoProducer_;
  std::string tdcRecInfoCollection_;
  std::string tdcRecInfoProducer_;
  std::string eventHeaderCollection_;
  std::string eventHeaderProducer_;
  
  // histos
  TH2F *h_xib,   *h_ampltdc, *h_Shape;
  TH1F *h_hodoX, *h_hodoY;
  TH1F *h_e1x1, *h_e3x3,  *h_e5x5;
  TH1F *h_e1e9, *h_e1e25, *h_e9e25;
  TH1F *h_e1x1_center, *h_e3x3_center,  *h_e5x5_center;
  TH2F *h_e1vsX,      *h_e1vsY;
  TH2F *h_e1e9vsX,    *h_e1e9vsY;
  TH2F *h_e1e25vsX,   *h_e1e25vsY;
  TH2F *h_e9e25vsX,   *h_e9e25vsY;
};



#endif
