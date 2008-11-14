#ifndef SiPixelRawToDigi_H
#define SiPixelRawToDigi_H

/** \class SiPixelRawToDigi_H
 *  Plug-in module that performs Raw data to digi conversion 
 *  for pixel subdetector
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class SiPixelFedCablingTree;
class TH1D;
class TFile;
class R2DTimerObserver;


class SiPixelRawToDigi : public edm::EDProducer {
public:

  /// ctor
  explicit SiPixelRawToDigi( const edm::ParameterSet& );

  /// dtor
  virtual ~SiPixelRawToDigi();

  /// initialisation. Retrieves cabling map from EventSetup. 
  virtual void beginJob( const edm::EventSetup& );

  /// dummy end of job 
  virtual void endJob() {}

  /// get data, convert to digis attach againe to Event
  virtual void produce( edm::Event&, const edm::EventSetup& );

private:

  edm::ParameterSet config_;
  const SiPixelFedCablingTree * cablingTree_;
  TH1D *hCPU, *hDigi;
  TFile * rootFile;
  R2DTimerObserver * theTimer;
  bool includeErrors;
  bool checkOrder;
};
#endif
