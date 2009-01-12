#ifndef SiStripElectronSeedProducer_h
#define SiStripElectronSeedProducer_h
  
  
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
 
#include "DataFormats/Common/interface/EDProduct.h"
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
  
class SiStripElectronSeedGenerator;
 
class SiStripElectronSeedProducer : public edm::EDProducer
{
 public:
  
  explicit SiStripElectronSeedProducer(const edm::ParameterSet& conf);
  
  virtual ~SiStripElectronSeedProducer();
  
  virtual void beginJob(edm::EventSetup const&iSetup);
  virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
 private:
  edm::InputTag superClusters_[2];
  const edm::ParameterSet conf_;
  SiStripElectronSeedGenerator *matcher_;
  };
  
#endif
