#ifndef FastSimulation_CaloRecHitsProducer_H
#define FastSimulation_CaloRecHitsProducer_H

//  F. Beaudette (LLR). Florian.Beaudette@cern.ch
//  Created 20/07/06 
//  The CaloRecHits producer.


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/EDProduct.h"


#include <string>

class HcalRecHitsMaker;
class EcalBarrelRecHitsMaker;
class EcalEndcapRecHitsMaker;
class EcalPreshowerRecHitsMaker;
class RandomEngine;
class CaloGeometryHelper;

class CaloRecHitsProducer : public edm::EDProducer
{

 public:

  explicit CaloRecHitsProducer(edm::ParameterSet const & p);
  virtual ~CaloRecHitsProducer();
  virtual void beginJob(const edm::EventSetup & c);
  virtual void endJob();
  virtual void produce(edm::Event & e, const edm::EventSetup & c);

 private:
  
  HcalRecHitsMaker * HcalRecHitsMaker_;
  EcalBarrelRecHitsMaker * EcalBarrelRecHitsMaker_;
  EcalEndcapRecHitsMaker * EcalEndcapRecHitsMaker_;
  EcalPreshowerRecHitsMaker * EcalPreshowerRecHitsMaker_;
  std::string EBrechitCollection_;
  std::string EErechitCollection_;
  std::string ESrechitCollection_;

   // The random engine
  RandomEngine* random;
  
  CaloGeometryHelper* myCaloGeometryHelper_ ;
};

#endif
