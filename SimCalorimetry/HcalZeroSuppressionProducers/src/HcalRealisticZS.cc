#include "SimCalorimetry/HcalZeroSuppressionProducers/src/HcalRealisticZS.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/Common/interface/EDCollection.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Selector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
using namespace std;

#include <iostream>

HcalRealisticZS::HcalRealisticZS(edm::ParameterSet const& conf):
  inputLabel_(conf.getParameter<edm::InputTag>("digiLabel"))
{
  int mode=conf.getParameter<int>("mode");
  HcalZeroSuppressionAlgo::ZSMode zmode;
  switch (mode) {
  case(0): zmode=HcalZeroSuppressionAlgo::zs_SingleChannel; break;
  case(1): zmode=HcalZeroSuppressionAlgo::zs_TriggerTowerOR; break;
  case(2): zmode=HcalZeroSuppressionAlgo::zs_AllDepthsOR; break;
  default:
    edm::LogWarning("Hcal") << "Unknown zero suppression mode " << mode << " for HBHE. Using single-channel mode.";
    zmode=HcalZeroSuppressionAlgo::zs_SingleChannel; 
  }

  algo_=std::auto_ptr<HcalZSAlgoRealistic>(new HcalZSAlgoRealistic(zmode,
							   conf.getParameter<int>("HBlevel"),
							   conf.getParameter<int>("HElevel"),
							   conf.getParameter<int>("HOlevel"),
								   conf.getParameter<int>("HFlevel")));

  produces<HBHEDigiCollection>();
  produces<HODigiCollection>();
  produces<HFDigiCollection>();
    
}
    
HcalRealisticZS::~HcalRealisticZS() {
}
    
void HcalRealisticZS::produce(edm::Event& e, const edm::EventSetup& eventSetup)
{
 
  edm::Handle<HBHEDigiCollection> hbhe;    
  edm::Handle<HODigiCollection> ho;    
  edm::Handle<HFDigiCollection> hf;    
  
  e.getByLabel(inputLabel_,hbhe);
  
  // create empty output
  std::auto_ptr<HBHEDigiCollection> zs_hbhe(new HBHEDigiCollection);
  
  e.getByLabel(inputLabel_,ho);
  
  // create empty output
  std::auto_ptr<HODigiCollection> zs_ho(new HODigiCollection);
  
  e.getByLabel(inputLabel_,hf);
  
  // create empty output
  std::auto_ptr<HFDigiCollection> zs_hf(new HFDigiCollection);
  
  // run the algorithm
  algo_->suppress(*(hbhe.product()),*zs_hbhe);
  algo_->suppress(*(ho.product()),*zs_ho);
  algo_->suppress(*(hf.product()),*zs_hf);
  
  edm::LogInfo("HcalZeroSuppression") << "Suppression (HBHE) input " << hbhe->size() << " digis, output " << zs_hbhe->size() << " digis" 
				      <<  " (HO) input " << ho->size() << " digis, output " << zs_ho->size() << " digis"
				      <<  " (HF) input " << hf->size() << " digis, output " << zs_hf->size() << " digis";
  

    // return result
    e.put(zs_hbhe);
    e.put(zs_ho);
    e.put(zs_hf);

}
