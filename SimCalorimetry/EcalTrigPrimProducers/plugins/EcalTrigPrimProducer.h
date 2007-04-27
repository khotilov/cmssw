#ifndef EcaltrigprimProducer_h
#define EcaltrigprimProducer_h
  
/** \class EcalTrigPrimProducer
 *
 * EcalTrigPrimProducer produces a EcalTrigPrimDigiCollection
 * Simulation as close as possible to hardware
 * Main algorithm is EcalTrigPrimFunctionalAlgo which is now
 * templated to take EBdataFrames/EEDataFrames as input
 *
 * \author Ursula Berthon, Stephanie Baffioni, Pascal Paganini,   LLR Palaiseau
 *
 * \version   1st Version may 2006
 * \version   2nd Version jul 2006
 * \version   3rd Version nov 2006
 * \version   4th Version apr 2007   full endcap
 *
 ************************************************************/

 
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
 
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
  
class TFile;
class TTree;
class EcalTrigPrimFunctionalAlgo;
 
class EcalTrigPrimProducer : public edm::EDProducer
{
 public:
  
  explicit EcalTrigPrimProducer(const edm::ParameterSet& conf);
  
  virtual ~EcalTrigPrimProducer();
  
  void beginJob(edm::EventSetup const& setup);

  virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
 private:
  EcalTrigPrimFunctionalAlgo *algo_;
  TFile *histfile_;
  TTree *valTree_;
  bool valid_;
  bool barrelOnly_;
  bool tcpFormat_;
  bool debug_;
  std::string label_;
  std::string instanceNameEB_;
  std::string instanceNameEE_;

  int binOfMaximum_;

   static const int nrSamples_; //nr samples to write, should not be changed, if not problems in EcalTriggerPrimitiveDigi class
  const edm::ParameterSet ps_;
};
  
#endif
 


