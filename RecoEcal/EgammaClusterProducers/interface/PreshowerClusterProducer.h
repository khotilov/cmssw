#ifndef RecoEcal_EgammaClusterProducers_PreshowerClusterProducer_h
#define RecoEcal_EgammaClusterProducers_PreshowerClusterProducer_h

#include <memory>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "RecoEcal/EgammaClusterAlgos/interface/PreshowerClusterAlgo.h"
#include "CondFormats/ESObjects/interface/ESGain.h"
#include "CondFormats/ESObjects/interface/ESMIPToGeVConstant.h"
#include "CondFormats/ESObjects/interface/ESEEIntercalibConstants.h"

class PreshowerClusterProducer : public edm::EDProducer {

 public:

  typedef math::XYZPoint Point;

  explicit PreshowerClusterProducer (const edm::ParameterSet& ps);

  ~PreshowerClusterProducer();

  virtual void produce( edm::Event& evt, const edm::EventSetup& es);
  void set(const edm::EventSetup& es);

 private:

  int nEvt_;         // internal counter of events

  //clustering parameters:
  edm::InputTag preshHitProducer_;         // name of module/plugin/producer producing hits
  edm::InputTag endcapSClusterProducer_;   // ditto SuperClusters

  // name out output collections
  std::string preshClusterCollectionX_;  
  std::string preshClusterCollectionY_;  


  int preshNclust_;
  float preshClustECut;
  double etThresh_;

  // association parameters:
  std::string assocSClusterCollection_;    // name of super cluster output collection

  edm::ESHandle<ESGain> esgain_;
  edm::ESHandle<ESMIPToGeVConstant> esMIPToGeV_;
  edm::ESHandle<ESEEIntercalibConstants> esEEInterCalib_;
  double mip_;
  double gamma_;
  double alpha_;

  PreshowerClusterAlgo * presh_algo; // algorithm doing the real work
   // The set of used DetID's
  //std::set<DetId> used_strips;

  PreshowerClusterAlgo::DebugLevel debugL;  

};
#endif

