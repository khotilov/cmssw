#ifndef PhotonIDProducer_h
#define PhotonIDProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "RecoEgamma/PhotonIdentification/interface/CutBasedPhotonIDAlgo.h"


class PhotonIDProducer : public edm::EDProducer
{
 public:

  explicit PhotonIDProducer(const edm::ParameterSet& conf);
  virtual ~PhotonIDProducer();

  //  virtual void beginJob(edm::EventSetup const& iSetup);
  virtual void produce(edm::Event& e, const edm::EventSetup& c);
   
 private:

  CutBasedPhotonIDAlgo* cutBasedAlgo_; 	   

  edm::ParameterSet conf_;

  std::string photonProducer_;
  std::string photonLabel_;
  std::string photonIDLabel_;
  std::string photonIDAssociation_; //association map

  bool doCutBased_;

};

#endif
