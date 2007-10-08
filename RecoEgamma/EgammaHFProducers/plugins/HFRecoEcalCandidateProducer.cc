/** \class HFRecoEcalCandidateProducers
 *
 *  \author Kevin Klapoetke (Minnesota)
 *
 * $Id:
 *
 */

#include <iostream>
#include <vector>
#include <memory>

// Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
//
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"

#include "RecoEgamma/EgammaHFProducers/plugins/HFRecoEcalCandidateProducer.h"


HFRecoEcalCandidateProducer::HFRecoEcalCandidateProducer(edm::ParameterSet const& conf):
  hfClusterShapes_(conf.getUntrackedParameter<std::string>("hfClusterShapes")),
  algo_(conf.getParameter<bool>("Correct"),conf.getParameter<double>("e9e25Cut"),conf.getParameter<double>("eCOREe9Cut"),conf.getParameter<double>("eSeLCut")) {

  produces<reco::RecoEcalCandidateCollection>();

} 

void HFRecoEcalCandidateProducer::produce(edm::Event & e, edm::EventSetup const& iSetup) {  
  
  //edm::Handle<reco::HFEMClusterShapeCollection> hf_clus;
  edm::Handle<reco::SuperClusterCollection> super_clus;
  edm::Handle<reco::HFEMClusterShapeAssociationCollection> hf_assoc;
  //e.getByLabel(hf_clus);
  e.getByLabel(hfClusterShapes_,super_clus);
  e.getByLabel(hfClusterShapes_,hf_assoc);

  
  
  // create return data
  std::auto_ptr<reco::RecoEcalCandidateCollection> retdata1(new reco::RecoEcalCandidateCollection());

  
  algo_.produce(super_clus,*hf_assoc,*retdata1);
 
  e.put(retdata1);

}












