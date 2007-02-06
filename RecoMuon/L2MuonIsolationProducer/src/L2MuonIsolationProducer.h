#ifndef RecoMuon_L2MuonIsolationProducer_H
#define RecoMuon_L2MuonIsolationProducer_H

/**  \class L2MuonIsolationProducer
 * 
 *   L2 HLT muon isolation producer
 *
 *   \author  J.Alcaraz
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"

#include "RecoMuon/MuonIsolation/interface/Cuts.h"

#include "RecoMuon/MuonIsolation/interface/CaloExtractor.h"

class L2MuonIsolationProducer : public edm::EDProducer {

 public:

  /// constructor with config
  L2MuonIsolationProducer(const edm::ParameterSet&);
  
  /// destructor
  virtual ~L2MuonIsolationProducer(); 
  
  /// Produce isolation maps
  virtual void produce(edm::Event&, const edm::EventSetup&);
  // ex virtual void reconstruct();

 private:
  
  // Muon track Collection Label
  std::string theSACollectionLabel;

  // Isolation cuts
  muonisolation::Cuts theCuts;

  // Option to write MuIsoDeposits into the event
  double optOutputIsoDeposits;

  // MuIsoExtractor
  muonisolation::CaloExtractor theCalExtractor;

};

#endif
