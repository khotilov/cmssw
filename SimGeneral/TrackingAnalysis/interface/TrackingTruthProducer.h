#ifndef TrackingAnalysis_TrackingTruthProducer_h
#define TrackingAnalysis_TrackingTruthProducer_h
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

class TrackingTruthProducer : public edm::EDProducer {

public:
  explicit TrackingTruthProducer( const edm::ParameterSet & );

private:
  void produce( edm::Event &, const edm::EventSetup & );

  edm::ParameterSet conf_;
  double distanceCut_;
  std::vector<std::string> dataLabels_;
  
};

#endif
