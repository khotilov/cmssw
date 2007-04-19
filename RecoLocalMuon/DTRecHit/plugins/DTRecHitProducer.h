#ifndef RecoLocalMuon_DTRecHitProducer_h
#define RecoLocalMuon_DTRecHitProducer_h

/** \class DTRecHitProducer
 *  Module for 1D DTRecHitPairs production. The concrete reconstruction algorithm
 *  is specified with the parameter "recAlgo" and must be configured with the
 *  "recAlgoConfig" parameter set.
 *
 *  $Date: 2006/03/14 13:06:15 $
 *  $Revision: 1.3 $
 *  \author G. Cerminara
 */

#include "FWCore/Framework/interface/EDProducer.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class DTRecHitBaseAlgo;

class DTRecHitProducer : public edm::EDProducer {
public:
  /// Constructor
  DTRecHitProducer(const edm::ParameterSet&);

  /// Destructor
  virtual ~DTRecHitProducer();

  /// The method which produces the rechits
  virtual void produce(edm::Event& event, const edm::EventSetup& setup);

private:
  // Switch on verbosity
  static bool debug;
  // The label to be used to retrieve DT digis from the event
  std::string theDTDigiLabel;
  // The reconstruction algorithm
  DTRecHitBaseAlgo *theAlgo;
//   static string theAlgoName;

};
#endif

