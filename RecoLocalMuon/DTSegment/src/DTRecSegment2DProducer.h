#ifndef DTRecHit_DTRecSegment2DProducer_h
#define DTRecHit_DTRecSegment2DProducer_h

/** \class DTRecSegment2DProducer
 *
 * Producer for DT segment in one projection.
 *  
 * $Date: 2006/04/20 07:45:24 $
 * $Revision: 1.4 $
 * \author Stefano Lacaprara - INFN Legnaro <stefano.lacaprara@pd.infn.it>
 * \author Riccardo Bellan - INFN TO <riccardo.bellan@cern.ch>
 *
 */

/* Base Class Headers */
#include "FWCore/Framework/interface/EDProducer.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class DTRecSegment2DBaseAlgo;

/* C++ Headers */

/* ====================================================================== */

/* Class DTRecSegment2DProducer Interface */

class DTRecSegment2DProducer : public edm::EDProducer {

 public:

  /// Constructor
  DTRecSegment2DProducer(const edm::ParameterSet&) ;

  /// Destructor
  virtual ~DTRecSegment2DProducer() ;
    
  // Operations

  /// The method which produces the 2D-segments
  virtual void produce(edm::Event& event, const edm::EventSetup& setup);

 protected:

 private:
  // Switch on verbosity
  bool debug;

  // The 2D-segments reconstruction algorithm
  DTRecSegment2DBaseAlgo* theAlgo;

  //static std::string theAlgoName;
  std::string theRecHits1DLabel;
};
#endif // DTRecHit_DTRecSegment2DProducer_h

