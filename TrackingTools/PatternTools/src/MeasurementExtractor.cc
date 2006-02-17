#include "TrackingTools/PatternTools/interface/MeasurementExtractor.h"

AlgebraicVector MeasurementExtractor::measuredParameters(const TransientTrackingRecHit& hit) {
  AlgebraicVector par5( theTSoS.localParameters().vector());
  AlgebraicMatrix H( hit.projectionMatrix());
  return H*par5;
}

AlgebraicSymMatrix MeasurementExtractor::measuredError(const TransientTrackingRecHit& hit) {
  
  AlgebraicSymMatrix err5( theTSoS.localError().matrix());
  AlgebraicMatrix H( hit.projectionMatrix());
  //  return AlgebraicSymMatrix( H * err5 * H.T());
  return err5.similarity(H);
}
