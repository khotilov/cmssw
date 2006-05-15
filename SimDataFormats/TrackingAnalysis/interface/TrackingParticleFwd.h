#ifndef TrackingAnalysis_TrackingParticle_h
#define TrackingAnalysis_TrackingParticle_h
#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefProd.h"

class TrackingParticle;
typedef std::vector<TrackingParticle> TrackingParticleCollection;
typedef edm::Ref<TrackingParticleCollection> TrackingParticleRef;
typedef edm::RefVector<TrackingParticleCollection> TrackingParticleRefVector;
typedef edm::RefProd<TrackingParticleCollection> TrackingParticleRefProd;

#endif

