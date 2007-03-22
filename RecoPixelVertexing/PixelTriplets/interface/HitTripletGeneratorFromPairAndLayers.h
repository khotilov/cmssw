#ifndef HitTripletGeneratorFromPairAndLayers_H
#define HitTripletGeneratorFromPairAndLayers_H

/** A HitTripletGenerator from HitPairGenerator and vector of
    Layers. The HitPairGenerator provides a set of hit pairs.
    For each pair the search for compatible hit(s) is done among
    provided Layers
 */

#include "RecoPixelVertexing/PixelTriplets/interface/HitTripletGenerator.h"
#include "RecoTracker/TkHitPairs/interface/HitPairGenerator.h"
#include <vector>
#include "RecoPixelVertexing/PixelTriplets/interface/CombinedHitTripletGenerator.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayer.h"

class HitTripletGeneratorFromPairAndLayers : public HitTripletGenerator {

public:
  typedef CombinedHitTripletGenerator::LayerCacheType       LayerCacheType;
  virtual ~HitTripletGeneratorFromPairAndLayers() {}
  virtual void init( const HitPairGenerator & pairs, 
    const std::vector<ctfseeding::SeedingLayer>& layers, LayerCacheType* layerCache) = 0; 
};
#endif


