#ifndef RecoTracker_TkSeedingLayers_HitExtractorPIX_H
#define RecoTracker_TkSeedingLayers_HitExtractorPIX_H

#include "RecoTracker/TkSeedingLayers/interface/SeedingLayer.h"
#include "HitExtractor.h"

#include <string>
#include <vector>


namespace ctfseeding {
class HitExtractorPIX : public HitExtractor {
public:
  HitExtractorPIX( SeedingLayer::Side & side, int idLayer, const std::string & hitProducer);
  HitExtractorPIX( SeedingLayer::Side & side, int idLayer, const std::string & hitProducer,
                   double hitErrorRPhi, double hitErrorRZ);
  virtual ~HitExtractorPIX(){}
  virtual std::vector<SeedingHit> hits(const edm::Event& , const edm::EventSetup& ) const;
  virtual HitExtractorPIX * clone() const { return new HitExtractorPIX(*this); }

private:
  SeedingLayer::Side theSide;
  int theIdLayer;
  std::string theHitProducer; 

  bool   theUseHitErrors; 
  double theHitErrorRPhi;
  double theHitErrorRZ;
};
}
#endif
