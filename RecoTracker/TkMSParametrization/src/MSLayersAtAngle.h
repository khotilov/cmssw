#ifndef MSLayersAtAngle_H
#define MSLayersAtAngle_H

/**
 *
 */

#include <vector>
#include <cmath>

#include "RecoTracker/TkMSParametrization/interface/MSLayer.h"
class PixelRecoLineRZ;

class MSLayersAtAngle {

public:
  MSLayersAtAngle() { }
  MSLayersAtAngle(const vector<MSLayer> & layers);
  void update(const MSLayer & layer);
  const MSLayer * findLayer(const MSLayer & layer) const;

  float sumX0D(const PixelRecoPointRZ & pointI,
               const PixelRecoPointRZ & pointO) const;
  float sumX0D(const PixelRecoPointRZ & pointI,
               const PixelRecoPointRZ & pointM,
               const PixelRecoPointRZ & pointO) const;

  int size() const { return theLayers.size(); }
  void print() const;

private:
  vector<MSLayer> theLayers;

private:
  typedef vector<MSLayer>::const_iterator LayerItr;
  LayerItr findLayer(const PixelRecoPointRZ & point,
                     LayerItr i1, LayerItr i2) const;
  float sum2RmRn(LayerItr i1, LayerItr i2,
                 float rTarget, const PixelRecoLineRZ & line) const;
};

#endif
