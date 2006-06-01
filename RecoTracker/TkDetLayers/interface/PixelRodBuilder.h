#ifndef TkDetLayers_PixelRodBuilder_h
#define TkDetLayers_PixelRodBuilder_h


#include "RecoTracker/TkDetLayers/interface/PixelRod.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"

/** A concrete builder for PixelRod 
 */

class PixelRodBuilder {  
 public:
  PixelRodBuilder(){};
  PixelRod* build(const GeometricDet* aRod,
		  const TrackerGeometry* theGeomDetGeometry);
  
};


#endif 
