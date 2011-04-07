#ifndef TkDetLayers_TIDRingBuilder_h
#define TkDetLayers_TIDRingBuilder_h


#include "TIDRing.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"

/** A concrete builder for TIDRing 
 */

class TIDRingBuilder {  
 public:
  TIDRingBuilder(){};
  TIDRing* build(const GeometricDet* aTIDRing,
		 const TrackerGeometry* theGeomDetGeometry);

  
};


#endif 
