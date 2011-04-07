#ifndef TkDetLayers_TIDLayerBuilder_h
#define TkDetLayers_TIDLayerBuilder_h


#include "TIDLayer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"

/** A concrete builder for TIDLayer 
 */

class TIDLayerBuilder {  
 public:
  TIDLayerBuilder(){};
  TIDLayer* build(const GeometricDet* aTIDLayer,
		  const TrackerGeometry* theGeomDetGeometry);

   
};


#endif 
