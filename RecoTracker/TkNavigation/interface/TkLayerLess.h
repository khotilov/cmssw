#ifndef TkNavigation_TkLayerLess_H
#define TkNavigation_TkLayerLess_H

#include "TrackingTools/GeomPropagators/interface/PropagationDirection.h"
#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"
#include "Geometry/Surface/interface/BoundCylinder.h"
#include "Geometry/Surface/interface/BoundDisk.h"
#include <functional>

/** Defines order of layers in the Tracker as seen by straight tracks
 *  coming from the interaction region.
 */

class TkLayerLess 
  : public binary_function< const DetLayer*,const DetLayer*,bool> {
public:

  TkLayerLess( PropagationDirection dir = alongMomentum) :
    theDir(dir) {}

  bool operator()( const DetLayer* a, const DetLayer* b) const {
    if (theDir == alongMomentum) return insideOutLess( a, b);
    else return insideOutLess( b, a);
  }

private:

  PropagationDirection theDir;

  bool insideOutLess( const DetLayer*,const DetLayer*) const;

  bool barrelForwardLess( const BarrelDetLayer* blb,
			  const ForwardDetLayer* fla) const;

};

#endif
