#ifndef RecoTracker_TkDetLayers_BoundDiskSector_h
#define RecoTracker_TkDetLayers_BoundDiskSector_h
 
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "DiskSectorBounds.h"

#pragma GCC visibility push(hidden)
class BoundDiskSector GCC11_FINAL : public BoundPlane {
 public:
 
 
  virtual ~BoundDiskSector() {}
 
  BoundDiskSector( const PositionType& pos, 
		   const RotationType& rot, 
		   Bounds* bounds) : Surface( pos,rot),
    BoundPlane( pos, rot, bounds) {}
  
  BoundDiskSector( const PositionType& pos, 
		   const RotationType& rot, 
		   const Bounds& bounds) : Surface( pos,rot),
    BoundPlane( pos, rot, bounds) {}
  
  float innerRadius() const { return bounds().innerRadius();}
  float outerRadius() const  { return bounds().outerRadius();}
  float phiExtension() const  { return bounds().phiExtension();}

  // hide
  DiskSectorBounds const & bounds() const {
    return static_cast<DiskSectorBounds const &>(BoundPlane::bounds());
  }

};
 
 
#pragma GCC visibility pop
#endif 

