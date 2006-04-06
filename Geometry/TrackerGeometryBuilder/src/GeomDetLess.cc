#include "Geometry/TrackerGeometryBuilder/interface/GeomDetLess.h"
//#include "CommonDet/DetUtilities/interface/DetExceptions.h"
// temporary solution
#include "Utilities/General/interface/CMSexception.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "Geometry/Surface/interface/Bounds.h"
#include "Geometry/Surface/interface/BoundPlane.h"

#include <iostream>

bool GeomDetLess::insideOutLess( const GeomDet* a, const GeomDet* b) const
{
  if (a == b) return false;

  DetId ida(a->geographicalId());
  DetId idb(b->geographicalId());

  if      ( (ida.subdetId() == StripSubdetector::TIB || ida.subdetId() == StripSubdetector::TOB || ida.subdetId() == PixelSubdetector::PixelBarrel) &&
	    (idb.subdetId() == StripSubdetector::TIB || idb.subdetId() == StripSubdetector::TOB || idb.subdetId() == PixelSubdetector::PixelBarrel)) {  // barrel with barrel
    return a->surface().position().perp() < b->surface().position().perp();
  }

  if      ( (ida.subdetId() == StripSubdetector::TEC || ida.subdetId() == StripSubdetector::TID || ida.subdetId() == PixelSubdetector::PixelEndcap) &&
	    (idb.subdetId() == StripSubdetector::TEC || idb.subdetId() == StripSubdetector::TID || idb.subdetId() == PixelSubdetector::PixelEndcap)) {  // fwd with fwd
    return fabs(a->surface().position().z()) < fabs(b->surface().position().z());
  }
  
  //
  //  here I have 1 barrel against one forward
  //

  if      ( (ida.subdetId() == StripSubdetector::TIB || ida.subdetId() == StripSubdetector::TOB || ida.subdetId() == PixelSubdetector::PixelBarrel) &&
	    (idb.subdetId() == StripSubdetector::TEC || idb.subdetId() == StripSubdetector::TID || idb.subdetId() == PixelSubdetector::PixelEndcap)) {  // barrel with barrel
    return barrelForwardLess( a, b);
  }else{
    return !barrelForwardLess( b, a);
  }
  
  //throw DetLogicError("GeomDetLess: arguments are not Barrel or Forward GeomDets");
  throw Genexception("GeomDetLess: arguments are not Ok");
  
}

bool GeomDetLess::barrelForwardLess( const GeomDet* bla, 
				     const GeomDet* flb) const
{
  
  BoundPlane s = bla->specificSurface();
  Bounds * b     = &(s.bounds());

  return (fabs(bla->surface().position().z())
	  + fabs(b->length()/2.)) 
	  < fabs( flb->position().z());
}


