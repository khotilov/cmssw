#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/ModifiedSurfaceGenerator.h"
#include "DataFormats/TrackingRecHit/interface/AlignmentPositionError.h"

GeomDet::GeomDet( BoundPlane* plane):
  thePlane(plane), theAlignmentPositionError(0) {}

GeomDet::GeomDet( const ReferenceCountingPointer<BoundPlane>& plane) :
  thePlane(plane), theAlignmentPositionError(0) {}

GeomDet::~GeomDet() 
{
  delete theAlignmentPositionError;
}

void GeomDet::move( const GlobalVector& displacement)
{
  //
  // Should recreate the surface like the set* methods ?
  //
  thePlane->move(displacement);
}

void GeomDet::rotate( const Surface::RotationType& rotation)
{
  //
  // Should recreate the surface like the set* methods ?
  //
  thePlane->rotate(rotation);
}

void GeomDet::setPosition( const Surface::PositionType& position, 
			   const Surface::RotationType& rotation)
{
  thePlane = ModifiedSurfaceGenerator<BoundPlane>(thePlane).atNewPosition(position,
									  rotation);
}

void GeomDet::setAlignmentPositionError (const AlignmentPositionError& ape) 
{
  if (!theAlignmentPositionError) {
    theAlignmentPositionError = new AlignmentPositionError(ape);
  } 
  else *theAlignmentPositionError = ape;

}
