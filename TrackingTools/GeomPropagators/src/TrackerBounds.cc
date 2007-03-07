#include "TrackingTools/GeomPropagators/interface/TrackerBounds.h"

#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"

//Ported from ORCA
//  $Date: 2007/03/06 14:34:33 $
//  $Revision: 1.2 $

void TrackerBounds::initialize() 
{
  const float epsilon = 0.001; // should not matter at all

  Surface::RotationType rot; // unit rotation matrix

  theCylinder = new BoundCylinder( Surface::PositionType(0,0,0), rot, 
				   SimpleCylinderBounds( radius()-epsilon, 
							 radius()+epsilon, 
							 -halfLength(), 
							 halfLength()));
  theNegativeDisk = 
    new BoundDisk( Surface::PositionType( 0, 0, -halfLength()), rot, 
		   SimpleDiskBounds( 0, radius(), -epsilon, epsilon));

  thePositiveDisk = 
    new BoundDisk( Surface::PositionType( 0, 0, halfLength()), rot, 
		   SimpleDiskBounds( 0, radius(), -epsilon, epsilon));


  theInit = true;
}

bool TrackerBounds::isInside(const GlobalPoint &point){
  return (point.perp() <= radius() &&
	  fabs(point.z()) <= halfLength());
}


// static initializers

ReferenceCountingPointer<BoundCylinder>  TrackerBounds::theCylinder = 0;
ReferenceCountingPointer<BoundDisk>      TrackerBounds::theNegativeDisk = 0;
ReferenceCountingPointer<BoundDisk>      TrackerBounds::thePositiveDisk = 0;

bool TrackerBounds::theInit = false;


