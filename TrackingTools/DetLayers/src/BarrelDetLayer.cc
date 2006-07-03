#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/Surface/interface/SimpleCylinderBounds.h"
#include "Geometry/Surface/interface/BoundingBox.h"
#include "Geometry/CommonDetUnit/interface/ModifiedSurfaceGenerator.h"

using namespace std;

BarrelDetLayer::~BarrelDetLayer() {}



//--- Extension of the interface
void BarrelDetLayer::setSurface( BoundCylinder* cp) { 
  theCylinder = cp;
}

bool BarrelDetLayer::contains(const Local3DPoint& p) const {
  return surface().bounds().inside(p);
}

void BarrelDetLayer::initialize() 
{
  setSurface( computeSurface());
}



//--- protected methods
BoundCylinder* BarrelDetLayer::computeSurface() {
  vector< const GeometricSearchDet*> comps = components();

  // Find extension in Z
  float theRmin = comps.front()->position().perp();   
  float theRmax = theRmin;
  float theZmin = comps.front()->position().z(); 
  float theZmax = theZmin;
  for ( vector< const GeometricSearchDet*>::const_iterator deti = comps.begin(); 
	deti != comps.end(); deti++) {
    vector<GlobalPoint> corners = 
      BoundingBox().corners( dynamic_cast<const BoundPlane&>((*deti)->surface()));
    for (vector<GlobalPoint>::const_iterator ic = corners.begin();
	 ic != corners.end(); ic++) {
      float r = ic->perp();
      float z = ic->z();
      theRmin = min( theRmin, r);
      theRmax = max( theRmax, r);
      theZmin = min( theZmin, z);
      theZmax = max( theZmax, z);
    }
    // in addition to the corners we have to check the middle of the 
    // det +/- thickness/2
    // , since the min  radius for some barrel dets is reached there
    float rdet = (**deti).position().perp();
    float thick = (**deti).surface().bounds().thickness();
    theRmin = min( theRmin, rdet-thick/2.F);
    theRmax = max( theRmax, rdet+thick/2.F);
  }
  
  // By default the barrel layers are positioned at the center of the 
  // global frame, and the axes of their local frame coincide with 
  // those of the global grame (z along the cylinder axis)
  PositionType pos(0.,0.,0.);
  RotationType rot;

  return new BoundCylinder( pos, rot, 
			    SimpleCylinderBounds( theRmin, theRmax, 
						  theZmin, theZmax));
}  


pair<bool, TrajectoryStateOnSurface>
BarrelDetLayer::compatible( const TrajectoryStateOnSurface& ts, 
			    const Propagator& prop, 
			    const MeasurementEstimator& est) const
{
  if(theCylinder == 0)  edm::LogError("DetLayers") 
    << "ERROR: BarrelDetLayer::compatible() is used before the layer surface is initialized" ;
  // throw an exception? which one?

  TrajectoryStateOnSurface myState = prop.propagate( ts, specificSurface());
  if ( !myState.isValid()) return make_pair( false, myState);

  // take into account the thickness of the layer
  float deltaZ = surface().bounds().thickness()/2. / 
    fabs( tan( myState.globalDirection().theta()));

  // take into account the error on the predicted state
  const float nSigma = 3.;
  if (myState.hasError()) {
    deltaZ += nSigma * sqrt( ts.cartesianError().position().czz());
  }
  //
  // check z assuming symmetric bounds around position().z()
  //
  deltaZ += surface().bounds().length()/2;
  return make_pair(fabs(myState.globalPosition().z()-surface().position().z())<deltaZ,
		   myState);  
}
