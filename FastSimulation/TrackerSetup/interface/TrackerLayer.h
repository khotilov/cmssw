#ifndef _TrackerLayer_H_
#define _TrackerLayer_H_

#include "DataFormats/GeometrySurface/interface/BoundSurface.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

#include <vector>

/** A class that gives some properties of the Tracker Layers in FAMOS
 */

class TrackerLayer {
public:
  
  /// constructor from private members
  TrackerLayer(BoundSurface* theSurface,
	       bool isForward,
	       unsigned int theLayerNumber,
	       double theModuleThickness = 0.,
	       double theResolutionAlongX = 0.,
	       double theResolutionAlongY = 0.,
	       double theHitEfficiency    = 1. ) :
    theSurface(theSurface), 
    isForward(isForward),
    theLayerNumber(theLayerNumber),
    theResolutionAlongX(theResolutionAlongX),
    theResolutionAlongY(theResolutionAlongY),
    theHitEfficiency(theHitEfficiency),
    theModuleThickness(theModuleThickness),
    theNumberOfFudgeFactors(0)
   { 
     isSensitive = (theLayerNumber<100);
     theFirstRing = 0;
     theLastRing = 0;
     if ( isForward ) { 
       theDisk = dynamic_cast<BoundDisk*>(theSurface);
       theDiskInnerRadius = theDisk->innerRadius();
       theDiskOuterRadius = theDisk->outerRadius();
       theCylinder = 0;
     } else {
       theCylinder = dynamic_cast<BoundCylinder*>(theSurface);
       theDisk = 0;
       theDiskInnerRadius = 0.;
       theDiskOuterRadius = 0.;
     }

   }

  TrackerLayer(BoundSurface* theSurface,
	       unsigned int theLayerNumber,
	       double theModuleThickness, 
	       unsigned int theFirstRing, 
	       unsigned int theLastRing) :
    theSurface(theSurface), 
    theLayerNumber(theLayerNumber),
    theFirstRing(theFirstRing),
    theLastRing(theLastRing),
    theModuleThickness(theModuleThickness),
    theNumberOfFudgeFactors(0)
   { 
     isSensitive = true;
     isForward = true;
     theResolutionAlongX = 0.;
     theResolutionAlongY = 0.;
     theHitEfficiency = 1.;
     theDisk = dynamic_cast<BoundDisk*>(theSurface);
     theDiskInnerRadius = theDisk->innerRadius();
     theDiskOuterRadius = theDisk->outerRadius();
     theCylinder = 0;
   }

  /// Is the layer sensitive ?
  inline bool sensitive() const { return isSensitive; }

  /// Is the layer forward ?
  inline bool forward() const { return isForward; }

  /// Returns the surface
  inline const BoundSurface& surface() const { return *theSurface; }

  /// Returns the cylinder
  inline BoundCylinder* cylinder() const { return theCylinder; }

  /// Returns the surface
  inline BoundDisk* disk() const { return theDisk; }

  /// Returns the layer number  
  inline unsigned int layerNumber() const { return theLayerNumber; }

  /// Returns the first ring  
  inline unsigned int firstRing() const { return theFirstRing; }

  /// Returns the lasst ring  
  inline unsigned int lastRing() const { return theLastRing; }

  /// Returns the resolution along x in cm (local coordinates)
  inline double resolutionAlongxInCm() const { return theResolutionAlongX; }

  /// Returns the resolution along y in cm(local coordinates)
  inline double resolutionAlongyInCm() const { return theResolutionAlongY; }

  /// Returns the hit reconstruction efficiency
  inline double hitEfficiency() const { return theHitEfficiency; }

  /// Returns the sensitive module thickness
  inline double moduleThickness() const { return theModuleThickness; }

  /// Returns the inner radius of a disk
  inline double diskInnerRadius() const { return theDiskInnerRadius; }
  /// Returns the outer radius of a disk
  inline double diskOuterRadius() const { return theDiskOuterRadius; }

  /// Set a fudge factor for material inhomogeneities in this layer
  void setFudgeFactor(double min, double max, double f) { 
    ++theNumberOfFudgeFactors;
    theDimensionMinValues.push_back(min);
    theDimensionMaxValues.push_back(max);
    theFudgeFactors.push_back(f);
  }

  /// Get the fudge factors back
  inline unsigned int fudgeNumber() const { return  theNumberOfFudgeFactors; }
  inline double fudgeMin(unsigned iFudge) const { 
    return (iFudge < theNumberOfFudgeFactors) ? theDimensionMinValues[iFudge] : 999.;
  }
  inline double fudgeMax(unsigned iFudge) const { 
    return (iFudge < theNumberOfFudgeFactors) ? theDimensionMaxValues[iFudge] : -999.;
  }
  inline double fudgeFactor(unsigned iFudge) const { 
    return (iFudge < theNumberOfFudgeFactors) ? theFudgeFactors[iFudge] : 0.;
  }

private:

  BoundSurface* theSurface;
  BoundDisk* theDisk;
  BoundCylinder* theCylinder;
  bool isForward;
  unsigned int theLayerNumber;
  unsigned int theFirstRing;
  unsigned int theLastRing;
  double theResolutionAlongX;
  double theResolutionAlongY;
  double theHitEfficiency;
  double theModuleThickness;
  bool isSensitive;
  double theDiskInnerRadius;
  double theDiskOuterRadius;

  /// These are fudges factors to account for the inhomogeneities of the material
  unsigned int  theNumberOfFudgeFactors;
  std::vector<double> theDimensionMinValues;
  std::vector<double> theDimensionMaxValues;
  std::vector<double> theFudgeFactors;
  

};
#endif

