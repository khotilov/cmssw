#ifndef Alignment_MuonAlignment_AlignableDTWheel_H
#define Alignment_MuonAlignment_AlignableDTWheel_H

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Alignment/CommonAlignment/interface/AlignableComposite.h"
#include "Alignment/CommonAlignment/interface/AlignableSurface.h"

#include "Alignment/MuonAlignment/interface/AlignableDTStation.h"

#include "DataFormats/TrackingRecHit/interface/AlignmentPositionError.h"
#include "Geometry/Vector/interface/Basic3DVector.h"

#include <vector>

class GeomDet;

/// Concrete class for muon DT Wheel alignable.
///
/// Misalignment can be de-/reactivated (forwarded to components).
///

class AlignableDTWheel : public AlignableComposite 
{

 public:

  AlignableDTWheel( const std::vector<AlignableDTStation*> dtStations );

  ~AlignableDTWheel();
  
  virtual std::vector<Alignable*> components() const 
  {

        std::vector<Alignable*> result;
        result.insert( result.end(), theDTStations.begin(), theDTStations.end() );
        return result;

  }
  
  typedef GlobalPoint           PositionType;
  typedef TkRotation<float>     RotationType;

  // gets the global position as the average over all positions of the layers
  PositionType computePosition() ;
  // get the global orientation
  RotationType computeOrientation() ; //see explanation for "theOrientation"
  // get the Surface
  AlignableSurface computeSurface() ;

  AlignableDTStation &station(int i);  
  
  /// Return alignable object identifier
  virtual int alignableObjectId() const { return AlignableObjectId::AlignableDTWheel; }

  /// Printout muon DT wheel information (not recursive)
  friend std::ostream& operator << ( std::ostream&, const AlignableDTWheel& );

  /// Recursive printout of the muon DT wheel structure
  void dump( void );


private:

  std::vector<AlignableDTStation*> theDTStations;


};

#endif




