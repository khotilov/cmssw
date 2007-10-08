#ifndef Alignment_MuonAlignment_AlignableDTStation_H
#define Alignment_MuonAlignment_AlignableDTStation_H

/** \class AlignableDTStation
 *  The alignable muon DT station.
 *
 *  $Date: 2007/03/02 15:46:12 $
 *  $Revision: 1.5 $
 *  \author Andre Sznajder - UERJ(Brazil)
 */


#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Alignment/CommonAlignment/interface/AlignableComposite.h"
#include "Alignment/CommonAlignment/interface/AlignableSurface.h"

#include "Alignment/MuonAlignment/interface/AlignableDTChamber.h"

#include "DataFormats/TrackingRecHit/interface/AlignmentPositionError.h"
#include "DataFormats/GeometryVector/interface/Basic3DVector.h"

#include <vector>

class GeomDet;

/// Concrete class for muon DT Station alignable.
///
/// Misalignment can be de-/reactivated (forwarded to components).
///

class AlignableDTStation : public AlignableComposite 
{

 public:

  AlignableDTStation( const std::vector<AlignableDTChamber*> dtChambers );

  ~AlignableDTStation();
  
  virtual std::vector<Alignable*> components() const 
  {

        std::vector<Alignable*> result;
        result.insert( result.end(), theDTChambers.begin(), theDTChambers.end() );
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

  AlignableDTChamber &chamber(int i);  
  
  /// Return alignable object identifier
  virtual StructureType alignableObjectId() const { return align::AlignableDTStation; }

  /// Printout muon DT Station information (not recursive)
  friend std::ostream& operator << ( std::ostream&, const AlignableDTStation& );

  /// Recursive printout of the muon DT Station structure
  void dump( void );


private:

  std::vector<AlignableDTChamber*> theDTChambers;


};

#endif




