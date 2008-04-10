#ifndef Alignment_MuonAlignment_AlignableDTBarrel_H
#define Alignment_MuonAlignment_AlignableDTBarrel_H

/** \class AlignableDTBarrel
 *  The alignable muon DT barrel.
 *
 *  $Date: 2008/02/14 09:39:19 $
 *  $Revision: 1.8 $
 *  \author Andre Sznajder - UERJ(Brazil)
 */


#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Alignment/CommonAlignment/interface/AlignableComposite.h"
#include "Alignment/CommonAlignment/interface/AlignableSurface.h"

#include "Alignment/MuonAlignment/interface/AlignableDTWheel.h"



#include <vector>

class GeomDet;

/// Concrete class for muon DT Barrel alignable.
///
/// Misalignment can be de-/reactivated (forwarded to components).
///

class AlignableDTBarrel : public AlignableComposite 
{

 public:

  AlignableDTBarrel( const std::vector<AlignableDTWheel*> dtWheels );

  ~AlignableDTBarrel();
  
  virtual std::vector<Alignable*> components() const 
  {

        std::vector<Alignable*> result;
        result.insert( result.end(), theDTWheels.begin(), theDTWheels.end() );
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

  AlignableDTWheel &wheel(int i);  
  
  /// Printout muon Barrel information (not recursive)
  friend std::ostream& operator << ( std::ostream&, const AlignableDTBarrel& );

  /// Recursive printout of the muon Barrel structure
  void dump( void );


  // Get alignments sorted by DetId
  Alignments* alignments() const;

  // Get alignment errors sorted by DetId
  AlignmentErrors* alignmentErrors() const;



private:

  std::vector<AlignableDTWheel*> theDTWheels;


};

#endif




