#ifndef Alignment_MuonAlignment_AlignableMuon_H
#define Alignment_MuonAlignment_AlignableMuon_H

/** \class AlignableMuon
 *  The alignable muon.
 *
 *  $Date: 2007/12/06 01:31:06 $
 *  $Revision: 1.15.4.1 $
 *  \author Andre Sznajder - UERJ(Brazil)
 */


#include <vector>


#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include <DataFormats/GeometryVector/interface/GlobalPoint.h>
#include <Geometry/CSCGeometry/interface/CSCLayer.h>

#include "Alignment/CommonAlignment/interface/AlignableComposite.h"


// Classes that will be used to construct the muon
class AlignableDTBarrel;
class AlignableDTWheel;
class AlignableDTStation;
class AlignableDTChamber;
class AlignableCSCEndcap;
class AlignableCSCStation;
class AlignableCSCRing;
class AlignableCSCChamber;




/// Constructor of the full muon geometry.
/// This object is stored to the EventSetup for further retrieval

class AlignableMuon: public AlignableComposite 
{

public:

  /// Constructor from geometries
  AlignableMuon( const DTGeometry* , const CSCGeometry* );

  /// Destructor
  ~AlignableMuon();
  
public:

  // Some typdefs to simplify notation
  typedef GlobalPoint           _PositionType;
  typedef TkRotation<float>     _RotationType;

  
  /// Return all components
  virtual std::vector<Alignable*> components() const { return theMuonComponents; }

  /// Alignable tracker has no mother
  virtual Alignable* mother() { return 0; }

  // Methods to return specific of components
  std::vector<Alignable*> DTSuperLayers();
  std::vector<Alignable*> DTChambers();
  std::vector<Alignable*> DTStations();
  std::vector<Alignable*> DTWheels();
  std::vector<Alignable*> DTBarrel();
  std::vector<Alignable*> CSCLayers();
  std::vector<Alignable*> CSCChambers();
  std::vector<Alignable*> CSCStations();
  std::vector<Alignable*> CSCRings();
  std::vector<Alignable*> CSCEndcaps();

  // Get DT alignments sorted by DetId
  Alignments* dtAlignments();

  // Get DT alignment errors sorted by DetId
  AlignmentErrors* dtAlignmentErrors();

  // Get CSC alignments sorted by DetId
  Alignments* cscAlignments();

  // Get CSC alignment errors sorted by DetId
  AlignmentErrors* cscAlignmentErrors();



private:
  
  // Get the position (centered at 0 by default)
  PositionType computePosition(); 

  // Get the global orientation (no rotation by default)
  RotationType computeOrientation();

  // Get the Surface
  AlignableSurface computeSurface();

  // Return alignable object identifier
  virtual StructureType alignableObjectId() const { return align::AlignableMuon; }

  // Get alignments sorted by DetId
  Alignments* alignments() const;

  // Get alignment errors sorted by DetId
  AlignmentErrors* alignmentErrors() const;



   // Sub-structure builders 

   // Build muon barrel
  void buildDTBarrel( const DTGeometry*  );

  // Build muon end caps
  void buildCSCEndcap( const CSCGeometry*  );

  /// Set mothers recursively
  void recursiveSetMothers( Alignable* alignable );


  // Containers of separate components

  std::vector<AlignableDTChamber*>   theDTChambers;
  std::vector<AlignableDTStation*>   theDTStations;
  std::vector<AlignableDTWheel*>     theDTWheels;
  std::vector<AlignableDTBarrel*>    theDTBarrel;
  
  std::vector<AlignableCSCChamber*>  theCSCChambers;
  std::vector<AlignableCSCStation*>  theCSCStations;
  std::vector<AlignableCSCRing*>     theCSCRings;
  std::vector<AlignableCSCEndcap*>   theCSCEndcaps;

  std::vector<Alignable*> theMuonComponents;

};

#endif //AlignableMuon_H




