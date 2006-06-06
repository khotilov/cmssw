#ifndef Alignment_CommonAlignment_Alignable_H
#define Alignment_CommonAlignment_Alignable_H

#include <iostream>

#include "Geometry/Surface/interface/Surface.h"
#include "Geometry/Vector/interface/LocalVector.h"
#include "DataFormats/TrackingRecHit/interface/AlignmentPositionError.h"
#include "Geometry/CommonDetUnit/interface/DetPositioner.h"
//#include "CondFormats/Alignment/interface/Alignments.h"
//#include "CondFormats/Alignment/interface/AlignmentErrors.h"

// Headers in the same package
#include "Alignment/CommonAlignment/interface/AlignableObjectId.h"
#include "Alignment/CommonAlignment/interface/AlignmentParameters.h"
#include "Alignment/CommonAlignment/interface/AlignableSurface.h"


///  Abstract base class for alignable entities.
///  Any Alignable object can be moved and rotated.
///  Also an alignment uncertainty can be set.
///  The class derives from DetPositioner, a friend class of
///  GeomDet, which allows to move the GeomDet. 

class Alignable : public DetPositioner
{  
  
public:
  
  typedef Surface::RotationType    RotationType;
  typedef Surface::PositionType    PositionType;
  typedef AlignableObjectId::AlignableObjectIdType AlignableObjectIdType;
  
  /// Default constructor
  Alignable() : theMisalignmentActive(true), 
				theAlignmentParameters(0) {}

  /// Destructor
  virtual ~Alignable() { delete theAlignmentParameters; }

  /// Set the AlignmentParameters
  void setAlignmentParameters( AlignmentParameters* dap );
  /// Get the AlignmentParameters
  AlignmentParameters* alignmentParameters() const;

  /// Return vector of all direct components
  virtual std::vector<Alignable*> components() const = 0;
  /// Return number of direct components
  inline const int size() const { return components().size(); }

  /// Movement with respect to the global reference frame
  virtual void move( const GlobalVector& displacement) = 0;

  /// Rotation interpreted in global reference frame
  virtual void rotateInGlobalFrame( const RotationType& rotation) = 0;
  
  /// Rotation intepreted in the local reference frame
  virtual void rotateInLocalFrame( const RotationType& rotation);
  
  /// Rotation around arbitratry global axis
  virtual void rotateAroundGlobalAxis( const GlobalVector axis, const float radians );

  /// Rotation around arbitratry local axis
  virtual void rotateAroundLocalAxis( const LocalVector axis, const float radians );

  /// Rotation around global x-axis
  virtual void rotateAroundGlobalX( const float radians );

  /// Rotation around local x-axis
  virtual void rotateAroundLocalX( const float radians );

  /// Rotation around global y-axis
  virtual void rotateAroundGlobalY( const float radians );

  /// Rotation around local y-axis
  virtual void rotateAroundLocalY( const float radians ); 

  /// Rotation around global z-axis
  virtual void rotateAroundGlobalZ( const float radians );

  /// Rotation around local z-axis
  virtual void rotateAroundLocalZ( const float radians);


  /// Return the global position of the object
  virtual const GlobalPoint globalPosition () const = 0;

  /// Return the global orientation of the object
  virtual const RotationType globalRotation () const = 0;

  /// Return the Surface (global position and orientation) of the object
  virtual const AlignableSurface& surface () const = 0;

  /// Return change of the global position since the creation of the object
  virtual const GlobalVector displacement() const { return theDisplacement; }

  /// Return change of orientation since the creation of the object 
  virtual const RotationType rotation() const { return theRotation; }

  /// Set the alignment position error
  virtual void 
  setAlignmentPositionError( const AlignmentPositionError& ape ) = 0;

  /// Add (or set if not already present) the AlignmentPositionError
  virtual void 
  addAlignmentPositionError( const AlignmentPositionError& ape ) = 0;

  /// add (or set if not already present) the AlignmentPositionError 
  /// which would result from a rotation (given in the GLOBAL frame
  /// of CMS) of the alignable object
  virtual void 
  addAlignmentPositionErrorFromRotation( const RotationType& rotation ) = 0;

  /// add (or set if not already present) the AlignmentPositionError 
  /// which would result from a rotation (given in the LOCAL frame
  /// of the Alignable)  of the alignable object
  virtual void 
  addAlignmentPositionErrorFromLocalRotation( const RotationType& rotation ) = 0;

  /// Restore original position
  virtual void deactivateMisalignment() = 0;

  /// Restore misaligned position
  virtual void reactivateMisalignment() = 0;

  /// Return true if the misalignment is active
  bool misalignmentActive() const { return theMisalignmentActive; }

  /// Return the alignable type identifier
  virtual int alignableObjectId() const = 0;

  /// Recursive printout of alignable information
  virtual void dump() const = 0;

  /// Return vector of alignment data
  //virtual Alignments* alignments() const = 0;
  
  /// Return vector of alignment errors
  //virtual AlignmentErrors* alignmentErrors() const = 0;

protected:

  virtual void addDisplacement( const GlobalVector& displacement );
  virtual void addRotation( const RotationType& rotation );

protected:
  bool theMisalignmentActive;           ///< (de)activation flag

private:

  GlobalVector theDisplacement;
  RotationType theRotation;

  AlignmentParameters* theAlignmentParameters;

};




#endif
