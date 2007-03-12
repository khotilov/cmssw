#ifndef Alignment_CommonAlignment_RigidBodyAlignmentParameters_h
#define Alignment_CommonAlignment_RigidBodyAlignmentParameters_h

#include "Alignment/CommonAlignment/interface/AlignmentParameters.h"

/// \class RigidBodyAlignmentParameters
///
/// Concrete class for alignment parameters and associated quantities 
/// [derived from AlignmentParameters]. The number of parameters
/// N_PARAM is fixed to 6 (3 translations + 3 rotations)
///
///  $Date: 2006/12/23 15:55:47 $
///  $Revision: 1.4 $
/// (last update by $Author: ewidl $)

class RigidBodyAlignmentParameters : public AlignmentParameters 
{

public:

  /// Give parameters a name
  enum AlignmentParameterName 
	{
	  dx = 0, dy, dz,
	  dalpha, dbeta, dgamma,
	  N_PARAM
	};

  /// Constructor creating parameters (no covariance) from current (mis-)placement of alignable
  explicit RigidBodyAlignmentParameters(Alignable* alignable);

  /// Constructor for full set of parameters
  RigidBodyAlignmentParameters( Alignable* alignable, 
				const AlgebraicVector& parameters, 
				const AlgebraicSymMatrix& covMatrix );

  /// Constructor for selection 
  RigidBodyAlignmentParameters( Alignable* alignable, const AlgebraicVector& parameters, 
				const AlgebraicSymMatrix& covMatrix, 
				const std::vector<bool>& selection );

  /// Destructor 
  virtual ~RigidBodyAlignmentParameters() {};

  /// Clone all parameters (for update of parameters)
  virtual RigidBodyAlignmentParameters* clone( const AlgebraicVector& parameters, 
					       const AlgebraicSymMatrix& covMatrix ) const;
 
  /// Clone selected parameters (for update of parameters)
    virtual RigidBodyAlignmentParameters* 
      cloneFromSelected( const AlgebraicVector& parameters, 
			 const AlgebraicSymMatrix& covMatrix ) const;
  
  /// Get all derivatives 
  virtual AlgebraicMatrix derivatives( const TrajectoryStateOnSurface& tsos,
				       AlignableDet* ) const;

  /// Get selected derivatives
  AlgebraicMatrix selectedDerivatives( const TrajectoryStateOnSurface& tsos, 
				       AlignableDet* ) const;

  /// Get translation parameters
  AlgebraicVector translation(void) const;

  /// Get rotation parameters
  AlgebraicVector rotation(void) const;

  /// calculate and return parameters in global frame 
  AlgebraicVector globalParameters(void) const;

  /// print parameters to screen 
  void print(void) const;

protected:

  AlgebraicVector displacementFromAlignable(Alignable* ali) const;

};

#endif

