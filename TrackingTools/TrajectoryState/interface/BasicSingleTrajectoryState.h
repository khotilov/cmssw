#ifndef BasicSingleTrajectoryState_H
#define BasicSingleTrajectoryState_H

#include "TrackingTools/TrajectoryState/interface/BasicTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryParametrization/interface/LocalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/LocalTrajectoryError.h"

#include "DataFormats/GeometrySurface/interface/ReferenceCounted.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/DeepCopyPointer.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "TrackingTools/TrajectoryParametrization/interface/TrajectoryStateExceptions.h"

class MagneticField;

/** Concrete implementation for the state of one trajectory on a surface.
 */

class BasicSingleTrajectoryState : public BasicTrajectoryState {
public:

  /** Constructor from FTS and surface. For surfaces with material
   *  the side of the surface should be specified explicitely.
   */
  BasicSingleTrajectoryState( const FreeTrajectoryState& fts,
			      const Surface& aSurface,
			      const SurfaceSide side = atCenterOfSurface);
  /** Constructor from global parameters and surface. For surfaces with material
   *  the side of the surface should be specified explicitely.
   */
  BasicSingleTrajectoryState( const GlobalTrajectoryParameters& par,
			      const Surface& aSurface,
			      const SurfaceSide side = atCenterOfSurface);
  /** Constructor from global parameters, errors and surface. For surfaces 
   *  with material the side of the surface should be specified explicitely.
   */
  BasicSingleTrajectoryState( const GlobalTrajectoryParameters& par,
			      const CartesianTrajectoryError& err,
			      const Surface& aSurface,
			      const SurfaceSide side = atCenterOfSurface);
  /** Constructor from global parameters, errors and surface. For surfaces 
   *  with material the side of the surface should be specified explicitely. 
   *  For multi-states the weight should be specified explicitely.
   */
  BasicSingleTrajectoryState( const GlobalTrajectoryParameters& par,
			      const CurvilinearTrajectoryError& err,
			      const Surface& aSurface,
			      const SurfaceSide side = atCenterOfSurface,
			      double weight = 1.);
  /** Constructor from global parameters, errors and surface. For multi-states the
   *  weight should be specified explicitely. For backward compatibility without
   *  specification of the side of the surface.
   */
  BasicSingleTrajectoryState( const GlobalTrajectoryParameters& par,
			      const CurvilinearTrajectoryError& err,
			      const Surface& aSurface,
			      double weight);
  /** Constructor from local parameters, errors and surface. For surfaces 
   *  with material the side of the surface should be specified explicitely.
   */
  BasicSingleTrajectoryState( const LocalTrajectoryParameters& par,
			      const Surface& aSurface,
			      const MagneticField* field,
			      const SurfaceSide side = atCenterOfSurface);
  /** Constructor from local parameters, errors and surface. For surfaces 
   *  with material the side of the surface should be specified explicitely. 
   *  For multi-states the weight should be specified explicitely.
   */
  BasicSingleTrajectoryState( const LocalTrajectoryParameters& par,
			      const LocalTrajectoryError& err,
			      const Surface& aSurface,
			      const MagneticField* field,
			      const SurfaceSide side = atCenterOfSurface,
			      double weight = 1.);
  /** Constructor from local parameters, errors and surface. For multi-states the
   *  weight should be specified explicitely. For backward compatibility without
   *  specification of the side of the surface.
   */
  BasicSingleTrajectoryState( const LocalTrajectoryParameters& par,
			      const LocalTrajectoryError& err,
			      const Surface& aSurface,
			      const MagneticField* field,
			      double weight);

/// construct invalid trajectory state (without parameters)
  BasicSingleTrajectoryState(const Surface& aSurface);

  virtual ~BasicSingleTrajectoryState();

  bool isValid() const {
    return theFreeState || theLocalParametersValid;
  }

  bool hasError() const {
    return (theFreeState && theFreeState->hasError()) || theLocalErrorValid;
  }

// access global parameters/errors
  const GlobalTrajectoryParameters& globalParameters() const {
    return freeTrajectoryState()->parameters();
  }
  GlobalPoint globalPosition() const {
    return freeTrajectoryState()->position();
  }
  GlobalVector globalMomentum() const {
    return freeTrajectoryState()->momentum();
  }
  GlobalVector globalDirection() const {
    return freeTrajectoryState()->momentum().unit();
  }  
  TrackCharge charge() const {
    return freeTrajectoryState()->charge();
  }
  double signedInverseMomentum() const {
    return freeTrajectoryState()->signedInverseMomentum();
  }
  double transverseCurvature() const {
    return freeTrajectoryState()->transverseCurvature();
  }
  const CartesianTrajectoryError& cartesianError() const {
    if(!hasError()) throw TrajectoryStateException(
     "TrajectoryStateOnSurface: attempt to access cartesian errors when none available");
    return freeTrajectoryState()->cartesianError();
  }
  const CurvilinearTrajectoryError& curvilinearError() const {
    if(!hasError()) throw TrajectoryStateException(
     "TrajectoryStateOnSurface: attempt to access curvilinear errors when none available");
    return freeTrajectoryState()->curvilinearError();
  }
  FreeTrajectoryState* freeTrajectoryState() const {
    if(!isValid()) throw TrajectoryStateException(
      "TrajectoryStateOnSurface is invalid and cannot return any parameters");
    checkGlobalParameters();
    if(hasError()) {
      checkCartesianError();
      checkCurvilinError();
    }
    return &(*theFreeState);
  }
  
// access local parameters/errors
  const LocalTrajectoryParameters& localParameters() const {
    if (!isValid()) throw TrajectoryStateException(
      "TrajectoryStateOnSurface is invalid and cannot return any parameters");
    if (!theLocalParametersValid)
      createLocalParameters();
    return theLocalParameters;
  }
  LocalPoint localPosition() const {
    return localParameters().position();
  }
  LocalVector localMomentum() const {
    return localParameters().momentum();
  }
  LocalVector localDirection() const {
    return localMomentum().unit();
  }
  const LocalTrajectoryError& localError() const {
    if (!hasError()) throw TrajectoryStateException(
      "TrajectoryStateOnSurface: attempt to access errors when none available");
    if (!theLocalErrorValid)
      createLocalError();
    return theLocalError;
  }

  const Surface& surface() const {
    return *theSurfaceP;
  }

  virtual double weight() const {return theWeight;} 

  void rescaleError(double factor) {
    if (!hasError()) throw TrajectoryStateException(
      "TrajectoryStateOnSurface: attempt to access errors when none available");
    if (theLocalErrorValid) theLocalError *= (factor*factor);
    if (theFreeState)
      theFreeState->rescaleError(factor);
  }

  BasicSingleTrajectoryState* clone() const {
    return new BasicSingleTrajectoryState(*this);
  }

  /// Position relative to material, defined relative to momentum vector.
  virtual SurfaceSide surfaceSide() const {
    return theSurfaceSide;
  }

private:

// create global parameters and errors from local
  void checkGlobalParameters() const;
  void checkCurvilinError() const;
  void checkCartesianError() const;

// create local parameters and errors from global
  void createLocalParameters() const;
  // create local errors from global
  void createLocalError() const;
  void createLocalErrorFromCartesianError() const;
  void createLocalErrorFromCurvilinearError() const;

private:

  mutable DeepCopyPointer<FreeTrajectoryState> theFreeState;

  mutable bool theGlobalParamsUp2Date;
  mutable bool theCartesianErrorUp2Date;
  mutable bool theCurvilinErrorUp2Date;

  mutable LocalTrajectoryParameters theLocalParameters;
  mutable LocalTrajectoryError      theLocalError;
  mutable bool                      theLocalParametersValid;
  mutable bool                      theLocalErrorValid;

  ConstReferenceCountingPointer<Surface> theSurfaceP;
  const SurfaceSide theSurfaceSide;
  double theWeight;
  const MagneticField* theField;

};

#endif
