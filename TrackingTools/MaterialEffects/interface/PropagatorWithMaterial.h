#ifndef _COMMONRECO_PROPAGATORWITHMATERIAL_H_
#define _COMMONRECO_PROPAGATORWITHMATERIAL_H_

#include "Geometry/CommonDetAlgo/interface/DeepCopyPointerByClone.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/MaterialEffects/interface/MaterialEffectsUpdator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

using namespace std;

/** Propagation including material effects.
 *
 *  Propagates using a specific for the geometrical part 
 *   and a MaterialEffectsUpdator to include multiple scattering and
 *   energy loss. By default material effects are included at the
 *   source in the case of forward propagation and at the destination
 *   for backward propagation. Material effects at the source can
 *   only be included when propagating from a TrajectoryStateOnSurface.
 */
class PropagatorWithMaterial : public Propagator {
  
public:
  /** Constructor with PropagationDirection and mass hypothesis.
   *  Uses AnalyticalPropagator and CombinedMaterialEffectsUpdator 
   *  with explicit mass hypothesis.
   */
  PropagatorWithMaterial (PropagationDirection dir, const float mass,  const MagneticField * mf=0);
  /** Constructor with explicit propagator and material effects objects.
   */
  PropagatorWithMaterial(const Propagator& Propagator,
			 const MaterialEffectsUpdator& MEUpdator,  const MagneticField * mf=0);

  ~PropagatorWithMaterial() {}

  virtual TrajectoryStateOnSurface propagate (const TrajectoryStateOnSurface& tsos, 
					      const Plane& plane) const
  {
    // should be implemented (in case underlying propagator has an independent
    // implementation)
    return propagateWithPath(tsos,plane).first;
  }

  virtual TrajectoryStateOnSurface propagate (const FreeTrajectoryState& fts, 
					      const Plane& plane) const
  {
    // should be implemented (in case underlying propagator has an independent
    // implementation)
    return propagateWithPath(fts,plane).first;
  }

  virtual pair<TrajectoryStateOnSurface,double> propagateWithPath (const TrajectoryStateOnSurface& tsos, 
								   const Plane& plane) const;

  virtual pair<TrajectoryStateOnSurface,double> propagateWithPath (const FreeTrajectoryState& fts, 
								   const Plane& plane) const;

  virtual TrajectoryStateOnSurface propagate (const TrajectoryStateOnSurface& tsos, 
					      const Cylinder& cylinder) const
  {
    // should be implemented (in case underlying propagator has an independent
    // implementation)
    return propagateWithPath(tsos,cylinder).first;
  }

  virtual TrajectoryStateOnSurface propagate (const FreeTrajectoryState& fts, 
					      const Cylinder& cylinder) const
  {
    // should be implemented (in case underlying propagator has an independent
    // implementation)
    return propagateWithPath(fts,cylinder).first;
  }

  virtual pair<TrajectoryStateOnSurface,double> propagateWithPath (const TrajectoryStateOnSurface& tsos, 
								   const Cylinder& cylinder) const;

  virtual pair<TrajectoryStateOnSurface,double> propagateWithPath (const FreeTrajectoryState& fts, 
								   const Cylinder& cylinder) const;

  /// Limit on change in azimuthal angle 
  virtual bool setMaxDirectionChange( float phiMax) { 
    return theGeometricalPropagator->setMaxDirectionChange(phiMax);
  }
  /// Propagation direction
  virtual void setPropagationDirection (PropagationDirection dir);

  enum MaterialLocation {atSource, atDestination, fromDirection};
  /** Choice of location for including material effects:
   *  fromDirection is equivalent to atSource for propagation alongMomentum
   *  and to atDestination for propagation oppositeToMomentum. 
   *  Inclusion of material effects at the source (either explicitely or
   *  implicitely) is not possible if propagating with anyDirection and
   *  will effectively disable material effects when propagating from
   *  a FreeTrajectoryState.
   */
  void setMaterialLocation (const MaterialLocation location) {
    theMaterialLocation = location;
  }
  /// Access to the geometrical propagator
  const Propagator& geometricalPropagator() const {
    return *theGeometricalPropagator;
  }
  /// Access to the MaterialEffectsUpdator
  const MaterialEffectsUpdator& materialEffectsUpdator() const {
    return *theMEUpdator;
  }

  virtual const MagneticField* magneticField() const {return field;}


  virtual PropagatorWithMaterial* clone() const
    {
      return new PropagatorWithMaterial(*theGeometricalPropagator,*theMEUpdator);
    }

private:
  /// Inclusion of material at the source?
  bool materialAtSource() const;

private:
  // Geometrical propagator
  DeepCopyPointerByClone<Propagator> theGeometricalPropagator;
  // Material effects
  DeepCopyPointerByClone<MaterialEffectsUpdator> theMEUpdator;
  typedef pair<TrajectoryStateOnSurface,double> TsosWP;
  // Use material at source?
  MaterialLocation theMaterialLocation;
  const MagneticField * field;
};

#endif


