#ifndef TRAJECTORYMANAGER_H
#define TRAJECTORYMANAGER_H

//Framework Headers
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//FAMOS Headers
#include "FastSimulation/MaterialEffects/interface/MaterialEffects.h"

/**
 * This class takes all the particles of a FSimEvent with no end vertex, 
 * and propagate them from their parent vertex through the tracker 
 * material layers, inside which Material Effects are simulated. The
 * propagation is stopped when one of the following configurations is
 * is met:
 *
 *   - the particle has reached the ECAL (possibly after several loops)
 *   - the particle does no longer pass the KineParticleFilter 
 *   - the propagation lasted for more than 100 loops 
 *   - the particle reached an end vertex (e.g., photon conversion)
 *
 * The FSimEvent is updated after each interaction, and the propagation 
 * continues with the updated or the newly created particles. The process
 * ends when there is no new particles to propagate.
 *
 * Charged particles now leave RecHit's, for later track fitting.
 *
 * \author: Florian Beaudette, Patrick Janot
 * $Date Last modification - 18-Jan-2004 
 */

class Pythia6Decays;
class TrackerGeometry;
class TrackerLayer;
class ParticlePropagator;
class FSimEvent;
class FSimTrack;
//class Histos;
//class FamosBasicRecHit;
//class RecHit;

class TrajectoryManager
{
 public:

  /// Default Constructor
  TrajectoryManager() {;}

  /// Constructor from a FSimEvent
  TrajectoryManager(FSimEvent* aSimEvent, 
		    const edm::ParameterSet& matEff,
		    bool activateDecays);

  /// Default Destructor
  ~TrajectoryManager();
  
  /// Does the real job
  void reconstruct();

/// Propagate the particle through the calorimeters
  void propagateToCalorimeters(ParticlePropagator& PP, 
			       int fsimi);


  /// Propagate a particle to a given tracker layer 
  /// (for electron pixel matching mostly)
  bool propagateToLayer(ParticlePropagator& PP,unsigned layer);

  /// Returns the pointer to geometry
  TrackerGeometry* theGeometry();

 private:

  /// Decay the particle and update the SimEvent with daughters 
  void updateWithDaughters(ParticlePropagator& PP, int fsimi);

 private:

  /// Add a RecHit
  //  FamosBasicRecHit* oneHit(const ParticlePropagator& PP, 
  //			   const TrackerLayer& layer,
  //			   unsigned ringNumber) const;

  FSimEvent* mySimEvent;

  TrackerGeometry* _theGeometry;
  
  MaterialEffects* theMaterialEffects;

  Pythia6Decays* myDecayEngine;

  //  Histos* myHistos;

};
#endif
