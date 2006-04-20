#ifndef FastSimulation_Event_FBaseSimEvent_H
#define FastSimulation_Event_FBaseSimEvent_H

//Framework Headers
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// CLHEP Headers
#include "CLHEP/Vector/LorentzVector.h"

#include <map>
#include <vector>

/** FSimEvent special features for FAMOS
 *
 * \author Patrick Janot, CERN
 * \date: 9-Dec-2003
 */

class FSimEvent;
class FSimTrack;
class FSimVertex;
class RawParticle;
class HepPDTable;
class KineParticleFilter;

class EmbdSimTrack;
class EmbdSimVertex;
class BaseEventVertexGenerator;
//class Histos;

namespace HepMC {
  class GenEvent;
  class GenParticle;
}

class FBaseSimEvent  
{

public:

  /// Default constructor
  FBaseSimEvent(const edm::ParameterSet& vtx,
		const edm::ParameterSet& kine);

  ///  usual virtual destructor
  ~FBaseSimEvent();

  /// fill the FBaseSimEvent from the current HepMC::GenEvent
  void fill(const HepMC::GenEvent& hev);
  
  /// print the original MCTruth event
  void printMCTruth(const HepMC::GenEvent& hev);

  /// Add the particles and their vertices to the list
  void addParticles(const HepMC::GenEvent& hev);

  /// print the FBaseSimEvent in an intelligible way
  void print() const;

  /// clear the FBaseSimEvent content before the next event
  void clear();

  /// Number of tracks
  unsigned int nTracks() const;
  /// Number of vertices
  unsigned int nVertices() const;
  /// Number of generator particles
  unsigned int nGenParts() const;

  /// Return track with given Id 
  FSimTrack& track(int id) const;
  /// Return vertex with given Id 
  FSimVertex& vertex(int id) const;
  /// return "reconstructed" charged tracks index.
  int chargedTrack(int id) const;

  /// return embedded track with given id
  const EmbdSimTrack & embdTrack(int i) const;
  /// return embedded vertex with given id
  const EmbdSimVertex & embdVertex(int i) const;
  /// return MC track with a given id
  const HepMC::GenParticle* embdGenpart(int i) const;

  /// The pointer to the vector of FSimTrack's 
  std::vector<FSimTrack>* tracks() const; 
  
  /// The pointer to the vector of FSimVertex's 
  std::vector<FSimVertex>* vertices() const;

  /// The pointer to the vector of GenParticle's 
  std::vector<HepMC::GenParticle*>* genparts() const;

  /// Add a new track to the Event and to the various lists
  //  int addSimTrack(HepMC::GenParticle* part, 
  //		  HepMC::GenVertex* originVertex, 
  //		  int ig=-1);
  int addSimTrack(const RawParticle* p, int iv, int ig=-1);

  /// Add a new vertex to the Event and to the various lists
  //  int addSimVertex(HepMC::GenVertex* decayVertex,int im=-1);
  int addSimVertex(const CLHEP::HepLorentzVector& decayVertex,int im=-1);

  const KineParticleFilter& filter() const { return *myFilter; } 

 private:

  std::vector<FSimTrack>* theSimTracks;
  std::vector<FSimVertex>* theSimVertices;
  std::vector<HepMC::GenParticle*>* theGenParticles;

  /// The particle filter
  KineParticleFilter* myFilter;

  double sigmaVerteX;
  double sigmaVerteY;
  double sigmaVerteZ;

  HepPDTable * tab;

  BaseEventVertexGenerator* theVertexGenerator;

  //  Histos* myHistos;

};

#endif // FBaseSimEvent_H
