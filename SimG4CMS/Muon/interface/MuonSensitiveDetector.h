#ifndef MuonSensitiveDetector_h
#define MuonSensitiveDetector_h

/** \class MuonSensitiveDetector
 *
 * implementation of SensitiveDetector for the muon detector;
 * a MuonSlaveSD handles the interfacing to the database;
 * numbering scheme are booked according
 * to the detector name
 * 
 * \author Arno Straessner, CERN <arno.straessner@cern.ch>
 *
 * Modification:
 * 19/05/03. P.Arce
 * Add SimTracks selection
 */

#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/SensitiveDetector/interface/SensitiveTkDetector.h"
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "Geometry/Vector/interface/LocalPoint.h"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"

#include <string>

class MuonSlaveSD;
class MuonSimHitNumberingScheme;
class UpdatablePSimHit;
class MuonSubDetector;
class MuonG4Numbering;
class SimHitPrinter;
class TrackInformation;
class G4Track;
class G4ProcessTypeEnumerator;
class G4TrackToParticleID;
class SimTrackManager;

class MuonSensitiveDetector : 
public SensitiveTkDetector,
public Observer<const BeginOfEvent*>
 {

 public:    
  MuonSensitiveDetector(std::string, const DDCompactView &,
			edm::ParameterSet const &,
			const SimTrackManager*);
  virtual ~MuonSensitiveDetector();
  virtual G4bool ProcessHits(G4Step *,G4TouchableHistory *);
  virtual uint32_t setDetUnitId(G4Step *);
  virtual void EndOfEvent(G4HCofThisEvent*);

  void fillHits(edm::PSimHitContainer&, std::string use);
  std::vector<std::string> getNames();
  std::string type();

  const MuonSlaveSD* GetSlaveMuon() const {
    return slaveMuon; }
  
 private:
  bool hitInChimney(G4Step * aStep);
  void update(const BeginOfEvent *);
  virtual void clearHits();

  Local3DPoint toOrcaUnits(Local3DPoint);
  Global3DPoint toOrcaUnits(Global3DPoint);
  
  TrackInformation* getOrCreateTrackInformation( const G4Track* theTrack );

 private:
  MuonSlaveSD* slaveMuon;
  MuonSimHitNumberingScheme* numbering;
  MuonSubDetector* detector;
  MuonG4Numbering* g4numbering;

  void storeVolumeAndTrack(G4Step *);
  bool newHit(G4Step *);
  void createHit(G4Step *);
  void updateHit(G4Step *);
  void saveHit();

  G4VPhysicalVolume * thePV;
  UpdatablePSimHit* theHit;
  uint32_t theDetUnitId; 
  unsigned int theTrackID;
 
  bool printHits;
  SimHitPrinter* thePrinter;
  Global3DPoint theGlobalEntry;

  //--- SimTracks cuts
  double STenergyPersistentCut;
  bool STallMuonsPersistent;

  G4ProcessTypeEnumerator* theG4ProcessTypeEnumerator;

  G4TrackToParticleID* myG4TrackToParticleID;

};

#endif // MuonSensitiveDetector_h
