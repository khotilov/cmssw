#ifndef SimG4CMS_TkAccumulatingSensitiveDetector_H
#define SimG4CMS_TkAccumulatingSensitiveDetector_H

/**
 * This is TkSensitiveDetector which accumulates all the SimHits coming from the same track in the
 * same volume, thus reducing the db size.
 */
#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/SensitiveDetector/interface/SensitiveTkDetector.h"
#include "SimG4CMS/Tracker/interface/TrackerG4SimHitNumberingScheme.h"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"

#include <string>

class EventAction;
class TrackingSlaveSD;
class FrameRotation;
class UpdatablePSimHit;
class G4ProcessTypeEnumerator;
class G4TrackToParticleID;

class TkAccumulatingSensitiveDetector : public SensitiveTkDetector//, private Observer<const BeginOfEvent*>
{ 
public:    
    TkAccumulatingSensitiveDetector(std::string, const DDCompactView &,
				    edm::ParameterSet const &);
    virtual ~TkAccumulatingSensitiveDetector();
    virtual bool ProcessHits(G4Step *,G4TouchableHistory *);
    virtual uint32_t setDetUnitId(G4Step*);
    virtual void EndOfEvent(G4HCofThisEvent*);

    void fillHits(edm::PSimHitContainer&, std::string use);
    std::vector<std::string> getNames();
    std::string type();

private:
    virtual void sendHit();
    virtual void updateHit(G4Step *);
    virtual bool newHit(G4Step *);
    virtual bool closeHit(G4Step *);
    virtual void createHit(G4Step *);
    void checkExitPoint(Local3DPoint);
    void update(const BeginOfEvent *);
    virtual void clearHits();
    Local3DPoint toOrcaRef(Local3DPoint ,G4VPhysicalVolume *);
    int tofBin(float);
    std::string myName;
    TrackingSlaveSD * slaveLowTof;
    TrackingSlaveSD * slaveHighTof;
    FrameRotation * myRotation;
    UpdatablePSimHit * mySimHit;
    std::string pname;
    Local3DPoint globalEntryPoint;
    Local3DPoint globalExitPoint;
    G4VPhysicalVolume * oldVolume;
    G4ProcessTypeEnumerator * theG4ProcessTypeEnumerator;
    TrackerG4SimHitNumberingScheme* tkG4SimHitNumberingScheme;
    double theSigma;
    uint32_t lastId;
    unsigned int lastTrack;
    int eventno;
    // cache stuff for debugging
    float px,py,pz;
    bool allowZeroEnergyLoss;
    bool printHits;
    bool neverAccumulate;
    G4TrackToParticleID * myG4TrackToParticleID;
};

#endif





