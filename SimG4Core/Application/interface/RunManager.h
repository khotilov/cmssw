#ifndef SimG4Core_RunManager_H
#define SimG4Core_RunManager_H

#include <memory>
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimG4Core/SensitiveDetector/interface/AttachSD.h"
#include "SimG4Core/SensitiveDetector/interface/SensitiveDetector.h"
#include "SimG4Core/SensitiveDetector/interface/SensitiveTkDetector.h"
#include "SimG4Core/SensitiveDetector/interface/SensitiveCaloDetector.h"

#include "SimG4Core/Notification/interface/SimActivityRegistry.h"

#include <memory>
#include "boost/shared_ptr.hpp"

namespace CLHEP {
  class HepJamesRandom;
}

namespace sim {
   class FieldBuilder;
}

class PrimaryTransformer;
class Generator;
// class BaseEventVertexGenerator;
class PhysicsList;
class SimWatcher;
class SimProducer;
class G4SimEvent;
class SimTrackManager;

class DDDWorld;

class G4RunManagerKernel;
class G4Run;
class G4Event;
class G4UserRunAction;

class RunManager
{
public:
    static RunManager * instance();
    static RunManager * init(edm::ParameterSet const & p); 
    virtual ~RunManager();
    void initG4(const edm::EventSetup & es);
    void initializeUserActions();
    void initializeRun();
    void terminateRun();
    void abortRun(bool softAbort=false);
    const G4Run * currentRun() const { return m_currentRun; }
    void produce(edm::Event& inpevt, const edm::EventSetup & es);
    void abortEvent();
    const Generator * generator() const { return m_generator; }
    const G4Event * currentEvent() const { return m_currentEvent; }
    G4SimEvent * simEvent() { return m_simEvent; }
    std::vector<SensitiveTkDetector*>& sensTkDetectors() { return m_sensTkDets; }
    std::vector<SensitiveCaloDetector*>& sensCaloDetectors() { return m_sensCaloDets; }
    void runRNDMstore(int run);
    void runRNDMrestore(int run);
    void eventRNDMstore(int run, int event);
    void eventRNDMrestore(int run, int event);

    std::vector<boost::shared_ptr<SimProducer> > producers() const {
       return m_producers;
    }
protected:
    G4Event * generateEvent(int evt, edm::Event& inpevt);
private:
    static RunManager * me;
    explicit RunManager(edm::ParameterSet const & p);
    G4RunManagerKernel * m_kernel;
    Generator * m_generator;
    bool m_nonBeam;
    // std::auto_ptr<BaseEventVertexGenerator> m_eventVertexGenerator;
    std::auto_ptr<PhysicsList> m_physicsList;
    PrimaryTransformer * m_primaryTransformer;
    CLHEP::HepJamesRandom * m_engine;
    bool m_managerInitialized;
    bool m_geometryInitialized;
    bool m_physicsInitialized;
    bool m_runInitialized;
    bool m_runTerminated;
    bool m_runAborted;
    bool m_pUseMagneticField;
    G4Run * m_currentRun;
    G4Event * m_currentEvent;
    G4SimEvent * m_simEvent;
    G4UserRunAction * m_userRunAction;
    bool m_rndmStore;
    bool m_rndmRestore;
    std::string m_PhysicsTablesDir;
    bool m_StorePhysicsTables;
    bool m_RestorePhysicsTables;
    int m_EvtMgrVerbosity;
    bool m_Override;
    int m_RunNumber;
    edm::ParameterSet m_pGeometry;
    edm::ParameterSet m_pField;
    edm::ParameterSet m_pGenerator;   
    edm::ParameterSet m_pVertexGenerator;
    edm::ParameterSet m_pPhysics; 
    edm::ParameterSet m_pRunAction;      
    edm::ParameterSet m_pEventAction;
    edm::ParameterSet m_pTrackingAction;
    edm::ParameterSet m_pSteppingAction;
    edm::ParameterSet m_p;
    AttachSD * m_attach;
    std::vector<SensitiveTkDetector*> m_sensTkDets;
    std::vector<SensitiveCaloDetector*> m_sensCaloDets;

    SimActivityRegistry m_registry;
    std::vector<boost::shared_ptr<SimWatcher> > m_watchers;
    std::vector<boost::shared_ptr<SimProducer> > m_producers;
    
    std::auto_ptr<SimTrackManager> m_trackManager;
    std::auto_ptr<sim::FieldBuilder> m_fieldBuilder;
};

#endif
