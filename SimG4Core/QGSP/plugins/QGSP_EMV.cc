#include "QGSP_EMV.hh"
#ifndef G4V7
#include "SimG4Core/QGSP/src/HadronPhysicsQGSP.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics71.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4QStoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh" 
 
#include "G4DataQuestionaire.hh"
 
QGSP_EMV::QGSP_EMV(G4LogicalVolumeToDDLogicalPartMap& map,
		   const edm::ParameterSet & p) : PhysicsList(map, p) {
    G4DataQuestionaire it(photon);
    std::cout << "You are using the simulation engine: QGSP_EMV 3.1" << std::endl;
  
    // EM Physics
    RegisterPhysics(new G4EmStandardPhysics71("standard EM v71"));
    // Synchroton Radiation & GN Physics
    RegisterPhysics(new G4EmExtraPhysics("extra EM"));
    // Decays
    RegisterPhysics(new G4DecayPhysics("decay"));
    // Hadron Elastic scattering
    RegisterPhysics(new G4HadronElasticPhysics("elastic")); 
    // Hadron Physics
    RegisterPhysics(new HadronPhysicsQGSP("hadron"));
    // Stopping Physics
    RegisterPhysics(new G4QStoppingPhysics("stopping"));
    // Ion Physics
    RegisterPhysics(new G4IonPhysics("ion"));
}
#endif
