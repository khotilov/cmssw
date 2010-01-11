//
// initial setup : E.Barberio & Joanna Weng 
// big changes : Soon Jun & Dongwook Jang
//
#include "SimG4Core/Application/interface/SteppingAction.h"
#include "SimG4Core/GFlash/interface/GflashEMShowerModel.h"

#include "SimGeneral/GFlash/interface/GflashEMShowerProfile.h"
#include "SimGeneral/GFlash/interface/GflashHit.h"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4VProcess.hh"
#include "G4VPhysicalVolume.hh" 
#include "G4LogicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4EventManager.hh"
#include "G4FastSimulationManager.hh"
#include "G4TouchableHandle.hh"
#include "G4VSensitiveDetector.hh"

GflashEMShowerModel::GflashEMShowerModel(G4String modelName, G4Envelope* envelope, edm::ParameterSet parSet)
  : G4VFastSimulationModel(modelName, envelope), theParSet(parSet) {

  theWatcherOn = parSet.getParameter<bool>("watcherOn");

  theProfile = new GflashEMShowerProfile(parSet);

  theGflashStep = new G4Step();
  theGflashTouchableHandle = new G4TouchableHistory();
  theGflashNavigator = new G4Navigator();

}

// -----------------------------------------------------------------------------------

GflashEMShowerModel::~GflashEMShowerModel() {

  if(theProfile) delete theProfile;
  if(theGflashStep) delete theGflashStep;
}

G4bool GflashEMShowerModel::IsApplicable(const G4ParticleDefinition& particleType) { 

  return ( &particleType == G4Electron::ElectronDefinition() ||
	   &particleType == G4Positron::PositronDefinition() ); 

}

// -----------------------------------------------------------------------------------
G4bool GflashEMShowerModel::ModelTrigger(const G4FastTrack & fastTrack ) {

  // Mininum energy cutoff to parameterize
  if(fastTrack.GetPrimaryTrack()->GetKineticEnergy() < 1.0*GeV) return false;
  if(excludeDetectorRegion(fastTrack)) return false;

  G4bool trigger = fastTrack.GetPrimaryTrack()->GetDefinition() == G4Electron::ElectronDefinition() || 
    fastTrack.GetPrimaryTrack()->GetDefinition() == G4Positron::PositronDefinition();

  if(!trigger) return false;

  // This will be changed accordingly when the way dealing with CaloRegion changes later.
  G4TouchableHistory* touch = (G4TouchableHistory*)(fastTrack.GetPrimaryTrack()->GetTouchable());
  G4VPhysicalVolume* pCurrentVolume = touch->GetVolume();
  if( pCurrentVolume == 0) return false;

  G4LogicalVolume* lv = pCurrentVolume->GetLogicalVolume();
  if(lv->GetRegion()->GetName() != "CaloRegion") return false;

  // The parameterization starts inside crystals
  std::size_t pos1 = lv->GetName().find("EBRY");
  std::size_t pos2 = lv->GetName().find("EFRY");
  /*
  std::size_t pos3 = lv->GetName().find("HVQ");
  std::size_t pos4 = lv->GetName().find("HF");
  if(pos1 == std::string::npos && pos2 == std::string::npos &&
     pos3 == std::string::npos && pos4 == std::string::npos) return false;
  */
  //@@@for now, HF is not a part of Gflash Envelopes
  if(pos1 == std::string::npos && pos2 == std::string::npos ) return false;

  return true;

}


// -----------------------------------------------------------------------------------
void GflashEMShowerModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) {

  // Kill the parameterised particle:
  fastStep.KillPrimaryTrack();
  fastStep.ProposePrimaryTrackPathLength(0.0);

  //input variables for GflashEMShowerProfile with showerType = 1 (shower starts inside crystals)
  G4int showerType = 1;
  G4double energy = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/GeV;
  G4double globalTime = fastTrack.GetPrimaryTrack()->GetStep()->GetPostStepPoint()->GetGlobalTime();
  G4double charge = fastTrack.GetPrimaryTrack()->GetStep()->GetPreStepPoint()->GetCharge();
  G4ThreeVector position = fastTrack.GetPrimaryTrack()->GetPosition() / cm;
  G4ThreeVector momentum = fastTrack.GetPrimaryTrack()->GetMomentum()/GeV;

  // Do actual parameterization. The result of parameterization is gflashHitList
  theProfile->initialize(showerType,energy,globalTime,charge,position,momentum);
  theProfile->parameterization();

  //make hits
  makeHits(fastTrack);
}

void GflashEMShowerModel::makeHits(const G4FastTrack& fastTrack) {

  std::vector<GflashHit>& gflashHitList = theProfile->getGflashHitList();

  theGflashStep->SetTrack(const_cast<G4Track*>(fastTrack.GetPrimaryTrack()));

  theGflashStep->GetPostStepPoint()->SetProcessDefinedStep(const_cast<G4VProcess*>
    (fastTrack.GetPrimaryTrack()->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()));
  theGflashNavigator->SetWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());

  std::vector<GflashHit>::const_iterator spotIter    = gflashHitList.begin();
  std::vector<GflashHit>::const_iterator spotIterEnd = gflashHitList.end();

  for( ; spotIter != spotIterEnd; spotIter++){

      //put touchable for each hit so that touchable history keeps track of each step.
      theGflashNavigator->LocateGlobalPointAndUpdateTouchableHandle(spotIter->getPosition(),G4ThreeVector(0,0,0),
								    theGflashTouchableHandle, false);
      updateGflashStep(spotIter->getPosition(),spotIter->getTime());

      // if there is a watcher defined in a job and the flag is turned on
      if(theWatcherOn) {
      	SteppingAction* userSteppingAction = (SteppingAction*) G4EventManager::GetEventManager()->GetUserSteppingAction();
      	userSteppingAction->m_g4StepSignal(theGflashStep);
      }

      // Send G4Step information to Hit/Digi if the volume is sensitive
      // Copied from G4SteppingManager.cc
    
      G4VPhysicalVolume* aCurrentVolume = theGflashStep->GetPreStepPoint()->GetPhysicalVolume();
      if( aCurrentVolume == 0 ) continue;

      G4LogicalVolume* lv = aCurrentVolume->GetLogicalVolume();
      if(lv->GetRegion()->GetName() != "CaloRegion") continue;

      theGflashStep->GetPreStepPoint()->SetSensitiveDetector(aCurrentVolume->GetLogicalVolume()->GetSensitiveDetector());
      G4VSensitiveDetector* aSensitive = theGflashStep->GetPreStepPoint()->GetSensitiveDetector();
      
      if( aSensitive == 0 ) continue;

      theGflashStep->SetTotalEnergyDeposit(spotIter->getEnergy());
      aSensitive->Hit(theGflashStep);
  }

}

void GflashEMShowerModel::updateGflashStep(G4ThreeVector spotPosition, G4double timeGlobal)
{
  theGflashStep->GetPostStepPoint()->SetGlobalTime(timeGlobal);
  theGflashStep->GetPreStepPoint()->SetPosition(spotPosition);
  theGflashStep->GetPostStepPoint()->SetPosition(spotPosition);
  theGflashStep->GetPreStepPoint()->SetTouchableHandle(theGflashTouchableHandle);
}

// -----------------------------------------------------------------------------------
G4bool GflashEMShowerModel::excludeDetectorRegion(const G4FastTrack& fastTrack) {

  G4bool isExcluded=false;

  //exclude regions where geometry are complicated
  G4double eta =   fastTrack.GetPrimaryTrack()->GetPosition().pseudoRapidity() ;
  if(fabs(eta) > 1.3 && fabs(eta) < 1.57) return true;

  return isExcluded;
}
