#include "SimG4Core/Application/interface/StackingAction.h"
#include "SimG4Core/Notification/interface/CurrentG4Track.h"
#include "SimG4Core/Notification/interface/NewTrackAction.h"
#include "SimG4Core/Notification/interface/TrackInformation.h"
#include "SimG4Core/Notification/interface/TrackInformationExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4VProcess.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RegionStore.hh"
 
StackingAction::StackingAction(const edm::ParameterSet & p): tracker(0),
							     beam(0), calo(0),
							     muon(0) {
  trackNeutrino  = p.getParameter<bool>("TrackNeutrino");
  killHeavy      = p.getParameter<bool>("KillHeavy");
  kmaxIon        = p.getParameter<double>("IonThreshold")*MeV;
  kmaxProton     = p.getParameter<double>("ProtonThreshold")*MeV;
  kmaxNeutron    = p.getParameter<double>("NeutronThreshold")*MeV;
  maxTrackTime   = p.getParameter<double>("MaxTrackTime")*ns;
  maxTrackTimes  = p.getParameter<std::vector<double> >("MaxTrackTimes");
  maxTimeNames   = p.getParameter<std::vector<std::string> >("MaxTimeNames");
  savePDandCinTracker = p.getUntrackedParameter<bool>("SavePrimaryDecayProductsAndConversionsInTracker",false);
  savePDandCinCalo    = p.getUntrackedParameter<bool>("SavePrimaryDecayProductsAndConversionsInCalo",false);
  savePDandCinMuon    = p.getUntrackedParameter<bool>("SavePrimaryDecayProductsAndConversionsInMuon",false);
  saveFirstSecondary  = p.getUntrackedParameter<bool>("SaveFirstLevelSecondary",false);

  edm::LogInfo("SimG4CoreApplication") << "StackingAction initiated with"
				       << " flag for saving decay products in "
				       << " Tracker: " << savePDandCinTracker
                                       << " in Calo: " << savePDandCinCalo
                                       << " in Muon: " << savePDandCinMuon
				       << "\n               saveFirstSecondary"
				       << ": " << saveFirstSecondary
				       << " Flag for tracking neutrino: "
				       << trackNeutrino << " Killing Flag "
				       << killHeavy << " protons below " 
				       << kmaxProton <<" MeV, neutrons below "
				       << kmaxNeutron << " MeV and ions"
				       << " below " << kmaxIon << " MeV\n"
				       << "               kill tracks with "
				       << "time larger than " << maxTrackTime
				       << " ns";
  for (unsigned int i=0; i<maxTrackTimes.size(); i++) {
    maxTrackTimes[i] *= ns;
    edm::LogInfo("SimG4CoreApplication") << "SteppingAction::MaxTrackTime for "
					 << maxTimeNames[i] << " is " 
					 << maxTrackTimes[i];
  }
  initPointer();
}

StackingAction::~StackingAction() {}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track * aTrack) {

  // G4 interface part
  G4ClassificationOfNewTrack classification = fUrgent;
  int flag = 0;

  NewTrackAction newTA;
  if (aTrack->GetCreatorProcess()==0 || aTrack->GetParentID()==0) {
    newTA.primary(aTrack);
  } else if (aTrack->GetTouchable() == 0) {
    classification = fKill;
  } else {
    const G4Track * mother = CurrentG4Track::track();
    if ((savePDandCinTracker && (isThisVolume(aTrack->GetTouchable(),tracker)||
				 isThisVolume(aTrack->GetTouchable(),beam))) ||
	(savePDandCinCalo && isThisVolume(aTrack->GetTouchable(),calo)) ||
	(savePDandCinMuon && isThisVolume(aTrack->GetTouchable(),muon)))
      flag = isItPrimaryDecayProductOrConversion(aTrack, *mother);
    if (saveFirstSecondary) flag = isItFromPrimary(*mother, flag);
    newTA.secondary(aTrack, *mother, flag);

    if (aTrack->GetTrackStatus() == fStopAndKill) classification = fKill;
    if (killHeavy) {
      int    pdg = aTrack->GetDefinition()->GetPDGEncoding();
      double ke  = aTrack->GetKineticEnergy()/MeV;
      if (((pdg/1000000000 == 1) && (((pdg/10000)%100) > 0) && 
	   (((pdg/10)%100) > 0) && (ke<kmaxIon)) || 
	  ((pdg == 2212) && (ke < kmaxProton)) ||
	  ((pdg == 2112) && (ke < kmaxNeutron))) classification = fKill;
    }
    if (!trackNeutrino) {
      int    pdg = std::abs(aTrack->GetDefinition()->GetPDGEncoding());
      if (pdg == 12 || pdg == 14 || pdg == 16 || pdg == 18) 
	classification = fKill;
    }
    if (isItLongLived(aTrack)) classification = fKill;
//  LogDebug("SimG4CoreApplication") << "StackingAction:Classify Track "
//                                   << aTrack->GetTrackID() << " Parent " 
//                                   << aTrack->GetParentID() << " Type "
//                                   << aTrack->GetDefinition()->GetParticleName() 
//                                   << " K.E. " << aTrack->GetKineticEnergy()/MeV
//                                   << " MeV as " << classification 
//                                   << " Flag " << flag;
  }
  return classification;
}

void StackingAction::NewStage() {}

void StackingAction::PrepareNewEvent() {}

void StackingAction::initPointer() {

  const G4PhysicalVolumeStore * pvs = G4PhysicalVolumeStore::GetInstance();
  if (pvs) {
    std::vector<G4VPhysicalVolume*>::const_iterator pvcite;
    for (pvcite = pvs->begin(); pvcite != pvs->end(); pvcite++) {
      if (savePDandCinTracker) {
        if ((*pvcite)->GetName() == "Tracker") tracker = (*pvcite);
        if ((*pvcite)->GetName() == "Beam")    beam    = (*pvcite);
      }
      if (savePDandCinCalo) {
        if ((*pvcite)->GetName() == "CALO")    calo    = (*pvcite);
      }
      if (savePDandCinMuon) {
        if ((*pvcite)->GetName() == "MUON")    muon    = (*pvcite);
      }
      if ( (!savePDandCinTracker || (tracker && beam)) && 
	   (!savePDandCinCalo || calo) && (!savePDandCinMuon || muon ) ) break;
    }
    edm::LogInfo("SimG4CoreApplication") << "Pointers for Tracker " << tracker
                                         << ", BeamPipe " << beam << ", Calo " 
					 << calo << ", Muon " << muon;
    if (tracker) edm::LogInfo("SimG4CoreApplication") << "Tracker vol name "
						      << tracker->GetName();
    if (beam)    edm::LogInfo("SimG4CoreApplication") << "BeamPipe vol name "
						      << beam->GetName();
    if (calo)    edm::LogInfo("SimG4CoreApplication")<< "Calorimeter vol name "
						     << calo->GetName();
    if (muon)    edm::LogInfo("SimG4CoreApplication") << "Muon vol name "
						      << muon->GetName();
  }

  const G4RegionStore * rs = G4RegionStore::GetInstance();
  unsigned int num = maxTimeNames.size();
  if (num > 0) {
    std::vector<double> tofs;
    if (rs) {
      std::vector<G4Region*>::const_iterator rcite;
      for (rcite = rs->begin(); rcite != rs->end(); rcite++) {
	for (unsigned int i=0; i<num; i++) {
	  if ((*rcite)->GetName() == (G4String)(maxTimeNames[i])) {
	    maxTimeRegions.push_back(*rcite);
	    tofs.push_back(maxTrackTimes[i]);
	    break;
	  }
	}
	if (tofs.size() == num) break;
      }
    }
    for (unsigned int i=0; i<tofs.size(); i++) {
      maxTrackTimes[i] = tofs[i];
      G4String name = "Unknown";
      if (maxTimeRegions[i]) name = maxTimeRegions[i]->GetName();
      edm::LogInfo("SimG4CoreApplication") << name << " with pointer " 
					   << maxTimeRegions[i]<<" KE cut off "
					   << maxTrackTimes[i];
    }
  }

}

bool StackingAction::isThisVolume(const G4VTouchable* touch, 
				  G4VPhysicalVolume* pv) const {

  bool flag = false;
  if (pv != 0 && touch !=0) {
    int level = ((touch->GetHistoryDepth())+1);
    if (level >= 3) {
      int  ii = level - 3;
      flag    = (touch->GetVolume(ii) == pv);
    }
  }
  return flag;
}

int StackingAction::isItPrimaryDecayProductOrConversion(const G4Track * aTrack,
							const G4Track & mother) const {

  int flag = 0;
  TrackInformationExtractor extractor;
  const TrackInformation & motherInfo(extractor(mother));
  // Check whether mother is a primary
  if (motherInfo.isPrimary()) {
    if (aTrack->GetCreatorProcess()->GetProcessType() == fDecay &&
	aTrack->GetCreatorProcess()->GetProcessName() == "Decay") flag = 1;
    else if (aTrack->GetCreatorProcess()->GetProcessType() == fElectromagnetic &&
	     aTrack->GetCreatorProcess()->GetProcessName() == "conv") flag = 2;
  }
  return flag;
}

int StackingAction::isItFromPrimary(const G4Track & mother, int flagIn) const {

  int flag = flagIn;
  if (flag != 1) {
    TrackInformationExtractor extractor;
    const TrackInformation & motherInfo(extractor(mother));
    if (motherInfo.isPrimary()) flag = 3;
  }
  return flag;
}

bool StackingAction::isItLongLived(const G4Track * aTrack) const {

  bool   flag = false;
  double time = (aTrack->GetGlobalTime())/nanosecond;
  double tofM = maxTrackTime;
  if (maxTimeRegions.size() > 0) {
    G4Region* reg = aTrack->GetTouchable()->GetVolume(0)->GetLogicalVolume()->GetRegion();
    for (unsigned int i=0; i<maxTimeRegions.size(); i++) {
      if (reg == maxTimeRegions[i]) {
	tofM = maxTrackTimes[i];
	break;
      }
    }
  }
  if (time > tofM) flag = true;
  return flag;
}
