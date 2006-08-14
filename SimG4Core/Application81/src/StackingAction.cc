#include "SimG4Core/Application81/interface/StackingAction.h"
#include "SimG4Core/Notification/interface/CurrentG4Track.h"
#include "SimG4Core/Notification/interface/NewTrackAction.h"
 
StackingAction::StackingAction() {}

StackingAction::~StackingAction() {}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
    NewTrackAction newTA;
    if (aTrack->GetCreatorProcess()==0 || aTrack->GetParentID()==0)
        newTA.primary(aTrack);
    else
    {
        const G4Track * mother = CurrentG4Track::track();
        newTA.secondary(aTrack, *mother);
    }
    // G4 interface part
    G4ClassificationOfNewTrack classification = fUrgent;
    return classification;
}

void StackingAction::NewStage() {}

void StackingAction::PrepareNewEvent() {}


