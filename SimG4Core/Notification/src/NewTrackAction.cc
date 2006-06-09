#include "SimG4Core/Notification/interface/NewTrackAction.h"
#include "SimG4Core/Notification/interface/TrackInformation.h"
#include "SimG4Core/Notification/interface/SimG4Exception.h"
#include "SimG4Core/Notification/interface/TrackInformationExtractor.h"

#include "G4Track.hh"

void NewTrackAction::primary(const G4Track * aTrack) const
{ primary(const_cast<G4Track *>(aTrack)); }

void NewTrackAction::primary(G4Track * aTrack) const
{ addUserInfoToPrimary(aTrack); }

void NewTrackAction::secondary(const G4Track * aSecondary,const G4Track & mother) const
{ secondary(const_cast<G4Track *>(aSecondary),mother); }

void NewTrackAction::secondary(G4Track * aSecondary,const G4Track & mother) const
{
    if (aSecondary->GetParentID() != mother.GetTrackID()) 
	throw SimG4Exception("NewTrackAction: secondary parent ID does not match mother id");
    TrackInformationExtractor extractor;
    const TrackInformation & motherInfo(extractor(mother));
    addUserInfoToSecondary(aSecondary,motherInfo);
}

void NewTrackAction::addUserInfoToPrimary(G4Track * aTrack) const
{
    TrackInformation * trkInfo = new TrackInformation();
    trkInfo->isPrimary(true);
    trkInfo->storeTrack(true);
    trkInfo->putInHistory();
    aTrack->SetUserInformation(trkInfo);  
}

void NewTrackAction::addUserInfoToSecondary(G4Track * aTrack,const TrackInformation & motherInfo) const
{
    TrackInformation * trkInfo = new TrackInformation();
    // transfer calo ID from mother (to be checked in TrackingAction)
    trkInfo->setIDonCaloSurface(motherInfo.getIDonCaloSurface());
    aTrack->SetUserInformation(trkInfo);  
}
