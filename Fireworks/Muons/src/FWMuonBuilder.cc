// -*- C++ -*-
//
// Package:     Muons
// Class  :     FWMuonBuilder
// $Id: FWMuonBuilder.cc,v 1.20 2010/04/20 20:49:44 amraktad Exp $
//

// system include files
#include "TEveTrackPropagator.h"
#include "TEveVSDStructs.h"
#include "TEveManager.h"
#include "TEveTrack.h"
#include "TEveStraightLineSet.h"
#include "TEveGeoNode.h"

// user include files
#include "Fireworks/Core/interface/FWEventItem.h"
#include "Fireworks/Core/interface/FWMagField.h"
#include "Fireworks/Core/interface/FWProxyBuilderBase.h"
#include "Fireworks/Core/interface/DetIdToMatrix.h"
#include "Fireworks/Candidates/interface/CandidateUtils.h"
#include "Fireworks/Tracks/interface/TrackUtils.h"
#include "Fireworks/Muons/interface/FWMuonBuilder.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"


namespace  {
std::vector<TEveVector> getRecoTrajectoryPoints( const reco::Muon* muon,
                                                 const FWEventItem* iItem )
{
   std::vector<TEveVector> points;
   const DetIdToMatrix* geom = iItem->getGeom();
   const std::vector<reco::MuonChamberMatch>& matches = muon->matches();
   Double_t localTrajectoryPoint[3];
   Double_t globalTrajectoryPoint[3];
   std::vector<reco::MuonChamberMatch>::const_iterator chamber = matches.begin();
   for ( ; chamber != matches.end(); ++chamber )
   {
      // expected track position
      localTrajectoryPoint[0] = chamber->x;
      localTrajectoryPoint[1] = chamber->y;
      localTrajectoryPoint[2] = 0;

      DetId id = chamber->id;
      const TGeoHMatrix* matrix = geom->getMatrix( chamber->id.rawId() );
      if ( matrix ) {
         matrix->LocalToMaster( localTrajectoryPoint, globalTrajectoryPoint );
         points.push_back(TEveVector(globalTrajectoryPoint[0],
                                     globalTrajectoryPoint[1],
                                     globalTrajectoryPoint[2]));
      }
   }
   return points;
}

//______________________________________________________________________________

void addMatchInformation( const reco::Muon* muon,
                          FWProxyBuilderBase* pb,
                          TEveElement* parentList,
                          bool showEndcap)
{
   std::set<unsigned int> ids;
   const DetIdToMatrix* geom = pb->context().getGeom();
   const std::vector<reco::MuonChamberMatch>& matches = muon->matches();
   //need to use auto_ptr since the segmentSet may not be passed to muonList
   std::auto_ptr<TEveStraightLineSet> segmentSet(new TEveStraightLineSet);
   segmentSet->SetLineWidth(4);
   std::vector<reco::MuonChamberMatch>::const_iterator chamber = matches.begin();
   for ( ; chamber != matches.end(); ++chamber )
   {
      DetId id = chamber->id;
      if ( ids.insert(id.rawId()).second &&  // ensure that we add same chamber only once
           ( id.subdetId() != MuonSubdetId::CSC || showEndcap ) ){
         TEveGeoShape* shape = geom->getShape( chamber->id.rawId() );
         if(0!=shape) {
            shape->RefMainTrans().Scale(0.999, 0.999, 0.999);
            shape->SetMainTransparency(65);
            pb->setupAddElement(shape, parentList);
         }
      }
      const TGeoHMatrix* matrix = geom->getMatrix( chamber->id.rawId() );
      if( matrix ) {
         // make muon segment 20 cm long along local z-axis
         // add segments
	 for( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->segmentMatches.begin(),
							       segmentEnd = chamber->segmentMatches.end();
	      segment != segmentEnd; ++segment )
         {
	    const double segmentLength = 15;

	    Double_t localSegmentInnerPoint[3];
	    Double_t localSegmentCenterPoint[3];
	    Double_t localSegmentOuterPoint[3];
	    Double_t globalSegmentInnerPoint[3];
	    Double_t globalSegmentCenterPoint[3];
	    Double_t globalSegmentOuterPoint[3];

	    localSegmentOuterPoint[0] = segment->x + segmentLength*segment->dXdZ;
	    localSegmentOuterPoint[1] = segment->y + segmentLength*segment->dYdZ;
	    localSegmentOuterPoint[2] = segmentLength;

	    localSegmentCenterPoint[0] = segment->x;
	    localSegmentCenterPoint[1] = segment->y;
	    localSegmentCenterPoint[2] = 0;

	    localSegmentInnerPoint[0] = segment->x - segmentLength*segment->dXdZ;
	    localSegmentInnerPoint[1] = segment->y - segmentLength*segment->dYdZ;
	    localSegmentInnerPoint[2] = -segmentLength;

	    matrix->LocalToMaster( localSegmentInnerPoint, globalSegmentInnerPoint );
	    matrix->LocalToMaster( localSegmentCenterPoint, globalSegmentCenterPoint );
	    matrix->LocalToMaster( localSegmentOuterPoint, globalSegmentOuterPoint );

	    if( globalSegmentInnerPoint[1]*globalSegmentOuterPoint[1] > 0 ) {
	       segmentSet->AddLine( globalSegmentInnerPoint[0], globalSegmentInnerPoint[1], globalSegmentInnerPoint[2],
				    globalSegmentOuterPoint[0], globalSegmentOuterPoint[1], globalSegmentOuterPoint[2] );
	    } else {
	       if( fabs(globalSegmentInnerPoint[1]) > fabs(globalSegmentOuterPoint[1]) )
		  segmentSet->AddLine( globalSegmentInnerPoint[0], globalSegmentInnerPoint[1], globalSegmentInnerPoint[2],
				       globalSegmentCenterPoint[0], globalSegmentCenterPoint[1], globalSegmentCenterPoint[2] );
	       else
		  segmentSet->AddLine( globalSegmentCenterPoint[0], globalSegmentCenterPoint[1], globalSegmentCenterPoint[2],
				       globalSegmentOuterPoint[0], globalSegmentOuterPoint[1], globalSegmentOuterPoint[2] );
	    }
         }
      }
   }
   if ( !matches.empty() ) pb->setupAddElement( segmentSet.release(), parentList );
}

//______________________________________________________________________________

bool
buggyMuon( const reco::Muon* muon,
           const DetIdToMatrix* geom )
{
   if (!muon->standAloneMuon().isAvailable() ||
       !muon->standAloneMuon()->extra().isAvailable() ) return false;
   const std::vector<reco::MuonChamberMatch>& matches = muon->matches();
   Double_t localTrajectoryPoint[3];
   Double_t globalTrajectoryPoint[3];
   std::vector<reco::MuonChamberMatch>::const_iterator chamber = matches.begin();
   for ( ; chamber != matches.end(); ++chamber )
   {
      // expected track position
      localTrajectoryPoint[0] = chamber->x;
      localTrajectoryPoint[1] = chamber->y;
      localTrajectoryPoint[2] = 0;

      DetId id = chamber->id;
      const TGeoHMatrix* matrix = geom->getMatrix( chamber->id.rawId() );
      if ( matrix ) {
         matrix->LocalToMaster( localTrajectoryPoint, globalTrajectoryPoint );
         double phi = atan2(globalTrajectoryPoint[1],globalTrajectoryPoint[0]);
         if ( cos( phi - muon->standAloneMuon()->innerPosition().phi()) < 0 )
            return true;
      }
   }
   return false;
}

}

//
// constructors and destructor
//
FWMuonBuilder::FWMuonBuilder()
{
   //NOTE: We call IncRefCount and IncDenyDestroy since TEveTrackPropagator actually has two reference counts being done on it
   // We only want the one using IncRefCount to actually cause the deletion which is why 'IncDenyDestroy' does not have a matching
   // DecDenyDestroy.  I'm still using a edm::FWEvePtr to hold the Propagator since I want to know if the propagator is deleted
   m_trackerPropagator.reset(new TEveTrackPropagator()); // propagate within tracker
   m_trackerPropagator->IncRefCount();
   m_trackerPropagator->IncDenyDestroy();
   m_trackerPropagator->SetMaxR( 850 );
   m_trackerPropagator->SetMaxZ( 1100 );
   m_trackerPropagator->SetDelta(0.05);
}

FWMuonBuilder::~FWMuonBuilder()
{
   m_trackerPropagator->DecRefCount();
}

//
// member functions
//
//______________________________________________________________________________

void
FWMuonBuilder::calculateField(const reco::Muon& iData, FWMagField* field)
{

   // if auto field estimation mode, do extra loop over muons.
   // use both inner and outer track if available
   if ( field->getAutodetect() ) {
     if ( fabs( iData.eta() ) > 2.0 || iData.pt() < 3 ) return;
     if ( iData.innerTrack().isAvailable() ){
       double estimate = fireworks::estimateField(*(iData.innerTrack()),true);
       if ( estimate >= 0 ) field->guessField( estimate );
	 
     }
     if ( iData.outerTrack().isAvailable() ){
       double estimate = fireworks::estimateField(*(iData.outerTrack()));
       if ( estimate >= 0 ) field->guessFieldIsOn( estimate > 0.5 );
     }
   }
}

//______________________________________________________________________________

void
FWMuonBuilder::buildMuon(FWProxyBuilderBase* pb,
                         const reco::Muon* muon,
                         TEveElement* tList,
                         bool showEndcap,
                         bool tracksOnly)
{
   calculateField(*muon, pb->context().getField());

   // workaround for missing GetFieldObj() in TEveTrackPropagator, default stepper is kHelix
   if (m_trackerPropagator->GetStepper() == TEveTrackPropagator::kHelix) {
      m_trackerPropagator->SetStepper(TEveTrackPropagator::kRungeKutta);
      m_trackerPropagator->SetMagFieldObj(pb->context().getField());
   }

   TEveRecTrack recTrack;
   recTrack.fBeta = 1.;

   // If we deal with a tracker muon we use the inner track and guide it
   // through the trajectory points from the reconstruction. Segments
   // represent hits. Matching between hits and the trajectory shows
   // how well the inner track matches with the muon hypothesis.
   //
   // In other cases we use a global muon track with a few states from 
   // the inner and outer tracks or just the outer track if it's the
   // only option

   if ( muon->isTrackerMuon() && 
	muon->innerTrack().isAvailable() &&
	muon->isMatchesValid() &&
	!buggyMuon( &*muon, pb->context().getGeom() ) )
   {
      TEveTrack* trk = fireworks::prepareTrack(*(muon->innerTrack()),
					       m_trackerPropagator.get(),
					       getRecoTrajectoryPoints(muon,pb->item()) );
      trk->MakeTrack();
      pb->setupAddElement(trk, tList);
      if ( ! tracksOnly )
	 addMatchInformation( &(*muon), pb, tList, showEndcap );
      return;
   } 

   if ( muon->isGlobalMuon() &&
	muon->globalTrack().isAvailable() )
   {
      std::vector<TEveVector> extraPoints;
      if ( muon->innerTrack().isAvailable() ){
	 extraPoints.push_back( TEveVector(muon->innerTrack()->innerPosition().x(),
					   muon->innerTrack()->innerPosition().y(),
					   muon->innerTrack()->innerPosition().z()) );
	 extraPoints.push_back( TEveVector(muon->innerTrack()->outerPosition().x(),
					   muon->innerTrack()->outerPosition().y(),
					   muon->innerTrack()->outerPosition().z()) );
      }
      if ( muon->outerTrack().isAvailable() ){
	 extraPoints.push_back( TEveVector(muon->outerTrack()->innerPosition().x(),
					   muon->outerTrack()->innerPosition().y(),
					   muon->outerTrack()->innerPosition().z()) );
	 extraPoints.push_back( TEveVector(muon->outerTrack()->outerPosition().x(),
					   muon->outerTrack()->outerPosition().y(),
					   muon->outerTrack()->outerPosition().z()) );
      }
      TEveTrack* trk = fireworks::prepareTrack(*(muon->globalTrack()),
                                               m_trackerPropagator.get(),
                                               extraPoints);
      trk->MakeTrack();
      pb->setupAddElement(trk, tList);
      return;
   }

   if ( muon->innerTrack().isAvailable() )
   {
      TEveTrack* trk = fireworks::prepareTrack(*(muon->innerTrack()),
                                               m_trackerPropagator.get());
      trk->MakeTrack();
      pb->setupAddElement(trk, tList);
      return;
   }

   if ( muon->outerTrack().isAvailable() )
   {
      TEveTrack* trk = fireworks::prepareTrack(*(muon->outerTrack()),
                                               m_trackerPropagator.get());
      trk->MakeTrack();
      pb->setupAddElement(trk, tList);
      return;
   }
   
   // if got that far it means we have nothing but a candidate
   // show it anyway.
   TEveTrack* trk = fireworks::prepareCandidate(*muon,
						m_trackerPropagator.get());
   trk->MakeTrack();
   pb->setupAddElement(trk, tList);
}
