// \class JetTracksAssociatorAtCaloFace JetTracksAssociatorAtCaloFace.cc 
//
// Original Author:  Andrea Rizzi
//         Created:  Wed Apr 12 11:12:49 CEST 2006
// Accommodated for Jet Package by: Fedor Ratnikov Jul. 30, 2007
// $Id: JetTracksAssociatorAtCaloFace.cc,v 1.3 2008/05/29 17:58:55 fedor Exp $
//
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"

#include "JetTracksAssociatorAtCaloFace.h"

JetTracksAssociatorAtCaloFace::JetTracksAssociatorAtCaloFace(const edm::ParameterSet& fConfig)
  : mJets (fConfig.getParameter<edm::InputTag> ("jets")),
    mTracks (fConfig.getParameter<edm::InputTag> ("tracks")),
    mAssociator (fConfig.getParameter<double> ("coneSize"))
{
  reco::TrackBase::TrackQuality trackQuality = 
    reco::TrackBase::qualityByName (fConfig.getParameter<std::string> ("trackQuality"));
  if (trackQuality == reco::TrackBase::undefQuality) { // we have a problem
    edm::LogError("JetTracksAssociatorAtCaloFace") << "Unknown trackQuality value '" 
						   << fConfig.getParameter<std::string> ("trackQuality")
						   << "'. See possible values in 'reco::TrackBase::qualityByName'";
  }
  mTrackQuality = int (trackQuality);
  produces<reco::JetTracksAssociation::Container> ();
}

JetTracksAssociatorAtCaloFace::~JetTracksAssociatorAtCaloFace() {}

void JetTracksAssociatorAtCaloFace::produce(edm::Event& fEvent, const edm::EventSetup& fSetup) {
  // get stuff from Event Setup
  edm::ESHandle<MagneticField> field_h;
  fSetup.get<IdealMagneticFieldRecord>().get(field_h);
  edm::ESHandle<Propagator> propagator_h;
  fSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator_h);

  // get stuff from Event
  edm::Handle <edm::View <reco::Jet> > jets_h;
  fEvent.getByLabel (mJets, jets_h);
  edm::Handle <reco::TrackCollection> tracks_h;
  fEvent.getByLabel (mTracks, tracks_h);
  
  std::auto_ptr<reco::JetTracksAssociation::Container> jetTracks (new reco::JetTracksAssociation::Container (reco::JetRefBaseProd(jets_h)));

  if ( 1 ) { // original code 

    // format inputs
    std::vector <edm::RefToBase<reco::Jet> > allJets;
    allJets.reserve (jets_h->size());
    for (unsigned i = 0; i < jets_h->size(); ++i) allJets.push_back (jets_h->refAt(i));
    std::vector <reco::TrackRef> allTracks;
    allTracks.reserve (tracks_h->size());
    reco::TrackBase::TrackQuality trackQuality = reco::TrackBase::TrackQuality (mTrackQuality); // convert back
    for (unsigned i = 0; i < tracks_h->size(); ++i) {
      if ((*tracks_h)[i].quality (trackQuality)) allTracks.push_back (reco::TrackRef (tracks_h, i));
    }
    // run algo
    mAssociator.produce (&*jetTracks, allJets, allTracks, *field_h, *propagator_h);

  } else { // new code
    
    reco::TrackBase::TrackQuality quality = static_cast<reco::TrackBase::TrackQuality>( mTrackQuality );
    mAssociator.produce( &*jetTracks, jets_h, tracks_h, quality, *field_h, *propagator_h );
    
  }

  // store output
  fEvent.put (jetTracks);
}
