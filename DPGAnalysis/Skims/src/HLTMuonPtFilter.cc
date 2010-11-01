/** \file
 *
 * $Date: 2009/02/13 15:37:12 $
 * $Revision: 1.1 $
 * \author Silvia Goy Lopez - CERN <silvia.goy.lopez@cern.ch>
 */

/* This Class Header */
#include "DPGAnalysis/Skims/interface/HLTMuonPtFilter.h"

/* Collaborating Class Header */
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"




/* C++ Headers */
using namespace std;
using namespace edm;

/* ====================================================================== */

/// Constructor
HLTMuonPtFilter::HLTMuonPtFilter(const edm::ParameterSet& pset) {

  // the name of the STA rec hits collection
  theSTAMuonLabel = pset.getParameter<std::string>("SALabel");

  theMinPt = pset.getParameter<double>("minPt"); // pt min (GeV)

  LogDebug("HLTMuonPt") << " SALabel : " << theSTAMuonLabel 
    << " Min Pt : " << theMinPt;
}

/// Destructor
HLTMuonPtFilter::~HLTMuonPtFilter() {
}

/* Operations */ 
bool HLTMuonPtFilter::filter(edm::Event& event, const edm::EventSetup& eventSetup) {
  bool accept = false;

  // Get the RecTrack collection from the event
  Handle<reco::TrackCollection> staTracks;
  event.getByLabel(theSTAMuonLabel, staTracks);
  
  reco::TrackCollection::const_iterator staTrack;
  
  for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack){
    
    if(staTrack->pt()>theMinPt){
      accept=true;
      return accept;
    }

  }

  return accept;


}

// define this as a plug-in
DEFINE_FWK_MODULE(HLTMuonPtFilter);
