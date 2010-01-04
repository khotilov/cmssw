/*
 *  See header file for a description of this class.
 *
 *  $Date: 2009/12/14 15:07:55 $
 *  $Revision: 1.4 $
 *  \author S. Bolognesi, Erik - CERN
 */

#include "DQM/Physics/src/BPhysicsOniaDQM.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h" 
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <vector>

#include <string>
#include <cmath>
using namespace std;
using namespace edm;
using namespace reco;



BPhysicsOniaDQM::BPhysicsOniaDQM(const ParameterSet& parameters) {

  // the services
  theDbe = NULL;
  
  // Muon Collection Label
  theMuonCollectionLabel = parameters.getParameter<InputTag>("MuonCollection");
  global_background = NULL;
  diMuonMass_global = NULL;
  tracker_background = NULL;
  diMuonMass_tracker = NULL;
  standalone_background = NULL;
  diMuonMass_standalone = NULL;
    
  
}

BPhysicsOniaDQM::~BPhysicsOniaDQM() { 
  

}


void BPhysicsOniaDQM::beginJob() {
 
  metname = "oniaAnalyzer";

  LogTrace(metname)<<"[BPhysicsOniaDQM] Parameters initialization";
  theDbe = Service<DQMStore>().operator->();
  if(theDbe!=NULL){
    theDbe->setCurrentFolder("Physics/BPhysics");  // Use folder with name of PAG
    global_background = theDbe->book1D("global_background", "Same-sign global-global dimuon mass", 750, 0, 15);
    diMuonMass_global = theDbe->book1D("diMuonMass_global", "Opposite-sign global-global dimuon mass", 750, 0, 15);
    tracker_background = theDbe->book1D("tracker_background", "Same-sign tracker-tracker (arbitrated) dimuon mass", 750, 0, 15);
    diMuonMass_tracker = theDbe->book1D("diMuonMass_tracker", "Opposite-sign tracker-tracker (arbitrated) dimuon mass", 750, 0, 15);
    standalone_background = theDbe->book1D("standalone_background", "Same-sign standalone-standalone dimuon mass", 500, 0, 15);
    diMuonMass_standalone = theDbe->book1D("diMuonMass_standalone", "Opposite-sign standalone-standalone dimuon mass", 500, 0, 15);
  }

}


void BPhysicsOniaDQM::analyze(const Event& iEvent, const EventSetup& iSetup) {

  LogTrace(metname)<<"[BPhysicsOniaDQM] Analysis of event # ";
  
  // Take the STA muon container
  Handle<MuonCollection> muons;
  iEvent.getByLabel(theMuonCollectionLabel,muons);

  if(muons.isValid()){
    for (MuonCollection::const_iterator recoMu1 = muons->begin(); recoMu1!=muons->end(); ++recoMu1){
      
      // only loop over the remaining muons if recoMu1 is one of the following
      if(recoMu1->isGlobalMuon() || recoMu1->isTrackerMuon() || recoMu1->isStandAloneMuon()){
	for (MuonCollection::const_iterator recoMu2 = recoMu1+1; recoMu2!=muons->end(); ++recoMu2){
	  
	  // fill the relevant histograms if recoMu2 satisfies one of the following
	  if (recoMu1->isGlobalMuon() && recoMu2->isGlobalMuon()){
	    math::XYZVector vec1 = recoMu1->globalTrack()->momentum();
	    math::XYZVector vec2 = recoMu2->globalTrack()->momentum();
	    float massJPsi = computeMass(vec1,vec2);
	    
	    
	    // if opposite charges, fill _global, else fill _background
	    if (((*recoMu1).charge()*(*recoMu2).charge())<0) {
	      if(diMuonMass_global!=NULL){
		diMuonMass_global->Fill(massJPsi);
	      }
	    } else {
	      if(global_background!=NULL){
		global_background->Fill (massJPsi);
	      }
	    }

	  }
	  if(recoMu1->isStandAloneMuon() && recoMu2->isStandAloneMuon() &&
	     fabs(recoMu1->outerTrack()->d0()) < 5 && fabs(recoMu1->outerTrack()->dz()) < 30 &&
	     fabs(recoMu2->outerTrack()->d0()) < 5 && fabs(recoMu2->outerTrack()->dz()) < 30){
	    math::XYZVector vec1 = recoMu1->outerTrack()->momentum();
	    math::XYZVector vec2 = recoMu2->outerTrack()->momentum();
	    float massJPsi = computeMass(vec1,vec2);
	    
	    // if opposite charges, fill _standalone, else fill _background
	    if (((*recoMu1).charge()*(*recoMu2).charge())<0) {
	      if(diMuonMass_standalone!=NULL){
		diMuonMass_standalone->Fill(massJPsi);
	      }
	    } else {
	      if(standalone_background!=NULL){
		standalone_background->Fill (massJPsi);
	      }
	    }
	  }
	  if(recoMu1->isTrackerMuon() && recoMu2->isTrackerMuon() &&
	     muon::isGoodMuon(*recoMu1, muon::TrackerMuonArbitrated) &&
	     muon::isGoodMuon(*recoMu2, muon::TrackerMuonArbitrated)){
	    math::XYZVector vec1 = recoMu1->innerTrack()->momentum();
	    math::XYZVector vec2 = recoMu2->innerTrack()->momentum();
	    float massJPsi = computeMass(vec1,vec2);
	    
	    // if opposite charges, fill _tracker, else fill _background
	    if (((*recoMu1).charge()*(*recoMu2).charge())<0) {
	      if(diMuonMass_tracker!=NULL){
		diMuonMass_tracker->Fill(massJPsi);
	      }
	    } else {
	      if(tracker_background!=NULL){
		tracker_background->Fill (massJPsi);
	      }
	    }
	  }
	}
      }
    } 
  }    
}




void BPhysicsOniaDQM::endJob(void) {
  LogTrace(metname)<<"[BPhysicsOniaDQM] EndJob";
}


float BPhysicsOniaDQM::computeMass(const math::XYZVector &vec1,const math::XYZVector &vec2){
  // mass of muon
  float massMu = 0.10566;
  float eMu1 = -999;
  if(massMu*massMu + vec1.Mag2()>0)
    eMu1 = sqrt(massMu*massMu + vec1.Mag2());
  float eMu2 = -999;
  if(massMu*massMu + vec2.Mag2()>0)
    eMu2 = sqrt(massMu*massMu + vec2.Mag2());

  float pJPsi = -999;
  if((vec1+vec2).Mag2()>0)
    pJPsi = sqrt((vec1+vec2).Mag2());
  float eJPsi = eMu1 + eMu2;

  float massJPsi = -999;
  if((eJPsi*eJPsi - pJPsi*pJPsi) > 0)
    massJPsi = sqrt(eJPsi*eJPsi - pJPsi*pJPsi);
 
 return massJPsi;
}



