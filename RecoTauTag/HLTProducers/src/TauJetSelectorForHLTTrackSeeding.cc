#include "RecoTauTag/HLTProducers/interface/TauJetSelectorForHLTTrackSeeding.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


TauJetSelectorForHLTTrackSeeding::TauJetSelectorForHLTTrackSeeding(const edm::ParameterSet& iConfig):
  inputTrackJetTag_(iConfig.getParameter< edm::InputTag > ("inputTrackJetTag")),
  inputCaloJetTag_(iConfig.getParameter< edm::InputTag > ("inputCaloJetTag")),
  inputTrackTag_(iConfig.getParameter< edm::InputTag > ("inputTrackTag")),
  ptMinCaloJet_(iConfig.getParameter< double > ("ptMinCaloJet")),
  etaMinCaloJet_(iConfig.getParameter< double > ("etaMinCaloJet")),
  etaMaxCaloJet_(iConfig.getParameter< double > ("etaMaxCaloJet")),
  tauConeSize_(iConfig.getParameter< double > ("tauConeSize")),
  isolationConeSize_(iConfig.getParameter< double > ("isolationConeSize")),
  fractionMinCaloInTauCone_(iConfig.getParameter< double > ("fractionMinCaloInTauCone")),
  fractionMaxChargedPUInCaloCone_(iConfig.getParameter< double > ("fractionMaxChargedPUInCaloCone")),
  ptTrkMaxInCaloCone_(iConfig.getParameter< double > ("ptTrkMaxInCaloCone")),
  nTrkMaxInCaloCone_(iConfig.getParameter< int > ("nTrkMaxInCaloCone"))
{
   //now do what ever initialization is needed
  produces<reco::TrackJetCollection>();
}


TauJetSelectorForHLTTrackSeeding::~TauJetSelectorForHLTTrackSeeding()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
void
TauJetSelectorForHLTTrackSeeding::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::auto_ptr< reco::TrackJetCollection > augmentedTrackJets (new reco::TrackJetCollection);

   edm::Handle<reco::TrackJetCollection> trackjets;
   iEvent.getByLabel(inputTrackJetTag_,trackjets);

   for (reco::TrackJetCollection::const_iterator trackjet = trackjets->begin();
	  trackjet != trackjets->end(); trackjet++) {
     augmentedTrackJets->push_back(*trackjet);
   }


   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(inputTrackTag_,tracks);

   edm::Handle<reco::CaloJetCollection> calojets;
   iEvent.getByLabel(inputCaloJetTag_,calojets);

   for (reco::CaloJetCollection::const_iterator calojet = calojets->begin();
	  calojet != calojets->end(); calojet++) {

     if ( calojet->pt() < ptMinCaloJet_ ) continue;
     double etaJet = calojet->eta();
     double phiJet = calojet->phi();
     if ( etaJet < etaMinCaloJet_ ) continue;
     if ( etaJet > etaMaxCaloJet_ ) continue;

     std::vector <CaloTowerPtr> theTowers = calojet->getCaloConstituents();
     double ptIn = 0.;
     double ptOut = 0.;
     for ( unsigned int itwr = 0; itwr < theTowers.size(); ++itwr ) { 
       double etaTwr = theTowers[itwr]->eta() - etaJet;
       double phiTwr = theTowers[itwr]->phi() - phiJet;
       double deltaR = sqrt( etaTwr*etaTwr + phiTwr*phiTwr );
       //std::cout << "Tower eta/phi/et : " << etaTwr << " " << phiTwr << " " << theTowers[itwr]->pt() << std::endl;
       if ( deltaR < tauConeSize_ ) { 
	 ptIn += theTowers[itwr]->pt(); 
       } else if ( deltaR < isolationConeSize_ ) { 
	 ptOut += theTowers[itwr]->pt(); 
       }
     }
     double ptTot = ptIn+ptOut;
     double fracIn = ptIn/ptTot;

     // We are looking for isolated tracks
     if ( fracIn < fractionMinCaloInTauCone_) continue;

     int ntrk = 0;
     double ptTrk = 0.;

     for (reco::TrackJetCollection::const_iterator trackjet = trackjets->begin();
	  trackjet != trackjets->end(); trackjet++) {
       for (unsigned itr=0; itr<trackjet->numberOfTracks(); ++itr) { 
	 edm::Ptr<reco::Track> track = trackjet->track(itr);
	 double trackEta = track->eta() - etaJet;
	 double trackPhi = track->phi() - phiJet;
	 double deltaR = sqrt( trackEta*trackEta + trackPhi*trackPhi );
	 if ( deltaR < isolationConeSize_ ) { 
	   ntrk++; 
	   ptTrk += track->pt();
	 }
       }
     }
     // We are looking for calojets without signal tracks already in
     if ( ntrk > nTrkMaxInCaloCone_ ) continue;
     if ( ptTrk > ptTrkMaxInCaloCone_ ) continue;

     int ntrk2 = 0;
     double ptTrk2 = 0.;

     for (reco::TrackCollection::const_iterator track = tracks->begin();
	  track != tracks->end(); track++) {
       double trackEta = track->eta() - etaJet;
       double trackPhi = track->phi() - phiJet;
       double deltaR = sqrt( trackEta*trackEta + trackPhi*trackPhi );
       if ( deltaR < isolationConeSize_ ) { 
	 ntrk2++; 
	 ptTrk2 += track->pt();
       }
     }
     // We are looking for signal jets, not PU jets
     double fractionChargedPU = ptTrk2/calojet->pt(); 
     if ( fractionChargedPU > fractionMaxChargedPUInCaloCone_ ) continue;
     /*
     std::cout << "Calo Jet " << calojet->pt() << " " << calojet->eta() 
	       << " " << ptIn << " " << ptOut << " " << fracIn
	       << " " << ptTrk << " " << ntrk 
	       << " " << fractionChargedPU
	       << std::endl;
     */
     math::XYZTLorentzVector p4(calojet->p4());
     math::XYZPoint vertex(calojet->vertex());
     augmentedTrackJets->push_back(reco::TrackJet(p4,vertex));
   }

   iEvent.put(augmentedTrackJets);
}

// ------------ method called once each job just before starting event loop  ------------
void 
TauJetSelectorForHLTTrackSeeding::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauJetSelectorForHLTTrackSeeding::endJob() {}

// ------------ method called when starting to processes a run  ------------
void 
TauJetSelectorForHLTTrackSeeding::beginRun(edm::Run&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void 
TauJetSelectorForHLTTrackSeeding::endRun(edm::Run&, edm::EventSetup const&) {}

