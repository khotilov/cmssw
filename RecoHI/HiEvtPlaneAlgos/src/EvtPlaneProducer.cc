// -*- C++ -*-
//
// Package:    EvtPlaneProducer
// Class:      EvtPlaneProducer
// 
/**\class EvtPlaneProducer EvtPlaneProducer.cc RecoHI/EvtPlaneProducer/src/EvtPlaneProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Sergey Petrushanko
//         Created:  Fri Jul 11 10:05:00 2008
// $Id: EvtPlaneProducer.cc,v 1.6 2009/09/08 10:52:01 edwenger Exp $
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Ref.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include <iostream>

using namespace std;

//
// class decleration
//

class EvtPlaneProducer : public edm::EDProducer {
   public:
      explicit EvtPlaneProducer(const edm::ParameterSet&);
      ~EvtPlaneProducer();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------

  bool useECAL_;
  bool useHCAL_;
  bool useTrackMidEta_;
  bool useTrackPosEta_;
  bool useTrackNegEta_;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
EvtPlaneProducer::EvtPlaneProducer(const edm::ParameterSet& iConfig)
{
   //register your products
  
  useECAL_ = iConfig.getUntrackedParameter<bool>("useECAL",true);
  useHCAL_ = iConfig.getUntrackedParameter<bool>("useHCAL",true);
  useTrackMidEta_ = iConfig.getUntrackedParameter<bool>("useTrackMidEta",true);
  useTrackPosEta_ = iConfig.getUntrackedParameter<bool>("useTrackPosEta",true);
  useTrackNegEta_ = iConfig.getUntrackedParameter<bool>("useTrackNegEta",true);




  produces<reco::EvtPlaneCollection>("recoLevel");

}


EvtPlaneProducer::~EvtPlaneProducer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
EvtPlaneProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace HepMC;

  //calorimetry part

      double ugol[9], ugol2[9];
      double tower_eta, tower_phi, tower_energy, tower_energy_e, tower_energy_h;
      double s1t, s2t, s1e, s2e, s1h, s2h;
      double s1[9], s2[9], TEnergy[144], TPhi[144];
      double pi = 3.14159;
      int numb;	

//      double planeA     =  0.;

//       cout << endl << "  Start of the event plane determination." << endl;

       for(int j=0;j<9;j++) {
        s1[j]  = 0.;
        s2[j]  = 0.;
       }
      
       for(int l=0;l<144;l++) {
        TEnergy[l]  = 0.;
        TPhi[l]  = 0.;
       }

      Handle<CaloTowerCollection> calotower;
      iEvent.getByLabel("towerMaker",calotower);
      
      if(!calotower.isValid()){
        cout << "Error! Can't get calotower product!" << endl;
       return ;
      }

	for (CaloTowerCollection::const_iterator j = calotower->begin();j !=calotower->end(); j++) {

//        cout << *j << std::endl;
//        cout << "ENERGY HAD " << j->hadEnergy()<< " ENERGY EM " <<j->emEnergy() 
//	  << " ETA " <<j->eta() << " PHI " <<j->phi() << std::endl;
	     
        tower_eta        = j->eta();
        tower_phi        = j->phi();
	tower_energy_e   = j->emEnergy();
	tower_energy_h   = j->hadEnergy();
        tower_energy     = tower_energy_e + tower_energy_h;
	
	s1t = tower_energy*sin(2.*tower_phi-pi);
	s2t = tower_energy*cos(2.*tower_phi-pi);
        s1e = tower_energy_e*sin(2.*tower_phi-pi);
	s2e = tower_energy_e*cos(2.*tower_phi-pi);
	s1h = tower_energy_h*sin(2.*tower_phi-pi);
	s2h = tower_energy_h*cos(2.*tower_phi-pi);

	 if (fabs(tower_eta)<3.){

	  numb = static_cast< int >(72.*(tower_phi/pi + 1.) - 0.5);
	  TEnergy[numb] += tower_energy;
	  TPhi[numb]     = tower_phi;

// barrel + endcap
  	  s1[0] +=  s1t;
	  s2[0] +=  s2t;
  	  s1[3] +=  s1h;
	  s2[3] +=  s2h;
  	  s1[6] +=  s1e;
	  s2[6] +=  s2e;
	  
// endcap
	 if (fabs(tower_eta)>1.5) {
  	  s1[2] +=  s1t;
	  s2[2] +=  s2t;
  	  s1[5] +=  s1h;
	  s2[5] +=  s2h;
  	  s1[8] +=  s1e;
	  s2[8] +=  s2e;
	 }
	 }
	 
// barrel
	 if (fabs(tower_eta)<1.5){
  	  s1[1] +=  s1t;
	  s2[1] +=  s2t;
  	  s1[4] +=  s1h;
	  s2[4] +=  s2h;
  	  s1[7] +=  s1e;
	  s2[7] +=  s2e;
	 }
	}
	
      for(int j1=0;j1<9;j1++) {
 
       if (s2[j1]==0.) {ugol[j1]=0.;}
       else {ugol[j1] = 0.5*atan(s1[j1]/s2[j1]);}
       
       if ( s2[j1] < 0 && s1[j1] <  0) ugol[j1] = ugol[j1] - pi/2.;
       if ( s2[j1] < 0 && s1[j1] >= 0) ugol[j1] = ugol[j1] + pi/2.;
       
       ugol2[j1] = ugol[j1] + pi/2.;
       if (ugol2[j1]>pi/2.) ugol2[j1] = ugol2[j1] - pi;
       
      }

/*       
       cout <<  endl << "   Azimuthal angle of reaction plane (with minimum)" << endl
       << "HCAL+ECAL (b+e)   " << ugol[0] << endl
       << "HCAL+ECAL (b)     " << ugol[1] << endl
       << "HCAL+ECAL (e)     " << ugol[2] << endl
       << "HCAL      (b+e)   " << ugol[3] << endl
       << "HCAL      (b)     " << ugol[4] << endl
       << "HCAL      (e)     " << ugol[5] << endl
       << "ECAL      (b+e)   " << ugol[6] << endl
       << "ECAL      (b)     " << ugol[7] << endl
       << "ECAL      (e)     " << ugol[8] << endl;

       cout <<  endl << "   Azimuthal angle of reaction plane (with maximum)" << endl
       << "HCAL+ECAL (b+e)   " << ugol2[0] << endl
       << "HCAL+ECAL (b)     " << ugol2[1] << endl
       << "HCAL+ECAL (e)     " << ugol2[2] << endl
       << "HCAL      (b+e)   " << ugol2[3] << endl
       << "HCAL      (b)     " << ugol2[4] << endl
       << "HCAL      (e)     " << ugol2[5] << endl
       << "ECAL      (b+e)   " << ugol2[6] << endl
       << "ECAL      (b)     " << ugol2[7] << endl
       << "ECAL      (e)     " << ugol2[8] << endl  << endl;

*/


//Tracking part

   double track_eta;
   double track_phi;
   double track_pt;

   double trackPsi_eta_mid;
   double trackPsi_eta_pos;
   double trackPsi_eta_neg;


   Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel("hiSelectedTracks", tracks);

   // cout << " TRACKS Size " << tracks->size() << endl;

   if(!tracks.isValid()){
     cout << "Error! Can't get selectTracks!" << endl;
     return ;
   }
   double trackSin_eta_mid = 0;
   double trackCos_eta_mid = 0;
   double trackSin_eta_pos = 0;
   double trackCos_eta_pos = 0;
   double trackSin_eta_neg = 0;
   double trackCos_eta_neg = 0;

   for(reco::TrackCollection::const_iterator j = tracks->begin(); j != tracks->end(); j++){
     track_eta = j->eta();
     track_phi = j->phi();
     track_pt = j->pt();
    
     if(fabs(track_eta)<0.75){
       trackSin_eta_mid+=sin(2*track_phi);
       trackCos_eta_mid+=cos(2*track_phi);
     }
   
     if((track_eta >= 0.75) && (track_eta < 2.0)){
       trackSin_eta_pos+=sin(2*track_phi);
       trackCos_eta_pos+=cos(2*track_phi);
     }
     if((track_eta <= -0.75) && (track_eta > -2.0)){
       trackSin_eta_neg+=sin(2*track_phi);
       trackCos_eta_neg+=cos(2*track_phi);
     }

    }


   trackPsi_eta_mid = 0.5*atan2(trackSin_eta_mid,trackCos_eta_mid);
   trackPsi_eta_pos = 0.5*atan2(trackSin_eta_pos,trackCos_eta_pos);
   trackPsi_eta_neg = 0.5*atan2(trackSin_eta_neg,trackCos_eta_neg);




   std::auto_ptr<EvtPlaneCollection> evtplaneOutput(new EvtPlaneCollection);
      

   EvtPlane ecalPlane(ugol2[6],s1[6],s2[6],"Ecal");
   EvtPlane hcalPlane(ugol2[3],s1[3],s2[3],"Hcal");
   EvtPlane caloPlane(ugol2[0],s1[0],s2[0],"Calo");

   EvtPlane EvtPlaneFromTracksMidEta(trackPsi_eta_mid,trackSin_eta_mid,trackCos_eta_mid,"EvtPlaneFromTracksMidEta");
   EvtPlane EvtPlaneFromTracksPosEta(trackPsi_eta_pos,trackSin_eta_pos,trackCos_eta_pos,"EvtPlaneFromTracksPosEta");
   EvtPlane EvtPlaneFromTracksNegEta(trackPsi_eta_neg,trackSin_eta_neg,trackCos_eta_neg,"EvtPlaneFromTracksNegEta");
  
   if(useTrackMidEta_) evtplaneOutput->push_back(EvtPlaneFromTracksMidEta);
   if(useTrackPosEta_) evtplaneOutput->push_back(EvtPlaneFromTracksPosEta);
   if(useTrackNegEta_) evtplaneOutput->push_back(EvtPlaneFromTracksNegEta);



   
   if(useECAL_) evtplaneOutput->push_back(ecalPlane);
   if(useHCAL_) evtplaneOutput->push_back(hcalPlane);
   if(useECAL_ && useHCAL_) evtplaneOutput->push_back(caloPlane);


   iEvent.put(evtplaneOutput, "recoLevel");


// cout << "  "<< planeA << endl;
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
EvtPlaneProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EvtPlaneProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EvtPlaneProducer);
