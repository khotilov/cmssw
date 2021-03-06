
// -*- C++ -*-
//
// Package:    TauAna
// Class:      TauAna
// 
/**\class TauAna

Description: To study MC Tau properties

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Dongwook Jang
//         Created:  Wed Oct 11 11:08:40 CDT 2006
// $Id: TauAna.cc,v 1.3 2008/02/16 06:32:50 dwjang Exp $
//
//

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TauReco/interface/CaloTauTagInfo.h"
#include "RecoTauTag/Pi0Tau/interface/Pi0.h"
#include "RecoTauTag/Pi0Tau/interface/Pi0Fwd.h"
#include "RecoTauTag/Pi0Tau/interface/Tau3D.h"
#include "RecoTauTag/Pi0Tau/interface/Tau3DFwd.h"
#include "RecoTauTag/Pi0Tau/interface/Tau3DAlgo.h"
#include "RecoTauTag/Pi0Tau/interface/TauVariables.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "RecoTauTag/Pi0Tau/test/TauAna.h"

// ROOT related includes
#include <TFile.h>
#include <TVector2.h>

#include <iostream>

using namespace std;
using namespace reco;

//
// constructors and destructor
//
TauAna::TauAna(const edm::ParameterSet& iConfig)
{

  trackCollectionName_ = iConfig.getParameter<string>("trackCollectionName");
  tauCollectionName_ = iConfig.getParameter<string>("tauCollectionName");
  pFCandidateProducerName_ = iConfig.getParameter<string>("pFCandidateProducerName");
  pFCandidateCollectionName_ = iConfig.getParameter<string>("pFCandidateCollectionName");

  histFileName_ = iConfig.getParameter<string>("histFileName");

}


TauAna::~TauAna() {}

//
// member functions
//

// ------------ method called to for each event  ------------
void TauAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  t_nSignalTracks    = 0;
  t_nSignalPi0s      = 0;
  t_nIsolationTracks = 0;
  t_nIsolationPi0s   = 0;
  t_tracksMomentum->SetXYZT(0,0,0,0);
  t_pi0sMomentum->SetXYZT(0,0,0,0);
  t_momentum->SetXYZT(0,0,0,0);
  t_seedTrackVertex->SetXYZ(0,0,0);

  edm::Handle<reco::TrackCollection> trackHandle;
  iEvent.getByLabel(trackCollectionName_,trackHandle);
  if(!trackHandle.isValid()) cout << "trackHandle is not valid" << endl;
  LogDebug("TauAna") << "size of default track collection : " << trackHandle->size() << endl;

  edm::Handle<reco::PFCandidateCollection> pFCandidateHandle;
  iEvent.getByLabel(pFCandidateProducerName_,pFCandidateCollectionName_,pFCandidateHandle);
  if(!pFCandidateHandle.isValid()) cout << "pFCandidateHandle is not valid" << endl;
  LogDebug("TauAna") << "size of PFCandidate collection : " << pFCandidateHandle->size() << endl;


  edm::Handle<reco::CaloTauTagInfoCollection> tauTagInfoHandle;
  iEvent.getByLabel(tauCollectionName_, tauTagInfoHandle);
  if(!tauTagInfoHandle.isValid()) cout << "CaloTauTagInfoCollection is not valid" << endl;
  LogDebug("TauAna") << "size of tauTag collection : " << tauTagInfoHandle->size() << endl;


  Tau3DAlgo tau3DAlgo(&trackHandle);
  tau3DAlgo.setUse3DAngle(true);
  tau3DAlgo.setSeedTrackThreshold(5.0);
  tau3DAlgo.setTauOuterConeSize(0.524);
  tau3DAlgo.fillTau3Ds(pFCandidateHandle);
  LogDebug("TauAna") << "size of tau3D collection : " << tau3DAlgo.tau3DCollection().size() << endl;

  vector<TauVariables> tauVars;
  for(reco::Tau3DCollection::const_iterator iter = tau3DAlgo.tau3DCollection().begin();
      iter != tau3DAlgo.tau3DCollection().end(); iter++){
    const Tau3D *tau3D = &*iter;

    TauVariables tauVar(tau3D,&tauTagInfoHandle);
    tauVar.setUse3DAngle(false);
    tauVar.setSignalConeSize(0.15);
    tauVar.setIsolationConeSize(0.5);
    tauVar.setUseVariableSignalCone(true);
    tauVar.setSignalConeFunction(5.0);
    tauVar.setUseVariableIsolationCone(false);
    tauVar.setIsolationConeFunction(5.0); // do not affect anything
    tauVar.setSeedTrackThreshold(5.0);
    tauVar.setShoulderTrackThreshold(1.0);
    tauVar.setPi0Threshold(1.0);
    tauVar.setDZTrackAssociation(2.0);
    tauVar.makeVariables(); // perform calculation
    tauVars.push_back(tauVar);

    LogDebug("TauAna") << "nSignalTracks, nSignalPi0s, nIsolationTracks, nIsolationPi0s : "
		       << tauVar.nSignalTracks() << ", "
		       << tauVar.nSignalPi0s() << ", "
		       << tauVar.nIsolationTracks() << ", "
		       << tauVar.nIsolationPi0s() << "\n";

    t_nSignalTracks    = tauVar.nSignalTracks(); 
    t_nSignalPi0s      = tauVar.nSignalPi0s(); 
    t_nIsolationTracks = tauVar.nIsolationTracks(); 
    t_nIsolationPi0s   = tauVar.nIsolationPi0s(); 
    t_tracksMomentum->SetXYZT(tauVar.tracksMomentum().X(), 
                              tauVar.tracksMomentum().Y(), 
                              tauVar.tracksMomentum().Z(), 
                              tauVar.tracksMomentum().E()); 
    t_pi0sMomentum->SetXYZT(tauVar.pi0sMomentum().X(), 
                            tauVar.pi0sMomentum().Y(), 
                            tauVar.pi0sMomentum().Z(), 
                            tauVar.pi0sMomentum().E()); 
    t_momentum->SetXYZT(tauVar.momentum().X(), 
                        tauVar.momentum().Y(), 
                        tauVar.momentum().Z(), 
                        tauVar.momentum().E()); 
    t_seedTrackVertex->SetXYZ(tauVar.seedTrack()->vertex().x(),
			      tauVar.seedTrack()->vertex().y(),
			      tauVar.seedTrack()->vertex().z());

 
    tree_->Fill(); 
  }

}// analyze


// ------------ method called once each job just before starting event loop  ------------
void TauAna::beginJob(const edm::EventSetup&)
{

  t_tracksMomentum  = new TLorentzVector(0.0,0.0,0.0,0.0);
  t_pi0sMomentum    = new TLorentzVector(0.0,0.0,0.0,0.0);
  t_momentum        = new TLorentzVector(0.0,0.0,0.0,0.0);
  t_seedTrackVertex = new TVector3(0.0,0.0,0.0);

  TFile *f = new TFile(histFileName_.c_str(),"RECREATE");
  if(f){
    tree_ = new TTree("tree","tree");
    tree_->SetAutoSave();
    tree_->Branch("nSignalTracks",&t_nSignalTracks,"nSignalTracks/I");
    tree_->Branch("nSignalPi0s",&t_nSignalPi0s,"nSignalPi0s/I");
    tree_->Branch("nIsolationTracks",&t_nIsolationTracks,"nIsolationTracks/I");
    tree_->Branch("nIsolationPi0s",&t_nIsolationPi0s,"nIsolationPi0s/I");
    //    tree_->Branch("seedTrackVertex","TVector3",&t_seedTrackVertex,32000,99);
    tree_->Branch("tracksMomentum.","TLorentzVector",&t_tracksMomentum,32000,99);
    tree_->Branch("pi0sMomentum.","TLorentzVector",&t_pi0sMomentum,32000,99);
    tree_->Branch("momentum.","TLorentzVector",&t_momentum,32000,99);
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void TauAna::endJob() {

  TFile *f = tree_->GetCurrentFile();
  if(f){
    f->cd();
    f->Write();
    f->Close();
  }

}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//define this as a plug-in
DEFINE_FWK_MODULE(TauAna);
