// -*- C++ -*-
//
// Package:    PrimaryVertexAnalyzer
// Class:      PrimaryVertexAnalyzer
// 
/**\class PrimaryVertexAnalyzer PrimaryVertexAnalyzer.cc Validation/RecoVertex/src/PrimaryVertexAnalyzer.cc

 Description: simple primary vertex analyzer

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wolfram Erdmann
//         Created:  Fri Jun  2 10:54:05 CEST 2006
// $Id: PrimaryVertexAnalyzer.h,v 1.8 2006/12/21 16:33:30 werdmann Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
 
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//generator level + CLHEP
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CLHEP/HepMC/GenEvent.h"
#include "CLHEP/HepMC/GenVertex.h"
#include "CLHEP/HepMC/GenParticle.h"
 
// vertex stuff
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// simulated vertices,..., add <use name=SimDataFormats/Vertex> and <../Track>
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

// Root
#include <TH1.h>
#include <TFile.h>




// class declaration
class PrimaryVertexAnalyzer : public edm::EDAnalyzer {


// auxiliary class holding simulated primary vertices
class simPrimaryVertex {
public:
  simPrimaryVertex(double x1,double y1,double z1):x(x1),y(y1),z(z1),ptsq(0),nGenTrk(0){};
  double x,y,z;
  HepLorentzVector ptot;
  double ptsq;
  int nGenTrk;
  std::vector<int> finalstateParticles;
  std::vector<int> simTrackIndex;
  std::vector<int> genVertex;
  const reco::Vertex *recVtx;
};



public:
  explicit PrimaryVertexAnalyzer(const edm::ParameterSet&);
  ~PrimaryVertexAnalyzer();
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob(edm::EventSetup const&);
  virtual void endJob();

private:

  bool matchVertex(const simPrimaryVertex  &vsim, 
		   const reco::Vertex       &vrec);
  bool isResonance(const HepMC::GenParticle * p);
  bool isFinalstateParticle(const HepMC::GenParticle * p);
  bool isCharged(const HepMC::GenParticle * p);
 
  void printRecVtxs(const edm::Handle<reco::VertexCollection> recVtxs);
  void printSimVtxs(const edm::Handle<edm::SimVertexContainer> simVtxs);
  void printSimTrks(const edm::Handle<edm::SimTrackContainer> simVtrks);
  std::vector<simPrimaryVertex> getSimPVs(const edm::Handle<edm::HepMCProduct> evtMC);
  std::vector<simPrimaryVertex> getSimPVs(const edm::Handle<edm::HepMCProduct> evt, 
					  const edm::Handle<edm::SimVertexContainer> simVtxs, 
					  const edm::Handle<edm::SimTrackContainer> simTrks);
  // ----------member data ---------------------------
  std::string recoTrackProducer_;
  std::string outputFile_;       // output file
  std::string vtxSample_;        // which vertices to analyze
  TFile*  rootFile_;             
  bool verbose_;
  edm::InputTag simG4_;
  double simUnit_;     
  edm::ESHandle < ParticleDataTable > pdt;
     

  std::map<std::string, TH1*> h;
};

