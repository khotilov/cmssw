#ifndef MuScleFit_H
#define MuScleFit_H

/** \class MuScleFit
 *  Analyzer of the Global muon tracks
 *
 *  $Date: 2009/11/10 11:15:29 $
 *  $Revision: 1.24 $
 *  \author C.Mariotti, S.Bolognesi - INFN Torino / T.Dorigo - INFN Padova
 */

// Base Class Headers
// ------------------
#include "FWCore/Framework/interface/EDLooper.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include <CLHEP/Vector/LorentzVector.h>
#include <vector>
// #include "MuonAnalysis/MomentumScaleCalibration/interface/Histograms.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "MuonAnalysis/MomentumScaleCalibration/interface/MuScleFitBase.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
// The following is required in CMSSW v2.0.x (was contained in Muon.h in 1.6.7)
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "DataFormats/PatCandidates/interface/Muon.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class TH1F;
class TH2F;
class TProfile;
class MuonServiceProxy;
class TTree;
class MuScleFitPlotter;

class MuScleFit: public edm::EDLooper, MuScleFitBase {

 public:
  // Constructor
  // -----------
  MuScleFit( const edm::ParameterSet& pset );

  // Destructor
  // ----------
  virtual ~MuScleFit();

  // Operations
  // ----------
  virtual void beginOfJob();
  virtual void endOfJob();

  virtual void startingNewLoop( unsigned int iLoop );

  virtual edm::EDLooper::Status endOfLoop( const edm::EventSetup& eventSetup, unsigned int iLoop );
  virtual void endOfFastLoop( const unsigned int iLoop );

  virtual edm::EDLooper::Status duringLoop( const edm::Event & event, const edm::EventSetup& eventSetup );
  /**
   * This method performs all needed operations on the muon pair. It reads the muons from SavedPair and uses the iev
   * counter to keep track of the event number. The iev is incremented internally and reset to 0 in startingNewLoop.
   */
  virtual edm::EDLooper::Status duringFastLoop();

  template<typename T>
  std::vector<reco::LeafCandidate> fillMuonCollection( const std::vector<T>& tracks );
 private:

 protected:
  /// Template method used to fill the track collection starting from reco::muons or pat::muons
  template<typename T>
  void takeSelectedMuonType(const T & muon, vector<reco::Track> & tracks);

  /// Check if two lorentzVector are near in deltaR
  bool checkDeltaR( reco::Particle::LorentzVector& genMu, reco::Particle::LorentzVector& recMu );
  /// Fill the reco vs gen and reco vs sim comparison histograms
  void fillComparisonHistograms( const reco::Particle::LorentzVector & genMu, const reco::Particle::LorentzVector & recoMu, const string & inputName, const int charge );

  /**
   * Simple method to check parameters consistency. It aborts the job if the parameters
   * are not consistent.
   */
  void checkParameters();

  MuonServiceProxy *theService;

  // Counters
  // --------
  int numberOfSimTracks;
  int numberOfSimMuons;
  int numberOfSimVertices;
  int numberOfEwkZ;

  bool ifHepMC;
  bool ifGenPart;

  // Constants
  // ---------
  double minResMass_hwindow[6];
  double maxResMass_hwindow[6];

  // Total number of loops
  // ---------------------
  unsigned int maxLoopNumber;
  unsigned int loopCounter;

  bool fastLoop;

  MuScleFitPlotter *plotter;

  // The reconstructed muon 4-momenta to be put in the tree
  // ------------------------------------------------------
  reco::Particle::LorentzVector recMu1, recMu2;
  int iev;
  int totalEvents_;

  bool compareToSimTracks_;
  edm::InputTag simTracksCollection_;
  bool PATmuons_;
  string genParticlesName_;
};

template<typename T>
std::vector<reco::LeafCandidate> MuScleFit::fillMuonCollection( const std::vector<T>& tracks )
{
  std::vector<reco::LeafCandidate> muons;
  typename std::vector<T>::const_iterator track;
  for( track = tracks.begin(); track != tracks.end(); ++track ) {
    reco::Particle::LorentzVector mu;
    mu = reco::Particle::LorentzVector(track->px(),track->py(),track->pz(),
                                       sqrt(track->p()*track->p() + MuScleFitUtils::mMu2));
    // Apply smearing if needed, and then bias
    // ---------------------------------------
    MuScleFitUtils::goodmuon++;
    if (debug_>0) cout <<setprecision(9)<< "Muon #" << MuScleFitUtils::goodmuon
                       << ": initial value   Pt = " << mu.Pt() << endl;
    if (MuScleFitUtils::SmearType>0) {
      mu = MuScleFitUtils::applySmearing( mu );
      if (debug_>0) cout << "Muon #" << MuScleFitUtils::goodmuon
                         << ": after smearing  Pt = " << mu.Pt() << endl;
    }
    if (MuScleFitUtils::BiasType>0) {
      mu = MuScleFitUtils::applyBias( mu, track->charge() );
      if (debug_>0) cout << "Muon #" << MuScleFitUtils::goodmuon
                         << ": after bias      Pt = " << mu.Pt() << endl;
    }
    reco::LeafCandidate muon(track->charge(),mu);
    // Store modified muon
    // -------------------
    muons.push_back (muon);
  }
  return muons;
}

template<typename T>
void MuScleFit::takeSelectedMuonType(const T & muon, vector<reco::Track> & tracks)
{
  // cout<<"muon "<<muon->isGlobalMuon()<<muon->isStandAloneMuon()<<muon->isTrackerMuon()<<endl;
  //NNBB: one muon can be of many kinds at once but with the theMuonType_ we are sure
  // to avoid double counting of the same muon
  if(muon->isGlobalMuon() && theMuonType_==1)
    tracks.push_back(*(muon->globalTrack()));
  else if(muon->isGlobalMuon() && theMuonType_==10) //particular case!!
    tracks.push_back(*(muon->innerTrack()));
  else if(muon->isStandAloneMuon() && theMuonType_==2)
    tracks.push_back(*(muon->outerTrack()));
  else if(muon->isTrackerMuon() && theMuonType_==3)
    tracks.push_back(*(muon->innerTrack()));
}

#endif
