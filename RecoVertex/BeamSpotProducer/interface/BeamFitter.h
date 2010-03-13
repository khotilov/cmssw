#ifndef BeamFitter_H
#define BeamFitter_H

/**_________________________________________________________________
   class:   BeamFitter.h
   package: RecoVertex/BeamSpotProducer
   


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
         Geng-Yuan Jeng, UC Riverside (Geng-Yuan.Jeng@cern.ch)
 
 version $Id: BeamFitter.h,v 1.18 2010/03/12 21:45:35 yumiceva Exp $

 ________________________________________________________________**/

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "RecoVertex/BeamSpotProducer/interface/BSTrkParameters.h"
#include "RecoVertex/BeamSpotProducer/interface/BSFitter.h"
#include "RecoVertex/BeamSpotProducer/interface/PVFitter.h"

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include <fstream>

class BeamFitter {
 public:
  BeamFitter() {}
  BeamFitter(const edm::ParameterSet& iConfig);
  virtual ~BeamFitter();

  void readEvent(const edm::Event& iEvent);

  bool runFitter(); 
  void runAllFitter();
  void resetTrkVector() { fBSvector.clear(); }
  void resetTotTrk() { ftotal_tracks=0; }
  void resetLSRange() { fbeginLumiOfFit=fendLumiOfFit=-1; }
  void dumpTxtFile(std::string &,bool);
  void write2DB();
  reco::BeamSpot getBeamSpot() { return fbeamspot; }
  std::vector<BSTrkParameters> getBSvector() { return fBSvector; }
  TH1F * getCutFlow() { return h1cutFlow; }
  void resetCutFlow() {
    h1cutFlow->Reset();
    for (unsigned int i=0; i<sizeof(countPass)/sizeof(countPass[0]); i++) countPass[i]=0;
  };
  int* getFitLSRange() {
    int *tmp=new int[2];
    tmp[0] = fbeginLumiOfFit;
    tmp[1] = fendLumiOfFit;
    return tmp;
  }
  
 private:

  std::vector<BSTrkParameters> fBSvector;
  reco::BeamSpot fbeamspot;
  BSFitter *fmyalgo;
  std::ofstream fasciiFile;
  std::ofstream fasciiFileBak;

  bool debug_;
  edm::InputTag tracksLabel_;
  bool writeTxt_;
  bool writeTxtBak_;
  std::string outputTxt_;
  std::string outputTxtBak_;
  double trk_MinpT_;
  double trk_MaxZ_;
  double trk_MaxEta_;
  double trk_MaxIP_;
  int trk_MinNTotLayers_;
  int trk_MinNPixLayers_;
  double trk_MaxNormChi2_;
  std::vector<std::string> trk_Algorithm_;
  std::vector<std::string> trk_Quality_;
  std::vector<reco::TrackBase::TrackQuality> quality_;
  std::vector<reco::TrackBase::TrackAlgorithm> algorithm_;
  double inputBeamWidth_;
  double convergence_;
  int ftotal_tracks;
  int min_Ntrks_;
  bool isMuon_;
  
  // ntuple
  TH1F* h1z;
  bool saveNtuple_;
  bool saveBeamFit_;
  std::string outputfilename_;
  TFile* file_;
  TTree* ftree_;
  double ftheta;
  double fpt;
  double feta;
  int    fcharge;
  double fnormchi2;
  double fphi0;
  double fd0;
  double fd0bs;
  double fsigmad0;
  double fz0;
  double fsigmaz0;
  int fnTotLayerMeas;
  int fnPixelLayerMeas;
  int fnStripLayerMeas;
  int fnTIBLayerMeas;
  int fnTIDLayerMeas;
  int fnTOBLayerMeas;
  int fnTECLayerMeas;
  int fnPXBLayerMeas;
  int fnPXFLayerMeas;
  double fd0phi_chi2;
  double fd0phi_d0;
  double fcov[7][7];
  double fvx;
  double fvy;
  int frun;
  int flumi;
  bool fquality;
  bool falgo;
  bool fpvValid;
  double fpvx, fpvy, fpvz;
  
  //beam fit results
  TTree* ftreeFit_;
  int frunFit;
  int fbeginLumiOfFit;
  int fendLumiOfFit;
  char fbeginTimeOfFit[30];
  char fendTimeOfFit[30];
  double fx;
  double fy;
  double fz;
  double fsigmaZ;
  double fdxdz;
  double fdydz;
  double fxErr;
  double fyErr;
  double fzErr;
  double fsigmaZErr;
  double fdxdzErr;
  double fdydzErr;
  double fwidthX;
  double fwidthY;
  double fwidthXErr;
  double fwidthYErr;

  TH1F *h1ntrks;
  TH1F *h1vz_event;
  TH1F *h1cutFlow;
  int countPass[9];

  PVFitter *MyPVFitter;

};

#endif
