#ifndef BeamSpotProducer_BeamSpotAnalyzer_h
#define BeamSpotProducer_BeamSpotAnalyzer_h

/**_________________________________________________________________
   class:   BeamSpotAnalyzer.h
   package: RecoVertex/BeamSpotProducer
   


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: BeamSpotAnalyzer.h,v 1.5 2008/05/14 15:31:57 yumiceva Exp $

________________________________________________________________**/


// C++ standard
#include <string>
// CMS
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoVertex/BeamSpotProducer/interface/BSTrkParameters.h"
#include "RecoVertex/BeamSpotProducer/interface/BeamFitter.h"

// ROOT
#include "TFile.h"
#include "TTree.h"

#include<fstream>

class BeamSpotAnalyzer : public edm::EDAnalyzer {
 public:
  explicit BeamSpotAnalyzer(const edm::ParameterSet&);
  ~BeamSpotAnalyzer();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  std::string outputfilename_;
  std::string fasciifilename;
  std::ofstream fasciiFile;
  std::string fasciiFileName;

  TFile* file_;
  TTree* ftree_;

  int    ftotalevents;
  double ftheta;
  double fpt;
  double feta;
  int    fcharge;
  double fchi2;
  double fndof;
  double fphi0;
  double fd0;
  double fsigmad0;
  double fz0;
  double fsigmaz0;
  unsigned int fnHit;
  unsigned int fnStripHit;
  unsigned int fnPixelHit;
  unsigned int fnTIBHit;
  unsigned int fnTOBHit;
  unsigned int fnTIDHit;
  unsigned int fnTECHit;
  unsigned int fnPXBHit;
  unsigned int fnPXFHit;
  double fd0phi_chi2;
  double fd0phi_d0;
  double fcov[7][7];
  
  std::vector< BSTrkParameters > fBSvector;
  
  bool write2DB_;
  bool runallfitters_;
  double inputBeamWidth_;
  int ftotal_tracks;

  BeamFitter * theBeamFitter;
};

#endif
