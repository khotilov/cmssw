#ifndef  HFDumpVertexH
#define  HFDumpVertexH
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TFile;
class TTree;
class TAna00Event;


// ----------------------------------------------------------------------
class HFDumpVertex : public edm::EDAnalyzer {
 public:
  explicit HFDumpVertex(const edm::ParameterSet&);
  ~HFDumpVertex();
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();


  int           fVerbose;
  std::string   fVertexLabel;
  std::string   fVertexTracksLabel;
  std::string   fSimVertexLabel;
  
  int nevt;

};

#endif
