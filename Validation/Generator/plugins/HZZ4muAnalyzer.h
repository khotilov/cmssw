#include "FWCore/Framework/interface/EDAnalyzer.h"

class TFile;
class TH1D;

class HZZ4muAnalyzer : public edm::EDAnalyzer
{
 public:
  //
  explicit HZZ4muAnalyzer( const edm::ParameterSet& ) ;
  virtual ~HZZ4muAnalyzer() {} // no need to delete ROOT stuff
  // as it'll be deleted upon closing TFile
  
  virtual void analyze( const edm::Event&, const edm::EventSetup& ) ;
  virtual void beginJob( const edm::EventSetup& ) ;
  virtual void endJob() ;
  
 private:
  //
  std::string fOutputFileName ;
  TFile*      fOutputFile ;
  TH1D*       fHist2muMass ;
  TH1D*       fHist4muMass ;
  TH1D*       fHistZZMass ;
  
};
