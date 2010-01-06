#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/WeightContainer.h"

#include "TH1D.h"
#include "TFile.h"

class MuMuAnalyzer : public edm::EDAnalyzer {
 public:
  explicit MuMuAnalyzer(const edm::ParameterSet&);
  ~MuMuAnalyzer();
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
 private:  
  std::string outputFilename;
  TH1D* weight_histo; int event;
  TH1F* MuMu_invmass_histo; TH1F* J1Eta_histo; TH1F* J2Eta_histo;
  TH1F* Pt_histo; TH1F* J1Pt_histo; TH1F* J2Pt_histo; TH1F* EJDelR_histo; 
  TH1F* E1Pt_histo; TH1F* E2Pt_histo; TH1F* ZPz_histo; TH1F* ZPt_histo;
  TH1F* JDelR_histo; TH1F* JDelPhi_histo; TH1F* J1Phi_histo; TH1F* J2Phi_histo;
};
