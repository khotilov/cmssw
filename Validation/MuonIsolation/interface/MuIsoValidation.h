// -*- C++ -*-
// MuIsoValidation.h
// Package:    Muon Isolation Validation
// Class:      MuIsoValidation
// 
/*

 Description: Muon Isolation Validation class

 NOTE: The static member variable declerations *should* include the key word "static", but 
	 I haven't found an elegant way to initalize the vectors.  Static primatives (e.g. int, 
	 float, ...) and simple static objects are easy to initialze.  Outside of the class 
	 decleration, you would write
	
 		int MuIsoValidation::CONST_INT = 5;
 		FooType MuIsoValidation::CONST_FOOT = Foo(constructor_argument);
	
	 but you can't do this if you want to, say, initalize a std::vector with a bunch of 
	 different values.  So, you can't make them static and you have to initialize them using 
	 a member method.  To keep it consistent, I've just initialized them all in the same 
	 method, even the simple types.
 
*/
//
// Original Author:  "C. Jess Riedel", UC Santa Barbara
//         Created:  Tue Jul 17 15:58:24 CDT 2007
//

//Base class
#include "FWCore/Framework/interface/EDAnalyzer.h"

//Member types
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//Other include files
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"

//----------------------------------------

//Forward declarations
class TH1;
class TH1I;
class TH1D;
class TH2;
class TProfile;

//------------------------------------------
//  Class Declaration: MuIsoValidation
//--------------------------------------
class MuIsoValidation : public edm::EDAnalyzer {
  //---------namespace and typedefs--------------
  typedef edm::View<reco::Muon>::const_iterator MuonIterator;
  typedef edm::RefToBase<reco::Muon> MuonBaseRef;
  typedef edm::Handle<reco::IsoDepositMap> MuIsoDepHandle;
  typedef const reco::IsoDeposit MuIsoDepRef;
  
public:
  //---------methods----------------------------
  explicit MuIsoValidation(const edm::ParameterSet&);
  ~MuIsoValidation();
  
  
private:
  //---------methods----------------------------
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void InitStatics();
  void RecordData(MuonIterator muon, MuIsoDepRef& tkDep,  
		  MuIsoDepRef& ecalDep, MuIsoDepRef& hcalDep, 
		  MuIsoDepRef& hoDep);//Fills Histograms with info from single muon
  void InitHistos();//adds title, bin information to member histograms
  void MakeLogBinsForProfile(Double_t* bin_edges, const double min, const double max);
  void FillHistos();//Fills histograms with data
  void NormalizeHistos(); //Normalize to number of muons
  TH1* GetTH1FromMonitorElement(MonitorElement* me);
  TH2* GetTH2FromMonitorElement(MonitorElement* me);
  TProfile* GetTProfileFromMonitorElement(MonitorElement* me);
  
  
  //----------Static Variables---------------
  
  //Collection labels
  edm::InputTag Muon_Tag;
  edm::InputTag tkIsoDeposit_Tag;
  edm::InputTag hcalIsoDeposit_Tag;
  edm::InputTag ecalIsoDeposit_Tag;
  edm::InputTag hoIsoDeposit_Tag;
  
  //root file name
  std::string rootfilename;
  // Directories within the rootfile
  std::string dirName;

  //Histogram parameters
  static const int NUM_VARS = 11;
  double L_BIN_WIDTH;//large bins
  double S_BIN_WIDTH;//small bins
  int LOG_BINNING_ENABLED;//pseudo log binning for profile plots
  int NUM_LOG_BINS;
  double LOG_BINNING_RATIO;
  bool requireCombinedMuon;
  
  std::string title_sam;
  std::string title_cone;
  std::string title_cd;
  
  std::vector<std::string> main_titles;//[NUM_VARS]
  std::vector<std::string> axis_titles;//[NUM_VARS]
  std::vector<std::string> names;//[NUM_VARS]
  std::vector< std::vector<double> > param;//[NUM_VARS][3]
  std::vector<int> isContinuous;//[NUM_VARs]
  
  //---------------Dynamic Variables---------------------
  
  //MonitorElement
  DQMStore* dbe;
  
  //The Data
  int theMuonData;//[number of muons]
  double theData[NUM_VARS];
  
  //Histograms
  MonitorElement* h_nMuons;
  std::vector<MonitorElement*> h_1D;//[NUM_VARS]
  std::vector<MonitorElement*> cd_plots;//[NUM_VARS]
  std::vector< std::vector<MonitorElement*> > h_2D;//[NUM_VARS][NUM_VARS]
  std::vector< std::vector<MonitorElement*> > p_2D;//[NUM_VARS][NUM_VARS]
  
  //Counters
  int nEvents;
  int nIncMuons;
  int nCombinedMuons;
  
  //enums for monitorElement
  enum {NOAXIS,XAXIS,YAXIS,ZAXIS};
};

