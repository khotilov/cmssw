// -*- C++ -*-
//
// Package:    JPsiAnalyzerPAT
// Class:      JPsiAnalyzerPAT
// 
/**\class JPsiAnalyzerPAT JPsiAnalyzerPAT.cc OctoberXTracking/JPsiAnalyzerPAT/src/JPsiAnalyzerPAT.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author: Roberto Covarelli 
//         Created:  Fri Oct  9 04:59:40 PDT 2009
// $Id: JPsiAnalyzerPAT.cc,v 1.50 2011/04/01 16:11:02 covarell Exp $
//
// based on: Onia2MuMu package V00-11-00
// changes done by: FT-HW

// system include files
#include <memory>
#include <fstream>
#include <ostream>
#include <iostream>
#include <math.h>

// ROOT/Roofit include files
#include <TStyle.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooCategory.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include <DataFormats/Math/interface/deltaR.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/BeamSpot/interface/BeamSpot.h>
#include <DataFormats/Math/interface/deltaR.h>
#include <MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h>

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace RooFit;

//
// class declaration
//
class JPsiAnalyzerPAT : public edm::EDAnalyzer {
   public:
      explicit JPsiAnalyzerPAT(const edm::ParameterSet&);
      ~JPsiAnalyzerPAT();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void makeCuts(int sign) ;
      pair< unsigned int, const pat::CompositeCandidate* > theBestQQ(int sign);
      void fillTreeAndDS(unsigned int theCat, const pat::CompositeCandidate* aCand, const edm::Event&);
      bool isMuonInAccept(const pat::Muon* aMuon);
      bool selGlobalMuon(const pat::Muon* aMuon);
      bool selTrackerMuon(const pat::Muon* aMuon);
      bool selCaloMuon(const pat::Muon* aMuon);
      bool selDimuon(const pat::CompositeCandidate* aCand);
      int getJpsiVarType(const double jpsivar, vector<double> vectbin);
      double CorrectMass(const reco::Muon& mu1,const reco::Muon& mu2, int mode);

      // additional functions by f
      void resetDSVariables();
      void analyzeGenerator(const edm::Handle<reco::GenParticleCollection>& genParticles);
      // void calcPol(TLorentzVector&, TLorentzVector&, std::vector< float >&, std::vector< float >& );
      void beginRun(const edm::Run &, const edm::EventSetup &);
      void hltReport(const edm::Event &iEvent ,const edm::EventSetup& iSetup);
      void matchMuonToHlt(const pat::Muon*, const pat::Muon*);
      void muonStationDistance (const pat::CompositeCandidate* aCand);

      // ROOT tree 
      TTree* tree_;//data; //*recoData;
      TFile* fOut_;

      // SMALL dataset and RooRealVars
      TFile* fOut2_;
      RooDataSet* data;
      RooRealVar* Jpsi_Mass;      
      RooRealVar* Jpsi_Pt;
      RooRealVar* Jpsi_Rap; 
      RooRealVar* Jpsi_ct;
      RooRealVar* Jpsi_ctErr;
      RooRealVar* Jpsi_ctTrue;      			
      RooCategory* Jpsi_Type;
      RooCategory* Jpsi_PtType;
      RooCategory* Jpsi_RapType;
      RooCategory* Jpsi_MatchType;
      RooCategory* Jpsi_Sign;   

      //1.) J/psi variables RECO
      // double JpsiMass, JpsiPt, JpsiRap;
      // double JpsiPx, JpsiPy, JpsiPz;
      TLorentzVector* JpsiP;
      double Jpsict, JpsictErr, JpsiVprob;
      int JpsiType,  JpsiCharge, MCType; //GG, GT and TT
      double JpsiDistM1, JpsiDphiM1, JpsiDrM1;
      double JpsiDistM2, JpsiDphiM2, JpsiDrM2;

      //2.) muon variables RECO
      // double muPosPx, muPosPy, muPosPz;
      TLorentzVector* muPosP;
      // double muNegPx, muNegPy, muNegPz;
      TLorentzVector* muNegP;

      //3.) J/psi variables GEN
      // double JpsiMass_Gen, JpsiPt_Gen, JpsiRap_Gen;
      // double JpsiPx_Gen,   JpsiPy_Gen, JpsiPz_Gen;
      TLorentzVector* JpsiP_Gen;
      double Jpsict_Gen;

      //4.)muon variables GEN
      // double muPosPx_Gen, muPosPy_Gen, muPosPz_Gen;
      TLorentzVector* muPosP_Gen;
      // double muNegPx_Gen, muNegPy_Gen, muNegPz_Gen;
      TLorentzVector* muNegP_Gen;

      //5.) Event related variables
      unsigned int eventNb, runNb, lumiBlock, nPriVtx;

      //6.) POL variables
      // std::vector<std::string> polVarNames_;
      // std::vector<std::string> polVarNamesGen_;
      // std::map<std::string, double> mapPolVarsToValue_;
      // std::map<std::string, double> mapPolVarsToValueGen_;

      //7.) TriggerNames Map
      std::map<std::string, int> mapTriggerNameToIntFired_;
      // USAGE:
      // ---> For all single or symmetric double triggers:
      // 0 : event not firing the corresponding trigger
      // 1 : event firing the corresponding trigger and all RECO muons matched to the required number of HLT objects (1 or 2)
      // 3 : event firing the corresponding trigger but at least one RECO muon is not matched to the HLT objects
      // 
      // ---> For the asymmetric triggers:
      // 0 : event not firing the corresponding trigger
      // 1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, POSITIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
      // -1 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, NEGATIVE-charge muon matched to the tighter HLT object (usually a L3 muon)
      // 2 : event firing the corresponding trigger, 2 RECO muons matched to 2 HLT objects, both matched to the tighter HLT object (usually a L3 muon)
      // 3 : event firing the corresponding trigger but at least one RECO muon is not matched to the HLT objects 
      std::map<std::string, int> mapTriggerNameToPrescaleFac_;

      Handle<pat::CompositeCandidateCollection > collAll;
      Handle<pat::CompositeCandidateCollection > collCalo;
      // Handle<TriggerResults> trigger;

      // data members
      InputTag       _patJpsi;
      InputTag       _patJpsiWithCalo;
      bool           _writeTree;
      string         _treefilename; 
      bool           _writeDataSet; 
      string         _datasetname;
      vector<string> _triggerForDataset;
      double         _massMin;
      double         _massMax;
      vector<double> _ptbinranges;
      vector<double> _etabinranges;
      bool           _onlythebest;
      bool           _applycuts;
      bool           _applyExpHitcuts;
      bool           _applyDiMuoncuts;
      bool           _useBS;
      bool           _useCalo;
      bool           _removeSignal;
      bool           _removeMuons;
      bool           _storeWs;
      bool           _writeOutCands;
      int            _MassCorr;
      // bool           _JSON;
      int            _oniaPDG;
      InputTag       _genParticles;
      bool           _isMC;
      bool           _storeAllMCEvents;
      bool           _isPromptMC;

      InputTag      _triggerresults;  

      vector<unsigned int>                     _thePassedCats[3];
      vector<const pat::CompositeCandidate*>   _thePassedCands[3];

      // number of events
      unsigned int nEvents;
      unsigned int passedTriggerResults_;
      unsigned int passedMuonSelectionCuts_;
      unsigned int passedTriggerMatch_;
      unsigned int passedTriggerResultsAnalyzer_;

      // limits 
      float JpsiMassMin;
      float JpsiMassMax;
      float JpsiCtMin;
      float JpsiCtMax;
      float JpsiPtMin;           // SET BY 
      float JpsiPtMax;           // DEFINITION
      float JpsiRapMin;          // OF BIN
      float JpsiRapMax;          // LIMITS

      math::XYZPoint RefVtx;
      ofstream* theTextFile;
      ofstream* JSON;
  
      int runtmp,lumitmp,count;
      int runmax,runmin;

      //muon distance
      PropagateToMuon prop1_, prop2_; 

      // Trigger Filter Studies
      edm::Handle< edm::TriggerResults> handleTriggerResults_;
      edm::InputTag tagTriggerResults_;
      HLTConfigProvider hltConfig_;
      bool hltConfigInit_;
      std::vector<std::string> HLTBitNames_;
      std::vector<std::string> HLTBitNames_SingleMu;
      std::vector<std::string> HLTLastFilterNames_SingleMu;
      std::vector<std::string> HLTBitNames_DoubleMu;
      std::vector<std::string> HLTLastFilterNames_DoubleMu;
      std::vector<std::string> HLTBitNames_MuL2Mu;
      std::vector<std::string> HLTLastFilterNames_MuL2Mu;
      std::vector<std::string> HLTBitNames_MuTrack;
      std::vector<std::string> HLTLastFilterNames_MuTrack;
      std::vector<std::string> HLTBitNames_MuTkMu;
      std::vector<std::string> HLTLastFilterNames_MuTkMu;
      std::map< std::string, std::string> mapTriggerToLastFilter_;
      
};    
      
// constants, enums and typedefs
enum {CS, HX, PHX, sGJ, GJ1, GJ2};
//

//
// constructors and destructor
//
JPsiAnalyzerPAT::JPsiAnalyzerPAT(const edm::ParameterSet& iConfig):
  _patJpsi(iConfig.getParameter<InputTag>("src")),
  _patJpsiWithCalo(iConfig.getParameter<InputTag>("srcWithCaloMuons")),
  _writeTree(iConfig.getParameter<bool>("writeTree")),
  _treefilename(iConfig.getParameter<string>("treeFileName")),	
  _writeDataSet(iConfig.getParameter<bool>("writeDataSet")),
  _datasetname(iConfig.getParameter<string>("dataSetName")),
  _triggerForDataset(iConfig.getParameter< vector<string> >("triggersForDataset")),
  _massMin(iConfig.getParameter<double>("massMin")),
  _massMax(iConfig.getParameter<double>("massMax")),
  _ptbinranges(iConfig.getParameter< vector<double> >("pTBinRanges")),	
  _etabinranges(iConfig.getParameter< vector<double> >("etaBinRanges")),
  _onlythebest(iConfig.getParameter<bool>("onlyTheBest")),		
  _applycuts(iConfig.getParameter<bool>("applyCuts")),
  _applyExpHitcuts(iConfig.getUntrackedParameter<bool>("applyExpHitCuts",false)),
  _applyDiMuoncuts(iConfig.getUntrackedParameter<bool>("applyDiMuonCuts",false)),
  _useBS(iConfig.getParameter<bool>("useBeamSpot")),
  _useCalo(iConfig.getUntrackedParameter<bool>("useCaloMuons",false)),
  _removeSignal(iConfig.getUntrackedParameter<bool>("removeSignalEvents",false)),
  _removeMuons(iConfig.getUntrackedParameter<bool>("removeTrueMuons",false)),
  _storeWs(iConfig.getUntrackedParameter<bool>("storeWrongSign",false)),
  _writeOutCands(iConfig.getUntrackedParameter<bool>("writeOutCandidates",false)),
  _MassCorr(iConfig.getParameter<int>("massCorrectionMode")),
  _oniaPDG(iConfig.getParameter<int>("oniaPDG")),
  _genParticles(iConfig.getParameter<InputTag>("genParticles")),
  _isMC(iConfig.getUntrackedParameter<bool>("isMC",false)),
  _storeAllMCEvents(iConfig.getUntrackedParameter<bool>("storeAllMCEvents",false)),
  _isPromptMC(iConfig.getUntrackedParameter<bool>("isPromptMC",false) ),
  prop1_(iConfig.getParameter<edm::ParameterSet>("propagatorStation1")),
  prop2_(iConfig.getParameter<edm::ParameterSet>("propagatorStation2")),
  tagTriggerResults_(iConfig.getParameter<InputTag>("triggerResultsLabel")),
  HLTBitNames_SingleMu(iConfig.getParameter< vector<string> >("HLTBitNames_SingleMu")),
  HLTLastFilterNames_SingleMu(iConfig.getParameter< vector<string> >("HLTLastFilterNames_SingleMu")),
  HLTBitNames_DoubleMu(iConfig.getParameter< vector<string> >("HLTBitNames_DoubleMu")),
  HLTLastFilterNames_DoubleMu(iConfig.getParameter< vector<string> >("HLTLastFilterNames_DoubleMu")),
  HLTBitNames_MuL2Mu(iConfig.getParameter< vector<string> >("HLTBitNames_MuL2Mu")),
  HLTLastFilterNames_MuL2Mu(iConfig.getParameter< vector<string> >("HLTLastFilterNames_MuL2Mu")),
  HLTBitNames_MuTrack(iConfig.getParameter< vector<string> >("HLTBitNames_MuTrack")),
  HLTLastFilterNames_MuTrack(iConfig.getParameter< vector<string> >("HLTLastFilterNames_MuTrack")),
  HLTBitNames_MuTkMu(iConfig.getParameter< vector<string> >("HLTBitNames_MuTkMu")),
  HLTLastFilterNames_MuTkMu(iConfig.getParameter< vector<string> >("HLTLastFilterNames_MuTkMu"))
{
   //now do what ever initialization is needed
  nEvents = 0; 
  // passedTriggerResults_=0;
  passedMuonSelectionCuts_=0;
  passedTriggerMatch_=0;
  // passedTriggerResultsAnalyzer_=0;

  JpsiMassMin = _massMin;
  JpsiMassMax = _massMax;
  JpsiCtMin = -2.0;
  JpsiCtMax = 3.5;

  if (_writeOutCands) theTextFile = new ofstream("passedCandidates.txt");
  
  if (HLTBitNames_SingleMu.size() != HLTLastFilterNames_SingleMu.size()) std::cout << "WARNING: Trigger names and last filters do not match in size!" << std::endl;
  std::vector<std::string>::iterator it2 = HLTLastFilterNames_SingleMu.begin();
  for(std::vector<std::string>::iterator it = HLTBitNames_SingleMu.begin(); it != HLTBitNames_SingleMu.end(); ++it){
    mapTriggerToLastFilter_[*it] = *it2;
    ++it2;   HLTBitNames_.push_back(*it);
  }
  if (HLTBitNames_DoubleMu.size() != HLTLastFilterNames_DoubleMu.size()) std::cout << "WARNING: Trigger names and last filters do not match in size!" << std::endl;
  it2 = HLTLastFilterNames_DoubleMu.begin();
  for(std::vector<std::string>::iterator it = HLTBitNames_DoubleMu.begin(); it != HLTBitNames_DoubleMu.end(); ++it){
    mapTriggerToLastFilter_[*it] = *it2;
    ++it2;   HLTBitNames_.push_back(*it);
  }
  if (2*HLTBitNames_MuL2Mu.size() != HLTLastFilterNames_MuL2Mu.size()) std::cout << "WARNING: Trigger names and last filters do not match in size!" << std::endl;
  it2 = HLTLastFilterNames_MuL2Mu.begin();
  for(std::vector<std::string>::iterator it = HLTBitNames_MuL2Mu.begin(); it != HLTBitNames_MuL2Mu.end(); ++it){
    std::string theName = *it;
    mapTriggerToLastFilter_[theName] = *it2;
    ++it2;   HLTBitNames_.push_back(theName);  
    theName += "_special";
    mapTriggerToLastFilter_[theName] = *it2;
    ++it2;    
  }
  if (HLTBitNames_MuTrack.size() != HLTLastFilterNames_MuTrack.size()) std::cout << "WARNING: Trigger names and last filters do not match in size!" << std::endl;
  it2 = HLTLastFilterNames_MuTrack.begin();
  for(std::vector<std::string>::iterator it = HLTBitNames_MuTrack.begin(); it != HLTBitNames_MuTrack.end(); ++it){
    mapTriggerToLastFilter_[*it] = *it2;
    ++it2;   HLTBitNames_.push_back(*it);
  }
  if (HLTBitNames_MuTkMu.size() != HLTLastFilterNames_MuTkMu.size()) std::cout << "WARNING: Trigger names and last filters do not match in size!" << std::endl;
  it2 = HLTLastFilterNames_MuTkMu.begin();
  for(std::vector<std::string>::iterator it = HLTBitNames_MuTkMu.begin(); it != HLTBitNames_MuTkMu.end(); ++it){
    mapTriggerToLastFilter_[*it] = *it2;
    ++it2;  HLTBitNames_.push_back(*it);
  }

  for(std::vector<std::string>::iterator it = HLTBitNames_.begin(); it != HLTBitNames_.end(); ++it){
      mapTriggerNameToIntFired_[*it] = -9999;
  }
}


JPsiAnalyzerPAT::~JPsiAnalyzerPAT()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if (_writeOutCands) theTextFile->close();
}

// ------------ method called once each job just before starting event loop  ------------
void
JPsiAnalyzerPAT::beginJob()
{
    //std::cout << "[JPsiAnalyzerPAT] --- beginJob " << std::endl;
  if (_writeTree) {
    fOut_ = new TFile(_treefilename.c_str(), "RECREATE");
    fOut_->cd();

    // TTree
    //load Branches
    tree_ = new TTree ("data", "CMSSW Quarkonia J/psi Polarization+Trigger Tree");

    JpsiP = new TLorentzVector();
    muPosP = new TLorentzVector();
    muNegP = new TLorentzVector();
    JpsiP_Gen = new TLorentzVector();
    muPosP_Gen = new TLorentzVector();
    muNegP_Gen = new TLorentzVector();

    // Event variables
    tree_->Branch("eventNb",             &eventNb,             "eventNb/I");
    tree_->Branch("runNb",               &runNb,               "runNb/I");
    tree_->Branch("lumiBlock",           &lumiBlock,           "lumiBlock/I");
    tree_->Branch("nPriVtx",             &nPriVtx,             "nPriVtx/I");

    // Jpsi Variables
    tree_->Branch("JpsiType",   &JpsiType,  "JpsiType/I");
    tree_->Branch("JpsiP",  "TLorentzVector", &JpsiP);
    // tree_->Branch("JpsiMass",   &JpsiMass,  "JpsiMass/D");
    // tree_->Branch("JpsiPt",     &JpsiPt,    "JpsiPt/D");
    // tree_->Branch("JpsiRap",    &JpsiRap,   "JpsiRap/D");
    tree_->Branch("JpsiCharge", &JpsiCharge,"JpsiCharge/I");
    // tree_->Branch("JpsiPx",     &JpsiPx,    "JpsiPx/D");
    // tree_->Branch("JpsiPy",     &JpsiPy,    "JpsiPy/D");
    // tree_->Branch("JpsiPz",     &JpsiPz,    "JpsiPz/D");
    tree_->Branch("Jpsict",     &Jpsict,    "Jpsict/D");
    tree_->Branch("JpsictErr",  &JpsictErr, "JpsictErr/D");
    tree_->Branch("JpsiVprob",  &JpsiVprob, "JpsiVprob/D");
    //muon distance 
    tree_->Branch("JpsiDistM1",   &JpsiDistM1,    "JpsiDistM1/D");
    tree_->Branch("JpsiDphiM1",   &JpsiDphiM1,    "JpsiDphiM1/D");
    tree_->Branch("JpsiDrM1",     &JpsiDrM1,      "JpsiDrM1/D");
    tree_->Branch("JpsiDistM2",   &JpsiDistM2,    "JpsiDistM2/D");
    tree_->Branch("JpsiDphiM2",   &JpsiDphiM2,    "JpsiDphiM2/D");
    tree_->Branch("JpsiDrM2",     &JpsiDrM2,      "JpsiDrM2/D");

    tree_->Branch("muPosP", "TLorentzVector", &muPosP);
    tree_->Branch("muNegP", "TLorentzVector", &muNegP);
    // tree_->Branch("muPosPx",    &muPosPx,   "muPosPx/D");
    // tree_->Branch("muPosPy",    &muPosPy,   "muPosPy/D");
    // tree_->Branch("muPosPz",    &muPosPz,   "muPosPz/D");
    // tree_->Branch("muNegPx",    &muNegPx,   "muNegPx/D");
    // tree_->Branch("muNegPy",    &muNegPy,   "muNegPy/D");
    // tree_->Branch("muNegPz",    &muNegPz,   "muNegPz/D");

    //add HLT Variables to TTree
    for(std::vector< std::string >:: iterator it = HLTBitNames_.begin(); it != HLTBitNames_.end(); ++it){
        std::string hlt_name= *it;
        tree_->Branch(hlt_name.c_str(), &(mapTriggerNameToIntFired_[*it]), (hlt_name + "/I").c_str());
	tree_->Branch((hlt_name+"_PreScale").c_str(), &(mapTriggerNameToPrescaleFac_[*it]), (hlt_name+"_PreScale" + "/I").c_str());
    }

    //add Generator Information
    if(_isMC){
        tree_->Branch("MCType",         &MCType,        "MCType/I");
	tree_->Branch("JpsiP_Gen",  "TLorentzVector", &JpsiP_Gen);
        // tree_->Branch("JpsiMass_Gen",   &JpsiMass_Gen,  "JpsiMass_Gen/D");
        // tree_->Branch("JpsiPt_Gen",     &JpsiPt_Gen,    "JpsiPt_Gen/D");
        // tree_->Branch("JpsiRap_Gen",    &JpsiRap_Gen,   "JpsiRap_Gen/D");
        // tree_->Branch("JpsiPx_Gen",     &JpsiPx_Gen,    "JpsiPx_Gen/D");
        // tree_->Branch("JpsiPy_Gen",     &JpsiPy_Gen,    "JpsiPy_Gen/D");
        // tree_->Branch("JpsiPz_Gen",     &JpsiPz_Gen,    "JpsiPz_Gen/D");
        tree_->Branch("Jpsict_Gen",     &Jpsict_Gen,    "Jpsict_Gen/D");
        tree_->Branch("muPosP_Gen",  "TLorentzVector", &muPosP_Gen);
        tree_->Branch("muNegP_Gen",  "TLorentzVector", &muNegP_Gen);
        // tree_->Branch("muPosPx_Gen",    &muPosPx_Gen,   "muPosPx_Gen/D");
        // tree_->Branch("muPosPy_Gen",    &muPosPy_Gen,   "muPosPy_Gen/D");
        // tree_->Branch("muPosPz_Gen",    &muPosPz_Gen,   "muPosPz_Gen/D");
        // tree_->Branch("muNegPx_Gen",    &muNegPx_Gen,   "muNegPx_Gen/D");
        // tree_->Branch("muNegPy_Gen",    &muNegPy_Gen,   "muNegPy_Gen/D");
        // tree_->Branch("muNegPz_Gen",    &muNegPz_Gen,   "muNegPz_Gen/D");
    }
  }
  if (_writeDataSet) {
    
     fOut2_ = new TFile(_datasetname.c_str(), "RECREATE");
     fOut2_->cd();

     Jpsi_PtType = new RooCategory("Jpsi_PtType","Category of Pt");
     Jpsi_RapType = new RooCategory("Jpsi_RapType","Category of Rap");

     JpsiPtMin = _ptbinranges[0];  cout << "Pt min = " << JpsiPtMin << endl;
     JpsiPtMax = _ptbinranges[_ptbinranges.size()-1];  cout << "Pt max = " << JpsiPtMax << endl;
     
     for(unsigned int i=0;i<_ptbinranges.size()-1;i++){
       char catname[100];
       sprintf(catname,"P%d",i+1);
       Jpsi_PtType->defineType(catname,i+1); 
       cout << "Pt bin " << i+1 << ": Min = " << _ptbinranges[i] << " Max = " << _ptbinranges[i+1] << endl;   
     }
     
     JpsiRapMin = _etabinranges[0];  cout << "Rap min = " << JpsiRapMin << endl;
     JpsiRapMax = _etabinranges[_etabinranges.size()-1];  cout << "Rap max = " << JpsiRapMax << endl;
     
     for(unsigned int i=0;i<_etabinranges.size()-1;i++){
       char catname[100];
       sprintf(catname,"E%d",i+1);
       Jpsi_RapType->defineType(catname,i+1); 
       cout << "Rap bin " << i+1 << ": Min = " << _etabinranges[i] << " Max = " << _etabinranges[i+1] << endl;   
     }
     
     Jpsi_Type = new RooCategory("Jpsi_Type","Category of Jpsi");
     Jpsi_MatchType = new RooCategory("Jpsi_MatchType","Category of matching");
     
     Jpsi_Type->defineType("GG",0);
     Jpsi_Type->defineType("GT",1);
     Jpsi_Type->defineType("TT",2);
     if (_useCalo) {
       Jpsi_Type->defineType("GC",3);  
       Jpsi_Type->defineType("TC",4);  
       Jpsi_Type->defineType("CC",5);  
     }
     
     Jpsi_MatchType->defineType("unmatched",0);
     Jpsi_MatchType->defineType("matched",1);
     
     Jpsi_Sign = new RooCategory("Jpsi_Sign","Sign of Jpsi dimuons");
     
     Jpsi_Sign->defineType("OS",0);
     Jpsi_Sign->defineType("SSP",1);
     Jpsi_Sign->defineType("SSM",2);

     Jpsi_Mass = new RooRealVar("Jpsi_Mass","J/psi mass",JpsiMassMin,JpsiMassMax,"GeV/c^{2}");
     Jpsi_Pt = new RooRealVar("Jpsi_Pt","J/psi pt",JpsiPtMin,JpsiPtMax,"GeV/c");
     Jpsi_Rap = new RooRealVar("Jpsi_Rap","J/psi eta",-JpsiRapMax,JpsiRapMax);
     Jpsi_ct = new RooRealVar("Jpsi_ct","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");
     Jpsi_ctErr = new RooRealVar("Jpsi_ctErr","J/psi ctau error",-1.,1.,"mm");
     Jpsi_ctTrue = new RooRealVar("Jpsi_ctTrue","J/psi ctau true",-100.,JpsiCtMax,"mm"); 		

     RooArgList varlist(*Jpsi_Mass,*Jpsi_ct,*Jpsi_Pt,*Jpsi_Rap,*Jpsi_Type,*Jpsi_MatchType);
     varlist.add(*Jpsi_ctTrue);   varlist.add(*Jpsi_PtType);
     varlist.add(*Jpsi_RapType);  varlist.add(*Jpsi_ctErr);
     varlist.add(*Jpsi_Sign);

     data = new RooDataSet("data","A sample",varlist);
  }
}


double JPsiAnalyzerPAT::CorrectMass(const reco::Muon& mu1,const reco::Muon& mu2, int mode){  
  double CMass=0;
  const double mumass=0.105658;
  double k1,k2;
  double pt1=mu1.innerTrack()->pt();
  double pt2=mu2.innerTrack()->pt();
  double eta1=mu1.innerTrack()->eta();
  double eta2=mu2.innerTrack()->eta();
  if (mode==1){
    k1=1.0009;//constant scale correction
    k2=1.0009;
  }
  if (mode==2){
    k1=1.0019-0.0004*pt1;
    k2=1.0019-0.0004*pt2; // pt dependent correction
  }
  if (mode==3){
      double a0=0.00038; //3.8 * pow(10,-4);
      double a1=0.0;
      double a2=0.0003; //3.0 * pow(10,-4);
      double a3=0.0;

      k1=1+a0+a1*fabs(eta1)+a2*eta1*eta1+a3*pt1;
      k2=1+a0+a1*fabs(eta2)+a2*eta2*eta2+a3*pt2;// pt and eta dependent
  }

  if (mode == 4){
      double a0=1.002;
      double a1=-0.002;
      double a2=0.001;
      double a3=-0.0001;

      k1=a0+a1*fabs(eta1)+a2*eta1*eta1+a3*pt1;
      k2=a0+a1*fabs(eta2)+a2*eta2*eta2+a3*pt2;// pt and eta dependent
  }

  math::XYZVector mom1=mu1.innerTrack()->momentum();
  math::XYZVector mom2=mu2.innerTrack()->momentum();
  mom1=k1*mom1; 
  mom2=k2*mom2;
  double E1=sqrt(mom1.mag2()+(mumass*mumass));
  double E2=sqrt(mom2.mag2()+(mumass*mumass));
  math::XYZVector momtot=mom1+mom2;
  CMass=sqrt((E1+E2)*(E1+E2)-momtot.mag2());
  return CMass;
}

// ------------ method called to for each event  ------------
void
JPsiAnalyzerPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
   nEvents++;

   // reset TTree Variables
   if (_writeTree) this->resetDSVariables();

   // check HLT TriggerReuslts
   this->hltReport(iEvent, iSetup);

   // Event related infos
   eventNb= iEvent.id().event() ;
   runNb=iEvent.id().run() ;
   lumiBlock= iEvent.luminosityBlock() ;

   Handle<reco::VertexCollection> privtxs;
   iEvent.getByLabel("offlinePrimaryVertices", privtxs);
   nPriVtx = privtxs->size();
   VertexCollection::const_iterator privtx;

   if ( privtxs->begin() != privtxs->end() ) {
     privtx=privtxs->begin();
     RefVtx = privtx->position();
   } else {
     RefVtx.SetXYZ(0.,0.,0.);
   }
   // }

   try {iEvent.getByLabel(_patJpsi,collAll);} 
   catch (...) {cout << "J/psi not present in event!" << endl;}

   if (_useCalo) {
     try {iEvent.getByLabel(_patJpsiWithCalo,collCalo);} 
     catch (...) {cout << "J/psi to calomuons not present in event!" << endl;}
   }

   _thePassedCats[0].clear();      _thePassedCands[0].clear();
   _thePassedCats[1].clear();      _thePassedCands[1].clear();
   _thePassedCats[2].clear();      _thePassedCands[2].clear();

   // APPLY CUTS
   int lastSign = 0;
   this->makeCuts(0);
   if (_storeWs) {
     this->makeCuts(1);
     this->makeCuts(2);
     lastSign = 2;
   }

   bool storeEvent = false;

   // BEST J/PSI? 
   if (_onlythebest) {  // yes, fill simply the best (possibly wrong-sign)

//      for (int iSign = 0; iSign <= lastSign; iSign++) {
     for (int iSign = lastSign; iSign >= 0; iSign--) {
       pair< unsigned int, const pat::CompositeCandidate* > theBest = theBestQQ(iSign);
       if (theBest.first < 10) {
           fillTreeAndDS(theBest.first, theBest.second, iEvent);
           passedMuonSelectionCuts_++;
	   if(iSign == 0) storeEvent = true;
       }
     }
   } else {   // no, fill all candidates passing cuts (possibly wrong-sign)

     for (int iSign = 0; iSign <= lastSign; iSign++) {
       for( unsigned int count = 0; count < _thePassedCands[iSign].size(); count++) { 
	 fillTreeAndDS(_thePassedCats[iSign].at(count), _thePassedCands[iSign].at(count),iEvent); 
       }
     }
   }

   //! FILL GENERATOR COLLECTION and store the event
  if ( _storeAllMCEvents || storeEvent ) {

     Handle<reco::GenParticleCollection> genParticles;
     iEvent.getByLabel( _genParticles, genParticles );
     if ( genParticles.isValid() )
       {
	 //std::cout << "------ analyze GENERATED JPsis:" << std::endl;
	 this->analyzeGenerator( genParticles );
       }
   
     // Write all Branches to the Tree ONLY 
     // - for the best candidate
     // - for the opposite sign
     if (_writeTree) tree_->Fill();
   }
}

void
JPsiAnalyzerPAT::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){

    //init HLTConfigProvider
    const std::string pro = tagTriggerResults_.process();
    bool changed = true;

    //bool init(const edm::Run& iRun, const edm::EventSetup& iSetup, const std::string& processName, bool& changed);
    hltConfigInit_ = false;
    if( hltConfig_.init(iRun, iSetup, pro, changed) ) hltConfigInit_ = true;

    prop1_.init(iSetup);
    prop2_.init(iSetup);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JPsiAnalyzerPAT::endJob() {
  
  cout << "Total number of events = " << nEvents << endl;
  // cout << "Analyzed runs from  " << runmin << "  to  " << runmax << endl; 
  cout << "============================================================" << endl;
  // cout << "Total number of passed candidates TRIGGER RESULTS ANALYZER   = " << passedTriggerResultsAnalyzer_ << endl;
  cout << "Total number of passed candidates MUON SELECTION CUTS        = " << passedMuonSelectionCuts_ << endl;
  cout << "Total number of passed candidates TRIGGER MATCH              = " << passedTriggerMatch_ << endl;

  /* if (_JSON){
    cout << "JSON file produced" << endl;
    *JSON << lumitmp <<"]]}";
    JSON->close();
    }*/

  // Write TTree to File
  if (_writeTree) {
    fOut_->cd();
    tree_->Write();
    fOut_->Close();
  }
  if (_writeDataSet) {
    fOut2_->cd();
    data->Write();
    fOut2_->Close();
  }
}

//! Fill the TTree with all RECO variables
void 
JPsiAnalyzerPAT::fillTreeAndDS(unsigned int theCat, const pat::CompositeCandidate* aCand, const edm::Event& iEvent){
  
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon2"));
  
  //1.) continue only if we have opposite sign muons <- REMOVED
  //
  if (muon1->charge()*muon2->charge() >= 0) {
    if(muon1->charge() == 0) {
      printf("pat::Muon with zero charge!\n");   return;
    }
    if(muon2->charge() == 0) {
      printf("pat::Muon with zero charge!\n");   return;
    }
    // return;
  }

  const pat::Muon *muonPos = 0, *muonNeg = 0;
  if(muon1->charge() > 0){ muonPos = muon1; muonNeg = muon2;}
  else if(muon1->charge() < 0){ muonPos = muon2; muonNeg = muon1;}

  //
  float theMass = aCand->mass();
  if (_MassCorr!=0){
    double CMass = CorrectMass(*muon1,*muon2,_MassCorr);
    if (CMass!=0.0) theMass = CMass;
  }

  float theRapidity = aCand->rapidity();
  // if (!_useRapidity) theRapidity = theRapidity;

  float theCtau; 
  if (_useBS) {theCtau = 10.*aCand->userFloat("ppdlBS");}
  else {theCtau = 10.*aCand->userFloat("ppdlPV");}

  float theCtauErr; 
  if (_useBS) {theCtauErr = 10.*aCand->userFloat("ppdlErrBS");}
  else {theCtauErr = 10.*aCand->userFloat("ppdlErrPV");}

  float theCharge = aCand->charge();

  // MC matching
  reco::GenParticleRef genJpsi = aCand->genParticleRef();
  bool isMatched = (genJpsi.isAvailable() && genJpsi->pdgId() == _oniaPDG);

  // Input DataSet Type: P, NP, BG J/psi
  if (isMatched && _isPromptMC) MCType= 0;
  if (isMatched && _isPromptMC == false) MCType=1;
  if (!isMatched && _isMC) MCType=2;

  if (isMatched && _removeSignal) return;

  reco::GenParticleRef genMu1 = muon1->genParticleRef();
  reco::GenParticleRef genMu2 = muon2->genParticleRef();
  bool isMuMatched = (genMu1.isAvailable() && genMu2.isAvailable() && 
		      genMu1->pdgId()*genMu2->pdgId() == -169 && 
		      genMu1->momentum().rho() > 2.5 && genMu2->momentum().rho() > 2.5);
  if (isMuMatched && _removeMuons) return;

  // PAT trigger match, 2 muons to the last filter used in the HLT path (new way)
  this->matchMuonToHlt(muonPos, muonNeg);

  // some counter
  passedTriggerMatch_++;
    
  // JpsiMass=theMass;

  if (_writeOutCands) *theTextFile << iEvent.id().run() << "\t" << iEvent.luminosityBlock() << "\t" << iEvent.id().event() << "\t" << theMass << "\n";

  // write out JPsi RECO information
  // JpsiPt=aCand->pt();
  // JpsiRap=theRapidity;
  JpsiCharge=theCharge;
  // std::cout << "[JPsiAnalyzerPAT::fillTreeAndDS] ----- JpsiCharge: " << theCharge << std::endl;
  // JpsiPx=aCand->px();
  // JpsiPy=aCand->py();
  // JpsiPz=aCand->pz();
  JpsiP->SetPxPyPzE(aCand->px(),aCand->py(),aCand->pz(),aCand->energy());
  Jpsict=theCtau;
  JpsictErr=theCtauErr;
  Jpsict_Gen=10.*aCand->userFloat("ppdlTrue");
  JpsiType=theCat;
  JpsiVprob=aCand->userFloat("vProb");
  this->muonStationDistance(aCand);
  
  // write out Muon RECO information
  float f_muPosPx, f_muPosPy, f_muPosPz;
  float f_muNegPx, f_muNegPy, f_muNegPz;
  f_muPosPx = muonPos->px();
  f_muPosPy = muonPos->py();
  f_muPosPz = muonPos->pz();
  f_muNegPx = muonNeg->px();
  f_muNegPy = muonNeg->py();
  f_muNegPz = muonNeg->pz();
  // muPosPx= f_muPosPx ;
  // muPosPy= f_muPosPy ;
  // muPosPz= f_muPosPz ;
  // muNegPx= f_muNegPx ;
  // muNegPy= f_muNegPy ;
  // muNegPz= f_muNegPz ;
  
  //write out Calculated Polarization variables
  Double_t muMass = 0.105658;
  
  Double_t enMuPos = sqrt(f_muPosPx*f_muPosPx + f_muPosPy*f_muPosPy + f_muPosPz*f_muPosPz + muMass*muMass);
  // TLorentzVector *muPosP = new TLorentzVector();
  muPosP->SetPxPyPzE(f_muPosPx, f_muPosPy, f_muPosPz, enMuPos);
  
  Double_t enMuNeg = sqrt(f_muNegPx*f_muNegPx + f_muNegPy*f_muNegPy + f_muNegPz*f_muNegPz + muMass*muMass);
  // TLorentzVector *muNegP = new TLorentzVector();
  muNegP->SetPxPyPzE(f_muNegPx, f_muNegPy, f_muNegPz, enMuNeg);
  
  //! Fill Polarization Variables;
  // std::vector< float > thisCosTh, thisPhi;
  // thisCosTh.resize(6); thisPhi.resize(6);
  // this->calcPol(*muPosP, *muNegP, thisCosTh, thisPhi);

  if (_writeDataSet) {
    
    bool trigOK = false;
    for (unsigned int iTrig = 0 ; iTrig < _triggerForDataset.size() ; iTrig++) {
      if (mapTriggerNameToIntFired_[_triggerForDataset.at(iTrig)] == 1 ||
	  mapTriggerNameToIntFired_[_triggerForDataset.at(iTrig)] == -1 ||
	  mapTriggerNameToIntFired_[_triggerForDataset.at(iTrig)] == 2 ) trigOK = true;
    }

    if (theMass > JpsiMassMin && theMass < JpsiMassMax && 
	theCtau > JpsiCtMin && theCtau < JpsiCtMax && 
	aCand->pt() > JpsiPtMin && aCand->pt() < JpsiPtMax && 
	fabs(theRapidity) > JpsiRapMin && fabs(theRapidity) < JpsiRapMax &&
	isMuonInAccept(muon1) && isMuonInAccept(muon2) &&
	trigOK) {

      int ss=999;
      if (muon1->charge() + muon2->charge() == 0) ss=0;
      if (muon1->charge() + muon2->charge() == 2) ss=1;
      if (muon1->charge() + muon2->charge() == -2) ss=2;

      Jpsi_Sign->setIndex(ss,kTRUE);
      
      Jpsi_Pt->setVal(aCand->pt()); 
      Jpsi_Rap->setVal(theRapidity); 
      Jpsi_Mass->setVal(theMass);
      Jpsi_ct->setVal(theCtau);
      Jpsi_ctErr->setVal(theCtauErr);
      // cout << "Type = " << theCat << " pt = " << aCand->pt() << " eta = " << theRapidity << endl;
      // cout << " PPDL = " << theCtau << " Mother = " << aCand->userInt("momPDGId") << " PPDL true = " << 10.*aCand->userFloat("ppdlTrue") << endl;
      Jpsi_Type->setIndex(theCat,kTRUE);
      Jpsi_MatchType->setIndex((int)isMatched,kTRUE);
      Jpsi_ctTrue->setVal(10.*aCand->userFloat("ppdlTrue"));
    
      Jpsi_PtType->setIndex(getJpsiVarType(aCand->pt(),_ptbinranges),kTRUE);
      Jpsi_RapType->setIndex(getJpsiVarType(fabs(theRapidity),_etabinranges),kTRUE);
      // Fill RooDataSet
      RooArgSet varlist_tmp(*Jpsi_Mass,*Jpsi_ct,*Jpsi_Pt,*Jpsi_Rap,*Jpsi_Type,*Jpsi_MatchType);   // temporarily remove tag-and-probe weights
      varlist_tmp.add(*Jpsi_ctTrue);   varlist_tmp.add(*Jpsi_PtType);
      varlist_tmp.add(*Jpsi_RapType);  varlist_tmp.add(*Jpsi_ctErr);
      varlist_tmp.add(*Jpsi_Sign);
      data->add(varlist_tmp);
    }
  }
}
        
void JPsiAnalyzerPAT::makeCuts(int sign) {

  if (collAll.isValid()) {

    for(vector<pat::CompositeCandidate>::const_iterator it=collAll->begin();
	it!=collAll->end();++it) {
      
      const pat::CompositeCandidate* cand = &(*it);	
      // cout << "Now checking candidate of type " << theJpsiCat << " with pt = " << cand->pt() << endl;
      const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(cand->daughter("muon1"));
      const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(cand->daughter("muon2"));
 
      bool thisSign = (sign == 0 && muon1->charge() + muon2->charge() == 0) || 
	(sign == 1 && muon1->charge() + muon2->charge() == 2) || 
	(sign == 2 && muon1->charge() + muon2->charge() == -2);

      if (thisSign) {	  
	  
        // global + global?
	if (muon1->isGlobalMuon() && muon2->isGlobalMuon() &&
	    muon1->isTrackerMuon() && muon2->isTrackerMuon()   ) {
	  if (!_applycuts || (selGlobalMuon(muon1) &&
			      selGlobalMuon(muon2) &&
			      selDimuon(cand) )) {
	    _thePassedCats[sign].push_back(0);  _thePassedCands[sign].push_back(cand);
            continue;
	  }
	}
	
        // global + tracker? (x2)    
	if (muon1->isGlobalMuon() && muon2->isTrackerMuon() &&
	    muon1->isTrackerMuon()    ) {
	  if (!_applycuts || (selGlobalMuon(muon1) &&
			      selTrackerMuon(muon2) &&
			      selDimuon(cand) )) {
	    _thePassedCats[sign].push_back(1);  _thePassedCands[sign].push_back(cand);
	    continue;
	  }
	}

        if (muon2->isGlobalMuon() && muon1->isTrackerMuon() &&
            muon2->isTrackerMuon() ) {
	  if (!_applycuts || (selGlobalMuon(muon2) &&
			      selTrackerMuon(muon1) &&
			      selDimuon(cand) )) {
	    _thePassedCats[sign].push_back(1);  _thePassedCands[sign].push_back(cand);
	    continue;
	  }
	}

        // tracker + tracker?  
        if (muon1->isTrackerMuon() && muon2->isTrackerMuon() ) {
	  if (!_applycuts || (selTrackerMuon(muon1) &&
			      selTrackerMuon(muon2) &&
			      selDimuon(cand) )) {
	    _thePassedCats[sign].push_back(2);  _thePassedCands[sign].push_back(cand);
	    continue;
	  }
	}
      }
    }
  }
  
  if (_useCalo && collCalo.isValid()) {

    for(vector<pat::CompositeCandidate>::const_iterator it=collCalo->begin();
	it!=collCalo->end();++it) {
      
      const pat::CompositeCandidate* cand = &(*it);
      
      const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(cand->daughter("muon1"));
      const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(cand->daughter("muon2"));

      bool thisSign = (sign == 0 && muon1->charge() + muon2->charge() == 0) || 
	(sign == 1 && muon1->charge() + muon2->charge() == 2) || 
	(sign == 2 && muon1->charge() + muon2->charge() == -2);

      if (thisSign && !(muon1->isTrackerMuon() && muon2->isTrackerMuon()) ) {

	// global + calo? (x2)
	if (muon1->isGlobalMuon() && muon2->isCaloMuon() ) {
	  if (!_applycuts || (selGlobalMuon(muon1) &&
			      selCaloMuon(muon2) &&
			      selDimuon(cand) )) {
	    _thePassedCats[sign].push_back(3);  _thePassedCands[sign].push_back(cand);
            continue;
	  }
	}

	if (muon2->isGlobalMuon() && muon1->isCaloMuon() ) {
	  if (!_applycuts || (selGlobalMuon(muon2) &&
			      selCaloMuon(muon1) &&
			      selDimuon(cand) )) {
	    _thePassedCats[sign].push_back(3);  _thePassedCands[sign].push_back(cand);
            continue;
	  }
	}
	
        // tracker + calo? (x2)    
	if (muon1->isTrackerMuon() && muon2->isCaloMuon() ) {
	  if (!_applycuts || (selTrackerMuon(muon1) &&
			      selCaloMuon(muon2) &&
			      selDimuon(cand) )) {
	    _thePassedCats[sign].push_back(4);  _thePassedCands[sign].push_back(cand);
	    continue;
	  }
	}

        if (muon2->isTrackerMuon() && muon1->isCaloMuon() ) {
	  if (!_applycuts || (selTrackerMuon(muon2) &&
			      selCaloMuon(muon1) &&
			      selDimuon(cand) )) {
	    _thePassedCats[sign].push_back(4);  _thePassedCands[sign].push_back(cand);
	    continue;
	  }
	}

        // calo + calo? 
        if (muon1->isCaloMuon() && muon2->isCaloMuon() ) {
	  if (!_applycuts || (selCaloMuon(muon1) &&
			      selCaloMuon(muon2) &&
			      selDimuon(cand) )) {
	    _thePassedCats[sign].push_back(5);  _thePassedCands[sign].push_back(cand);
	    continue;
	  }
	}
      }
    }
  }

  return;
}

pair< unsigned int, const pat::CompositeCandidate* > 
JPsiAnalyzerPAT::theBestQQ(int sign) {

  unsigned int theBestCat = 99;
  const pat::CompositeCandidate* theBestCand = new pat::CompositeCandidate();

  for( unsigned int i = 0; i < _thePassedCands[sign].size(); i++) { 
    if (_thePassedCats[sign].at(i) < theBestCat) {
      theBestCat = _thePassedCats[sign].at(i);
      theBestCand = _thePassedCands[sign].at(i);
    }
  }

  pair< unsigned int, const pat::CompositeCandidate* > result = make_pair(theBestCat, theBestCand );
  return result;

}

bool
JPsiAnalyzerPAT::isMuonInAccept(const pat::Muon* aMuon) {
   // *USE* muon kinematical cuts (eta dependent momentum / pT cuts )
   return (fabs(aMuon->eta()) < 2.4 &&
           ((fabs(aMuon->eta()) < 1.3 && aMuon->pt() > 3.3) ||
           (fabs(aMuon->eta()) > 1.3 && fabs(aMuon->eta()) < 2.2 && aMuon->p() > 2.9) ||
           (fabs(aMuon->eta()) > 2.2 && aMuon->pt() > 0.8)));

   // *REMOVE* muon kinematical cuts (eta dependent momentum / pT cuts )
   // by just returning TRUE
   //  return true;
}

bool
JPsiAnalyzerPAT::selGlobalMuon(const pat::Muon* aMuon) {

  TrackRef iTrack = aMuon->innerTrack();
  const reco::HitPattern& p = iTrack->hitPattern();
  const reco::HitPattern& ei = iTrack->trackerExpectedHitsInner();
  const reco::HitPattern& eo = iTrack->trackerExpectedHitsOuter();

  TrackRef gTrack = aMuon->globalTrack();
  const reco::HitPattern& q = gTrack->hitPattern();

  bool trackOK = false;
  // cooler way of cutting on tracks
  if (_applyExpHitcuts) {
    float fHits = iTrack->found() / (iTrack->found() + iTrack->lost() + ei.numberOfHits() + eo.numberOfHits());
    trackOK = (fHits >= 0.8 && (p.hasValidHitInFirstPixelBarrel() || p.hasValidHitInFirstPixelEndcap() ));
  // old way of cutting on tracks  
  } else trackOK = (iTrack->found() > 10);

  return (// isMuonInAccept(aMuon) &&
	  trackOK &&
	  gTrack->chi2()/gTrack->ndof() < 20.0 &&
          q.numberOfValidMuonHits() > 0 &&
	  iTrack->chi2()/iTrack->ndof() < 1.8 &&
	  aMuon->muonID("TrackerMuonArbitrated") &&
	  aMuon->muonID("TMOneStationTight") &&
          p.pixelLayersWithMeasurement() > 1 &&
	  fabs(iTrack->dxy(RefVtx)) < 3.0 &&
          fabs(iTrack->dz(RefVtx)) < 15.0 );
}

bool 
JPsiAnalyzerPAT::selTrackerMuon(const pat::Muon* aMuon) {
  
  TrackRef iTrack = aMuon->innerTrack();
  const reco::HitPattern& p = iTrack->hitPattern();
  const reco::HitPattern& ei = iTrack->trackerExpectedHitsInner();
  const reco::HitPattern& eo = iTrack->trackerExpectedHitsOuter();

  bool trackOK = false;
  // cooler way of cutting on tracks
  if (_applyExpHitcuts) {
    float fHits = iTrack->found() / (iTrack->found() + iTrack->lost() + ei.numberOfHits() + eo.numberOfHits());
    trackOK = (fHits >= 0.8 && (p.hasValidHitInFirstPixelBarrel() || p.hasValidHitInFirstPixelEndcap() ));
  // old way of cutting on tracks  
  } else trackOK = (iTrack->found() > 10);

  return (// isMuonInAccept(aMuon) &&
	  trackOK &&
 	  iTrack->chi2()/iTrack->ndof() < 1.8 &&
	  aMuon->muonID("TrackerMuonArbitrated") &&
	  aMuon->muonID("TMOneStationTight") &&
          p.pixelLayersWithMeasurement() > 1 &&
	  fabs(iTrack->dxy(RefVtx)) < 3.0 &&
          fabs(iTrack->dz(RefVtx)) < 15.0 );
}

bool 
JPsiAnalyzerPAT::selCaloMuon(const pat::Muon* aMuon) {
  
  TrackRef iTrack = aMuon->innerTrack();
  const reco::HitPattern& p = iTrack->hitPattern();
  const reco::HitPattern& ei = iTrack->trackerExpectedHitsInner();
  const reco::HitPattern& eo = iTrack->trackerExpectedHitsOuter();
 
  bool trackOK = false;
  // cooler way of cutting on tracks
  if (_applyExpHitcuts) {
    float fHits = iTrack->found() / (iTrack->found() + iTrack->lost() + ei.numberOfHits() + eo.numberOfHits());
    trackOK = (fHits >= 0.8 && (p.hasValidHitInFirstPixelBarrel() || p.hasValidHitInFirstPixelEndcap() ));
  // old way of cutting on tracks  
  } else trackOK = (iTrack->found() > 10);

  return (// isMuonInAccept(aMuon) &&
	  aMuon->caloCompatibility() > 0.89 &&
	  trackOK &&
	  iTrack->chi2()/iTrack->ndof() < 1.8 &&
          p.pixelLayersWithMeasurement() > 1 &&
	  fabs(iTrack->dxy(RefVtx)) < 3.0 &&
          fabs(iTrack->dz(RefVtx)) < 15.0 );
}

bool 
JPsiAnalyzerPAT::selDimuon(const pat::CompositeCandidate* aCand) {
  
  if (!_applyDiMuoncuts) return true;
  return ( aCand->userFloat("vProb") > 0.01 );
}

int 
JPsiAnalyzerPAT::getJpsiVarType(const double jpsivar, vector<double> vectbin) {

  for(unsigned int i=0;i<vectbin.size()-1;i++) {
    if(jpsivar > vectbin[i] && jpsivar < vectbin[i+1]) return i+1;
  }

  return -999;
}

// reset the global DataSet variables
void
JPsiAnalyzerPAT::resetDSVariables(){

    //reset J/psi RECO variables
    // JpsiMass=-9999.;
    // JpsiPt=-9999.;
    // JpsiRap=-9999.;
    JpsiCharge=-9999;
    // JpsiPx=-9999.;
    // JpsiPy=-9999.;
    // JpsiPz=-9999.;
    Jpsict=-9999.;
    JpsictErr=-9999.;
    Jpsict_Gen=-9999.;
    JpsiVprob=-9999.;
    JpsiDistM1=-9999.;
    JpsiDphiM1=-9999.; 
    JpsiDrM1=-9999.;
    JpsiDistM2=-9999.;
    JpsiDphiM2=-9999.;
    JpsiDrM2=-9999.;

    JpsiType=-1;

    //reset MUON RECO variables
    /* muPosPx=-9999.;
    muPosPy=-9999.;
    muPosPz=-9999.;
    muNegPx=-9999.;
    muNegPy=-9999.;
    muNegPz=-9999.;*/
    if(_isMC){
        MCType=-1;

        //reset J/psi GEN variables
        /* JpsiMass_Gen=-9999.;
        JpsiPt_Gen=-9999.;
        JpsiRap_Gen=-9999.;
        JpsiPx_Gen=-9999.;
        JpsiPy_Gen=-9999.;
        JpsiPz_Gen=-9999.; */

        //reset MUON GEN variables
        /* muPosPx_Gen=-9999.;
        muPosPy_Gen=-9999.;
        muPosPz_Gen=-9999.;
        muNegPx_Gen=-9999.;
        muNegPy_Gen=-9999.;
        muNegPz_Gen=-9999.; */
    }

    //reset EVENT information
    eventNb= 0 ;
    runNb= 0 ;
    nPriVtx= 0 ;
    lumiBlock= 0 ;

    //reset Trigger Variables
    for(std::map< std::string, int >::iterator clearIt= mapTriggerNameToIntFired_.begin(); clearIt != mapTriggerNameToIntFired_.end(); clearIt++){
        clearIt->second=0;
    }
    for(std::map< std::string, int >::iterator clearIt= mapTriggerNameToPrescaleFac_.begin(); clearIt != mapTriggerNameToPrescaleFac_.end(); clearIt++){
        clearIt->second=1;
    }

    JpsiP->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    muPosP->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    muNegP->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    JpsiP_Gen->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    muPosP_Gen->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    muNegP_Gen->SetPtEtaPhiM(-999.,-999.,-999., 999.);
    
}

//! fill Generator Information
void
JPsiAnalyzerPAT::analyzeGenerator(const edm::Handle<reco::GenParticleCollection>& genParticles)
{
    using namespace trigger;

    std::vector < const reco::Candidate* > genMuons;
    //bool genjpsi= false;
    reco::Candidate::size_type nrD;

    //int count= 0;
    for( size_t i = 0; i < genParticles->size(); ++ i )
    {
        // std::cout << "analyzeGenerator: " << i << std::endl;
        const reco::Candidate & cand = (*genParticles)[ i ];
        int Mc_particleID = cand.pdgId();
        if (abs(Mc_particleID) == _oniaPDG && cand.status()==2 )//&& cand.pt() >= 1)
        {
//          std::cout << "------::analyzeGenerator:: gen JPsi's: pt=" << cand.pt() << "; eta=" << cand.eta() << "; phi=" << cand.phi() << std::endl;

            //Fill global TTree variables
            // JpsiMass_Gen=cand.mass();
            // JpsiPt_Gen=cand.pt();
            // JpsiRap_Gen=cand.rapidity();
            // JpsiPx_Gen=cand.px();
            // JpsiPy_Gen=cand.py();
            // JpsiPz_Gen=cand.pz();
          Double_t enGen = sqrt(cand.px()*cand.px() + cand.py()*cand.py() + cand.pz()*cand.pz() + cand.mass()*cand.mass());
	  JpsiP_Gen->SetPxPyPzE(cand.px(),cand.py(),cand.pz(),enGen);

            //Jpsict_Gen=10.*cand.userFloat("ppdlTrue"));

            nrD= cand.numberOfDaughters();
            int count_muon=0;
            for(reco::Candidate::size_type t=0; t < nrD; t++){
                const reco::Candidate* muon= cand.daughter(t);
                int pID = muon->pdgId();
//              std::cout << "------::analyzeGenerator:: gen JPsi's daughter pdgId: " << pID << std::endl;

                if (abs(pID) == 13 && cand.daughter(t)->status()==1)
                {
                    genMuons.push_back(muon);
//                  std::cout << "------::analyzeGenerator:: gen JPsi's daughter #: " << count_muon << std::endl;
//                  std::cout << " muon" << count_muon << " pt=     " << muon->pt() << std::endl;
//                  std::cout << " muon" << count_muon << " eta=     " << muon->eta() << std::endl;
//                  std::cout << " moun" << count_muon << " phi=     " << muon->phi() << std::endl;
                    count_muon++;
                }
            }


            if ( genMuons.empty() ) break;

            const reco::Candidate* muon1= genMuons.front();
            const reco::Candidate* muon2= genMuons.back();

            // look for opposite charge gen muon pair
            if (muon1->charge()*muon2->charge() <= 0){
                const reco::Candidate *muonPos = 0, *muonNeg = 0;

                if(muon1->charge() > 0){ muonPos = muon1; muonNeg = muon2;}
                else if(muon1->charge() < 0){ muonPos = muon2; muonNeg = muon1;}

                float f_muPosPx, f_muPosPy, f_muPosPz;
                float f_muNegPx, f_muNegPy, f_muNegPz;

                f_muPosPx = muonPos->px();
                f_muPosPy = muonPos->py();
                f_muPosPz = muonPos->pz();

                f_muNegPx = muonNeg->px();
                f_muNegPy = muonNeg->py();
                f_muNegPz = muonNeg->pz();

                // fill global TTree variables - gen muon
                // muPosPx_Gen=muonPos->px();
                // muPosPy_Gen=muonPos->py();
                // muPosPz_Gen=muonPos->pz();

                // muNegPx_Gen=muonNeg->px();
                // muNegPy_Gen=muonNeg->py();
                // muNegPz_Gen=muonNeg->pz();

                // fill Polarization variables - gen muons
                Double_t muMass = 0.105658;

                Double_t enMuPos = sqrt(f_muPosPx*f_muPosPx + f_muPosPy*f_muPosPy + f_muPosPz*f_muPosPz + muMass*muMass);
                // TLorentzVector *muPos = new TLorentzVector();
                muPosP_Gen->SetPxPyPzE(f_muPosPx, f_muPosPy, f_muPosPz, enMuPos);

                Double_t enMuNeg = sqrt(f_muNegPx*f_muNegPx + f_muNegPy*f_muNegPy + f_muNegPz*f_muNegPz + muMass*muMass);
                // TLorentzVector *muNeg = new TLorentzVector();
                muNegP_Gen->SetPxPyPzE(f_muNegPx, f_muNegPy, f_muNegPz, enMuNeg);

                //! Fill Polarization Variables;
                // std::vector< float > thisCosTh, thisPhi;
                // thisCosTh.resize(6); thisPhi.resize(6);
                // this->calcPol(*muPosP_Gen, *muNegP_Gen, thisCosTh, thisPhi);
     
            }
        } // end loop over genParticles
    }
}

void
JPsiAnalyzerPAT::hltReport(const edm::Event &iEvent ,const edm::EventSetup& iSetup)
{

    std::map<std::string, bool> mapTriggernameToTriggerFired;
    std::map<std::string, unsigned int> mapTriggernameToHLTbit;
    // std::map<std::string, unsigned int> mapTriggerNameToPrescaleFac;

    for(std::vector<std::string>::const_iterator it= HLTBitNames_.begin(); it !=HLTBitNames_.end(); ++it){
        mapTriggernameToTriggerFired[*it]=false;
        mapTriggernameToHLTbit[*it]=1000;
        // mapTriggerNameToPrescaleFac[*it]=0;
    }

    // HLTConfigProvider
    if ( hltConfigInit_ ) {
        
        //! Use HLTConfigProvider
      const unsigned int n= hltConfig_.size();
      for (std::map<std::string, unsigned int>::iterator it = mapTriggernameToHLTbit.begin(); it != mapTriggernameToHLTbit.end(); it++) {
	unsigned int triggerIndex= hltConfig_.triggerIndex( it->first );
	if (triggerIndex >= n) {
	  //std::cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerName " << it->first << " not available in config!" << std::endl;
	}
	else {
	  it->second= triggerIndex;
	  //std::cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerName " << it->first << " available in config!" << std::endl;
	}
      }
    }
    
    // Get Trigger Results
    try {
      iEvent.getByLabel( tagTriggerResults_, handleTriggerResults_ );
      //cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResult is present in current event" << endl;
    }
    catch(...) {
      //cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResults NOT present in current event" << endl;
    }
    if ( handleTriggerResults_.isValid() ){
      //cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResults IS valid in current event" << endl;
      
      // loop over Trigger Results to check if paths was fired
      for(std::vector< std::string >::iterator itHLTNames= HLTBitNames_.begin(); itHLTNames != HLTBitNames_.end(); itHLTNames++){
	const std::string triggerPathName =  *itHLTNames;
	//std::cout << "[FloJPsiAnalyzer::hltReport] --- TriggerName --- TriggerName LOOP" << std::endl;
	
	if ( mapTriggernameToHLTbit[triggerPathName] < 1000 ) {
	  if (handleTriggerResults_->accept( mapTriggernameToHLTbit[triggerPathName] ) ){
          //std::cout << "[FloJPsiAnalyzer::hltReport] --- TriggerName " << triggerPathName << " fired!" << std::endl;
	    mapTriggerNameToIntFired_[triggerPathName] = 3;
	  }

	 //-------prescale factor------------
	  if (!_isMC) {
	    const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,triggerPathName));
	    //std::cout << "[FloJPsiAnalyzer::prescalvalues] --- TriggerName"<<triggerPathName<<" prescales first "<< prescales.first <<" prescales second "<< prescales.second <<std::endl;
	    mapTriggerNameToPrescaleFac_[triggerPathName] = prescales.first * prescales.second;
	  }
	}
      }
    } else cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerResults NOT valid in current event" << endl;
}

void
JPsiAnalyzerPAT::matchMuonToHlt(const pat::Muon* muon1, const pat::Muon* muon2)
{

    std::string HLTL3MuCollName = "hltL3MuonCandidates::" + tagTriggerResults_.process();
    std::string HLTL2MuCollName = "hltL2MuonCandidates::" + tagTriggerResults_.process();
    std::string HLTTrackCollName = "hltMuTrackJpsiCtfTrackCands::" + tagTriggerResults_.process();
    std::string HLTTkMuCollName = "hltMuTkMuJpsiTrackerMuonCands::" + tagTriggerResults_.process();
    
    //! Loop over Trigger Paths and match muons to last Filter/collection
    for ( std::map<std::string, int>::iterator it = mapTriggerNameToIntFired_.begin(); it != mapTriggerNameToIntFired_.end(); it ++ ) {

        std::string triggerName = it->first;

        //! just use Triggers which are in TriggerResults; value == 3
        if ( it->second != 3 ) continue;

        std::string hltLastFilterName = mapTriggerToLastFilter_[triggerName];

        const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter( hltLastFilterName );
        const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter( hltLastFilterName );
        bool pass1 = mu1HLTMatches.size() > 0;
        bool pass2 = mu2HLTMatches.size() > 0; 

        // treat "MuX_TrackX" Trigger separately: Match by Tracker collection: hltMuTrackJpsiCtfTrackCands
	std::vector<std::string>::iterator theCheck = std::find(HLTBitNames_MuTrack.begin(),HLTBitNames_MuTrack.end(),triggerName);
	if (theCheck != HLTBitNames_MuTrack.end()) {
                bool matchedMu3[2] = {false, false}, matchedTrack[2] = {false, false};
                for (unsigned k = 0; k < mu1HLTMatches.size(); ++k) {
		    if (mu1HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[0] = true;	     
                    if (mu1HLTMatches[k].collection() == HLTTrackCollName ) matchedTrack[0] = true;
                }
                for (unsigned k = 0; k < mu2HLTMatches.size(); ++k) {
                    if (mu2HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[1] = true;
                    if (mu2HLTMatches[k].collection() == HLTTrackCollName ) matchedTrack[1] = true;
                }
                if( matchedMu3[0] && matchedTrack[1] )
		  mapTriggerNameToIntFired_[triggerName] = 1;
		else if( matchedMu3[1] && matchedTrack[0] )
		  mapTriggerNameToIntFired_[triggerName] = -1;
		if( matchedMu3[0] && matchedTrack[1] && matchedMu3[1] && matchedTrack[0] )
		  mapTriggerNameToIntFired_[triggerName] = 2;
        }

//         // treat "MuX_TkMuX" Trigger separately: Match by Tracker collection:hltMuTkMuJpsiTrackerMuonCands
	theCheck = std::find(HLTBitNames_MuTkMu.begin(),HLTBitNames_MuTkMu.end(),triggerName);
	if (theCheck != HLTBitNames_MuTkMu.end()) {

	  bool matchedMu3[2] = {false, false}, matchedTrack[2] = {false, false};
	  for (unsigned k = 0; k < mu1HLTMatches.size(); ++k) {
	    if (mu1HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[0] = true;
	    if (mu1HLTMatches[k].collection() == HLTTkMuCollName ) matchedTrack[0] = true;
	  }
	  for (unsigned k = 0; k < mu2HLTMatches.size(); ++k) {
	    if (mu2HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[1] = true;
	    if (mu2HLTMatches[k].collection() == HLTTkMuCollName ) matchedTrack[1] = true;
	  }
	  if( matchedMu3[0] && matchedTrack[1] )
	    mapTriggerNameToIntFired_[triggerName] = 1;
	  else if( matchedMu3[1] && matchedTrack[0] )
	    mapTriggerNameToIntFired_[triggerName] = -1;
	  if( matchedMu3[0] && matchedTrack[1] && matchedMu3[1] && matchedTrack[0] )
	    mapTriggerNameToIntFired_[triggerName] = 2;
        }
	
	// treat "MuX_L2MuX" Trigger separately: Match by L2 collection: hltL2MuonCandidates and on a different SaveTag'ed filter
        theCheck = std::find(HLTBitNames_MuL2Mu.begin(),HLTBitNames_MuL2Mu.end(),triggerName);
	if (theCheck != HLTBitNames_MuL2Mu.end()) {
        
	        std::string triggerNameSpecial = triggerName + "_special";
		std::string hltLastFilterName2 = mapTriggerToLastFilter_[triggerNameSpecial];
		
		const pat::TriggerObjectStandAloneCollection mu1Level2Matches = muon1->triggerObjectMatchesByFilter( hltLastFilterName2 );
		const pat::TriggerObjectStandAloneCollection mu2Level2Matches = muon2->triggerObjectMatchesByFilter( hltLastFilterName2 );
    
                bool matchedMu3[2] = {false, false}, matchedMu2[2] = {false, false};
                for (unsigned k = 0; k < mu1HLTMatches.size(); ++k) {
                    if (mu1HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[0] = true;
		}
		for (unsigned k = 0; k < mu1Level2Matches.size(); ++k) {
                    if (mu1Level2Matches[k].collection() == HLTL2MuCollName ) matchedMu2[0] = true;                }
                for (unsigned k = 0; k < mu2HLTMatches.size(); ++k) {
                    if (mu2HLTMatches[k].collection() == HLTL3MuCollName ) matchedMu3[1] = true;
		}
		for (unsigned k = 0; k < mu2Level2Matches.size(); ++k) {
                    if (mu2Level2Matches[k].collection() == HLTL2MuCollName ) matchedMu2[1] = true;
                }
                if( matchedMu3[0] && matchedMu2[1] )
		  mapTriggerNameToIntFired_[triggerName] = 1;
		else if( matchedMu3[1] && matchedMu2[0] )
		  mapTriggerNameToIntFired_[triggerName] = -1;
		if( matchedMu3[0] && matchedMu2[1] && matchedMu3[1] && matchedMu2[0] )
		  mapTriggerNameToIntFired_[triggerName] = 2;
        }

        // All the other Paths match by last filter:
	theCheck = std::find(HLTBitNames_DoubleMu.begin(),HLTBitNames_DoubleMu.end(),triggerName);
	if (theCheck != HLTBitNames_DoubleMu.end() && pass1 == true && pass2 == true ) mapTriggerNameToIntFired_[triggerName] = 1;
	theCheck = std::find(HLTBitNames_SingleMu.begin(),HLTBitNames_SingleMu.end(),triggerName);
	if (theCheck != HLTBitNames_SingleMu.end() && (pass1 == true || pass2 == true) ) mapTriggerNameToIntFired_[triggerName] = 1;
    }
	
}

void 
JPsiAnalyzerPAT::muonStationDistance (const pat::CompositeCandidate* aCand) 
{
  
   const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon1"));
   const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon2"));
   const reco::Candidate &d1 = *muon1; 
   const reco::Candidate &d2 = *muon2;
   const reco::RecoCandidate *mu1 = dynamic_cast<const reco::RecoCandidate *>(&d1);
   const reco::RecoCandidate *mu2 = dynamic_cast<const reco::RecoCandidate *>(&d2);
    
   // Propagate to station 1
   TrajectoryStateOnSurface prop1_Stat1 = prop1_.extrapolate(*mu1);
   TrajectoryStateOnSurface prop2_Stat1 = prop1_.extrapolate(*mu2);
   if (prop1_Stat1.isValid() && prop2_Stat1.isValid()) {
     JpsiDphiM1 = deltaPhi<float>(prop1_Stat1.globalPosition().phi(), prop2_Stat1.globalPosition().phi());
     JpsiDrM1   = hypot(JpsiDphiM1, std::abs<float>(prop1_Stat1.globalPosition().eta() - prop2_Stat1.globalPosition().eta()));
     JpsiDistM1 = (prop1_Stat1.globalPosition()-prop2_Stat1.globalPosition()).mag();
   }
   // Propagate to station 2
   TrajectoryStateOnSurface prop1_Stat2 = prop2_.extrapolate(*mu1);
   TrajectoryStateOnSurface prop2_Stat2 = prop2_.extrapolate(*mu2);
   if (prop1_Stat2.isValid() && prop2_Stat2.isValid()) {
     JpsiDphiM2 = deltaPhi<float>(prop1_Stat2.globalPosition().phi(), prop2_Stat2.globalPosition().phi());
     JpsiDrM2   = hypot(JpsiDphiM2, std::abs<float>(prop1_Stat2.globalPosition().eta() - prop2_Stat2.globalPosition().eta()));
     JpsiDistM2 = (prop1_Stat2.globalPosition()-prop2_Stat2.globalPosition()).mag();
   }

}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiAnalyzerPAT);
