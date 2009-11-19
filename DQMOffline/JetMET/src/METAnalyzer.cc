/*
 *  See header file for a description of this class.
 *
 *  $Date: 2009/11/12 17:29:36 $
 *  $Revision: 1.2 $
 *  \author A.Apresyan - Caltech
 */

#include "DQMOffline/JetMET/interface/METAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/MET.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include <string>
using namespace std;
using namespace edm;
using namespace reco;
using namespace math;

// ***********************************************************
METAnalyzer::METAnalyzer(const edm::ParameterSet& pSet) {

  parameters = pSet;

}

// ***********************************************************
METAnalyzer::~METAnalyzer() { }

void METAnalyzer::beginJob(edm::EventSetup const& iSetup,DQMStore * dbe) {

  evtCounter = 0;
  metname = "METAnalyzer";

  // trigger information
  HLTPathsJetMBByName_ = parameters.getParameter<std::vector<std::string > >("HLTPathsJetMB");

  _hlt_HighPtJet = parameters.getParameter<std::string>("HLT_HighPtJet");
  _hlt_LowPtJet  = parameters.getParameter<std::string>("HLT_LowPtJet");
  _hlt_HighMET   = parameters.getParameter<std::string>("HLT_HighMET");
  _hlt_LowMET    = parameters.getParameter<std::string>("HLT_LowMET");
  _hlt_Ele       = parameters.getParameter<std::string>("HLT_Ele");
  _hlt_Muon      = parameters.getParameter<std::string>("HLT_Muon");

  // MET information
  theMETCollectionLabel       = parameters.getParameter<edm::InputTag>("METCollectionLabel");
  _source                       = parameters.getParameter<std::string>("Source");

  // Other data collections
  HcalNoiseRBXCollectionTag   = parameters.getParameter<edm::InputTag>("HcalNoiseRBXCollection");
  HcalNoiseSummaryTag         = parameters.getParameter<edm::InputTag>("HcalNoiseSummary");
  theJetCollectionLabel       = parameters.getParameter<edm::InputTag>("JetCollectionLabel");

  // misc
  _verbose     = parameters.getParameter<int>("verbose");
  _etThreshold = parameters.getParameter<double>("etThreshold"); // MET threshold
  _allhist     = parameters.getParameter<bool>("allHist");       // Full set of monitoring histograms
  _allSelection= parameters.getParameter<bool>("allSelection");  // Plot with all sets of event selection

  _highPtJetThreshold = parameters.getParameter<double>("HighPtJetThreshold"); // High Pt Jet threshold
  _lowPtJetThreshold = parameters.getParameter<double>("LowPtJetThreshold");   // Low Pt Jet threshold
  _highMETThreshold = parameters.getParameter<double>("HighMETThreshold");     // High MET threshold
  _lowMETThreshold = parameters.getParameter<double>("LowMETThreshold");       // Low MET threshold

  //
  jetID = new reco::helper::JetIDHelper(parameters.getParameter<ParameterSet>("JetIDParams"));

  // DQStore stuff
  LogTrace(metname)<<"[METAnalyzer] Parameters initialization";
  std::string DirName = "JetMET/MET/"+_source;
  dbe->setCurrentFolder(DirName);

  metME = dbe->book1D("metReco", "metReco", 4, 1, 5);
  metME->setBinLabel(2,"MET",1);

  _dbe = dbe;

  _FolderNames.push_back("All");
  _FolderNames.push_back("Cleanup");
  _FolderNames.push_back("HcalNoiseFilter");
  _FolderNames.push_back("HcalNoiseFilterTight");
  _FolderNames.push_back("JetID");
  _FolderNames.push_back("JetIDTight");

  for (std::vector<std::string>::const_iterator ic = _FolderNames.begin(); 
       ic != _FolderNames.end(); ic++){
    if (*ic=="All")                  bookMESet(DirName+"/"+*ic);
    if (*ic=="Cleanup")              bookMESet(DirName+"/"+*ic);
    if (_allSelection){
    if (*ic=="HcalNoiseFilter")      bookMESet(DirName+"/"+*ic);
    if (*ic=="HcalNoiseFilterTight") bookMESet(DirName+"/"+*ic);
    if (*ic=="JetID")                bookMESet(DirName+"/"+*ic);
    if (*ic=="JetIDTight")           bookMESet(DirName+"/"+*ic);
    }
  }
}

// ***********************************************************
void METAnalyzer::endJob() {

  delete jetID;

}

// ***********************************************************
void METAnalyzer::bookMESet(std::string DirName)
{

  bool bLumiSecPlot=false;
  if (DirName.find("All")!=std::string::npos) bLumiSecPlot=true;

  bookMonitorElement(DirName,bLumiSecPlot);

  if (_hlt_HighPtJet.size()){
    bookMonitorElement(DirName+"/"+"HighPtJet",false);
    meTriggerName_HighPtJet = _dbe->bookString("triggerName_HighPtJet", _hlt_HighPtJet);
  }  

  if (_hlt_LowPtJet.size()){
    bookMonitorElement(DirName+"/"+"LowPtJet",false);
    meTriggerName_LowPtJet = _dbe->bookString("triggerName_LowPtJet", _hlt_LowPtJet);
  }

  if (_hlt_HighMET.size()){
    bookMonitorElement(DirName+"/"+"HighMET",false);
    meTriggerName_HighMET = _dbe->bookString("triggerName_HighMET", _hlt_HighMET);
  }

  if (_hlt_LowMET.size()){
    bookMonitorElement(DirName+"/"+"LowMET",false);
    meTriggerName_LowMET = _dbe->bookString("triggerName_LowMET", _hlt_LowMET);
  }

  if (_hlt_Ele.size()){
    bookMonitorElement(DirName+"/"+"Ele",false);
    meTriggerName_Ele = _dbe->bookString("triggerName_Ele", _hlt_Ele);
  }

  if (_hlt_Muon.size()){
    bookMonitorElement(DirName+"/"+"Muon",false);
    meTriggerName_Muon = _dbe->bookString("triggerName_Muon", _hlt_Muon);
  }

}

// ***********************************************************
void METAnalyzer::bookMonitorElement(std::string DirName, bool bLumiSecPlot=false)
{

  if (_verbose) std::cout << "booMonitorElement " << DirName << std::endl;
  _dbe->setCurrentFolder(DirName);
 
  meNevents              = _dbe->book1D("METTask_Nevents", "METTask_Nevents"   ,1,0,1);
  meMEx                = _dbe->book1D("METTask_MEx",   "METTask_MEx"   ,500,-500,500);
  meMEy                = _dbe->book1D("METTask_MEy",   "METTask_MEy"   ,500,-500,500);
  meEz                 = _dbe->book1D("METTask_Ez",    "METTask_Ez"    ,500,-500,500);
  meMETSig             = _dbe->book1D("METTask_METSig","METTask_METSig",51,0,51);
  meMET                = _dbe->book1D("METTask_MET",   "METTask_MET"   ,500,0,1000);
  meMETPhi             = _dbe->book1D("METTask_METPhi","METTask_METPhi",80,-TMath::Pi(),TMath::Pi());
  meSumET              = _dbe->book1D("METTask_SumET", "METTask_SumET" ,500,0,2000);

  meNeutralEMFraction  = _dbe->book1D("METTask_NeutralEMFraction", "METTask_NeutralEMFraction" ,50,0.,1.);
  meNeutralHadFraction = _dbe->book1D("METTask_NeutralHadFraction","METTask_NeutralHadFraction",50,0.,1.);
  meChargedEMFraction  = _dbe->book1D("METTask_ChargedEMFraction", "METTask_ChargedEMFraction" ,50,0.,1.);
  meChargedHadFraction = _dbe->book1D("METTask_ChargedHadFraction","METTask_ChargedHadFraction",50,0.,1.);
  meMuonFraction       = _dbe->book1D("METTask_MuonFraction",      "METTask_MuonFraction"      ,50,0.,1.);

  meMETIonFeedbck      = _dbe->book1D("METTask_METIonFeedbck", "METTask_METIonFeedbck" ,500,0,1000);
  meMETHPDNoise        = _dbe->book1D("METTask_METHPDNoise",   "METTask_METHPDNoise"   ,500,0,1000);
  meMETRBXNoise        = _dbe->book1D("METTask_METRBXNoise",   "METTask_METRBXNoise"   ,500,0,1000);

  if (_allhist){
    if (bLumiSecPlot){
      meMExLS              = _dbe->book2D("METTask_MEx_LS","METTask_MEx_LS",200,-200,200,50,0.,500.);
      meMEyLS              = _dbe->book2D("METTask_MEy_LS","METTask_MEy_LS",200,-200,200,50,0.,500.);
    }
  }
}

// ***********************************************************
void METAnalyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{

  //
  //--- htlConfig_
  
//   hltConfig_.init(processname_);
//   if (!hltConfig_.init(processname_)) {
//     processname_ = "FU";
//     if (!hltConfig_.init(processname_)){
//       LogDebug("METAnalyzer") << "HLTConfigProvider failed to initialize.";
//     }
//   }

//   if (_verbose) std::cout << hltConfig_.triggerIndex(_hlt_HighPtJet) << std::endl;
//   if (_verbose) std::cout << hltConfig_.triggerIndex(_hlt_LowPtJet)  << std::endl;
//   if (_verbose) std::cout << hltConfig_.triggerIndex(_hlt_HighMET)   << std::endl;
//   if (_verbose) std::cout << hltConfig_.triggerIndex(_hlt_LowMET)    << std::endl;
//   if (_verbose) std::cout << hltConfig_.triggerIndex(_hlt_Ele)       << std::endl;
//   if (_verbose) std::cout << hltConfig_.triggerIndex(_hlt_Muon)      << std::endl;

}

// ***********************************************************
void METAnalyzer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup, DQMStore * dbe)
{
  
  //
  //--- Check the time length of the Run from the lumi section plots

  std::string dirName = "JetMET/MET/"+_source+"/";
  _dbe->setCurrentFolder(dirName);

  TH1F* tlumisec;

  MonitorElement *meLumiSec = _dbe->get("aaa");
  meLumiSec = _dbe->get("JetMET/lumisec");

  int totlsec=0;
  double totltime=0.;
  if ( meLumiSec->getRootObject() ) {
    tlumisec = meLumiSec->getTH1F();
    for (int i=0; i<500; i++){
      if (tlumisec->GetBinContent(i+1)) totlsec++;
    }
    totltime = double(totlsec*90); // one lumi sec ~ 90 (sec)
  }

  if (totltime==0.) totltime=1.; 

  //
  //--- Make the integrated plots with rate (Hz)

  for (std::vector<std::string>::const_iterator ic = _FolderNames.begin(); ic != _FolderNames.end(); ic++)
    {

      std::string DirName;
      DirName = dirName+*ic;

      makeRatePlot(DirName,totltime);
      if (_hlt_HighPtJet.size()) makeRatePlot(DirName+"/"+_hlt_HighPtJet,totltime);
      if (_hlt_LowPtJet.size())  makeRatePlot(DirName+"/"+_hlt_LowPtJet,totltime);
      if (_hlt_HighMET.size())   makeRatePlot(DirName+"/"+_hlt_HighMET,totltime);
      if (_hlt_LowMET.size())    makeRatePlot(DirName+"/"+_hlt_LowMET,totltime);
      if (_hlt_Ele.size())       makeRatePlot(DirName+"/"+_hlt_Ele,totltime);
      if (_hlt_Muon.size())      makeRatePlot(DirName+"/"+_hlt_Muon,totltime);

    }
}

// ***********************************************************
void METAnalyzer::makeRatePlot(std::string DirName, double totltime)
{

  _dbe->setCurrentFolder(DirName);
  MonitorElement *meMET = _dbe->get(DirName+"/"+"METTask_MET");

  TH1F* tMET;
  TH1F* tMETRate;

  if ( meMET )
    if ( meMET->getRootObject() ) {
      tMET     = meMET->getTH1F();
      
      // Integral plot & convert number of events to rate (hz)
      tMETRate = (TH1F*) tMET->Clone("METTask_METRate");
      for (int i = tMETRate->GetNbinsX()-1; i>=0; i--){
	tMETRate->SetBinContent(i+1,tMETRate->GetBinContent(i+2)+tMET->GetBinContent(i+1));
      }
      for (int i = 0; i<tMETRate->GetNbinsX(); i++){
	tMETRate->SetBinContent(i+1,tMETRate->GetBinContent(i+1)/double(totltime));
      }      

      meMETRate      = _dbe->book1D("METTask_METRate",tMETRate);
      
    }

}

// ***********************************************************
void METAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
			    const edm::TriggerResults& triggerResults) {

  if (_verbose) std::cout << "METAnalyzer analyze" << std::endl;

  LogTrace(metname)<<"[METAnalyzer] Analyze MET";

  metME->Fill(2);

  // ==========================================================  
  // Trigger information 
  //
  _trig_JetMB=0;
  _trig_HighPtJet=0;
  _trig_LowPtJet=0;
  _trig_HighMET=0;
  _trig_LowMET=0;
  if(&triggerResults) {   

    /////////// Analyzing HLT Trigger Results (TriggerResults) //////////

    //
    //
    // Check how many HLT triggers are in triggerResults 
    int ntrigs = triggerResults.size();
    if (_verbose) std::cout << "ntrigs=" << ntrigs << std::endl;
    
    //
    //
    // If index=ntrigs, this HLT trigger doesn't exist in the HLT table for this data.
    edm::TriggerNames triggerNames; // TriggerNames class
    triggerNames.init(triggerResults);
    
    //
    //
    // count number of requested Jet or MB HLT paths which have fired
    for (unsigned int i=0; i!=HLTPathsJetMBByName_.size(); i++) {
      unsigned int triggerIndex = triggerNames.triggerIndex(HLTPathsJetMBByName_[i]);
      if (triggerIndex<triggerResults.size()) {
	if (triggerResults.accept(triggerIndex)) {
	  _trig_JetMB++;
	}
      }
    }
    // for empty input vectors (n==0), take all HLT triggers!
    if (HLTPathsJetMBByName_.size()==0) _trig_JetMB=triggerResults.size()-1;

    //
    if (_verbose) std::cout << "triggerNames size" << " " << triggerNames.size() << std::endl;
    if (_verbose) std::cout << _hlt_HighPtJet << " " << triggerNames.triggerIndex(_hlt_HighPtJet) << std::endl;
    if (_verbose) std::cout << _hlt_LowPtJet  << " " << triggerNames.triggerIndex(_hlt_LowPtJet)  << std::endl;
    if (_verbose) std::cout << _hlt_HighMET   << " " << triggerNames.triggerIndex(_hlt_HighMET)   << std::endl;
    if (_verbose) std::cout << _hlt_LowMET    << " " << triggerNames.triggerIndex(_hlt_LowMET)    << std::endl;
    if (_verbose) std::cout << _hlt_Ele       << " " << triggerNames.triggerIndex(_hlt_Ele)       << std::endl;
    if (_verbose) std::cout << _hlt_Muon      << " " << triggerNames.triggerIndex(_hlt_Muon)      << std::endl;

    if (triggerNames.triggerIndex(_hlt_HighPtJet) != triggerNames.size() &&
	triggerResults.accept(triggerNames.triggerIndex(_hlt_HighPtJet))) _trig_HighPtJet=1;

    if (triggerNames.triggerIndex(_hlt_LowPtJet)  != triggerNames.size() &&
	triggerResults.accept(triggerNames.triggerIndex(_hlt_LowPtJet)))  _trig_LowPtJet=1;

    if (triggerNames.triggerIndex(_hlt_HighMET)   != triggerNames.size() &&
        triggerResults.accept(triggerNames.triggerIndex(_hlt_HighMET)))   _trig_HighMET=1;

    if (triggerNames.triggerIndex(_hlt_LowMET)    != triggerNames.size() &&
        triggerResults.accept(triggerNames.triggerIndex(_hlt_LowMET)))    _trig_LowMET=1;

    if (triggerNames.triggerIndex(_hlt_Ele)       != triggerNames.size() &&
        triggerResults.accept(triggerNames.triggerIndex(_hlt_Ele)))       _trig_Ele=1;

    if (triggerNames.triggerIndex(_hlt_Muon)      != triggerNames.size() &&
        triggerResults.accept(triggerNames.triggerIndex(_hlt_Muon)))      _trig_Muon=1;
    
  } else {

    edm::LogInfo("MetAnalyzer") << "TriggerResults::HLT not found, "
      "automatically select events"; 

    // TriggerResults object not found. Look at all events.    
    _trig_JetMB=1;
  }

  // ==========================================================
  // MET information
  
  // **** Get the MET container  
  edm::Handle<reco::METCollection> tcmetcoll;
  iEvent.getByLabel(theMETCollectionLabel, tcmetcoll);
  
  if(!tcmetcoll.isValid()) return;

  const METCollection *tcmetcol = tcmetcoll.product();
  const MET *tcmet;
  tcmet = &(tcmetcol->front());
    
  LogTrace(metname)<<"[METAnalyzer] Call to the MET analyzer";

  // ==========================================================
  //
  edm::Handle<HcalNoiseRBXCollection> HRBXCollection;
  iEvent.getByLabel(HcalNoiseRBXCollectionTag,HRBXCollection);
  if (!HRBXCollection.isValid()) {
    LogDebug("") << "METAnalyzer: Could not find HcalNoiseRBX Collection" << std::endl;
    if (_verbose) std::cout << "METAnalyzer: Could not find HcalNoiseRBX Collection" << std::endl;
  }
  
  edm::Handle<HcalNoiseSummary> HNoiseSummary;
  iEvent.getByLabel(HcalNoiseSummaryTag,HNoiseSummary);
  if (!HNoiseSummary.isValid()) {
    LogDebug("") << "METAnalyzer: Could not find Hcal NoiseSummary product" << std::endl;
    if (_verbose) std::cout << "METAnalyzer: Could not find Hcal NoiseSummary product" << std::endl;
  }

  edm::Handle<reco::CaloJetCollection> caloJets;
  iEvent.getByLabel(theJetCollectionLabel, caloJets);
  if (!caloJets.isValid()) {
    LogDebug("") << "METAnalyzer: Could not find jet product" << std::endl;
    if (_verbose) std::cout << "METAnalyzer: Could not find jet product" << std::endl;
  }

  // ==========================================================
  // MET sanity check

  //   if (_source=="MET") validateMET(*tcmet, tcCandidates);
  
  // ==========================================================
  // JetID 

  if (_verbose) std::cout << "JetID starts" << std::endl;
  
  //
  // --- Loose cuts, not  specific for now!
  //
  bool bJetID=true;
  for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); 
       cal!=caloJets->end(); ++cal){ 
    jetID->calculate(iEvent, *cal);
    if (_verbose) std::cout << jetID->n90Hits() << " " 
			    << jetID->restrictedEMF() << " "
			    << cal->pt() << std::endl;
    if (cal->pt()>10.){
      //
      // for all regions
      if (jetID->n90Hits()<2)  bJetID=false; 
      if (jetID->fHPD()>=0.98) bJetID=false; 
      //if (jetID->restrictedEMF()<0.01) bJetID=false; 
      //
      // for non-forward
      if (fabs(cal->eta())<2.55){
	if (cal->emEnergyFraction()<=0.01) bJetID=false; 
      }
      // for forward
      else {
	if (cal->emEnergyFraction()<=-0.9) bJetID=false; 
	if (cal->pt()>80.){
	if (cal->emEnergyFraction()>= 1.0) bJetID=false; 
	}
      } // forward vs non-forward
    }   // pt>10 GeV/c
  }     // calor-jets loop
 
  //
  // --- Tight cuts
  //
  bool bJetIDTight=true;
  bJetIDTight=bJetID;
  for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); 
       cal!=caloJets->end(); ++cal){
    jetID->calculate(iEvent, *cal);
    if (cal->pt()>25.){
      //
      // for all regions
      if (jetID->fHPD()>=0.95) bJetIDTight=false; 
      //
      // for 1.0<|eta|<1.75
      if (fabs(cal->eta())>=1.00 && fabs(cal->eta())<1.75){
	if (cal->pt()>80. && cal->emEnergyFraction()>=1.) bJetIDTight=false; 
      }
      //
      // for 1.75<|eta|<2.55
      else if (fabs(cal->eta())>=1.75 && fabs(cal->eta())<2.55){
	if (cal->pt()>80. && cal->emEnergyFraction()>=1.) bJetIDTight=false; 
      }
      //
      // for 2.55<|eta|<3.25
      else if (fabs(cal->eta())>=2.55 && fabs(cal->eta())<3.25){
	if (cal->pt()< 50.                   && cal->emEnergyFraction()<=-0.3) bJetIDTight=false; 
	if (cal->pt()>=50. && cal->pt()< 80. && cal->emEnergyFraction()<=-0.2) bJetIDTight=false; 
	if (cal->pt()>=80. && cal->pt()<340. && cal->emEnergyFraction()<=-0.1) bJetIDTight=false; 
	if (cal->pt()>=340.                  && cal->emEnergyFraction()<=-0.1 
                                             && cal->emEnergyFraction()>=0.95) bJetIDTight=false; 
      }
      //
      // for 3.25<|eta|
      else if (fabs(cal->eta())>=3.25){
	if (cal->pt()< 50.                   && cal->emEnergyFraction()<=-0.3
                                             && cal->emEnergyFraction()>=0.90) bJetIDTight=false; 
	if (cal->pt()>=50. && cal->pt()<130. && cal->emEnergyFraction()<=-0.2
                                             && cal->emEnergyFraction()>=0.80) bJetIDTight=false; 
	if (cal->pt()>=130.                  && cal->emEnergyFraction()<=-0.1 
                                             && cal->emEnergyFraction()>=0.70) bJetIDTight=false; 
      }
    }   // pt>10 GeV/c
  }     // calor-jets loop
  
  if (_verbose) std::cout << "JetID ends" << std::endl;


  // ==========================================================
  // HCAL Noise filter
  
  bool bHcalNoiseFilter      = HNoiseSummary->passLooseNoiseFilter();
  bool bHcalNoiseFilterTight = HNoiseSummary->passTightNoiseFilter();

  // ==========================================================
  // Reconstructed MET Information - fill MonitorElements
  
  std::string DirName = "JetMET/MET/"+_source;
  
  for (std::vector<std::string>::const_iterator ic = _FolderNames.begin(); 
       ic != _FolderNames.end(); ic++){
    if (*ic=="All")                                   fillMESet(iEvent, DirName+"/"+*ic, *tcmet);
    if (*ic=="Cleanup" && bHcalNoiseFilter && bJetID) fillMESet(iEvent, DirName+"/"+*ic, *tcmet);
    if (_allSelection) {
    if (*ic=="HcalNoiseFilter"      && bHcalNoiseFilter )       fillMESet(iEvent, DirName+"/"+*ic, *tcmet);
    if (*ic=="HcalNoiseFilterTight" && bHcalNoiseFilterTight )  fillMESet(iEvent, DirName+"/"+*ic, *tcmet);
    if (*ic=="JetID"      && bJetID)                            fillMESet(iEvent, DirName+"/"+*ic, *tcmet);
    if (*ic=="JetIDTight" && bJetIDTight)                       fillMESet(iEvent, DirName+"/"+*ic, *tcmet);
    }
  }
}

// ***********************************************************
void METAnalyzer::fillMESet(const edm::Event& iEvent, std::string DirName, 
			      const reco::MET& tcmet)
{

  _dbe->setCurrentFolder(DirName);

  bool bLumiSecPlot=false;
  if (DirName.find("All")) bLumiSecPlot=true;

  if (_trig_JetMB) fillMonitorElement(iEvent,DirName,"",tcmet, bLumiSecPlot);
  if (_hlt_HighPtJet.size() && _trig_HighPtJet) fillMonitorElement(iEvent,DirName,"HighPtJet",tcmet,false);
  if (_hlt_LowPtJet.size() && _trig_LowPtJet) fillMonitorElement(iEvent,DirName,"LowPtJet",tcmet,false);
  if (_hlt_HighMET.size() && _trig_HighMET) fillMonitorElement(iEvent,DirName,"HighMET",tcmet,false);
  if (_hlt_LowMET.size() && _trig_LowMET) fillMonitorElement(iEvent,DirName,"LowMET",tcmet,false);
  if (_hlt_Ele.size() && _trig_Ele) fillMonitorElement(iEvent,DirName,"Ele",tcmet,false);
  if (_hlt_Muon.size() && _trig_Muon) fillMonitorElement(iEvent,DirName,"Muon",tcmet,false);
}

// ***********************************************************
void METAnalyzer::fillMonitorElement(const edm::Event& iEvent, std::string DirName, 
					 std::string TriggerTypeName, 
					 const reco::MET& tcmet, bool bLumiSecPlot)
{

  if (TriggerTypeName=="HighPtJet") {
    if (!selectHighPtJetEvent(iEvent)) return;
  }
  else if (TriggerTypeName=="LowPtJet") {
    if (!selectLowPtJetEvent(iEvent)) return;
  }
  else if (TriggerTypeName=="HighMET") {
    if (tcmet.pt()<_highMETThreshold) return;
  }
  else if (TriggerTypeName=="LowMET") {
    if (tcmet.pt()<_lowMETThreshold) return;
  }
  else if (TriggerTypeName=="Ele") {
    if (!selectWElectronEvent(iEvent)) return;
  }
  else if (TriggerTypeName=="Muon") {
    if (!selectWMuonEvent(iEvent)) return;
  }
  
// Reconstructed MET Information
  double tcSumET  = tcmet.sumEt();
  double tcMETSig = tcmet.mEtSig();
  double tcEz     = tcmet.e_longitudinal();
  double tcMET    = tcmet.pt();
  double tcMEx    = tcmet.px();
  double tcMEy    = tcmet.py();
  double tcMETPhi = tcmet.phi();

  //
  int myLuminosityBlock;
  //  myLuminosityBlock = (evtCounter++)/1000;
  myLuminosityBlock = iEvent.luminosityBlock();
  //

  if (TriggerTypeName!="") DirName = DirName +"/"+TriggerTypeName;

  if (_verbose) std::cout << "_etThreshold = " << _etThreshold << std::endl;
  if (tcMET>_etThreshold){
    
    meMEx    = _dbe->get(DirName+"/"+"METTask_MEx");    if (meMEx    && meMEx->getRootObject())    meMEx->Fill(tcMEx);
    meMEy    = _dbe->get(DirName+"/"+"METTask_MEy");    if (meMEy    && meMEy->getRootObject())    meMEy->Fill(tcMEy);
    meMET    = _dbe->get(DirName+"/"+"METTask_MET");    if (meMET    && meMET->getRootObject())    meMET->Fill(tcMET);
    meMETPhi = _dbe->get(DirName+"/"+"METTask_METPhi"); if (meMETPhi && meMETPhi->getRootObject()) meMETPhi->Fill(tcMETPhi);
    meSumET  = _dbe->get(DirName+"/"+"METTask_SumET");  if (meSumET  && meSumET->getRootObject())  meSumET->Fill(tcSumET);
    meMETSig = _dbe->get(DirName+"/"+"METTask_METSig"); if (meMETSig && meMETSig->getRootObject()) meMETSig->Fill(tcMETSig);
    meEz     = _dbe->get(DirName+"/"+"METTask_Ez");     if (meEz     && meEz->getRootObject())     meEz->Fill(tcEz);

    meMETIonFeedbck = _dbe->get(DirName+"/"+"METTask_METIonFeedbck");  if (meMETIonFeedbck && meMETIonFeedbck->getRootObject()) meMETIonFeedbck->Fill(tcMET);
    meMETHPDNoise   = _dbe->get(DirName+"/"+"METTask_METHPDNoise");    if (meMETHPDNoise   && meMETHPDNoise->getRootObject())   meMETHPDNoise->Fill(tcMET);
    meMETRBXNoise   = _dbe->get(DirName+"/"+"METTask_METRBXNoise");    if (meMETRBXNoise   && meMETRBXNoise->getRootObject())   meMETRBXNoise->Fill(tcMET);
        
    if (_allhist){
      if (bLumiSecPlot){
	meMExLS = _dbe->get(DirName+"/"+"METTask_MExLS"); if (meMExLS && meMExLS->getRootObject()) meMExLS->Fill(tcMEx,myLuminosityBlock);
	meMEyLS = _dbe->get(DirName+"/"+"METTask_MEyLS"); if (meMEyLS && meMEyLS->getRootObject()) meMEyLS->Fill(tcMEy,myLuminosityBlock);
      }
    } // _allhist
  } // et threshold cut
}

// ***********************************************************
bool METAnalyzer::selectHighPtJetEvent(const edm::Event& iEvent){

  bool return_value=false;

  edm::Handle<reco::CaloJetCollection> caloJets;
  iEvent.getByLabel(theJetCollectionLabel, caloJets);
  if (!caloJets.isValid()) {
    LogDebug("") << "METAnalyzer: Could not find jet product" << std::endl;
    if (_verbose) std::cout << "METAnalyzer: Could not find jet product" << std::endl;
  }

  for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); 
       cal!=caloJets->end(); ++cal){
    if (cal->pt()>_highPtJetThreshold){
      return_value=true;
    }
  }
  
  return return_value;
}

// // ***********************************************************
bool METAnalyzer::selectLowPtJetEvent(const edm::Event& iEvent){

  bool return_value=false;

  edm::Handle<reco::CaloJetCollection> caloJets;
  iEvent.getByLabel(theJetCollectionLabel, caloJets);
  if (!caloJets.isValid()) {
    LogDebug("") << "METAnalyzer: Could not find jet product" << std::endl;
    if (_verbose) std::cout << "METAnalyzer: Could not find jet product" << std::endl;
  }

  for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); 
       cal!=caloJets->end(); ++cal){
    if (cal->pt()>_lowPtJetThreshold){
      return_value=true;
    }
  }

  return return_value;

}

// ***********************************************************
bool METAnalyzer::selectWElectronEvent(const edm::Event& iEvent){

  bool return_value=false;

  /*
    W-electron event selection comes here
   */

  return return_value;

}

// ***********************************************************
bool METAnalyzer::selectWMuonEvent(const edm::Event& iEvent){

  bool return_value=false;

  /*
    W-muon event selection comes here
   */

  return return_value;

}
