////////////////////////////////////////////////////////////////////////////////
//
// JetResponseAnalyzer
// -------------------
//
//            07/04/2008 Kostas Kousouris       <kkousour@fnal.gov>
//                       Roger Wolf             <roger.wolf@cern.ch>
//                       Philipp Schieferdecker <philipp.schieferdecker@cern.ch>
////////////////////////////////////////////////////////////////////////////////


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
 
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "PhysicsTools/Utilities/interface/deltaPhi.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

#include <memory>
#include <vector>
#include <sstream>
#include <cmath>


using namespace std;


////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////

class JetResponseAnalyzer : public edm::EDAnalyzer
{
public:
  // construction/destruction
  explicit JetResponseAnalyzer(const edm::ParameterSet& iConfig);
  virtual ~JetResponseAnalyzer();

private:
  // member functions
  void beginJob(const edm::EventSetup& iSetup);
  void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);
  void endJob(){;}

private:
  // member data
  bool            doHistos_;
  bool            doTree_;
  bool            doFlavor_;
  bool            doJetPt_;
  bool            doRefPt_;
  bool            doRelRsp_;
  bool            doAbsRsp_;
  
  bool            doBalancing_;
  double          deltaRMax_;
  double          deltaPhiMin_;

  edm::InputTag   srcRef_;
  edm::InputTag   srcRefToJetMap_;
  edm::InputTag   srcRefToPartonMap_;
  bool            getFlavorFromMap_;
  double          deltaRPartonMax_;
  
  unsigned int    nRefMax_;
  unsigned int    nBinsRelRsp_;
  unsigned int    nBinsAbsRsp_;
  unsigned int    nBinsPt_;
  unsigned int    nBinsEta_;
  unsigned int    nBinsPhi_;

  double          etaBarrelMin_;
  double          etaBarrelMax_;

  double          relRspMin_;
  double          relRspMax_;
  double          absRspMin_;
  double          absRspMax_;

  vector<double>  binsPt_;
  vector<double>  binsEta_;
  vector<double>  binsPhi_;
  
  vector<TH1F**>  jetPtVsJetPt_;
  vector<TH1F**>  refPtVsRefPt_;
  vector<TH1F**>  jetPtVsRefPt_;
  vector<TH1F**>  refPtVsRefPtBarrel_;
  vector<TH1F**>  jetPtVsRefPtBarrel_;
  vector<TH1F**>  jetEtaVsJetEta_;
  vector<TH1F**>  jetPhiVsJetPhi_;
  vector<TH1F***> jetPtVsJetEtaJetPt_;
  vector<TH1F***> refPtVsJetEtaRefPt_;
  vector<TH1F***> jetPtVsJetEtaRefPt_;
  
  vector<TH1F**>  relRspVsJetPt_;
  vector<TH1F**>  relRspVsRefPt_;
  vector<TH1F**>  relRspVsRefPtBarrel_;
  vector<TH1F**>  relRspVsJetEta_;
  vector<TH1F**>  relRspVsJetPhi_;
  vector<TH1F***> relRspVsJetEtaJetPt_;
  vector<TH1F***> relRspVsJetEtaRefPt_;

  vector<TH1F**>  absRspVsJetPt_;
  vector<TH1F**>  absRspVsRefPt_;
  vector<TH1F**>  absRspVsRefPtBarrel_;
  vector<TH1F**>  absRspVsJetEta_;
  vector<TH1F**>  absRspVsJetPhi_;
  vector<TH1F***> absRspVsJetEtaJetPt_;
  vector<TH1F***> absRspVsJetEtaRefPt_;

  TTree*          tree_;
  unsigned char   nref_;
  int             refpdgid_[100];
  float           refpt_[100];
  float           refeta_[100];
  float           refphi_[100];
  float           jtpt_[100];
  float           jteta_[100];
  float           jtphi_[100];
  float           refdrjt_[100];
  float           refdphijt_[100];

};


////////////////////////////////////////////////////////////////////////////////
// define local methods
////////////////////////////////////////////////////////////////////////////////

/// get the suffix for the histogram name, e.g. JetPt100to150
string getSuffix(const string& varname,int ibin,const vector<double>& bins);

/// get the index of the histogram corresponding to x
int getIndex(double x,const vector<double>& binsx);

/// fill the appropriate histogram (histos), based on x and binsx
void fillHisto(double value,double x,
	       const vector<double>& binsx,const vector<TH1F**>& histos);

/// fill the appropriate histogram (histos), based on pdgid, x and binsx
void fillHisto(int pdgid,double value,double x,
	       const vector<double>& binsx,const vector<TH1F**>& histos);

/// fill the appropriate histogram (histos), based on x, y, binsx, and binsy
void fillHisto(double value,double x,double y,
	       const vector<double>& binsx,const vector<double>& binsy,
	       const vector<TH1F***>& histos);

/// fill the appropriate histogram (histos), based on pdgid, x, y, binsx, and binsy
void fillHisto(int pdgid,double value,double x,double y,
	       const vector<double>& binsx,const vector<double>& binsy,
	       const vector<TH1F***>& histos);



////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
JetResponseAnalyzer::JetResponseAnalyzer(const edm::ParameterSet& iConfig)
  : doHistos_(iConfig.getParameter<bool>("doHistos"))
  , doTree_(iConfig.getParameter<bool>  ("doTree"))
  , doFlavor_(iConfig.getParameter<bool>("doFlavor"))
  , doJetPt_(iConfig.getParameter<bool> ("doJetPt"))
  , doRefPt_(iConfig.getParameter<bool> ("doRefPt"))
  , doBalancing_(false)
  , deltaRMax_(0.0)
  , deltaPhiMin_(3.141)
  , srcRef_(iConfig.getParameter<edm::InputTag>        ("srcRef"))
  , srcRefToJetMap_(iConfig.getParameter<edm::InputTag>("srcRefToJetMap"))
  , getFlavorFromMap_(false)
  , deltaRPartonMax_(0.0)
  , nRefMax_(iConfig.getParameter<unsigned int>("nRefMax"))
{
  if (iConfig.exists("deltaRMax")) {
    doBalancing_=false;
    deltaRMax_=iConfig.getParameter<double>("deltaRMax");
  }
  else if (iConfig.exists("deltaPhiMin")) {
    doBalancing_=true;
    deltaPhiMin_=iConfig.getParameter<double>("deltaPhiMin");
  }
  else
    throw cms::Exception("MissingParameter")<<"Set *either* deltaRMax (matching)"
					    <<" *or* deltaPhiMin (balancing)";
  
  
  if (doFlavor_&&iConfig.exists("srcRefToPartonMap")) {
    srcRefToPartonMap_=iConfig.getParameter<edm::InputTag>("srcRefToPartonMap");
    deltaRPartonMax_  =iConfig.getParameter<double>       ("deltaRPartonMax");
    getFlavorFromMap_=true;
  }
  
  
  if (doHistos_) {

    nBinsPt_     =iConfig.getParameter<unsigned int>    ("nBinsPt");
    nBinsEta_    =iConfig.getParameter<unsigned int>    ("nBinsEta");
    nBinsPhi_    =iConfig.getParameter<unsigned int>    ("nBinsPhi");
    etaBarrelMin_=iConfig.getParameter<double>          ("etaBarrelMin");
    etaBarrelMax_=iConfig.getParameter<double>          ("etaBarrelMax");
    binsPt_      =iConfig.getParameter< vector<double> >("binsPt");
    binsEta_     =iConfig.getParameter< vector<double> >("binsEta");
    binsPhi_     =iConfig.getParameter< vector<double> >("binsPhi");
    
    nBinsRelRsp_ =iConfig.getParameter<unsigned int>    ("nBinsRelRsp");
    nBinsAbsRsp_ =iConfig.getParameter<unsigned int>    ("nBinsAbsRsp");
    
    doRelRsp_=(nBinsRelRsp_>0);
    if (doRelRsp_) {
      relRspMin_=iConfig.getParameter<double>("relRspMin");
      relRspMax_=iConfig.getParameter<double>("relRspMax");
    }
    doAbsRsp_=(nBinsAbsRsp_>0);
    if (doAbsRsp_) {
      absRspMin_=iConfig.getParameter<double>("absRspMin");
      absRspMax_=iConfig.getParameter<double>("absRspMax");
    }
  }
}


//______________________________________________________________________________
JetResponseAnalyzer::~JetResponseAnalyzer()
{

}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void JetResponseAnalyzer::beginJob(const edm::EventSetup& iSetup)
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration,
				"TFileService missing from configuration!");
  
  if (doHistos_) {

    // define flavors
    vector<string> flavor;
    flavor.push_back("");
    if (doFlavor_) {
      flavor.push_back("uds_");
      flavor.push_back("c_");
      flavor.push_back("b_");
      flavor.push_back("g_");
    }

    // book pT histograms
    if (binsPt_.size()>=2) {
      for (unsigned int iPt=0;iPt<binsPt_.size()-1;++iPt) {

	string hname; double ptMin=binsPt_[iPt]; double ptMax=binsPt_[iPt+1];
	
	if (doJetPt_) {
	  jetPtVsJetPt_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"JetPt_"+getSuffix("JetPt",iPt,binsPt_);
	    jetPtVsJetPt_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
						      ";p_{T} [GeV]",
						      nBinsPt_,ptMin,ptMax);
	  }
	}
	
	if (doRefPt_) {
	  refPtVsRefPt_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"RefPt_"+getSuffix("RefPt",iPt,binsPt_);
	    refPtVsRefPt_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
						      ";p_{T}^{ref} [GeV]",
						      nBinsPt_,ptMin,ptMax);
	  }
	}
	
	if (doRefPt_) {
	  jetPtVsRefPt_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"JetPt_"+getSuffix("RefPt",iPt,binsPt_);
	    jetPtVsRefPt_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
						      ";p_{T} [GeV]",
						      2*nBinsPt_,
						      (ptMin>100.)*0.25*ptMin,
						      1.25*ptMax);
	  }
	}
	
	if (doRefPt_) {
	  refPtVsRefPtBarrel_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"RefPt_Barrel_"+getSuffix("RefPt",iPt,binsPt_);
	    refPtVsRefPtBarrel_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
							    ";p_{T}^{ref} [GeV]",
							    nBinsPt_,ptMin,ptMax);
	  }
	}
	
	if (doRefPt_) {
	  jetPtVsRefPtBarrel_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"JetPt_Barrel_"+getSuffix("RefPt",iPt,binsPt_);
	    jetPtVsRefPtBarrel_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
							    ";p_{T} [GeV]",
							    2*nBinsPt_,
							    (ptMin>100)*0.25*ptMin,
							    1.25*ptMax);
	  }
	}
	
	if (doRelRsp_&&doJetPt_) {
	  relRspVsJetPt_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"RelRsp_"+getSuffix("JetPt",iPt,binsPt_);
	    relRspVsJetPt_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
						       ";p_{T}/p_{T}^{ref}",
						       nBinsRelRsp_,
						       relRspMin_,relRspMax_);
	  }
	}
	
	if (doRelRsp_&&doRefPt_) {
	  relRspVsRefPt_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"RelRsp_"+getSuffix("RefPt",iPt,binsPt_);
	    relRspVsRefPt_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
						       ";p_{T}/p_{T}^{ref}",
						       nBinsRelRsp_,
						       relRspMin_,relRspMax_);
	  }
	}
	
	if (doRelRsp_&&doRefPt_) {
	  relRspVsRefPtBarrel_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"RelRsp_Barrel_"+getSuffix("RefPt",iPt,binsPt_);
	    relRspVsRefPtBarrel_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
							     ";p_{T}/p_{T}^{ref}",
							     nBinsRelRsp_,
							     relRspMin_,
							     relRspMax_);
	  }
	}
	
	if (doAbsRsp_&&doJetPt_) {
	  absRspVsJetPt_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"AbsRsp_"+getSuffix("JetPt",iPt,binsPt_);
	    absRspVsJetPt_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
						       ";p_{T}-p_{T}^{ref} [GeV]",
						       nBinsAbsRsp_,
						       absRspMin_,absRspMax_);
	  }
	}
	
	if (doAbsRsp_&&doRefPt_) {
	  absRspVsRefPt_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"AbsRsp_"+getSuffix("RefPt",iPt,binsPt_);
	    absRspVsRefPt_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
						       ";p_{T}-p_{T}^{ref} [GeV]",
						       nBinsAbsRsp_,
						       absRspMin_,absRspMax_);
	  }
	}
	
	if (doAbsRsp_&&doRefPt_) {
	  absRspVsRefPtBarrel_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"AbsRsp_Barrel_"+getSuffix("RefPt",iPt,binsPt_);
	    absRspVsRefPtBarrel_.back()[iFlv]=
	      fs->make<TH1F>(hname.c_str(),";p_{T}-p_{T}^{ref} [GeV]",
			     nBinsAbsRsp_,absRspMin_,absRspMax_);
	  }
	}
	
      }
    }
    
    // book eta histograms
    if (binsEta_.size()>=2) {
      for (unsigned int iEta=0;iEta<binsEta_.size()-1;++iEta) {
	
	string hname; double etaMin=binsEta_[iEta]; double etaMax=binsEta_[iEta+1];

	if (1) {
	  jetEtaVsJetEta_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"JetEta_"+getSuffix("JetEta",iEta,binsEta_);
	    jetEtaVsJetEta_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),";#eta",
							nBinsEta_,etaMin,etaMax);
	  }
	}
	
	if (doRelRsp_) {
	  relRspVsJetEta_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"RelRsp_"+getSuffix("JetEta",iEta,binsEta_);
	    relRspVsJetEta_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
							";p_{T}/p_{T}^{ref}",
							nBinsRelRsp_,
							relRspMin_,relRspMax_);
	  }
	}
	
	if (doAbsRsp_) {
	  absRspVsJetEta_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"AbsRsp_"+getSuffix("JetEta",iEta,binsEta_);
	    absRspVsJetEta_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
							";p_{T}-p_{T}^{ref} [GeV]",
							nBinsAbsRsp_,
							absRspMin_,absRspMax_);
	  }
	}
	
      }
    }
    
    // book phi histograms
    if (binsPhi_.size()>=2) {
      for (unsigned int iPhi=0;iPhi<binsPhi_.size()-1;++iPhi) {
	
	string hname; double phiMin=binsPhi_[iPhi]; double phiMax=binsPhi_[iPhi+1];

	if (1) {
	  jetPhiVsJetPhi_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"JetPhi_"+getSuffix("JetPhi",iPhi,binsPhi_);
	    jetPhiVsJetPhi_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),";#phi",
							nBinsPhi_,phiMin,phiMax);
	  }
	}
	
	if (doRelRsp_) {
	  relRspVsJetPhi_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"RelRsp_"+getSuffix("JetPhi",iPhi,binsPhi_);
	    relRspVsJetPhi_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
							";p_{T}/p_{T}^{ref}",
							nBinsRelRsp_,
							relRspMin_,relRspMax_);
	  }
	}
	
	if (doAbsRsp_) {
	  absRspVsJetPhi_.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {
	    hname=flavor[iFlv]+"AbsRsp_"+getSuffix("JetPhi",iPhi,binsPhi_);
	    absRspVsJetPhi_.back()[iFlv]=fs->make<TH1F>(hname.c_str(),
							";p_{T}-p_{T}^{ref} [GeV]",
							nBinsAbsRsp_,
							absRspMin_,absRspMax_);
	  }
	}
	
      }
    }
    
    // book eta/pT histograms
    if (binsPt_.size()>=2&&binsEta_.size()>=2) {
      for (unsigned int iEta=0;iEta<binsEta_.size()-1;++iEta) {

	TH1F*** jetPtJetPt(0);
	TH1F*** refPtRefPt(0);
	TH1F*** jetPtRefPt(0);
	TH1F*** relRspJetPt(0);
	TH1F*** relRspRefPt(0);
	TH1F*** absRspJetPt(0);
	TH1F*** absRspRefPt(0);
	
	if (doJetPt_) {
	  jetPtJetPt=new TH1F**[binsPt_.size()-1];
	  for (unsigned int iPt=0;iPt<binsPt_.size()-1;iPt++)
	    jetPtJetPt[iPt]=new TH1F*[flavor.size()];
	}
	
	if (doRefPt_) {
	  refPtRefPt =new TH1F**[binsPt_.size()];
	  for (unsigned int iPt=0;iPt<binsPt_.size()-1;iPt++)
	    refPtRefPt[iPt]=new TH1F*[flavor.size()];
	}
	
	if (doRefPt_) {
	  jetPtRefPt =new TH1F**[binsPt_.size()];
	  for (unsigned int iPt=0;iPt<binsPt_.size()-1;iPt++)
	    jetPtRefPt[iPt]=new TH1F*[flavor.size()];
	}
	
	if (doRelRsp_&&doJetPt_) {
	  relRspJetPt=new TH1F**[binsPt_.size()-1];
	  for (unsigned int iPt=0;iPt<binsPt_.size()-1;iPt++)
	    relRspJetPt[iPt]=new TH1F*[flavor.size()];
	}
	
	if (doRelRsp_&&doRefPt_) {
	  relRspRefPt=new TH1F**[binsPt_.size()-1];
	  for (unsigned int iPt=0;iPt<binsPt_.size()-1;iPt++)
	    relRspRefPt[iPt]=new TH1F*[flavor.size()];
	}
	
	if (doAbsRsp_&&doJetPt_) {
	  absRspJetPt=new TH1F**[binsPt_.size()-1];
	  for (unsigned int iPt=0;iPt<binsPt_.size()-1;iPt++)
	    absRspJetPt[iPt]=new TH1F*[flavor.size()];
	}
	
	if (doAbsRsp_&&doRefPt_) {
	  absRspRefPt=new TH1F**[binsPt_.size()-1];
	  for (unsigned int iPt=0;iPt<binsPt_.size()-1;iPt++)
	    absRspRefPt[iPt]=new TH1F*[flavor.size()];
	}
	
	string jetEtaSuffix=getSuffix("JetEta",iEta,binsEta_);
	
	for (unsigned int iPt=0;iPt<binsPt_.size()-1;iPt++) {	

	  string hname; double ptMin=binsPt_[iPt]; double ptMax=binsPt_[iPt+1];
	  
	  string jetPtSuffix=getSuffix("JetPt",iPt,binsPt_);
	  string refPtSuffix=getSuffix("RefPt",iPt,binsPt_);
	  
	  for (unsigned int iFlv=0;iFlv<flavor.size();iFlv++) {

	    if (doJetPt_) {
	      hname=flavor[iFlv]+"JetPt_"+jetEtaSuffix+"_"+jetPtSuffix;
	      jetPtJetPt[iPt][iFlv]=fs->make<TH1F>(hname.c_str(),";p_{T}",
						   nBinsPt_,ptMin,ptMax);
	    }
	    
	    if (doRefPt_) {
	      hname=flavor[iFlv]+"RefPt_"+jetEtaSuffix+"_"+refPtSuffix;
	      refPtRefPt[iPt][iFlv]=fs->make<TH1F>(hname.c_str(),";p_{T}^{ref}",
						   nBinsPt_,ptMin,ptMax);
	    }
	    
	    if (doRefPt_) {
	      hname=flavor[iFlv]+"JetPt_"+jetEtaSuffix+"_"+refPtSuffix;
	      jetPtRefPt[iPt][iFlv]=fs->make<TH1F>(hname.c_str(),";p_{T}",
						   2*nBinsPt_,
						   (ptMin>100.)*0.25*ptMin,
						   1.25*ptMax);
	    }
	    
	    if (doRelRsp_&&doJetPt_) {
	      hname=flavor[iFlv]+"RelRsp_"+jetEtaSuffix+"_"+jetPtSuffix;
	      relRspJetPt[iPt][iFlv]=fs->make<TH1F>(hname.c_str(),
						    ";p_{T}/p_{T}^{ref}",
						    nBinsRelRsp_,
						    relRspMin_,relRspMax_);
	    }
	    
	    if (doRelRsp_&&doRefPt_) {
	      hname=flavor[iFlv]+"RelRsp_"+jetEtaSuffix+"_"+refPtSuffix;
	      relRspRefPt[iPt][iFlv]=fs->make<TH1F>(hname.c_str(),
						    ";p_{T}/p_{T}^{ref}",
						    nBinsRelRsp_,
						    relRspMin_,relRspMax_);
	    }
	    
	    if (doAbsRsp_&&doJetPt_) {
	      hname=flavor[iFlv]+"AbsRsp_"+jetEtaSuffix+"_"+jetPtSuffix;
	      absRspJetPt[iPt][iFlv]=fs->make<TH1F>(hname.c_str(),
						    ";p_{T}-p_{T}^{ref} [GeV]",
						    nBinsAbsRsp_,
						    absRspMin_,absRspMax_);
	    }
	    
	    if (doAbsRsp_&&doRefPt_) {
	      hname=flavor[iFlv]+"AbsRsp_"+jetEtaSuffix+"_"+refPtSuffix;
	      absRspRefPt[iPt][iFlv]=fs->make<TH1F>(hname.c_str(),
						    ";p_{T}-p_{T}^{ref} [GeV]",
						    nBinsAbsRsp_,
						    absRspMin_,absRspMax_);
	    }
	    
	  }

	}
	
	if (doJetPt_)            jetPtVsJetEtaJetPt_ .push_back(jetPtJetPt);
	if (doRefPt_)            refPtVsJetEtaRefPt_ .push_back(refPtRefPt);
	if (doRefPt_)            jetPtVsJetEtaRefPt_ .push_back(jetPtRefPt);
	if (doRelRsp_&&doJetPt_) relRspVsJetEtaJetPt_.push_back(relRspJetPt);
	if (doRelRsp_&&doRefPt_) relRspVsJetEtaRefPt_.push_back(relRspRefPt);
	if (doAbsRsp_&&doJetPt_) absRspVsJetEtaJetPt_.push_back(absRspJetPt);
	if (doAbsRsp_&&doRefPt_) absRspVsJetEtaRefPt_.push_back(absRspRefPt);
      }
    }
  } // doHistos_
  
  if (doTree_) {
    tree_=fs->make<TTree>("t","t");
    tree_->Branch("nref",  &nref_,  "nref/b");
    if (doFlavor_) tree_->Branch("refpdgid",refpdgid_,"refpdgid[nref]/I");
    tree_->Branch("refpt",  refpt_, "refpt[nref]/F");
    tree_->Branch("refeta", refeta_,"refeta[nref]/F");
    tree_->Branch("refphi", refphi_,"refphi[nref]/F");
    tree_->Branch("jtpt",   jtpt_,  "jtpt[nref]/F");
    tree_->Branch("jteta",  jteta_, "jteta[nref]/F");
    tree_->Branch("jtphi",  jtphi_, "jtphi[nref]/F");
    if (doBalancing_) tree_->Branch("refdphijt",refdphijt_,"refdphijt[nref]/F");
    else    	      tree_->Branch("refdrjt",  refdrjt_,  "refdrjt[nref]/F");
  }
  
}


//______________________________________________________________________________
void JetResponseAnalyzer::analyze(const edm::Event&      iEvent,
				  const edm::EventSetup& iSetup)
{
  nref_=0;
  
  edm::Handle<reco::CandidateView>    refs;
  edm::Handle<reco::CandViewMatchMap> refToJetMap;
  edm::Handle<reco::CandViewMatchMap> refToPartonMap;
  
  iEvent.getByLabel(srcRef_,        refs);
  iEvent.getByLabel(srcRefToJetMap_,refToJetMap);
  if (getFlavorFromMap_) iEvent.getByLabel(srcRefToPartonMap_,refToPartonMap);
  
  if (doBalancing_&&refToJetMap->size()!=1) return;
  
  unsigned int nRef=(nRefMax_==0) ? refs->size() : std::min(nRefMax_,refs->size());
  
  for (unsigned int iRef=0;iRef<nRef;iRef++) {

    reco::CandidateBaseRef ref=refs->refAt(iRef);
    reco::CandViewMatchMap::const_iterator itMatch=refToJetMap->find(ref);
    if (itMatch==refToJetMap->end()) continue;
    reco::CandidateBaseRef jet=itMatch->val;
    
    refdrjt_[nref_]  =reco::deltaR(*jet,*ref);
    refdphijt_[nref_]=reco::deltaPhi(*jet,*ref);
    
    if ((!doBalancing_&&refdrjt_[nref_]>deltaRMax_)||
	(doBalancing_&&std::abs(refdphijt_[nref_])<deltaPhiMin_)) continue;
    
    refpdgid_[nref_]=0;
    if (getFlavorFromMap_) {
      itMatch=refToPartonMap->find(ref);
      if (itMatch!=refToPartonMap->end()) {
	double refdrparton=reco::deltaR(*itMatch->key,*itMatch->val);
	if (refdrparton<deltaRPartonMax_)
	  refpdgid_[nref_]=itMatch->val->pdgId();
      }
    }
    else {
      refpdgid_[nref_]=ref->pdgId();
    }
    
    refpt_[nref_]   =ref->pt();
    refeta_[nref_]  =ref->eta();
    refphi_[nref_]  =ref->phi();
    jtpt_[nref_]    =jet->pt();
    jteta_[nref_]   =jet->eta();
    jtphi_[nref_]   =jet->phi();
    nref_++;
    
    double absRsp=jet->pt()-ref->pt();
    double relRsp=jet->pt()/ref->pt();
    
    if (doHistos_) {
      if (jet->eta()>=etaBarrelMin_&&jet->eta()<=etaBarrelMax_) {
	if (doRefPt_) {
	  fillHisto(ref->pt(),ref->pt(),binsPt_,refPtVsRefPtBarrel_);
	  fillHisto(jet->pt(),ref->pt(),binsPt_,jetPtVsRefPtBarrel_);
	  if (doFlavor_) {
	    fillHisto(refpdgid_[nref_],
		      ref->pt(),ref->pt(),binsPt_,refPtVsRefPtBarrel_);
	    fillHisto(refpdgid_[nref_],
		      jet->pt(),ref->pt(),binsPt_,jetPtVsRefPtBarrel_);
	  }
	}
	if (doRelRsp_&&doRefPt_) {
	  fillHisto(relRsp,ref->pt(),binsPt_,relRspVsRefPtBarrel_);
	  if (doFlavor_)
	    fillHisto(refpdgid_[nref_],
		      relRsp,ref->pt(),binsPt_,relRspVsRefPtBarrel_);
	}
	if (doAbsRsp_&&doRefPt_) {
	  fillHisto(absRsp,ref->pt(),binsPt_,absRspVsRefPtBarrel_);
	  if (doFlavor_)
	    fillHisto(refpdgid_[nref_],
		      absRsp,ref->pt(),binsPt_,absRspVsRefPtBarrel_);
	}
      }
      
      if (doJetPt_) {
	fillHisto(jet->pt(), jet->pt(), binsPt_, jetPtVsJetPt_);
	if (doFlavor_) fillHisto(refpdgid_[nref_],
				 jet->pt(), jet->pt(), binsPt_, jetPtVsJetPt_);
      }
      if (doRefPt_) {
	fillHisto(ref->pt(), ref->pt(), binsPt_, refPtVsRefPt_);
	fillHisto(jet->pt(), ref->pt(), binsPt_, jetPtVsRefPt_);
	if (doFlavor_) {
	  fillHisto(refpdgid_[nref_],ref->pt(), ref->pt(), binsPt_, refPtVsRefPt_);
	  fillHisto(refpdgid_[nref_],jet->pt(), ref->pt(), binsPt_, jetPtVsRefPt_);
	}
      }
      
      if (binsEta_.size()>=2) {
	fillHisto(jet->eta(),jet->eta(),binsEta_,jetEtaVsJetEta_);
	if (doFlavor_) fillHisto(refpdgid_[nref_],
				 jet->eta(),jet->eta(),binsEta_,jetEtaVsJetEta_);
      }
      
      if (binsPhi_.size()>=2) {
	fillHisto(jet->phi(),jet->phi(),binsPhi_,jetPhiVsJetPhi_);
	if (doFlavor_) fillHisto(refpdgid_[nref_],
				 jet->phi(),jet->phi(),binsPhi_,jetPhiVsJetPhi_);
      }
      
      if (doJetPt_&&binsEta_.size()>=2) {
	fillHisto(jet->pt(), jet->eta(),jet->pt(),
		  binsEta_,binsPt_,jetPtVsJetEtaJetPt_);
	if (doFlavor_) fillHisto(refpdgid_[nref_],
				 jet->pt(), jet->eta(),jet->pt(),
				 binsEta_,binsPt_,jetPtVsJetEtaJetPt_);
      }
      if (doRefPt_&&binsEta_.size()>=2) {
	fillHisto(ref->pt(), jet->eta(),ref->pt(),
		  binsEta_,binsPt_,refPtVsJetEtaRefPt_);
	fillHisto(jet->pt(), jet->eta(),ref->pt(),
		  binsEta_,binsPt_,jetPtVsJetEtaRefPt_);
	if (doFlavor_) {
	  fillHisto(refpdgid_[nref_],ref->pt(), jet->eta(),ref->pt(),
		    binsEta_,binsPt_,refPtVsJetEtaRefPt_);
	  fillHisto(refpdgid_[nref_],jet->pt(), jet->eta(),ref->pt(),
		    binsEta_,binsPt_,jetPtVsJetEtaRefPt_);
	}
      }
      
      if (doRelRsp_) {
	if (doJetPt_) {
	  fillHisto(relRsp,jet->pt(), binsPt_, relRspVsJetPt_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],
				   relRsp,jet->pt(), binsPt_, relRspVsJetPt_);
	}
	if (doRefPt_) {
	  fillHisto(relRsp,ref->pt(), binsPt_, relRspVsRefPt_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],
				   relRsp,ref->pt(), binsPt_, relRspVsRefPt_);
	}
	
	if (binsEta_.size()>=2) {
	  fillHisto(relRsp,jet->eta(),binsEta_,relRspVsJetEta_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],
				   relRsp,jet->eta(),binsEta_,relRspVsJetEta_);
	}
	
	if (binsPhi_.size()>=2) {
	  fillHisto(relRsp,jet->phi(),binsPhi_,relRspVsJetPhi_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],
				   relRsp,jet->phi(),binsPhi_,relRspVsJetPhi_);
	}
	
	if (doJetPt_&&binsEta_.size()>=2) {
	  fillHisto(relRsp,jet->eta(),jet->pt(),
		    binsEta_,binsPt_,relRspVsJetEtaJetPt_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],relRsp,jet->eta(),jet->pt(),
				   binsEta_,binsPt_,relRspVsJetEtaJetPt_);
	}
	if (doRefPt_&&binsEta_.size()>=2) {
	  fillHisto(relRsp,jet->eta(),ref->pt(),
		    binsEta_,binsPt_,relRspVsJetEtaRefPt_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],relRsp,jet->eta(),ref->pt(),
				   binsEta_,binsPt_,relRspVsJetEtaRefPt_);
	}
      }
      
      if (doAbsRsp_) {
	if (doJetPt_) {
	  fillHisto(absRsp,jet->pt(), binsPt_, absRspVsJetPt_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],
				   absRsp,jet->pt(), binsPt_, absRspVsJetPt_);
	}
	if (doRefPt_) {
	  fillHisto(absRsp,ref->pt(), binsPt_, absRspVsRefPt_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],
				   absRsp,ref->pt(), binsPt_, absRspVsRefPt_);
	}
	
	if (binsEta_.size()>=2) {
	  fillHisto(absRsp,jet->eta(),binsEta_,absRspVsJetEta_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],
				   absRsp,jet->eta(),binsEta_,absRspVsJetEta_);
	}
	
	if (binsPhi_.size()>=2) {
	  fillHisto(absRsp,jet->phi(),binsPhi_,absRspVsJetPhi_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],
				   absRsp,jet->phi(),binsPhi_,absRspVsJetPhi_);
	}
	
	if (doJetPt_&&binsEta_.size()>=2) {
	  fillHisto(absRsp,jet->eta(),jet->pt(),
		    binsEta_,binsPt_,absRspVsJetEtaJetPt_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],absRsp,jet->eta(),jet->pt(),
				   binsEta_,binsPt_,absRspVsJetEtaJetPt_);
	}
	if (doRefPt_&&binsEta_.size()>=2) {
	  fillHisto(absRsp,jet->eta(),ref->pt(),
		    binsEta_,binsPt_,absRspVsJetEtaRefPt_);
	  if (doFlavor_) fillHisto(refpdgid_[nref_],absRsp,jet->eta(),ref->pt(),
				   binsEta_,binsPt_,absRspVsJetEtaRefPt_);
	}
      }
      
    }
  }
  
  if (doTree_) tree_->Fill();
}


//______________________________________________________________________________
string getSuffix(const string& varname,int ibin,const vector<double>& bins)
{
  stringstream ss; ss<<varname<<bins[ibin]<<"to"<<bins[ibin+1];
  return ss.str();
}

//______________________________________________________________________________
int getIndex(double x,const vector<double>& binsx)
{
  for (unsigned int ix=0;ix<binsx.size()-1;ix++)
    if (x>=binsx[ix]&&x<binsx[ix+1]) return ix;
  return -1;
}

//______________________________________________________________________________
void fillHisto(double value,double x,
	       const vector<double>& binsx,const vector<TH1F**>& histos)
{
  int ix=getIndex(x,binsx);
  if (ix>=0) histos[ix][0]->Fill(value);
}

//______________________________________________________________________________
void fillHisto(int pdgid,double value,double x,
	       const vector<double>& binsx,const vector<TH1F**>& histos)
{
  int abspdgid=std::abs(pdgid);
  int iflv(-1);
  if (abspdgid>=1&&abspdgid<=3) iflv=1;
  else if (abspdgid== 4)        iflv=2;
  else if (abspdgid== 5)        iflv=3;
  else if (abspdgid==21)        iflv=4;
  else return;

  int ix=getIndex(x,binsx);
  if (ix>=0) histos[ix][iflv]->Fill(value);
}


//______________________________________________________________________________
void fillHisto(double value,double x,double y,
	       const vector<double>& binsx,const vector<double>& binsy,
	       const vector<TH1F***>& histos)
{
  int ix=getIndex(x,binsx);
  int iy=getIndex(y,binsy);
  if (ix>=0&&iy>=0) histos[ix][iy][0]->Fill(value);
}

//______________________________________________________________________________
void fillHisto(int pdgid,double value,double x,double y,
	       const vector<double>& binsx,const vector<double>& binsy,
	       const vector<TH1F***>& histos)
{
  int abspdgid=std::abs(pdgid);
  int iflv(-1);
  if (abspdgid>=1&&abspdgid<=3) iflv=1;
  else if (abspdgid== 4)        iflv=2;
  else if (abspdgid== 5)        iflv=3;
  else if (abspdgid==21)        iflv=4;
  else return;
  
  int ix=getIndex(x,binsx);
  int iy=getIndex(y,binsy);
  if (ix>=0&&iy>=0) histos[ix][iy][iflv]->Fill(value);
}



////////////////////////////////////////////////////////////////////////////////
// define JetEfficiencyAnalyzer as a plugin
////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(JetResponseAnalyzer);
