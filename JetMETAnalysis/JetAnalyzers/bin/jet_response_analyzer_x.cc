////////////////////////////////////////////////////////////////////////////////
//
// jet_response_analyzer_x
// -----------------------
//
// DESCRIPTION: jet_response_analyzer_x takes a file written with
// the CMSSW fwk analyzer JetResponseAnalyzer, and turns the trees
// in each directory (for each algorithm) into histograms. The exact
// same result can be achieved with the fwk analyzer directly! But if
// you decide to operate the fwk analyzer with doTree=true,
// doHistos=false, you can fit the output e.g. on a local disk and
// decide on the cuts and binning now. A lot more flexibility, for a
// lot more disk space ...
//
//            07/23/2008 Kostas Kousouris                    <kkousour@fnal.gov>
//                       Roger Wolf                         <roger.wolf@cern.ch>
//                       Philipp Schieferdecker <philipp.schieferdecker@cern.ch>
////////////////////////////////////////////////////////////////////////////////


#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TKey.h>

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>


using namespace std;


////////////////////////////////////////////////////////////////////////////////
// define local functions
////////////////////////////////////////////////////////////////////////////////

/// get the suffix for the histogram name, e.g. JetPt100to150
string get_suffix(const string& varname,int ibin,const vector<float>& bins);

/// get the index of the histogram corresponding to x
int get_index(float x,const vector<float>& binsx);

/// fill the appropriate histogram (histos), based on x and binsx
void fill_histo(float value,float x,
		const vector<float>& binsx,const vector<TH1F**>& histos);

/// fill the appropriate histogram (histos), based on pdgid, x and binsx
void fill_histo(int pdgid,float value,float x,
		const vector<float>& binsx,const vector<TH1F**>& histos);

/// fill the appropriate histogram (histos), based on x, y, binsx, and binsy
void fill_histo(float value,float x,float y,
		const vector<float>& binsx,const vector<float>& binsy,
		const vector<TH1F***>& histos);

/// fill the appropriate histogram (histos), based on pdgid, x, y, binsx, and binsy
void fill_histo(int pdgid,float value,float x,float y,
		const vector<float>& binsx,const vector<float>& binsy,
		const vector<TH1F***>& histos);

/// check if a vector of strings contains a certain element
bool contains(const vector<string>& collection,const string& element);


//______________________________________________________________________________
int main(int argc,char**argv)
{
  //
  // evaluate command-line / configuration file options
  //
  CommandLine cl;
  if (!cl.parse(argc,argv)) return 0;

  string         input        = cl.getValue<string> ("input");
  vector<float>  binspt       = cl.getVector<float> ("binspt",        "");
  vector<float>  binseta      = cl.getVector<float> ("binseta",       "");
  vector<float>  binsphi      = cl.getVector<float> ("binsphi",       "");
  string         treename     = cl.getValue<string> ("treename",     "t");
  string         output       = cl.getValue<string> ("output","jra.root");
  int            nbinspt      = cl.getValue<int>    ("nbinspt",       50);
  int            nbinseta     = cl.getValue<int>    ("nbinseta",      25);
  int            nbinsphi     = cl.getValue<int>    ("nbinsphi",      25);
  float          etabarrelmin = cl.getValue<float>  ("etabarrelmin",-1.0);
  float          etabarrelmax = cl.getValue<float>  ("etabarrelmax",+1.0);
  bool           dobalance    = cl.getValue<bool>   ("dobalance",  false);
  bool           doflavor     = cl.getValue<bool>   ("doflavor",   false);
  float          drmax        = cl.getValue<float>  ("drmax",        0.3);
  float          dphimin      = cl.getValue<float>  ("dphimin",      2.7);
  bool           dojetpt      = cl.getValue<bool>   ("dojetpt",     true);
  bool           dorefpt      = cl.getValue<bool>   ("dorefpt",     true);
  int            nbinsrelrsp  = cl.getValue<int>    ("nbinsrelrsp",  100);
  float          relrspmin    = cl.getValue<float>  ("relrspmin",    0.0);
  float          relrspmax    = cl.getValue<float>  ("relrspmax",    2.0);
  int            nbinsabsrsp  = cl.getValue<int>    ("nbinsabsrsp",  600);
  float          absrspmin    = cl.getValue<float>  ("absrspmin",-1000.0);
  float          absrspmax    = cl.getValue<float>  ("absrspmax",  200.0);
  vector<string> algs         = cl.getVector<string>("algs",          "");
  
  if (!cl.check()) return 0;
  cl.print();

  bool dorelrsp=(nbinsrelrsp>0);
  bool doabsrsp=(nbinsabsrsp>0);

  
  //
  // open input/output files and loop over input directories/trees (=algorithms!)
  //
  TFile* ifile = new TFile(input.c_str(),"READ");
  if (!ifile->IsOpen()) {  cout<<"Can't open "<<input<<endl; return 0; }
  
  TFile* ofile = new TFile(output.c_str(),"RECREATE");
  if (!ofile->IsOpen()) { cout<<"Can't create "<<output<<endl; return 0; }

  TIter next(ifile->GetListOfKeys());
  TKey* key(0);
  while ((key=(TKey*)next())) {
    if (strcmp(key->GetClassName(),"TDirectoryFile")!=0) continue;
    
    TDirectoryFile* idir = (TDirectoryFile*)key->ReadObj();
    string alg(idir->GetName());
    if (algs.size()>0&&!contains(algs,alg)) continue;
    
    cout<<alg<<" ... "<<flush;

    TTree* tree = (TTree*)idir->Get("t");
    if (0==tree) { cout<<"no tree found."<<endl; continue; }
    
    
    //
    // setup the tree for reading
    //
    unsigned char nref;
    int   refpdgid[100];
    float refpt[100];
    float refeta[100];
    float refphi[100];
    float jtpt[100];
    float jteta[100];
    float jtphi[100];
    float refdrjt[100];
    float refdphijt[100];
    
    tree->SetBranchAddress("nref",   &nref);
    if (doflavor) tree->SetBranchAddress("refpdgid",refpdgid);
    tree->SetBranchAddress("refpt",   refpt);
    tree->SetBranchAddress("refeta",  refeta);
    tree->SetBranchAddress("refphi",  refphi);
    tree->SetBranchAddress("jtpt",    jtpt);
    tree->SetBranchAddress("jteta",   jteta);
    tree->SetBranchAddress("jtphi",   jtphi);
    
    if (dobalance) {
      if (0==tree->GetBranch("refdphijt")) {
	cout<<"dobalance, but no branch 'refdphijt' in tree, skip!"<<endl;
	continue;
      }
      else tree->SetBranchAddress("refdphijt",refdphijt);
    }
    else {
      if (0==tree->GetBranch("refdrjt")) {
	cout<<"!dobalance, but no branch 'refdrjt' in tree, skip!"<<endl;
	continue;
      }
      else tree->SetBranchAddress("refdrjt",refdrjt);
    }
    
    
    //
    // create directory in output file and book histograms
    //
    TDirectoryFile* odir = (TDirectoryFile*)ofile->mkdir(alg.c_str());
    if (0==odir) { cout<<"failed to create directory."<<endl; continue; }
    odir->cd();
    
    // declare histograms
    vector<TH1F**>  jetPtVsJetPt;
    vector<TH1F**>  refPtVsRefPt;
    vector<TH1F**>  refPtVsRefPtBarrel;
    vector<TH1F**>  jetEtaVsJetEta;
    vector<TH1F**>  jetPhiVsJetPhi;
    vector<TH1F***> jetPtVsJetEtaJetPt;
    vector<TH1F***> refPtVsJetEtaRefPt;
    
    vector<TH1F**>  relRspVsJetPt;
    vector<TH1F**>  relRspVsRefPt;
    vector<TH1F**>  relRspVsRefPtBarrel;
    vector<TH1F**>  relRspVsJetEta;
    vector<TH1F**>  relRspVsJetPhi;
    vector<TH1F***> relRspVsJetEtaJetPt;
    vector<TH1F***> relRspVsJetEtaRefPt;
    
    vector<TH1F**>  absRspVsJetPt;
    vector<TH1F**>  absRspVsRefPt;
    vector<TH1F**>  absRspVsRefPtBarrel;
    vector<TH1F**>  absRspVsJetEta;
    vector<TH1F**>  absRspVsJetPhi;
    vector<TH1F***> absRspVsJetEtaJetPt;
    vector<TH1F***> absRspVsJetEtaRefPt;

    // define flavors
    vector<string> flavor;
    flavor.push_back("");
    if (doflavor) {
      flavor.push_back("uds_");
      flavor.push_back("c_");
      flavor.push_back("b_");
      flavor.push_back("g_");
    }
    
    // book pT histograms
    if (binspt.size()>=2) {
      for (unsigned int ipt=0;ipt<binspt.size()-1;++ipt) {
	
	string hname; float ptmin=binspt[ipt]; float ptmax=binspt[ipt+1];
	
	if (dojetpt) {
	  jetPtVsJetPt.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"JetPt_"+get_suffix("JetPt",ipt,binspt);
	    jetPtVsJetPt.back()[iflv]=new TH1F(hname.c_str(),";p_{T} [GeV]",
					       nbinspt,ptmin,ptmax);
	  }
	}
	
	if (dorefpt) {
	  refPtVsRefPt.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"RefPt_"+get_suffix("RefPt",ipt,binspt);
	    refPtVsRefPt.back()[iflv]=new TH1F(hname.c_str(),";p_{T}^{ref} [GeV]",
					       nbinspt,ptmin,ptmax);
	  }
	}
	
	if (dorefpt) {
	  refPtVsRefPtBarrel.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"RefPt_Barrel_"+get_suffix("RefPt",ipt,binspt);
	    refPtVsRefPtBarrel.back()[iflv]=new TH1F(hname.c_str(),
						     ";p_{T}^{ref} [GeV]",
						     nbinspt,ptmin,ptmax);
	  }
	}
	
	if (dorelrsp&&dojetpt) {
	  relRspVsJetPt.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"RelRsp_"+get_suffix("JetPt",ipt,binspt);
	    relRspVsJetPt.back()[iflv]=new TH1F(hname.c_str(),";p_{T}/p_{T}^{ref}",
						nbinsrelrsp,relrspmin,relrspmax);
	  }
	}

	if (dorelrsp&&dorefpt) {
	  relRspVsRefPt.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"RelRsp_"+get_suffix("RefPt",ipt,binspt);
	    relRspVsRefPt.back()[iflv]=new TH1F(hname.c_str(),";p_{T}/p_{T}^{ref}",
						nbinsrelrsp,relrspmin,relrspmax);
	  }
	}

	if (dorelrsp&&dorefpt) {
	  relRspVsRefPtBarrel.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"RelRsp_Barrel_"+get_suffix("RefPt",ipt,binspt);
	    relRspVsRefPtBarrel.back()[iflv]=new TH1F(hname.c_str(),
						      ";p_{T}/p_{T}^{ref}",
						      nbinsrelrsp,
						      relrspmin,relrspmax);
	  }
	}
	
	if (doabsrsp&&dojetpt) {
	  absRspVsJetPt.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"AbsRsp_"+get_suffix("JetPt",ipt,binspt);
	    absRspVsJetPt.back()[iflv]=new TH1F(hname.c_str(),
						";p_{T}-p_{T}^{ref} [GeV]",
						nbinsabsrsp,absrspmin,absrspmax);
	  }
	}
	
	if (doabsrsp&&dorefpt) {
	  absRspVsRefPt.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"AbsRsp_"+get_suffix("RefPt",ipt,binspt);
	    absRspVsRefPt.back()[iflv]=new TH1F(hname.c_str(),
						";p_{T}-p_{T}^{ref} [GeV]",
						nbinsabsrsp,absrspmin,absrspmax);
	  }
	}
	
	if (doabsrsp&&dorefpt) {
	  absRspVsRefPtBarrel.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"AbsRsp_Barrel_"+get_suffix("RefPt",ipt,binspt);
	    absRspVsRefPtBarrel.back()[iflv]=new TH1F(hname.c_str(),
						      ";p_{T}-p_{T}^{ref} [GeV]",
						      nbinsabsrsp,
						      absrspmin,absrspmax);
	  }
	}
	
      }
    }
    
    // book eta histograms
    if (binseta.size()>=2) {
      for (unsigned int ieta=0;ieta<binseta.size()-1;++ieta) {

	string hname; float etamin=binseta[ieta]; float etamax=binseta[ieta+1];

	if (1) {
	  jetEtaVsJetEta.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"JetEta_"+get_suffix("JetEta",ieta,binseta);
	    jetEtaVsJetEta.back()[iflv]=new TH1F(hname.c_str(),";#eta",
						 nbinseta,etamin,etamax);
	  }
	}
	
	if (dorelrsp) {
	  relRspVsJetEta.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"RelRsp_"+get_suffix("JetEta",ieta,binseta);
	    relRspVsJetEta.back()[iflv]=new TH1F(hname.c_str(),";p_{T}/p_{T}^{ref}",
						 nbinsrelrsp,relrspmin,relrspmax);
	  }
	}
	
	if (doabsrsp) {
	  absRspVsJetEta.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"AbsRsp_"+get_suffix("JetEta",ieta,binseta);
	    absRspVsJetEta.back()[iflv]=new TH1F(hname.c_str(),
						 ";p_{T}-p_{T}^{ref} [GeV]",
						 nbinsabsrsp,absrspmin,absrspmax);
	  }
	}

      }
    }
    
    // book phi histograms
    if (binsphi.size()>=2) {
      for (unsigned int iphi=0;iphi<binsphi.size()-1;++iphi) {

	string hname; float phimin=binsphi[iphi]; float phimax=binsphi[iphi+1];

	if (1) {
	  jetPhiVsJetPhi.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"JetPhi_"+get_suffix("JetPhi",iphi,binsphi);
	    jetPhiVsJetPhi.back()[iflv]=new TH1F(hname.c_str(),";#phi",
						 nbinsphi,phimin,phimax);
	  }
	}
	
	if (dorelrsp) {
	  relRspVsJetPhi.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"RelRsp_"+get_suffix("JetPhi",iphi,binsphi);
	    relRspVsJetPhi.back()[iflv]=new TH1F(hname.c_str(),";p_{T}/p_{T}^{ref}",
						 nbinsrelrsp,relrspmin,relrspmax);
	  }
	}
	
	if (doabsrsp) {
	  absRspVsJetPhi.push_back(new TH1F*[flavor.size()]);
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    hname=flavor[iflv]+"AbsRsp_"+get_suffix("JetPhi",iphi,binsphi);
	    absRspVsJetPhi.back()[iflv]=new TH1F(hname.c_str(),
						 ";p_{T}-p_{T}^{ref} [GeV]",
						 nbinsabsrsp,absrspmin,absrspmax);
	  }
	}
	
      }
    }
    
    // book eta/pT histograms
    if (binspt.size()>=2&&binseta.size()>=2) {
      for (unsigned int ieta=0;ieta<binseta.size()-1;++ieta) {
	
	TH1F*** jetPtJetPt(0);
	TH1F*** refPtRefPt(0);
	TH1F*** relRspJetPt(0);
	TH1F*** relRspRefPt(0);
	TH1F*** absRspJetPt(0);
	TH1F*** absRspRefPt(0);
	
	if (dojetpt) {
	  jetPtJetPt=new TH1F**[binspt.size()-1];
	  for (unsigned int ipt=0;ipt<binspt.size()-1;ipt++)
	    jetPtJetPt[ipt]=new TH1F*[flavor.size()];
	}
	
	if (dorefpt) {
	  refPtRefPt =new TH1F**[binspt.size()];
	  for (unsigned int ipt=0;ipt<binspt.size()-1;ipt++)
	    refPtRefPt[ipt]=new TH1F*[flavor.size()];
	}
	
	if (dorelrsp&&dojetpt) {
	  relRspJetPt=new TH1F**[binspt.size()-1];
	  for (unsigned int ipt=0;ipt<binspt.size()-1;ipt++)
	    relRspJetPt[ipt]=new TH1F*[flavor.size()];
	}
	
	if (dorelrsp&&dorefpt) {
	  relRspRefPt=new TH1F**[binspt.size()-1];
	  for (unsigned int ipt=0;ipt<binspt.size()-1;ipt++)
	    relRspRefPt[ipt]=new TH1F*[flavor.size()];
	}
	
	if (doabsrsp&&dojetpt) {
	  absRspJetPt=new TH1F**[binspt.size()-1];
	  for (unsigned int ipt=0;ipt<binspt.size()-1;ipt++)
	    absRspJetPt[ipt]=new TH1F*[flavor.size()];
	}
	
	if (doabsrsp&&dorefpt) {
	  absRspRefPt=new TH1F**[binspt.size()-1];
	  for (unsigned int ipt=0;ipt<binspt.size()-1;ipt++)
	    absRspRefPt[ipt]=new TH1F*[flavor.size()];
	}
	
	string jetEtaSuffix=get_suffix("JetEta",ieta,binseta);
	
	for (unsigned int ipt=0;ipt<binspt.size()-1;ipt++) {

	  string hname; float ptmin=binspt[ipt]; float ptmax=binspt[ipt+1];

	  string jetPtSuffix=get_suffix("JetPt",ipt,binspt);
	  string refPtSuffix=get_suffix("RefPt",ipt,binspt);
	  
	  for (unsigned int iflv=0;iflv<flavor.size();iflv++) {
	    
	    if (dojetpt) {
	      hname=flavor[iflv]+"JetPt_"+jetEtaSuffix+"_"+jetPtSuffix;
	      jetPtJetPt[ipt][iflv]=new TH1F(hname.c_str(),";p_{T}",
					     nbinspt,ptmin,ptmax);
	    }
	    
	    if (dorefpt) {
	      hname=flavor[iflv]+"RefPt_"+jetEtaSuffix+"_"+refPtSuffix;
	      refPtRefPt[ipt][iflv]=new TH1F(hname.c_str(),";p_{T}^{ref}",
					     nbinspt,ptmin,ptmax);
	    }
	    
	    if (dorelrsp&&dojetpt) {
	      hname=flavor[iflv]+"RelRsp_"+jetEtaSuffix+"_"+jetPtSuffix;
	      relRspJetPt[ipt][iflv]=new TH1F(hname.c_str(),";p_{T}/p_{T}^{ref}",
					      nbinsrelrsp,relrspmin,relrspmax);
	    }
	    
	    if (dorelrsp&&dorefpt) {
	      hname=flavor[iflv]+"RelRsp_"+jetEtaSuffix+"_"+refPtSuffix;
	      relRspRefPt[ipt][iflv]=new TH1F(hname.c_str(),";p_{T}/p_{T}^{ref}",
					      nbinsrelrsp,relrspmin,relrspmax);
	    }
	    
	    if (doabsrsp&&dojetpt) {
	      hname=flavor[iflv]+"AbsRsp_"+jetEtaSuffix+"_"+jetPtSuffix;
	      absRspJetPt[ipt][iflv]=new TH1F(hname.c_str(),
					      ";p_{T}-p_{T}^{ref} [GeV]",
					      nbinsabsrsp,absrspmin,absrspmax);
	    }
	    
	    if (doabsrsp&&dorefpt) {
	      hname=flavor[iflv]+"AbsRsp_"+jetEtaSuffix+"_"+refPtSuffix;
	      absRspRefPt[ipt][iflv]=new TH1F(hname.c_str(),
					      ";p_{T}-p_{T}^{ref} [GeV]",
					      nbinsabsrsp,absrspmin,absrspmax);
	    }
	  }
	}
	if (dojetpt)            jetPtVsJetEtaJetPt .push_back(jetPtJetPt);
	if (dorefpt)            refPtVsJetEtaRefPt .push_back(refPtRefPt);
	if (dorelrsp&&dojetpt) relRspVsJetEtaJetPt.push_back(relRspJetPt);
	if (dorelrsp&&dorefpt) relRspVsJetEtaRefPt.push_back(relRspRefPt);
	if (doabsrsp&&dojetpt) absRspVsJetEtaJetPt.push_back(absRspJetPt);
	if (doabsrsp&&dorefpt) absRspVsJetEtaRefPt.push_back(absRspRefPt);
      }
    }
    
    //
    // fill histograms
    //
    unsigned int nevt = (unsigned int)tree->GetEntries();
    for (unsigned int ievt=0;ievt<nevt;ievt++) {
      tree->GetEntry(ievt);
      for (unsigned char iref=0;iref<nref;iref++) {
	if (( dobalance&&refdphijt[iref]<dphimin)||
	    (!dobalance&&refdrjt[iref]>drmax)) continue;
	
	float absrsp = jtpt[iref]-refpt[iref];
	float relrsp = jtpt[iref]/refpt[iref];
	
	if (jteta[iref]>=etabarrelmin&&jteta[iref]<=etabarrelmax) {
	  if (dorefpt) {
	    fill_histo(refpt[iref],refpt[iref],binspt,refPtVsRefPtBarrel);
	    if (doflavor)
	      fill_histo(refpdgid[iref],
			 refpt[iref],refpt[iref],binspt,refPtVsRefPtBarrel);
	  }
	  if (dorelrsp&&dorefpt) {
	    fill_histo(relrsp,refpt[iref],binspt,relRspVsRefPtBarrel);
	    if (doflavor) fill_histo(refpdgid[iref],
				     relrsp,refpt[iref],binspt,relRspVsRefPtBarrel);
	  }
	  if (doabsrsp&&dorefpt) {
	    fill_histo(absrsp,refpt[iref],binspt,absRspVsRefPtBarrel);
	    if (doflavor) fill_histo(refpdgid[iref],
				     absrsp,refpt[iref],binspt,absRspVsRefPtBarrel);
	  }
	}
	
	if (dojetpt) {
	  fill_histo(jtpt[iref], jtpt[iref], binspt,jetPtVsJetPt);
	  if (doflavor) fill_histo(refpdgid[iref],
				   jtpt[iref], jtpt[iref], binspt,jetPtVsJetPt);
	}
	if (dorefpt) {
	  fill_histo(refpt[iref],refpt[iref],binspt,refPtVsRefPt);
	  if (doflavor) fill_histo(refpdgid[iref],
				   refpt[iref],refpt[iref],binspt,refPtVsRefPt);
	}
	
	fill_histo(jteta[iref],jteta[iref],binseta,jetEtaVsJetEta);
	if (doflavor) fill_histo(refpdgid[iref],
				 jteta[iref],jteta[iref],binseta,jetEtaVsJetEta);
	
	fill_histo(jtphi[iref],jtphi[iref],binsphi,jetPhiVsJetPhi);
	if (doflavor) fill_histo(refpdgid[iref],
				 jtphi[iref],jtphi[iref],binsphi,jetPhiVsJetPhi);
 
	if (dojetpt) {
	  fill_histo(jtpt[iref], jteta[iref],jtpt[iref],
		     binseta,binspt,jetPtVsJetEtaJetPt);
	  if (doflavor)
	    fill_histo(refpdgid[iref],jtpt[iref], jteta[iref],jtpt[iref],
		       binseta,binspt,jetPtVsJetEtaJetPt);
	}
	if (dorefpt) {
	  fill_histo(refpt[iref],jteta[iref],refpt[iref],
		     binseta,binspt,refPtVsJetEtaRefPt);
	  if (doflavor)
	    fill_histo(refpdgid[iref],refpt[iref],jteta[iref],refpt[iref],
		       binseta,binspt,refPtVsJetEtaRefPt);
	}
	
	if (dorelrsp) {
	  if (dojetpt) {
	    fill_histo(relrsp,jtpt[iref], binspt,relRspVsJetPt);
	    if (doflavor) fill_histo(refpdgid[iref],
				     relrsp,jtpt[iref], binspt,relRspVsJetPt);
	  }
	  if (dorefpt) {
	    fill_histo(relrsp,refpt[iref],binspt,relRspVsRefPt);
	    if (doflavor) fill_histo(refpdgid[iref],
				     relrsp,refpt[iref],binspt,relRspVsRefPt);
	  }

	  fill_histo(relrsp,jteta[iref],binseta,relRspVsJetEta);
	  if (doflavor) fill_histo(refpdgid[iref],
				   relrsp,jteta[iref],binseta,relRspVsJetEta);
	  
	  fill_histo(relrsp,jtphi[iref],binsphi,relRspVsJetPhi);
	  if (doflavor) fill_histo(refpdgid[iref],
				   relrsp,jtphi[iref],binsphi,relRspVsJetPhi);
	  
	  if (dojetpt) {
	    fill_histo(relrsp,jteta[iref],jtpt[iref],
		       binseta,binspt,relRspVsJetEtaJetPt);
	    if (doflavor) fill_histo(refpdgid[iref],relrsp,jteta[iref],jtpt[iref],
				     binseta,binspt,relRspVsJetEtaJetPt);
	  }
	  if (dorefpt) {
	    fill_histo(relrsp,jteta[iref],refpt[iref],
		       binseta,binspt,relRspVsJetEtaRefPt);
	    if (doflavor) fill_histo(refpdgid[iref],relrsp,jteta[iref],refpt[iref],
				     binseta,binspt,relRspVsJetEtaRefPt);
	  }
	}
	
	if (doabsrsp) {
	  if (dojetpt) {
	    fill_histo(absrsp,jtpt[iref], binspt,absRspVsJetPt);
	    if (doflavor) fill_histo(refpdgid[iref],
				     absrsp,jtpt[iref], binspt,absRspVsJetPt);
	  }
	  if (dorefpt) {
	    fill_histo(absrsp,refpt[iref],binspt,absRspVsRefPt);
	    if (doflavor) fill_histo(refpdgid[iref],
				     absrsp,refpt[iref],binspt,absRspVsRefPt);
	  }

	  fill_histo(absrsp,jteta[iref],binseta,absRspVsJetEta);
	  if (doflavor) fill_histo(refpdgid[iref],
				   absrsp,jteta[iref],binseta,absRspVsJetEta);

	  fill_histo(absrsp,jtphi[iref],binsphi,absRspVsJetPhi);
	  if (doflavor) fill_histo(refpdgid[iref],
				   absrsp,jtphi[iref],binsphi,absRspVsJetPhi);
	  
	  if (dojetpt) {
	    fill_histo(absrsp,jteta[iref],jtpt[iref],
		       binseta,binspt,absRspVsJetEtaJetPt);
	    if (doflavor) fill_histo(refpdgid[iref],absrsp,jteta[iref],jtpt[iref],
				     binseta,binspt,absRspVsJetEtaJetPt);
	  }
	  if (dorefpt) {
	    fill_histo(absrsp,jteta[iref],refpt[iref],
		       binseta,binspt,absRspVsJetEtaRefPt);
	    if (doflavor) fill_histo(refpdgid[iref],absrsp,jteta[iref],refpt[iref],
				     binseta,binspt,absRspVsJetEtaRefPt);
	  }
	}
	
      }
    }
    
    cout<<" DONE."<<endl;
  }


  //
  // close files
  //
  cout<<"close output file "<<output<<" ... "<<flush;
  ofile->Write();
  gROOT->GetListOfFiles()->Remove(ofile);
  ofile->Close();
  delete ofile;
  cout<<"DONE."<<endl;
  
  ifile->Close();
  delete ifile;

  return 0;
}



////////////////////////////////////////////////////////////////////////////////
// implement local functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
string get_suffix(const string& varname,int ibin,const vector<float>& bins)
{
  stringstream ss; ss<<varname<<bins[ibin]<<"to"<<bins[ibin+1];
  return ss.str();
}


//______________________________________________________________________________
int get_index(float x,const vector<float>& binsx)
{
  for (unsigned int ix=0;ix<binsx.size()-1;ix++)
    if (x>=binsx[ix]&&x<binsx[ix+1]) return ix;
  return -1;
}


//______________________________________________________________________________
void fill_histo(float value,float x,
		const vector<float>& binsx,const vector<TH1F**>& histos)
{
  int ix=get_index(x,binsx);
  if (ix>=0) histos[ix][0]->Fill(value);
}


//______________________________________________________________________________
void fill_histo(int pdgid,float value,float x,
		const vector<float>& binsx,const vector<TH1F**>& histos)
{
  int abspdgid=std::abs(pdgid);
  int iflv(-1);
  if (abspdgid>=1&&abspdgid<=3) iflv=1;
  else if (abspdgid== 4)        iflv=2;
  else if (abspdgid== 5)        iflv=3;
  else if (abspdgid==21)        iflv=4;
  else return;

  int ix=get_index(x,binsx);
  if (ix>=0) histos[ix][iflv]->Fill(value);
}


//______________________________________________________________________________
void fill_histo(float value,float x,float y,
		const vector<float>& binsx,const vector<float>& binsy,
		const vector<TH1F***>& histos)
{
  int ix=get_index(x,binsx);
  int iy=get_index(y,binsy);
  if (ix>=0&&iy>=0) histos[ix][iy][0]->Fill(value);
}


//______________________________________________________________________________
void fill_histo(int pdgid,float value,float x,float y,
		const vector<float>& binsx,const vector<float>& binsy,
		const vector<TH1F***>& histos)
{
  int abspdgid=std::abs(pdgid);
  int iflv(-1);
  if (abspdgid>=1&&abspdgid<=3) iflv=1;
  else if (abspdgid== 4)        iflv=2;
  else if (abspdgid== 5)        iflv=3;
  else if (abspdgid==21)        iflv=4;
  else return;

  int ix=get_index(x,binsx);
  int iy=get_index(y,binsy);
  if (ix>=0&&iy>=0) histos[ix][iy][iflv]->Fill(value);
}


//______________________________________________________________________________
bool contains(const vector<string>& collection,const string& element)
{
  vector<string>::const_iterator it;
  for (it=collection.begin();it!=collection.end();++it)
    if ((*it)==element) return true;
  return false;
}
