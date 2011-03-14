#include "anaXS.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/hpl.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/functions.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"

#include "TF1.h"
#include "TKey.h"
#include "TTree.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TVirtualPad.h"  // access to gPad
#include "TCanvas.h"

#include <iomanip>

using std::cout;
using std::endl;
using std::string;
using std::vector;

// ----------------------------------------------------------------------
// Usage: 
// -----
// 
// lxplus308>r
// Loading libPhysics.so
// root [0] gSystem->Load("lib/libUtil.so");                     
// root [1] gSystem->Load("lib/libAnaClasses.so");
// root [2] anaXS a("/afs/cern.ch/user/b/bora/scratch0/root/xs", 40)            
// --> Loading rootfiles in /afs/cern.ch/user/u/ursl/public/root/ for 40
// root [3] a.makeAll(3)                                          
//
// NOTE: (1) this works for ME, YOU have to copy the files and directories in ~ursl/public/root
//       (2) If it does not work the first time, .q, and try again... (don't ask, I don't know)
// ----------------------------------------------------------------------


ClassImp(anaXS)

// ----------------------------------------------------------------------
anaXS::anaXS(const char *dir, int i) { 
  init(dir, i);
  if ( gRandom ) delete gRandom;
  gRandom = new TRandom3();
}

// ----------------------------------------------------------------------
void anaXS::init(const char *dir, int i) {

  fFont = 42; 
  fMode = i;  

  c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (c0 == 0) {
    cout << "TCanvas c0 not found. Creating my own version." << endl;
    c0 = new TCanvas("c0","--c0--",356,0,656,700);
  }
  tl  = new TLatex();
  tl->SetNDC(kTRUE); 
  tl->SetTextSize(0.07);
  tl->SetTextFont(fFont);

  pl  = new TLine();
  pa  = new TArrow();
  box = new TBox();

  f0 = new TF1("f0", f_p1, 0., 12., 2); 

  f1 = new TF1("f1", f_p1aG, 0., 12., 5);
  f1->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 

  f2 = new TF1("f2", f_Upsilons, 0., 12., 9);  // signal + bg
  f2->SetParNames("Area 1", "Peak 1", "Sigma 1", 
		  "Area 2", "Peak 2", "Sigma 2", 
		  "Area 3", "Offset", "Slope"); 

  f3 = new TF1("f3", f_argus, 5.2, 5.29, 2); // for bg 
  f4 = new TF1("f4", f_argus, 5.2, 5.29, 2); // for bg variation
  f5 = new TF1("f5", f_expo, -0.4, 0.4, 2);
  
  f6 = new TF1("f6", f_p1acb, 0., 12., 7);
  f6->SetParNames("Mean", "Sigma", "Alpha", "n", "N","Offset", "Slope");
  
  f7 = new TF1("f7", f_fnov, 0., 1., 4);
  f8 = new TF1("f8", f_fnov, 0., 1., 4);

  f9 = new TF1("f9", f_p1ag, -0.2, 0.2, 5);
  
  f10 = new TF1("f10", f_cb, 0., 12., 5);
  f11 = new TF1("f11", f_cbaG, 0., 12., 8);
  f12 = new TF1("f12", f_p1acbaG, 0., 12., 10);
  f12->SetParNames("Mean", "Sigma_CB", "Alpha", "n", "N", "Area", "Peak", "Sigma_G", "Offset", "Slope");
  
  f13 = new TF1("f13", f_p1a3cb, 0., 12., 11);
  f13->SetParNames("Mean1", "Sigma1", "Alpha1", "n1", "N1", "Mean2", "Sigma2", "N2", "N3", "Offset", "Slope");
  
  fDirectory = string(dir); 
  fPtDirectory = fDirectory + string(Form("/pt-%02d", i)); 

  cout << "--> Loading rootfiles in " << fDirectory <<  " for " << i << endl;
  system(Form("/bin/rm -f %s/anaXS.tex", fDirectory.c_str()));
  fNumbersFileName = fDirectory + string("/anaXS.tex");
  loadFiles(dir, i);

  fS1VectorPos.clear();
  fS1VectorNeg.clear();
  fS2VectorPos.clear();
  fS2VectorNeg.clear();
  fS3VectorPos.clear();
  fS3VectorNeg.clear();

  fS1VectorMcpPos.clear();
  fS1VectorMcpNeg.clear();
  fS2VectorMcpPos.clear();
  fS2VectorMcpNeg.clear();
  
  /////////////////////////
  
  fS1Vector.clear();
  fS2Vector.clear();
  fS3Vector.clear();
  fS4Vector.clear();
  if (fS1Yield) fS1Yield->Reset();
  if (fS2Yield) fS2Yield->Reset();
  if (fS3Yield) fS3Yield->Reset();
  if (fS1YieldCorrected) fS1YieldCorrected->Reset();
  if (fS1YieldComparison) fS1YieldComparison->Reset();
  if (fS1YieldPt) fS1YieldPt->Reset();
  if (fS1YieldEta) fS1YieldEta->Reset();
  if (fAnaEff) fAnaEff->Reset();
  if (fAnaEff_2S) fAnaEff_2S->Reset();
  if (fAnaEff_3S) fAnaEff_3S->Reset();
  if (fAllGenRes) fAllGenRes->Reset();
  if (fAllGenRes_2S) fAllGenRes_2S->Reset();
  if (fAllGenRes_3S) fAllGenRes_3S->Reset();
  if (fAllGenResPt) fAllGenResPt->Reset();
  if (fAllGenResEta) fAllGenResEta->Reset();
  if (fAcceptance) fAcceptance->Reset();
  if (fPreSelEff) fPreSelEff->Reset();
  if (fPreSelEff_2S) fPreSelEff_2S->Reset();
  if (fPreSelEff_3S) fPreSelEff_3S->Reset();
  if (fPreSelAfter) fPreSelAfter->Reset();
  if (fPreSelAfter_2S) fPreSelAfter_2S->Reset();
  if (fPreSelAfter_3S) fPreSelAfter_3S->Reset();
  if (fPreSelBefore) fPreSelBefore->Reset();
  if (fPreSelBefore_2S) fPreSelBefore_2S->Reset();
  if (fPreSelBefore_3S) fPreSelBefore_3S->Reset();
  if (fRecoGenRes) fRecoGenRes->Reset();
  if (fRecoGenRes_2S) fRecoGenRes_2S->Reset();
  if (fRecoGenRes_3S) fRecoGenRes_3S->Reset();
  if (fAcceptanceProjPt) fAcceptanceProjPt->Reset();

  /////////////////////////
  
  if (fS1YieldPos) fS1YieldPos->Reset();
  if (fS1YieldNeg) fS1YieldNeg->Reset();
  if (fS2YieldPos) fS2YieldPos->Reset();
  if (fS2YieldNeg) fS2YieldNeg->Reset();
  if (fS3YieldPos) fS3YieldPos->Reset();
  if (fS3YieldNeg) fS3YieldNeg->Reset();
  
  if (fS1MctPos) fS1MctPos->Reset();
  if (fS1MctNeg) fS1MctNeg->Reset();
  if (fS2MctPos) fS2MctPos->Reset();
  if (fS2MctNeg) fS2MctNeg->Reset();
  if (fS1McpPos) fS1McpPos->Reset();
  if (fS1McpNeg) fS1McpNeg->Reset();
  if (fS2McpPos) fS2McpPos->Reset();
  if (fS2McpNeg) fS2McpNeg->Reset();

  fPtTnpNeg = new PidTable(1); 
  fPtTnpPos = new PidTable(1); 
  fPtMctNeg = new PidTable(1); 
  fPtMctPos = new PidTable(1); 
  fPtMcpNeg = new PidTable(1); 
  fPtMcpPos = new PidTable(1); 
  fPtMmbNeg = new PidTable(1); 
  fPtMmbPos = new PidTable(1); 
  
  //fPtTrigCorr = new PidTable("PidTables/DATA/Upsilon/PtTrigCorr.dat");
  //fPtMuidCorr = new PidTable("PidTables/DATA/Upsilon/PtMuidCorr.dat");
  //fPtTrigCorr = new PidTable("PidTables/DATA/Upsilon/PtTrigCorr_1Sbin.dat");
  //fPtMuidCorr = new PidTable("PidTables/DATA/Upsilon/PtMuidCorr_1Sbin.dat");
  fPtTrigCorr = new PidTable("PidTables/DATA/Upsilon/PtTrigCorr_2Sbin.dat");
  fPtMuidCorr = new PidTable("PidTables/DATA/Upsilon/PtMuidCorr_2Sbin.dat");  
  //fPtTrigCorr = new PidTable("PidTables/DATA/Upsilon/PtTrigCorr_3Sbin.dat");
  //fPtMuidCorr = new PidTable("PidTables/DATA/Upsilon/PtMuidCorr_3Sbin.dat");    
  
}


// ----------------------------------------------------------------------
void anaXS::loadFiles(const char *dir, int i) {

  // -- Upsilon merging
  if (0 == i) {
    string ufile = fDirectory + string("/") + string("upsilon/101201.fl10.mm.ups1s.xsReader_3Sbin.default.root");    

    fM[0] = new TFile(ufile.c_str()); lM[0] = 1.;
    ufile = fDirectory + string("/") + string("upsilon/101201.fl10.mm.ups2s.xsReader_3Sbin.default.root");
    fM[1] = new TFile(ufile.c_str()); lM[1] = 1.66; 
    ufile = fDirectory + string("/") + string("upsilon/101201.fl10.mm.ups3s.xsReader_3Sbin.default.root");
    fM[2] = new TFile(ufile.c_str()); lM[2] = 3.43; 
    cout << "Got the Files for Merging" << endl;
  }

  // -- all lumis in /pb
  if (0 != i) {
    //    string ufile = fDirectory + string("/upsilon/upsilon_hemispheres.root");
    string ufile;
    string jfile; 
    if (1 == i) {
      jfile = fDirectory + string("/jpsi/jpsi.root");
    } else if (2 == i) {
      jfile = fDirectory + string("/jpsi/JPsi1_64_new.root");          // if (tp > 4)
    } else if (3 == i) {
      jfile = fDirectory + string("/jpsi/Jpsi_sel1.root");             // if (TMath::Abs(te - pe) < 0.1)
    } else if (4 == i) {
      jfile = fDirectory + string("/jpsi/Jpsi_sel2.root");             // if (TMath::Abs(te) - TMath::Abs(pe) < 0.2 && tp>4.)
    } else if (5 == i) {
      jfile = fDirectory + string("/jpsi/Jpsi_sel1_newbinning.root");  // if (TMath::Abs(te - pe) < 0.1)
    } else if (6 == i) {
      jfile = fDirectory + string("/jpsi/Jpsi_sel2_newbinning.root");  // if (TMath::Abs(te) - TMath::Abs(pe) < 0.2 && tp>4.)
    } else if (7 == i) {
      jfile = fDirectory + string("/jpsi/Jpsi_0010.root");             // |te - pe| < 0.2 and tp>4, L ~ 0.75 pb^-1
    } else if (8 == i) {
      jfile = fDirectory + string("/jpsi/JPsi_ptcut0.4_etacut0.2.root");// |te - pe| < 0.2 and tp>4, L ~ 0.75 pb^-1. new version!
    } else if (9 == i) {
      jfile = fDirectory + string("/jpsi/JPsi_ptcut6.root");            // pt > 6 , L ~ 1.5 pb^-1
    } else if (10 == i) {
      jfile = fDirectory + string("/jpsi/JPsi_pt8.root");               // pt > 8, l ~ 1.5 pb^-1
    } else if (11 == i) {
      jfile = fDirectory + string("/jpsi/JPsi_upsbin_pt6.root");        // pt > 6, l ~ 1.5 pb^-1
    } else if (12 == i) {
      jfile = fDirectory + string("/jpsi/JPsi_upsbin_pt8.root");        // pt > 8, l ~ 1.5 pb^-1
    } else if (13 == i) {
      jfile = fDirectory + string("/jpsi/JPsi_upsbin_pt10.root");       // pt > 10, l ~ 1.5 pb^-1
    } else if (14 == i) {
      jfile = fDirectory + string("/jpsi/JPsi_upsbin_eta0.7.root");     // eta>0.7, l ~ 1.5 pb^-1 ??
    } else if (15 == i) {
      jfile = fDirectory + string("/jpsi/JPsi_upsbin_eta0.7_pt6.root");  // eta>0.7, pT>6 l ~ 1.5 pb^-1 ??
    } else if (20 == i) {
      ufile = fDirectory + string("/upsilon/summer09_upsilon_eta0.5.root");
      ufile = fDirectory + string("/upsilon/upsilon_hemispheres.root");     
      jfile = fDirectory + string("/jpsi/JPsi_upsbin_pt6.root");  
    } else if (30 == i) {
      ufile = fDirectory + string("/upsilon/UpsTagandprobe_10TeV_sameHem.root");
      jfile = fDirectory + string("/jpsi/JpsiTagandprobe_10TeV_pt6.root");  
    } else if (31 == i) {
      ufile = fDirectory + string("/upsilon/UpsTagandprobe_10TeV_nocut.root");
      jfile = fDirectory + string("/jpsi/JpsiTagandprobe_10TeV_nocut.root");  
    } else if (40 == i) {
      //ufile = fDirectory + string("/upsilon/101201.fl10.mm.ups1s.xsReader.trueyield_100.default.root");
      //ufile = fDirectory + string("/upsilon/101201.fl10.mm.ups1s.xsReader_1Sbin.default.root");
      ufile = fDirectory + string("/upsilon/101201.fl10.mm.COMBINED.xsReader_1Sbin.default.root");
      //jfile = fDirectory + string("/jpsi/111112.dimuons.xsReader.default.root");
      jfile = fDirectory + string("/upsilon/130211.nov4rereco_v2.dimuons.xsReader_1Sbin.default.root");
     
    } else {
      cout << "Don't know which J/psi file to open for i = " << i << ". Specify this in anaXS::loadfiles()" << endl;
      return;
    }
    fM[0] = new TFile(ufile.c_str());  lM[0] = 1.0; lM[2] = 3.43;
    fM[1] = new TFile(jfile.c_str());  lM[1] = 1.66;
  }

}


// ----------------------------------------------------------------------
void anaXS::combineUpsilons() {
  
  string ufile = fDirectory + string("/upsilon/101201.fl10.mm.COMBINED.xsReader_3Sbin.default.root");
  TFile *f = new TFile(ufile.c_str(), "RECREATE"); 

  fM[0]->cd();

  TH1D *h0, *h1, *h2; 
  TH2D *H0, *H1, *H2; 

  int ptbin(-1); int ybin(-1);
  double PTbin[20]; double Ybin[10];
  
  TObject *obj;
  TKey    *key;
  TIter next(gFile->GetListOfKeys());

  while ((key = (TKey*)next())) {
    obj = key->ReadObj();
    if (obj->InheritsFrom(TH1D::Class())) {
      h0 = (TH1D*)obj;
      cout << "TH1D: " << h0->GetName() << endl;
      
      for ( int i=1; i <= ptbin; i++ ){
	for ( int j=1; j <= ybin; j++ ){      
      
	  char *HistName = Form("AnaEff_%.1dS,rapidity%.1f_%.1f,pt%.1f_%.1f", 1, Ybin[j], Ybin[j+1], PTbin[i], PTbin[i+1]);
	  if ( !strcmp(h0->GetName(),HistName) ){
	    cout << HistName  << endl;
	    h1 = (TH1D*)(fM[1]->Get(Form("AnaEff_%.1dS,rapidity%.1f_%.1f,pt%.1f_%.1f", 2, Ybin[j], Ybin[j+1], PTbin[i], PTbin[i+1])));
	    h2 = (TH1D*)(fM[2]->Get(Form("AnaEff_%.1dS,rapidity%.1f_%.1f,pt%.1f_%.1f", 3, Ybin[j], Ybin[j+1], PTbin[i], PTbin[i+1])));
	    h1->SetDirectory(f);
	    h2->SetDirectory(f);
	    goto end;
	  }
	  
	  char *histName = Form("MuIDEff_%.1dS,rapidity%.1f_%.1f,pt%.1f_%.1f", 1, Ybin[j], Ybin[j+1], PTbin[i], PTbin[i+1]);
	  if ( !strcmp(h0->GetName(),histName) ){
	    cout << histName  << endl;
	    h1 = (TH1D*)(fM[1]->Get(Form("MuIDEff_%.1dS,rapidity%.1f_%.1f,pt%.1f_%.1f", 2, Ybin[j], Ybin[j+1], PTbin[i], PTbin[i+1])));
	    h2 = (TH1D*)(fM[2]->Get(Form("MuIDEff_%.1dS,rapidity%.1f_%.1f,pt%.1f_%.1f", 3, Ybin[j], Ybin[j+1], PTbin[i], PTbin[i+1])));
	    h1->SetDirectory(f);
	    h2->SetDirectory(f);
	    goto end;
	  }
	  
	  char *histname = Form("TrigEff_%.1dS,rapidity%.1f_%.1f,pt%.1f_%.1f", 1, Ybin[j], Ybin[j+1], PTbin[i], PTbin[i+1]);
	  if ( !strcmp(h0->GetName(),histname) ){
	    cout << histname  << endl;
	    h1 = (TH1D*)(fM[1]->Get(Form("TrigEff_%.1dS,rapidity%.1f_%.1f,pt%.1f_%.1f", 2, Ybin[j], Ybin[j+1], PTbin[i], PTbin[i+1])));
	    h2 = (TH1D*)(fM[2]->Get(Form("TrigEff_%.1dS,rapidity%.1f_%.1f,pt%.1f_%.1f", 3, Ybin[j], Ybin[j+1], PTbin[i], PTbin[i+1])));
	    h1->SetDirectory(f);
	    h2->SetDirectory(f);
	    goto end;
	  }	  

	}
      }
      
      if ( !strcmp(h0->GetName(),"MuIDEff_1S,OverAll") ){
	h1 = (TH1D*)(fM[1]->Get("MuIDEff_2S,OverAll"));
	h2 = (TH1D*)(fM[2]->Get("MuIDEff_3S,OverAll"));
	h1->SetDirectory(f);
	h2->SetDirectory(f);
	goto end;
      }
      
      if ( !strcmp(h0->GetName(),"TrigEff_1S,OverAll") ){
	h1 = (TH1D*)(fM[1]->Get("TrigEff_2S,OverAll"));
	h2 = (TH1D*)(fM[2]->Get("TrigEff_3S,OverAll"));
	h1->SetDirectory(f);
	h2->SetDirectory(f);
	goto end;
      }     
      
      if ( !strcmp(h0->GetName(),"AnaEff_1S,OverAll") ){
	h1 = (TH1D*)(fM[1]->Get("AnaEff_2S,OverAll"));
	h2 = (TH1D*)(fM[2]->Get("AnaEff_3S,OverAll"));
	h1->SetDirectory(f);
	h2->SetDirectory(f);
	goto end;
      }
      
      
      h1 = (TH1D*)(fM[1]->Get(h0->GetName()));
      h2 = (TH1D*)(fM[2]->Get(h0->GetName()));

      if (0) {
	zone(2,2);
	h0->DrawCopy();
	c0->cd(2);
	h1->Draw();
	c0->cd(3);
	h2->Draw();
      }

      h1->Scale(lM[0]/lM[1]);
      h2->Scale(lM[0]/lM[2]);
      h0->Add(h1);
      h0->Add(h2);
      
     // h0->Scale(5.);
    end:
      integerEntries(h0);
      
      h0->SetDirectory(f);

      if (0) {
	c0->cd(4);
	h0->Draw();
	
	c0->Update(); 
	c0->Draw();
      }
    } else if (obj->InheritsFrom(TH2D::Class())) {
      H0 = (TH2D*)obj;
      cout << "TH2D: " << H0->GetName() << endl;
      
      if ( !strcmp(H0->GetName(),"mt,pt-eta") ){
	
	ptbin = H0->GetXaxis()->GetNbins();
	ybin = H0->GetYaxis()->GetNbins();
	//cout << ptbin << " " << ybin << endl;
	
	for ( int i=1; i <= ptbin; i++ ){
	  PTbin[i]=H0->GetXaxis()->GetBinLowEdge(i);
	  if ( i == ptbin  ){
	    PTbin[i+1]=H0->GetXaxis()->GetBinUpEdge(i);
	  }
	}
	
	for ( int i=1; i <= ybin; i++ ){
	  Ybin[i]=H0->GetYaxis()->GetBinLowEdge(i);
	  if ( i == ybin  ){
	    Ybin[i+1]=H0->GetYaxis()->GetBinUpEdge(i);
	  }
	}
	
      }
      
      for ( int i=1; i <= ptbin+1; i++ ){
	for ( int j=1; j <= ybin+1; j++ ){
	  //cout<< PTbin[i] << "  " << Ybin[j]  << endl;
	}
      } 
      
      if ( !strcmp(H0->GetName(),"AllGenRes_1S") ){
	H1 = (TH2D*)(fM[1]->Get("AllGenRes_2S"));
	H2 = (TH2D*)(fM[2]->Get("AllGenRes_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }
      
      if ( !strcmp(H0->GetName(),"RecoGenRes_1S") ){
	H1 = (TH2D*)(fM[1]->Get("RecoGenRes_2S"));
	H2 = (TH2D*)(fM[2]->Get("RecoGenRes_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }
      
      if ( !strcmp(H0->GetName(),"PreSel_afterVtx_1S") ){
	H1 = (TH2D*)(fM[1]->Get("PreSel_afterVtx_2S"));
	H2 = (TH2D*)(fM[2]->Get("PreSel_afterVtx_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }      
      
      if ( !strcmp(H0->GetName(),"PreSel_beforeVtx_1S") ){
	H1 = (TH2D*)(fM[1]->Get("PreSel_beforeVtx_2S"));
	H2 = (TH2D*)(fM[2]->Get("PreSel_beforeVtx_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }         
      
      if ( !strcmp(H0->GetName(),"MuIDCheck_after_1S") ){
	H1 = (TH2D*)(fM[1]->Get("MuIDCheck_after_2S"));
	H2 = (TH2D*)(fM[2]->Get("MuIDCheck_after_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }            
      
      if ( !strcmp(H0->GetName(),"MuIDCheck_before_1S") ){
	H1 = (TH2D*)(fM[1]->Get("MuIDCheck_before_2S"));
	H2 = (TH2D*)(fM[2]->Get("MuIDCheck_before_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }         
      
      if ( !strcmp(H0->GetName(),"MuIDCheck_Numa_1S") ){
	H1 = (TH2D*)(fM[1]->Get("MuIDCheck_Numa_2S"));
	H2 = (TH2D*)(fM[2]->Get("MuIDCheck_Numa_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }            
      
      if ( !strcmp(H0->GetName(),"MuIDCheck_Deno_1S") ){
	H1 = (TH2D*)(fM[1]->Get("MuIDCheck_Deno_2S"));
	H2 = (TH2D*)(fM[2]->Get("MuIDCheck_Deno_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }       
      
      if ( !strcmp(H0->GetName(),"TrigCheck_after_1S") ){
	H1 = (TH2D*)(fM[1]->Get("TrigCheck_after_2S"));
	H2 = (TH2D*)(fM[2]->Get("TrigCheck_after_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }            
      
      if ( !strcmp(H0->GetName(),"TrigCheck_before_1S") ){
	H1 = (TH2D*)(fM[1]->Get("TrigCheck_before_2S"));
	H2 = (TH2D*)(fM[2]->Get("TrigCheck_before_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }         
      
      if ( !strcmp(H0->GetName(),"TrigCheck_Numa_1S") ){
	H1 = (TH2D*)(fM[1]->Get("TrigCheck_Numa_2S"));
	H2 = (TH2D*)(fM[2]->Get("TrigCheck_Numa_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }            
      
      if ( !strcmp(H0->GetName(),"TrigCheck_Deno_1S") ){
	H1 = (TH2D*)(fM[1]->Get("TrigCheck_Deno_2S"));
	H2 = (TH2D*)(fM[2]->Get("TrigCheck_Deno_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }        
      
      if ( !strcmp(H0->GetName(),"TrueYield_1S") ){
	H1 = (TH2D*)(fM[1]->Get("TrueYield_2S"));
	H2 = (TH2D*)(fM[2]->Get("TrueYield_3S"));
	H0->SetDirectory(f);
	H1->SetDirectory(f);
	H2->SetDirectory(f);
	continue;
      }       
      
      H1 = (TH2D*)(fM[1]->Get(H0->GetName()));
      H2 = (TH2D*)(fM[2]->Get(H0->GetName()));
      
      H1->Scale(lM[0]/lM[1]);
      H2->Scale(lM[0]/lM[2]);
      H0->Add(H1);
      H0->Add(H2);
      
     // H0->Scale(5.);
      H0->SetDirectory(f);
    }
  }
  
  f->Write();
  f->Close();
}

void anaXS::makeAllMC(int channel) {
  
  // -- Upsilon
  if (channel & 1) {
    init(fDirectory.c_str(), fMode);
    // -- fill histograms
    fSample = string("upsilon");
    ReadHistograms(fM[0], "UpsilonMass", "AnaEff_1S", "AnaEff_2S", "AnaEff_3S", "MuIDEff_1S", "MuIDEff_2S", "MuIDEff_3S", "TrigEff_1S", "TrigEff_2S", "TrigEff_3S", "mt,pt-eta");
    
    // -- add backgrounds
    addBackground(fS1Vector, 0.3);
    
    
    Pull(1);
    Pull(2);
    //FITUpsilon(1);
    //GetAnaEff();
    //GetPreSelEff();
    //GetMuIDEff(1);
    //GetTrigEff(1);
    //CorrectedYields(1);   // 1- FOR MC, 2 FOR DATA
    //PlotProjections(1);   // 1- FOR MC, 2 FOR DATA
        
  }

  // -- Trigger Eff Study
  if (channel & 2) {
    init(fDirectory.c_str(), fMode);
    // -- fill histograms
    fSample = string("jpsi");
    ReadHistograms(fM[0], "UpsilonMass", "AnaEff_1S", "AnaEff_2S", "AnaEff_3S", "MuIDEff_1S", "MuIDEff_2S", "MuIDEff_3S", "TrigEff_1S", "TrigEff_2S", "TrigEff_3S", "mt,pt-eta");
    
    GetTrigEff(1);
        
  }  
  
  
  
}

void anaXS::makeAllDATA(int channel) {
  
  // -- Upsilon
  if (channel & 1) {
    init(fDirectory.c_str(), fMode);
    // -- fill histograms
    fSample = string("upsilon");
    
    ReadHistogramsDATA0(fM[0], "AnaEff_1S", "AnaEff_2S", "AnaEff_3S", "mt,pt-eta");
    ReadHistogramsDATA1(fM[1], "UpsilonMass", "MuIDEff_1S",  "TrigEff_1S",  "mt,pt-eta");
    
    
    // -- add backgrounds
    //addBackground(fS1Vector, 0.3);
    
    
    FITUpsilon(1);
    GetAnaEff();
    GetPreSelEff();
    GetMuIDEff(2);
    GetTrigEff(2);
    CorrectedYields(2);   // 1- FOR MC, 2 FOR DATA
    PlotProjections(2);   // 1- FOR MC, 2 FOR DATA
    
  }

  
  
  
}


// ----------------------------------------------------------------------
void anaXS::makeAll(int channel) {

  // -- bias plots
  if (0) {
    //    biasPlots("/data/ursl/muonID-rootfiles/biased/upsilon/ups_biasstudy.root", "upsilon-all.png", 15);
    biasPlots("/data/ursl/muonID-rootfiles/biased/jpsi/jpsi_biasstudy.root", "jpsi-all.eps", 15);
    return;
  }

  // -- total mass
  if (0) {
    // -- fill histograms
    fSample = string("upsilon");
    readHistograms(fM[0], "mm", "mt", "mmbar", "Matched", "mt,pt-eta");
    
    // -- add backgrounds
    addBackground(fS1VectorPos, 0.3);
    addBackground(fS2VectorPos, 1.);
    
    addBackground(fS1VectorNeg, 0.3);
    addBackground(fS2VectorNeg, 1.);
    
    totalMass();

    init(fDirectory.c_str(), fMode);
    fSample = string("jpsi");
    readHistograms(fM[1], "mm", "mt", "mmbar", "Matched", "mt,pt-eta");
    McpYields();
    
    // -- add backgrounds
    addBackground(fS1VectorPos, 0.3);
    addBackground(fS2VectorPos, 1.);
    
    addBackground(fS1VectorNeg, 0.3);
    addBackground(fS2VectorNeg, 1.);

    totalMass();

  }


  // -- Upsilon
  if (channel & 1) {
    init(fDirectory.c_str(), fMode);
    // -- fill histograms
    fSample = string("upsilon");
    readHistograms(fM[0], "mm", "mt", "mmbar", "Matched", "mt,pt-eta");
    McpYields();
    
    // -- add backgrounds
    addBackground(fS1VectorPos, 0.3);
    addBackground(fS2VectorPos, 1.);
    addBackground(fS3VectorPos, 1.);
    
    addBackground(fS1VectorNeg, 0.3);
    addBackground(fS2VectorNeg, 1.);
    addBackground(fS3VectorNeg, 1.);

    fitUpsilon(2);
    fillPidTables(); 
    //  validation();
    projections(); 
    allDifferences(1); 

  }


  // -- J/psi
  if (channel & 2) {    
    init(fDirectory.c_str(), fMode);
    // -- fill histograms
    fSample = string("jpsi");
    readHistograms(fM[1], "mm", "mt", "mmbar", "Matched", "mt,pt-eta");
    McpYields();
    
    // -- add backgrounds
    addBackground(fS1VectorPos, 0.3);
    addBackground(fS2VectorPos, 1.);
    addBackground(fS3VectorPos, 1.);
    
    addBackground(fS1VectorNeg, 0.3);
    addBackground(fS2VectorNeg, 1.);
    addBackground(fS3VectorNeg, 1.);
    
    fitJpsi(1);
    fillPidTables(); 
    //  validation();
    projections();    
   // TFile *file = new TFile("ptDifference.root","RECREATE"); 
    allDifferences(2); 
  }

  if (1) {
    init(fDirectory.c_str(), fMode);
    fSample = string("upsilon");
    readPidTables("upsilon"); 
    projections(); 
    //    allDifferences(1); 

    init(fDirectory.c_str(), fMode);
    fSample = string("jpsi");
    readPidTables("jpsi"); 
    projections(); 
    //    allDifferences(2); 

  }
  
  if ((channel&1) && (channel&2)) allDifferences(3); 
  


}

void anaXS::HistoFlip2D(TH2D *H){
  
  TH2D *hflip;
  int bin(-1);
  if ( !strcmp(H->GetName(),"AllGenRes_1S")  ){
    hflip = new TH2D("AllGenRes_1S", "AllGenRes_1S", 
		     H->GetNbinsY(), H->GetYaxis()->GetXbins()->GetArray(),
		     H->GetNbinsX(), H->GetXaxis()->GetXbins()->GetArray()
		     );
  }

  if ( !strcmp(H->GetName(),"RecoGenRes_1S")  ){
    hflip = new TH2D("RecoGenRes_1S", "RecoGenRes_1S", 
		     H->GetNbinsY(), H->GetYaxis()->GetXbins()->GetArray(),
		     H->GetNbinsX(), H->GetXaxis()->GetXbins()->GetArray()
		     );
  }
  //makeCanvas(2);
  //c2->Divide(3,1);
  //c2->cd(1);
  //H->Draw("colz");
  //cout << "H->GetNbinsX() = " << H->GetNbinsX() <<" H->GetNbinsY() = " << H->GetNbinsY()<< endl;
  for ( int ipt = 1; ipt <= H->GetNbinsX(); ++ipt ){
    for ( int iy = 1; iy <= H->GetNbinsY(); ++iy ){
      
       bin = H->GetBinContent(ipt,iy);
       //cout << ipt <<"  "<< iy <<"  "<<  bin << endl;
       hflip->SetBinContent(iy,ipt,bin); 
      
    }
  }
  //c2->cd(2);
  //hflip->Draw("colz");
  if ( !strcmp(H->GetName(),"AllGenRes_1S")  ) fAllGenRes = hflip;
  if ( !strcmp(H->GetName(),"RecoGenRes_1S")  ) fRecoGenRes = hflip;
    //c2->cd(3);
  //H->Draw("colz"); 
  
}

void anaXS::plotAcceptance(){
  
  double deno(-1.); double numa(-1.); double acceptance(-1);
  
  for ( int iy = 1; iy <= fAllGenRes->GetNbinsX(); ++iy ){
    for ( int ipt = 1; ipt <= fAllGenRes->GetNbinsY(); ++ipt ){
      deno = fAllGenRes->GetBinContent(iy,ipt);
      numa = fRecoGenRes->GetBinContent(iy,ipt);
      acceptance = numa/deno;
      fAcceptance->SetBinContent(iy,ipt,acceptance);
    }
  }
  
  for ( int iy = 1; iy <= fAllGenRes_2S->GetNbinsX(); ++iy ){
    for ( int ipt = 1; ipt <= fAllGenRes_2S->GetNbinsY(); ++ipt ){
      deno = fAllGenRes_2S->GetBinContent(iy,ipt);
      numa = fRecoGenRes_2S->GetBinContent(iy,ipt);
      acceptance = numa/deno;
      fAcceptance_2S->SetBinContent(iy,ipt,acceptance);
    }
  }  
  
  for ( int iy = 1; iy <= fAllGenRes_3S->GetNbinsX(); ++iy ){
    for ( int ipt = 1; ipt <= fAllGenRes_3S->GetNbinsY(); ++ipt ){
      deno = fAllGenRes_3S->GetBinContent(iy,ipt);
      numa = fRecoGenRes_3S->GetBinContent(iy,ipt);
      acceptance = numa/deno;
      fAcceptance_3S->SetBinContent(iy,ipt,acceptance);
    }
  }  
  
  makeCanvas(1);
  c1->Divide(3,1);
  c1->cd(1);
  fAcceptance->Draw("colz");
  c1->cd(2);
  fAcceptance_2S->Draw("colz");
  c1->cd(3);
  fAcceptance_3S->Draw("colz");  
  
}

void anaXS::CorrectedYields(int mode){
  
  makeCanvas(2);
  c2->Divide(2,2);
  //c2->cd(1);
  //fS1Yield->Draw("colz");
  //c2->cd(2);
  //fAnaEff->Draw("colz");
  //c2->cd(3);
  //fAllGenRes->Draw("colz");
  //c2->cd(4);
  //fRecoGenRes->Draw("colz");
  
  plotAcceptance();
  c2->cd(1);
  fAcceptance->Draw("colz");
  c2->cd(2);
  fAnaEff->Draw("colz");
  c2->cd(3);
  fMuIDEff->Draw("colz");  
  c2->cd(4);
  fTrigEff->Draw("colz");   
  

  if ( mode == 1  ){
    double binErr(0);
    double bin(0.); double bin_ratio(-1); double yield(0);
    
    
    for ( int iy = 1; iy <= fS1Yield->GetNbinsX(); ++iy ){
      for ( int ipt = 1; ipt <= fS1Yield->GetNbinsY(); ++ipt ){
	
	bin = fS1Yield->GetBinContent(iy,ipt) / ( fMuIDEff->GetBinContent(iy,ipt) * fTrigEff->GetBinContent(iy,ipt) * fAcceptance->GetBinContent(iy,ipt) * fAnaEff->GetBinContent(iy,ipt) * fPreSelEff->GetBinContent(iy,ipt) );
	fS1YieldCorrected->SetBinContent(iy,ipt,bin);
	fS1YieldCorrected->SetBinError(iy,ipt,(bin/fS1Yield->GetBinContent(iy,ipt))*fS1Yield->GetBinError(iy,ipt) );
	cout << "fAnaEff->GetBinContent(iy,ipt) = "<< fAnaEff->GetBinContent(iy,ipt) << endl;
	cout << "fMuIDEff->GetBinContent(iy,ipt) = "<< fMuIDEff->GetBinContent(iy,ipt) << endl;
	cout << "fTrigEff->GetBinContent(iy,ipt) = "<< fTrigEff->GetBinContent(iy,ipt) << endl;
	cout << "fAllGenRes->GetBinContent(iy,ipt) = "<< fAllGenRes->GetBinContent(iy,ipt) <<endl;
	cout << "fRecoGenRes->GetBinContent(iy,ipt) = "<< fRecoGenRes->GetBinContent(iy,ipt) <<endl;
	cout << "fAcceptance->GetBinContent(iy,ipt) = "<< fAcceptance->GetBinContent(iy,ipt) << endl;
	cout << "fPreSelEff->GetBinContent(iy,ipt) = "<< fPreSelEff->GetBinContent(iy,ipt) << endl;
	cout << "fS1Yield->GetBinContent(iy,ipt) = "<< fS1Yield->GetBinContent(iy,ipt) <<endl;
	cout << "fS1YieldCorrected->GetBinContent(iy,ipt) = "<< fS1YieldCorrected->GetBinContent(iy,ipt) <<endl;
	bin_ratio = (fS1YieldCorrected->GetBinContent(iy,ipt) - fAllGenRes->GetBinContent(iy,ipt))/ fAllGenRes->GetBinContent(iy,ipt);
	cout << "bin_ratio = "  << bin_ratio << endl;
	fS1YieldComparison->SetBinContent(iy,ipt,bin_ratio);
	
      }
    }
    
    double Normalization_2S(1.66);

    for ( int iy = 1; iy <= fS2Yield->GetNbinsX(); ++iy ){
      for ( int ipt = 1; ipt <= fS2Yield->GetNbinsY(); ++ipt ){
	
	bin = fS2Yield->GetBinContent(iy,ipt) / ( fMuIDEff_2->GetBinContent(iy,ipt) * fTrigEff_2->GetBinContent(iy,ipt) * fAcceptance_2S->GetBinContent(iy,ipt) * fAnaEff_2S->GetBinContent(iy,ipt) * fPreSelEff_2S->GetBinContent(iy,ipt) );
	fS2YieldCorrected->SetBinContent(iy,ipt,bin);
	fS2YieldCorrected->SetBinError(iy,ipt,(bin/fS2Yield->GetBinContent(iy,ipt))*fS2Yield->GetBinError(iy,ipt) );
	cout << "fAnaEff_2S->GetBinContent(iy,ipt) = "<< fAnaEff_2S->GetBinContent(iy,ipt) << endl;
	cout << "fMuIDEff_2S->GetBinContent(iy,ipt) = "<< fMuIDEff_2->GetBinContent(iy,ipt) << endl;
	cout << "fTrigEff_2S->GetBinContent(iy,ipt) = "<< fTrigEff_2->GetBinContent(iy,ipt) << endl;
	cout << "fAllGenRes_2S->GetBinContent(iy,ipt) = "<< fAllGenRes_2S->GetBinContent(iy,ipt) <<endl;
	cout << "fRecoGenRes_2S->GetBinContent(iy,ipt) = "<< fRecoGenRes_2S->GetBinContent(iy,ipt) <<endl;
	cout << "fAcceptance_2S->GetBinContent(iy,ipt) = "<< fAcceptance_2S->GetBinContent(iy,ipt) << endl;
	cout << "fPreSelEff_2S->GetBinContent(iy,ipt) = "<< fPreSelEff_2S->GetBinContent(iy,ipt) << endl;
	cout << "fS2Yield->GetBinContent(iy,ipt) = "<< fS2Yield->GetBinContent(iy,ipt) <<endl;
	cout << "fS2YieldCorrected->GetBinContent(iy,ipt) = "<< fS2YieldCorrected->GetBinContent(iy,ipt) <<endl;
	bin_ratio = ((fS2YieldCorrected->GetBinContent(iy,ipt) * lM[1] ) - fAllGenRes_2S->GetBinContent(iy,ipt)) / (fAllGenRes_2S->GetBinContent(iy,ipt));
	cout << "bin_ratio = "  << bin_ratio << endl;
	cout << "Normalization_2S = " << lM[1]  << endl;
	fS2YieldComparison->SetBinContent(iy,ipt,bin_ratio);
	
      }
    }
    
    double Normalization_3S(3.43);

    for ( int iy = 1; iy <= fS3Yield->GetNbinsX(); ++iy ){
      for ( int ipt = 1; ipt <= fS3Yield->GetNbinsY(); ++ipt ){
	
	bin = fS3Yield->GetBinContent(iy,ipt) / ( fMuIDEff_3->GetBinContent(iy,ipt) * fTrigEff_3->GetBinContent(iy,ipt) * fAcceptance_3S->GetBinContent(iy,ipt) * fAnaEff_3S->GetBinContent(iy,ipt) * fPreSelEff_3S->GetBinContent(iy,ipt) );
	fS3YieldCorrected->SetBinContent(iy,ipt,bin);
	fS3YieldCorrected->SetBinError(iy,ipt,(bin/fS3Yield->GetBinContent(iy,ipt))*fS3Yield->GetBinError(iy,ipt) );
	cout << "fAnaEff_2S->GetBinContent(iy,ipt) = "<< fAnaEff_3S->GetBinContent(iy,ipt) << endl;
	cout << "fMuIDEff_2S->GetBinContent(iy,ipt) = "<< fMuIDEff_3->GetBinContent(iy,ipt) << endl;
	cout << "fTrigEff_2S->GetBinContent(iy,ipt) = "<< fTrigEff_3->GetBinContent(iy,ipt) << endl;
	cout << "fAllGenRes_2S->GetBinContent(iy,ipt) = "<< fAllGenRes_3S->GetBinContent(iy,ipt) <<endl;
	cout << "fRecoGenRes_2S->GetBinContent(iy,ipt) = "<< fRecoGenRes_3S->GetBinContent(iy,ipt) <<endl;
	cout << "fAcceptance_2S->GetBinContent(iy,ipt) = "<< fAcceptance_3S->GetBinContent(iy,ipt) << endl;
	cout << "fPreSelEff_2S->GetBinContent(iy,ipt) = "<< fPreSelEff_3S->GetBinContent(iy,ipt) << endl;
	cout << "fS2Yield->GetBinContent(iy,ipt) = "<< fS3Yield->GetBinContent(iy,ipt) <<endl;
	cout << "fS2YieldCorrected->GetBinContent(iy,ipt) = "<< fS3YieldCorrected->GetBinContent(iy,ipt) <<endl;
	bin_ratio = ((fS3YieldCorrected->GetBinContent(iy,ipt) * lM[2] ) - fAllGenRes_3S->GetBinContent(iy,ipt)) / (fAllGenRes_3S->GetBinContent(iy,ipt));
	cout << "bin_ratio = "  << bin_ratio << endl;
	cout << "Normalization_3S = " << lM[2]  << endl;
	fS3YieldComparison->SetBinContent(iy,ipt,bin_ratio);
	
      }
    }
    
    makeCanvas(1);
    c1->Divide(3,1);
    c1->cd(1);
    fS1YieldComparison->Draw("colz");
    c1->cd(2);
    fS2YieldComparison->Draw("colz");
    c1->cd(3);
    fS3YieldComparison->Draw("colz");
    
  }
  
  if ( mode == 2 ){
    double binErr(0);
    double bin(0.); double bin_ratio(-1);
    
    for ( int iy = 1; iy <= fS1Yield->GetNbinsX(); ++iy ){
      for ( int ipt = 1; ipt <= fS1Yield->GetNbinsY(); ++ipt ){
	
	bin = fS1Yield->GetBinContent(iy,ipt) / ( fMuIDEff->GetBinContent(iy,ipt) * fTrigEff->GetBinContent(iy,ipt) * fAcceptance->GetBinContent(iy,ipt) * fAnaEff->GetBinContent(iy,ipt) * fPreSelEff->GetBinContent(iy,ipt)  );
	fS1YieldCorrected->SetBinContent(iy,ipt,bin);
	fS1YieldCorrected->SetBinError(iy,ipt,(bin/fS1Yield->GetBinContent(iy,ipt))*fS1Yield->GetBinError(iy,ipt) );
	cout << "fAnaEff->GetBinContent(iy,ipt) = "<< fAnaEff->GetBinContent(iy,ipt) << endl;
	cout << "fMuIDEff->GetBinContent(iy,ipt) = "<< fMuIDEff->GetBinContent(iy,ipt) << endl;
	cout << "fTrigEff->GetBinContent(iy,ipt) = "<< fTrigEff->GetBinContent(iy,ipt) << endl;
	cout << "fAcceptance->GetBinContent(iy,ipt) = "<< fAcceptance->GetBinContent(iy,ipt) << endl;
	cout << "fPreSelEff->GetBinContent(iy,ipt) = "<< fPreSelEff->GetBinContent(iy,ipt) << endl;
	cout << "fS1Yield->GetBinContent(iy,ipt) = "<< fS1Yield->GetBinContent(iy,ipt) <<endl;
	cout << "fS1YieldCoreected->GetBinContent(iy,ipt) = "<< fS1YieldCorrected->GetBinContent(iy,ipt) <<endl;
      }
    }
    
    for ( int iy = 1; iy <= fS2Yield->GetNbinsX(); ++iy ){
      for ( int ipt = 1; ipt <= fS2Yield->GetNbinsY(); ++ipt ){
	
	bin = fS2Yield->GetBinContent(iy,ipt) / ( fMuIDEff->GetBinContent(iy,ipt) * fTrigEff->GetBinContent(iy,ipt) * fAcceptance_2S->GetBinContent(iy,ipt) * fAnaEff_2S->GetBinContent(iy,ipt) * fPreSelEff_2S->GetBinContent(iy,ipt)  );
	fS2YieldCorrected->SetBinContent(iy,ipt,bin);
	fS2YieldCorrected->SetBinError(iy,ipt,(bin/fS2Yield->GetBinContent(iy,ipt))*fS2Yield->GetBinError(iy,ipt) );
	cout << "fAnaEff_2S->GetBinContent(iy,ipt) = "<< fAnaEff_2S->GetBinContent(iy,ipt) << endl;
	cout << "fMuIDEff->GetBinContent(iy,ipt) = "<< fMuIDEff->GetBinContent(iy,ipt) << endl;
	cout << "fTrigEff->GetBinContent(iy,ipt) = "<< fTrigEff->GetBinContent(iy,ipt) << endl;
	cout << "fAcceptance_2S->GetBinContent(iy,ipt) = "<< fAcceptance_2S->GetBinContent(iy,ipt) << endl;
	cout << "fPreSelEff_2S->GetBinContent(iy,ipt) = "<< fPreSelEff_2S->GetBinContent(iy,ipt) << endl;
	cout << "fS2Yield->GetBinContent(iy,ipt) = "<< fS2Yield->GetBinContent(iy,ipt) <<endl;
	cout << "fS2YieldCoreected->GetBinContent(iy,ipt) = "<< fS2YieldCorrected->GetBinContent(iy,ipt) <<endl;
      }
    }
    
    for ( int iy = 1; iy <= fS3Yield->GetNbinsX(); ++iy ){
      for ( int ipt = 1; ipt <= fS3Yield->GetNbinsY(); ++ipt ){
	
	bin = fS3Yield->GetBinContent(iy,ipt) / ( fMuIDEff->GetBinContent(iy,ipt) * fTrigEff->GetBinContent(iy,ipt) * fAcceptance_3S->GetBinContent(iy,ipt) * fAnaEff_3S->GetBinContent(iy,ipt) * fPreSelEff_3S->GetBinContent(iy,ipt)  );
	fS3YieldCorrected->SetBinContent(iy,ipt,bin);
	fS3YieldCorrected->SetBinError(iy,ipt,(bin/fS3Yield->GetBinContent(iy,ipt))*fS3Yield->GetBinError(iy,ipt) );
	cout << "fAnaEff_3S->GetBinContent(iy,ipt) = "<< fAnaEff_3S->GetBinContent(iy,ipt) << endl;
	cout << "fMuIDEff->GetBinContent(iy,ipt) = "<< fMuIDEff->GetBinContent(iy,ipt) << endl;
	cout << "fTrigEff->GetBinContent(iy,ipt) = "<< fTrigEff->GetBinContent(iy,ipt) << endl;
	cout << "fAcceptance_3S->GetBinContent(iy,ipt) = "<< fAcceptance_3S->GetBinContent(iy,ipt) << endl;
	cout << "fPreSelEff_3S->GetBinContent(iy,ipt) = "<< fPreSelEff_3S->GetBinContent(iy,ipt) << endl;
	cout << "fS3Yield->GetBinContent(iy,ipt) = "<< fS3Yield->GetBinContent(iy,ipt) <<endl;
	cout << "fS3YieldCoreected->GetBinContent(iy,ipt) = "<< fS3YieldCorrected->GetBinContent(iy,ipt) <<endl;
      }
    }
   
    
    makeCanvas(1);
    c1->Divide(3,1);
    c1->cd(1);
    fS1YieldCorrected->Draw("colz");
    c1->cd(2);
    fS2YieldCorrected->Draw("colz");
    c1->cd(3);
    fS3YieldCorrected->Draw("colz");
    
  }
  
}
// ----------------------------------------------------------------------
void anaXS::fillPidTables() {

  // -- MM/MT
  fPtTnpNeg->readFromHist(gDirectory, "fS1YieldNeg", "fS2YieldNeg");
  fPtTnpNeg->setHistName("TNP negative muons");                     
  fPtTnpNeg->dumpToFile((fPtDirectory+string("/PtTnpNeg-") + fSample + string(".dat")).c_str());

  fPtTnpPos->readFromHist(gDirectory, "fS1YieldPos", "fS2YieldPos");
  fPtTnpPos->setHistName("TNP positive muons");                     
  fPtTnpPos->dumpToFile((fPtDirectory+string("/PtTnpPos-") + fSample + string(".dat")).c_str());

  // -- MC all 
  fPtMctNeg->readFromHist(gDirectory, "fS1MctNeg", "fS2MctNeg");
  fPtMctNeg->setHistName("MC truth negative muons");                     
  fPtMctNeg->dumpToFile((fPtDirectory+string("/PtMctNeg-") + fSample + string(".dat")).c_str());

  fPtMctPos->readFromHist(gDirectory, "fS1MctPos", "fS2MctPos");
  fPtMctPos->setHistName("MC truth positive muons");                     
  fPtMctPos->dumpToFile((fPtDirectory+string("/PtMctPos-") + fSample + string(".dat")).c_str());

  // -- MC probe
  fPtMcpNeg->readFromHist(gDirectory, "fS1McpNeg", "fS2McpNeg");
  fPtMcpNeg->setHistName("MC probe negative muons");                     
  fPtMcpNeg->dumpToFile((fPtDirectory+string("/PtMcpNeg-") + fSample + string(".dat")).c_str());

  fPtMcpPos->readFromHist(gDirectory, "fS1McpPos", "fS2McpPos");
  fPtMcpPos->setHistName("MC probe positive muons");                     
  fPtMcpPos->dumpToFile((fPtDirectory+string("/PtMcpPos-") + fSample + string(".dat")).c_str());

  // -- MM/(MM + MMbar)
  fPtMmbNeg->readFromHist(gDirectory, "fS1YieldNeg", "fS3YieldNeg");
  fPtMmbNeg->setHistName("TNP2 negative muons");                     
  fPtMmbNeg->dumpToFile((fPtDirectory+string("/PtMmbNeg-") + fSample + string(".dat")).c_str());

  fPtMmbPos->readFromHist(gDirectory, "fS1YieldPos", "fS3YieldPos");
  fPtMmbPos->setHistName("TNP2 positive muons");                     
  fPtMmbPos->dumpToFile((fPtDirectory+string("/PtMmbPos-") + fSample + string(".dat")).c_str());


}


// ----------------------------------------------------------------------
void anaXS::readPidTables(const char *sample) {
  fSample = string(sample); 
  string tname = fPtDirectory + string("/PtTnpNeg-") + fSample + string(".dat"); 
  fPtTnpNeg->readFromFile(tname.c_str());
  fPtTnpNeg->setHistName("TNP negative muons");                     

  tname = fPtDirectory + string("/PtTnpPos-") + fSample + string(".dat"); 
  fPtTnpPos->readFromFile(tname.c_str());
  fPtTnpPos->setHistName("TNP positive muons");                     
  
  tname = fPtDirectory + string("/PtMcpNeg-") + fSample + string(".dat"); 
  fPtMcpNeg->readFromFile(tname.c_str());
  fPtMcpNeg->setHistName("MC probe negative muons");                     

  tname = fPtDirectory + string("/PtMcpPos-") + fSample + string(".dat"); 
  fPtMcpPos->readFromFile(tname.c_str());
  fPtMcpPos->setHistName("MC probe positive muons");                     
  
  tname = fPtDirectory + string("/PtMctNeg-") + fSample + string(".dat"); 
  fPtMctNeg->readFromFile(tname.c_str());
  fPtMctNeg->setHistName("MC truth negative muons");                     

  tname = fPtDirectory + string("/PtMctPos-") + fSample + string(".dat"); 
  fPtMctPos->readFromFile(tname.c_str());
  fPtMctPos->setHistName("MC truth positive muons");                     

  tname = fPtDirectory + string("/PtMmbNeg-") + fSample + string(".dat"); 
  fPtMmbNeg->readFromFile(tname.c_str());
  fPtMmbNeg->setHistName("TNP 2 negative muons");                     

  tname = fPtDirectory + string("/PtMmbPos-") + fSample + string(".dat"); 
  fPtMmbPos->readFromFile(tname.c_str());
  fPtMmbPos->setHistName("TNP 2 positive muons");                     

  
}


// ----------------------------------------------------------------------
void anaXS::allDifferences(int sample) {
  // sample = 0x1 Upsilon
  //          0x2 J/psi
  
  double MIN(-0.1), MAX(0.1);    
 
  // -- mmbar vs mm/mt
  if (sample & 0x1) {
    //MIN = -0.05; MAX = 0.05;	
    ptDifference(Form("%s/PtTnpPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMmbPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 MIN, MAX, "upsilon-tnp-mmb-pos.eps");

    ptDifference(Form("%s/PtTnpNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMmbNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 MIN, MAX, "upsilon-tnp-mmb-neg.eps");
  }
  
  if (sample & 0x2) {
    ptDifference(Form("%s/PtTnpPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 Form("%s/PtMmbPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "jpsi-tnp-mmb-pos.eps");
    
    ptDifference(Form("%s/PtTnpNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 Form("%s/PtMmbNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "jpsi-tnp-mmb-neg.eps");
  }

  // -- fit bias
  if (sample & 0x1) {
    //MIN = -0.05; MAX = 0.05;  
    ptDifference(Form("%s/PtTnpPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMcpPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 MIN, MAX, "upsilon-tnp-mcp-pos.eps");
  
    ptDifference(Form("%s/PtTnpNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMcpNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 MIN, MAX, "upsilon-tnp-mcp-neg.eps");
  }  

  if (sample & 0x2) {
    ptDifference(Form("%s/PtTnpPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 Form("%s/PtMcpPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "jpsi-tnp-mcp-pos.eps");
    
    ptDifference(Form("%s/PtTnpNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 Form("%s/PtMcpNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "jpsi-tnp-mcp-neg.eps");
  }

  // -- selection bias
  if (sample & 0x1) {
    //MIN = -0.05; MAX = 0.05;
    ptDifference(Form("%s/PtMctPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMcpPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 MIN, MAX, "upsilon-mct-mcp-pos.eps");
    
    ptDifference(Form("%s/PtMctNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMcpNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 MIN, MAX, "upsilon-mct-mcp-neg.eps");
  }


  if (sample & 0x2) {
    ptDifference(Form("%s/PtMctPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 Form("%s/PtMcpPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "jpsi-mct-mcp-pos.eps");
    
    ptDifference(Form("%s/PtMctNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 Form("%s/PtMcpNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "jpsi-mct-mcp-neg.eps");
  }

  // -- The final difference
  if (sample & 0x1) {
    //MIN = -0.05; MAX = 0.05;
    ptDifference(Form("%s/PtTnpPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMctPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 MIN, MAX, "upsilon-tnp-mct-pos.eps");
    
    ptDifference(Form("%s/PtTnpNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMctNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 MIN, MAX, "upsilon-tnp-mct-neg.eps");
  }

  if (sample & 0x2) {
    ptDifference(Form("%s/PtTnpPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 Form("%s/PtMctPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "jpsi-tnp-mct-pos.eps");
    
    ptDifference(Form("%s/PtTnpNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 Form("%s/PtMctNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "jpsi-tnp-mct-neg.eps");
  }

  // -- NEG vs POS
  if (sample & 0x1) {
   // MIN = -0.05; MAX = 0.05;
    ptDifference(Form("%s/PtTnpNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtTnpPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 MIN, MAX, "upsilon-neg-pos.eps");
  }

  if (sample & 0x2) {
    ptDifference(Form("%s/PtTnpNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 Form("%s/PtTnpPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "jpsi-neg-pos.eps");
  }

  // -- Upsilon vs. J/psi
  if (sample & 0x2 && sample & 0x1) {
    ptDifference(Form("%s/PtTnpPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtTnpPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "upsilon-jpsi-tnp-pos.eps");
    
    ptDifference(Form("%s/PtTnpNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtTnpNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "upsilon-jpsi-tnp-neg.eps");
    
    ptDifference(Form("%s/PtMcpPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMcpPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "upsilon-jpsi-mcp-pos.eps");
    
    ptDifference(Form("%s/PtMcpNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMcpNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "upsilon-jpsi-mcp-neg.eps");
    
    ptDifference(Form("%s/PtMctPos-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMctPos-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "upsilon-jpsi-mct-pos.eps");
    
    ptDifference(Form("%s/PtMctNeg-%s.dat", fPtDirectory.c_str(), "upsilon"), 
		 Form("%s/PtMctNeg-%s.dat", fPtDirectory.c_str(), "jpsi"), 
		 MIN, MAX, "upsilon-jpsi-mct-neg.eps");
  }

}

// ----------------------------------------------------------------------
void anaXS::ptDifference(const char *a, const char *b, double MIN, double MAX, const char *fname) {

  
  PidTable td(a); 
  PidTable tb(b);

  td.setHistMinMax(0.8, 1.1);
  TH1D *hd = new TH1D("Hd", "", 40, MIN, MAX); 
  setTitles(hd, "difference", "# bins", 0.08, 0.9, 0.8, 0.08);
  setFilledHist(hd, kBlack, kYellow, 1000, 2); 

  td.subtract(tb);
  TH2D *h1 = td.get2dHist("H1", ""); setTitles(h1, "#eta", "p_{T} [GeV]", 0.08, 0.9, 0.8, 0.08); 
  h1->SetMinimum(MIN);
  h1->SetMaximum(MAX);
  
  td.eff2d(h1); 
  td.effHist(hd);

  gStyle->SetOptStat(0); 
  makeCanvas(1); 
  c1->Clear(); 
  c1->Divide(2,1); 

  c1->cd(1);   shrinkPad(0.15, 0.15, 0.15); 
  h1->Draw("colz");

  c1->cd(2);   shrinkPad(0.15, 0.15); 
  hd->Draw();
  tl->SetTextSize(0.08); 
  tl->DrawLatex(0.12, 0.92, Form("mean: %4.3f", hd->GetMean())); 
  tl->DrawLatex(0.50, 0.92, Form("RMS:  %4.3f", hd->GetRMS())); 
  c1->SaveAs(Form("%s/bias-%s", fPtDirectory.c_str(), fname)); 
   
  
  
  if ( fname == "jpsi-mct-mcp-pos.eps" ){
  	TFile file("BiasPos.root","recreate");
  	TH2D *BiasPos = (TH2D*)h1->Clone();
	BiasPos->SetName("BiasPos");
	BiasPos->Write();
  	//file.Write();
	file.Close();
  }
  
  if ( fname == "jpsi-mct-mcp-neg.eps" ){
  	TFile file1("BiasNeg.root","recreate");
  	TH2D *BiasNeg = (TH2D*)h1->Clone();
	BiasNeg->SetName("BiasNeg");
	BiasNeg->Write();
  	//file.Write();
	file1.Close();
  }  
  
}

void anaXS::PlotProjections(int mode) {

  
  if ( mode == 1 ){
    double bin_contentAll(0); double bin_contentReco(0); double bin_ratioA(0);
    double bin_contentYield(0); double bin_ratio(0);
    double bin_contentYieldErr(0); double yieldErr(0);
        
    for (int j = 1; j <= fAcceptance->GetNbinsY(); ++j){
      for (int i = 1; i <= fAcceptance->GetNbinsX(); ++i) {
    	bin_contentAll += fAllGenRes->GetCellContent(i,j);
	bin_contentReco += fRecoGenRes->GetCellContent(i,j);
	bin_contentYield += fS1YieldCorrected->GetCellContent(i,j);
	bin_contentYieldErr += fS1YieldCorrected->GetCellError(i,j)*fS1YieldCorrected->GetCellError(i,j);
       }
      bin_ratio = bin_contentReco/bin_contentAll;
      //cout << "bin_ratio = "  << bin_ratio << endl;
      fAcceptanceProjPt->SetBinContent(j,bin_ratio);
      cout << bin_contentYield << endl;
      fS1YieldPt->SetBinContent(j,bin_contentYield/fS1YieldPt->GetBinWidth(j));
      yieldErr = TMath::Sqrt(bin_contentYieldErr/fS1YieldPt->GetBinWidth(j));
      cout << yieldErr  << endl;
      fS1YieldPt->SetBinError(j,yieldErr);
      fAllGenResPt->SetBinContent(j,bin_contentAll/fAllGenResPt->GetBinWidth(j));
      bin_ratio=0;bin_contentAll=0;bin_contentReco=0;
      bin_contentYield=0;
    }
    
    for (int j = 1; j <= fAcceptance_2S->GetNbinsY(); ++j){
      for (int i = 1; i <= fAcceptance_2S->GetNbinsX(); ++i) {
    	bin_contentAll += fAllGenRes_2S->GetCellContent(i,j);
	bin_contentReco += fRecoGenRes_2S->GetCellContent(i,j);
	bin_contentYield += fS2YieldCorrected->GetCellContent(i,j);
	bin_contentYieldErr += fS2YieldCorrected->GetCellError(i,j)*fS2YieldCorrected->GetCellError(i,j);
       }
      bin_ratio = bin_contentReco/bin_contentAll;
      //cout << "bin_ratio = "  << bin_ratio << endl;
      fAcceptanceProjPt_2S->SetBinContent(j,bin_ratio);
      cout << bin_contentYield << endl;
      bin_contentYield *= lM[1];
      fS2YieldPt->SetBinContent(j,bin_contentYield/fS2YieldPt->GetBinWidth(j));
      yieldErr = TMath::Sqrt(bin_contentYieldErr/fS2YieldPt->GetBinWidth(j));
      cout << yieldErr  << endl;
      fS2YieldPt->SetBinError(j,yieldErr);
      fAllGenResPt_2S->SetBinContent(j,bin_contentAll/fAllGenResPt_2S->GetBinWidth(j));
      bin_ratio=0;bin_contentAll=0;bin_contentReco=0;
      bin_contentYield=0;
    }
    
    for (int j = 1; j <= fAcceptance_3S->GetNbinsY(); ++j){
      for (int i = 1; i <= fAcceptance_3S->GetNbinsX(); ++i) {
    	bin_contentAll += fAllGenRes_3S->GetCellContent(i,j);
	bin_contentReco += fRecoGenRes_3S->GetCellContent(i,j);
	bin_contentYield += fS3YieldCorrected->GetCellContent(i,j);
	bin_contentYieldErr += fS3YieldCorrected->GetCellError(i,j)*fS3YieldCorrected->GetCellError(i,j);
       }
      bin_ratio = bin_contentReco/bin_contentAll;
      //cout << "bin_ratio = "  << bin_ratio << endl;
      fAcceptanceProjPt_3S->SetBinContent(j,bin_ratio);
      cout << bin_contentYield << endl;
      bin_contentYield *= lM[2];
      fS3YieldPt->SetBinContent(j,bin_contentYield/fS3YieldPt->GetBinWidth(j));
      yieldErr = TMath::Sqrt(bin_contentYieldErr/fS3YieldPt->GetBinWidth(j));
      cout << yieldErr  << endl;
      fS3YieldPt->SetBinError(j,yieldErr);
      fAllGenResPt_3S->SetBinContent(j,bin_contentAll/fAllGenResPt_3S->GetBinWidth(j));
      bin_ratio=0;bin_contentAll=0;bin_contentReco=0;
      bin_contentYield=0;
    }
    makeCanvas(1);
    fAcceptanceProjPt->GetXaxis()->SetTitle("P_{T}");
    fAcceptanceProjPt->SetMinimum(0.2);
    fAcceptanceProjPt->SetMaximum(1.05);
    fAcceptanceProjPt->SetMarkerStyle(22);
    fAcceptanceProjPt->SetMarkerColor(2);
    fAcceptanceProjPt->SetLineColor(2);
    fAcceptanceProjPt->Draw("p");
    fAcceptanceProjPt_2S->SetMarkerStyle(23);
    fAcceptanceProjPt_2S->SetMarkerColor(3);
    fAcceptanceProjPt_2S->SetLineColor(3);
    fAcceptanceProjPt_2S->Draw("psame");
    fAcceptanceProjPt_3S->SetMarkerStyle(24);
    fAcceptanceProjPt_3S->SetMarkerColor(4);
    fAcceptanceProjPt_3S->SetLineColor(4);
    fAcceptanceProjPt_3S->Draw("psame");
    legg = new TLegend(0.6,0.6,0.8,0.8);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.06); legg->SetTextFont(132); 
    legg->SetHeader("Acceptance");
    legge = legg->AddEntry(fAcceptanceProjPt,  "Y(1S) ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(fAcceptanceProjPt_2S,  "Y(2S) ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(fAcceptanceProjPt_3S,  "Y(3S) ","p"); legge->SetTextColor(kBlack);
    legg->Draw();
    makeCanvas(1);
    c1->Divide(3,1);
    c1->cd(1);
    fS1YieldPt->SetTitle("Yield Comparison Ups(1S)");
    fS1YieldPt->GetXaxis()->SetTitle("P_{T}");
    fS1YieldPt->SetMinimum(0.);
    fS1YieldPt->SetMaximum(270000.);
    fS1YieldPt->SetMarkerStyle(21);
    fS1YieldPt->SetMarkerColor(3);
    fS1YieldPt->SetLineColor(3);
    fAllGenResPt->SetMarkerStyle(20);
    fAllGenResPt->SetMarkerColor(4);
    fAllGenResPt->SetLineColor(4);
    fS1YieldPt->Draw("p");
    fAllGenResPt->Draw("psame");
    legg = new TLegend(0.6,0.6,0.8,0.8);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Yield Comparison Ups(1S)");
    legge = legg->AddEntry(fS1YieldPt,  "Reconstructed Yield ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(fAllGenResPt,  "True Yield","p"); legge->SetTextColor(kBlack);
    legg->Draw();
    
    c1->cd(2);
    fS2YieldPt->SetTitle("Yield Comparison Ups(2S)");
    fS2YieldPt->GetXaxis()->SetTitle("P_{T}");
    fS2YieldPt->SetMinimum(0.);
    fS2YieldPt->SetMaximum(120000.);
    fS2YieldPt->SetMarkerStyle(21);
    fS2YieldPt->SetMarkerColor(3);
    fS2YieldPt->SetLineColor(3);
    fAllGenResPt_2S->SetMarkerStyle(20);
    fAllGenResPt_2S->SetMarkerColor(4);
    fAllGenResPt_2S->SetLineColor(4);
    fS2YieldPt->Draw("p");
    fAllGenResPt_2S->Draw("psame");
    legg = new TLegend(0.6,0.6,0.8,0.8);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Yield Comparison Ups(2S)");
    legge = legg->AddEntry(fS2YieldPt,  "Reconstructed Yield ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(fAllGenResPt_2S,  "True Yield","p"); legge->SetTextColor(kBlack);
    legg->Draw();
    
    c1->cd(3);
    fS3YieldPt->SetTitle("Yield Comparison Ups(3S)");
    fS3YieldPt->GetXaxis()->SetTitle("P_{T}");
    fS3YieldPt->SetMinimum(0.);
    fS3YieldPt->SetMaximum(100000.);
    fS3YieldPt->SetMarkerStyle(21);
    fS3YieldPt->SetMarkerColor(3);
    fS3YieldPt->SetLineColor(3);
    fAllGenResPt_3S->SetMarkerStyle(20);
    fAllGenResPt_3S->SetMarkerColor(4);
    fAllGenResPt_3S->SetLineColor(4);
    fS3YieldPt->Draw("p");
    fAllGenResPt_3S->Draw("psame");
    legg = new TLegend(0.6,0.6,0.8,0.8);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Yield Comparison Ups(3S)");
    legge = legg->AddEntry(fS3YieldPt,  "Reconstructed Yield ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(fAllGenResPt_3S,  "True Yield","p"); legge->SetTextColor(kBlack);
    legg->Draw();
    
    
    bin_contentAll=0; bin_contentYield=0; bin_contentYieldErr=0; yieldErr=0;
    for (int i = 1; i <= fAcceptance->GetNbinsX(); ++i) {
      for (int j = 1; j <= fAcceptance->GetNbinsY(); ++j){	
	
	bin_contentAll += fAllGenRes->GetCellContent(i,j);
	bin_contentYield += fS1YieldCorrected->GetCellContent(i,j);
	bin_contentYieldErr += fS1YieldCorrected->GetCellError(i,j)*fS1YieldCorrected->GetCellError(i,j);
      }
      
      fS1YieldEta->SetBinContent(i,bin_contentYield);
      fAllGenResEta->SetBinContent(i,bin_contentAll);
      yieldErr = TMath::Sqrt(bin_contentYieldErr);
      fS1YieldEta->SetBinError(i,yieldErr);
      bin_contentAll=0;
      bin_contentYield=0;
      bin_contentYieldErr=0;
    }
    
    for (int i = 1; i <= fAcceptance_2S->GetNbinsX(); ++i) {
      for (int j = 1; j <= fAcceptance_2S->GetNbinsY(); ++j){	
	
	bin_contentAll += fAllGenRes_2S->GetCellContent(i,j);
	bin_contentYield += fS2YieldCorrected->GetCellContent(i,j);
	bin_contentYieldErr += fS2YieldCorrected->GetCellError(i,j)*fS2YieldCorrected->GetCellError(i,j);
      }
      bin_contentYield *= lM[1];
      fS2YieldEta->SetBinContent(i,bin_contentYield);
      fAllGenResEta_2S->SetBinContent(i,bin_contentAll);
      yieldErr = TMath::Sqrt(bin_contentYieldErr);
      fS2YieldEta->SetBinError(i,yieldErr);
      bin_contentAll=0;
      bin_contentYield=0;
      bin_contentYieldErr=0;
    }
    
    for (int i = 1; i <= fAcceptance_3S->GetNbinsX(); ++i) {
      for (int j = 1; j <= fAcceptance_3S->GetNbinsY(); ++j){	
	
	bin_contentAll += fAllGenRes_3S->GetCellContent(i,j);
	bin_contentYield += fS3YieldCorrected->GetCellContent(i,j);
	bin_contentYieldErr += fS3YieldCorrected->GetCellError(i,j)*fS3YieldCorrected->GetCellError(i,j);
      }
      bin_contentYield *= lM[2];
      fS3YieldEta->SetBinContent(i,bin_contentYield);
      fAllGenResEta_3S->SetBinContent(i,bin_contentAll);
      yieldErr = TMath::Sqrt(bin_contentYieldErr);
      fS3YieldEta->SetBinError(i,yieldErr);
      bin_contentAll=0;
      bin_contentYield=0;
      bin_contentYieldErr=0;
    }
     
    //makeCanvas(1);
    //c1->Divide(3,1);
    //c1->cd(1);
    fS1YieldEta->SetTitle("Yield Comparison Ups(1S)");
    fS1YieldEta->GetXaxis()->SetTitle("Rapidity");
    fS1YieldEta->SetMinimum(300000.);
    fS1YieldEta->SetMaximum(650000.);
    fS1YieldEta->SetMarkerStyle(21);
    fS1YieldEta->SetMarkerColor(3);
    fS1YieldEta->SetLineColor(3);
    fAllGenResEta->SetMarkerStyle(20);
    fAllGenResEta->SetMarkerColor(4);
    fAllGenResEta->SetLineColor(4);
    //fS1YieldEta->Draw("p");
    //fAllGenResEta->Draw("psame");
    legg = new TLegend(0.6,0.7,0.8,0.9);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Yield Comparison");
    legge = legg->AddEntry(fS1YieldEta,  "Reconstructed Yield  ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(fAllGenResEta,  "True Yield","p"); legge->SetTextColor(kBlack);
    //legg->Draw();
    //c1->cd(2);
    fS2YieldEta->SetTitle("Yield Comparison Ups(2S)");
    fS2YieldEta->GetXaxis()->SetTitle("Rapidity");
    fS2YieldEta->SetMinimum(100000.);
    fS2YieldEta->SetMaximum(350000.);
    fS2YieldEta->SetMarkerStyle(21);
    fS2YieldEta->SetMarkerColor(3);
    fS2YieldEta->SetLineColor(3);
    fAllGenResEta_2S->SetMarkerStyle(20);
    fAllGenResEta_2S->SetMarkerColor(4);
    fAllGenResEta_2S->SetLineColor(4);
    //fS2YieldEta->Draw("p");
    //fAllGenResEta_2S->Draw("psame");
    legg = new TLegend(0.6,0.7,0.8,0.9);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Yield Comparison");
    legge = legg->AddEntry(fS2YieldEta,  "Reconstructed Yield  ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(fAllGenResEta_2S,  "True Yield","p"); legge->SetTextColor(kBlack);
    //legg->Draw();
    //c1->cd(3);
    fS3YieldEta->SetTitle("Yield Comparison Ups(3S)");
    fS3YieldEta->GetXaxis()->SetTitle("Rapidity");
    fS3YieldEta->SetMinimum(500000.);
    fS3YieldEta->SetMaximum(250000.);
    fS3YieldEta->SetMarkerStyle(21);
    fS3YieldEta->SetMarkerColor(3);
    fS3YieldEta->SetLineColor(3);
    fAllGenResEta_3S->SetMarkerStyle(20);
    fAllGenResEta_3S->SetMarkerColor(4);
    fAllGenResEta_3S->SetLineColor(4);
    //fS3YieldEta->Draw("p");
    //fAllGenResEta_3S->Draw("psame");
    legg = new TLegend(0.6,0.7,0.8,0.9);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Yield Comparison");
    legge = legg->AddEntry(fS3YieldEta,  "Reconstructed Yield  ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(fAllGenResEta_3S,  "True Yield","p"); legge->SetTextColor(kBlack);
    //legg->Draw();
     
    
    
  }
  
  if ( mode == 2 ){
    TH1D *hICHEP; TH1D *hICHEP_ETA;
    double bin_contentAll(0); double bin_contentReco(0); double bin_ratioA(0);
    double bin_contentYield(0); double bin_ratio(0);
    double lumi(31339);
    double xsection(0);
    double bin_contentYieldErr(0); double xsectionErr(0);
    hICHEP = new TH1D("hICHEP", "hICHEP", 
			  fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			  );
    
    hICHEP_ETA = new TH1D("hICHEP_ETA", "hICHEP_ETA", 
		      fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
		      );    
    
    plotAcceptance();
    for (int j = 1; j <= fAcceptance->GetNbinsY(); ++j){
      for (int i = 1; i <= fAcceptance->GetNbinsX(); ++i) {
	bin_contentAll += fAllGenRes->GetCellContent(i,j);
	bin_contentReco += fRecoGenRes->GetCellContent(i,j);
	bin_contentYield += fS1YieldCorrected->GetCellContent(i,j);
	bin_contentYieldErr += fS1YieldCorrected->GetCellError(i,j)*fS1YieldCorrected->GetCellError(i,j);
      }
      bin_ratio = bin_contentReco/bin_contentAll;
      //cout << "bin_ratio = "  << bin_ratio << endl;
      fAcceptanceProjPt->SetBinContent(j,bin_ratio);
      cout << bin_contentYield << endl;
      xsection = bin_contentYield/lumi;
      fS1YieldPt->SetBinContent(j,xsection/fS1YieldPt->GetBinWidth(j));
      xsectionErr = TMath::Sqrt(bin_contentYieldErr)/lumi;
      fS1YieldPt->SetBinError(j,xsectionErr/fS1YieldPt->GetBinWidth(j));
      cout << " Ups(1S) Results "  << endl;
      cout << "xsection = " << fS1YieldPt->GetBinContent(j) << "+/-" << fS1YieldPt->GetBinError(j) << "(%" << (fS1YieldPt->GetBinError(j)/fS1YieldPt->GetBinContent(j))*100  << ")" << endl;      
      bin_ratio=0;bin_contentAll=0;bin_contentReco=0;
      bin_contentYield=0;xsection=0;xsectionErr=0;bin_contentYieldErr=0;
    }
    
    for (int j = 1; j <= fAcceptance_2S->GetNbinsY(); ++j){
      for (int i = 1; i <= fAcceptance_2S->GetNbinsX(); ++i) {
	bin_contentAll += fAllGenRes_2S->GetCellContent(i,j);
	bin_contentReco += fRecoGenRes_2S->GetCellContent(i,j);
	bin_contentYield += fS2YieldCorrected->GetCellContent(i,j);
	bin_contentYieldErr += fS2YieldCorrected->GetCellError(i,j)*fS2YieldCorrected->GetCellError(i,j);
      }
      bin_ratio = bin_contentReco/bin_contentAll;
      //cout << "bin_ratio = "  << bin_ratio << endl;
      fAcceptanceProjPt_2S->SetBinContent(j,bin_ratio);
      cout << bin_contentYield << endl;
      xsection = bin_contentYield/lumi;
      fS2YieldPt->SetBinContent(j,xsection/fS2YieldPt->GetBinWidth(j));
      xsectionErr = TMath::Sqrt(bin_contentYieldErr)/lumi;
      fS2YieldPt->SetBinError(j,xsectionErr/fS2YieldPt->GetBinWidth(j));
      cout << " Ups(2S) Results "  << endl;
      cout << "xsection = " << fS2YieldPt->GetBinContent(j) << "+/-" << fS2YieldPt->GetBinError(j) << "(%" << (fS2YieldPt->GetBinError(j)/fS2YieldPt->GetBinContent(j))*100  << ")" << endl;   
      bin_ratio=0;bin_contentAll=0;bin_contentReco=0;
      bin_contentYield=0;xsection=0;xsectionErr=0;bin_contentYieldErr=0;
    }
    
    for (int j = 1; j <= fAcceptance_3S->GetNbinsY(); ++j){
      for (int i = 1; i <= fAcceptance_3S->GetNbinsX(); ++i) {
	bin_contentAll += fAllGenRes_3S->GetCellContent(i,j);
	bin_contentReco += fRecoGenRes_3S->GetCellContent(i,j);
	bin_contentYield += fS3YieldCorrected->GetCellContent(i,j);
	bin_contentYieldErr += fS3YieldCorrected->GetCellError(i,j)*fS3YieldCorrected->GetCellError(i,j);
      }
      bin_ratio = bin_contentReco/bin_contentAll;
      //cout << "bin_ratio = "  << bin_ratio << endl;
      fAcceptanceProjPt_3S->SetBinContent(j,bin_ratio);
      cout << bin_contentYield << endl;
      xsection = bin_contentYield/lumi;
      fS3YieldPt->SetBinContent(j,xsection/fS3YieldPt->GetBinWidth(j));
      xsectionErr = TMath::Sqrt(bin_contentYieldErr)/lumi;
      fS3YieldPt->SetBinError(j,xsectionErr/fS3YieldPt->GetBinWidth(j));
      cout << " Ups(3S) Results "  << endl;
      cout << "xsection = " << fS3YieldPt->GetBinContent(j) << "+/-" << fS3YieldPt->GetBinError(j) << "(%" << (fS3YieldPt->GetBinError(j)/fS3YieldPt->GetBinContent(j))*100  << ")"<< endl;     
      bin_ratio=0;bin_contentAll=0;bin_contentReco=0;
      bin_contentYield=0;xsection=0;xsectionErr=0;bin_contentYieldErr=0;
    }
     
    
    ////// 1s Bin
    /*    
    hICHEP->SetBinContent(1,0.30); hICHEP->SetBinError(1,0.05);
    hICHEP->SetBinContent(2,0.90); hICHEP->SetBinError(2,0.15);
    hICHEP->SetBinContent(3,1.04); hICHEP->SetBinError(3,0.17); 
    hICHEP->SetBinContent(4,0.88); hICHEP->SetBinError(4,0.14);    
    hICHEP->SetBinContent(5,0.90); hICHEP->SetBinError(5,0.15);  
    hICHEP->SetBinContent(6,0.82); hICHEP->SetBinError(6,0.14);
    hICHEP->SetBinContent(7,0.64); hICHEP->SetBinError(7,0.11);
    hICHEP->SetBinContent(8,0.51); hICHEP->SetBinError(8,0.09);
    hICHEP->SetBinContent(9,0.33); hICHEP->SetBinError(9,0.06); 
    hICHEP->SetBinContent(10,0.25); hICHEP->SetBinError(10,0.05);    
    hICHEP->SetBinContent(11,0.36/hICHEP->GetBinWidth(11)); hICHEP->SetBinError(11,0.06/hICHEP->GetBinWidth(11));  
    hICHEP->SetBinContent(12,0.18/hICHEP->GetBinWidth(12)); hICHEP->SetBinError(12,0.03/hICHEP->GetBinWidth(12));    
    hICHEP->SetBinContent(13,0.14/hICHEP->GetBinWidth(13)); hICHEP->SetBinError(13,0.03/hICHEP->GetBinWidth(13));    
    hICHEP->SetBinContent(14,0.06/hICHEP->GetBinWidth(14)); hICHEP->SetBinError(14,0.01/hICHEP->GetBinWidth(14));  
    hICHEP->SetBinContent(15,0.06/hICHEP->GetBinWidth(15)); hICHEP->SetBinError(15,0.01/hICHEP->GetBinWidth(15));    
    */		  
    
    
    ////// 2S Bin
    
    hICHEP->SetBinContent(1,0.25/hICHEP->GetBinWidth(1)); hICHEP->SetBinError(1,0.05/hICHEP->GetBinWidth(1));
    hICHEP->SetBinContent(2,0.48/hICHEP->GetBinWidth(2)); hICHEP->SetBinError(2,0.09/hICHEP->GetBinWidth(2));
    hICHEP->SetBinContent(3,0.41/hICHEP->GetBinWidth(3)); hICHEP->SetBinError(3,0.08/hICHEP->GetBinWidth(3));
    hICHEP->SetBinContent(4,0.41/hICHEP->GetBinWidth(4)); hICHEP->SetBinError(4,0.08/hICHEP->GetBinWidth(4));
    hICHEP->SetBinContent(5,0.21/hICHEP->GetBinWidth(5)); hICHEP->SetBinError(5,0.04/hICHEP->GetBinWidth(5));
    hICHEP->SetBinContent(6,0.09/hICHEP->GetBinWidth(6)); hICHEP->SetBinError(6,0.02/hICHEP->GetBinWidth(6));
    hICHEP->SetBinContent(7,0.04/hICHEP->GetBinWidth(7)); hICHEP->SetBinError(7,0.01/hICHEP->GetBinWidth(7));
    hICHEP->SetBinContent(8,0.02/hICHEP->GetBinWidth(8)); hICHEP->SetBinError(8,0.007/hICHEP->GetBinWidth(8));
    
    
    ////// 3S Bin
    /*
    hICHEP->SetBinContent(1,0.26/hICHEP->GetBinWidth(1)); hICHEP->SetBinError(1,0.05/hICHEP->GetBinWidth(1));
    hICHEP->SetBinContent(2,0.29/hICHEP->GetBinWidth(2)); hICHEP->SetBinError(2,0.07/hICHEP->GetBinWidth(2));
    hICHEP->SetBinContent(3,0.24/hICHEP->GetBinWidth(3)); hICHEP->SetBinError(3,0.05/hICHEP->GetBinWidth(3));
    hICHEP->SetBinContent(4,0.16/hICHEP->GetBinWidth(4)); hICHEP->SetBinError(4,0.03/hICHEP->GetBinWidth(4));
    hICHEP->SetBinContent(5,0.05/hICHEP->GetBinWidth(5)); hICHEP->SetBinError(5,0.013/hICHEP->GetBinWidth(5));
    hICHEP->SetBinContent(6,0.03/hICHEP->GetBinWidth(6)); hICHEP->SetBinError(6,0.075/hICHEP->GetBinWidth(6));
    */
    
    makeCanvas(1);
    //c1->Divide(3,1);
    c1->cd(1);
    fAcceptanceProjPt->GetXaxis()->SetTitle("P_{T}");
    fAcceptanceProjPt->SetMinimum(0.2);
    fAcceptanceProjPt->SetMaximum(1.05);
    fAcceptanceProjPt->SetMarkerStyle(22);
    fAcceptanceProjPt->SetMarkerColor(2);
    fAcceptanceProjPt->SetLineColor(2);
    //fAcceptanceProjPt->Draw("p");
    legg = new TLegend(0.6,0.6,0.8,0.8);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Acceptance");
    legge = legg->AddEntry(fAcceptanceProjPt,  "Y(1S) Acceptance ","p"); legge->SetTextColor(kBlack);
    //legg->Draw();
    //makeCanvas(1);
    fS1YieldPt->SetTitle("Differential XSection");
    fS1YieldPt->GetXaxis()->SetTitle("P_{T}");
    fS1YieldPt->SetMinimum(0.);
    fS1YieldPt->SetMaximum(1.5);
    fS1YieldPt->SetMarkerStyle(21);
    fS1YieldPt->SetMarkerColor(3);
    fS1YieldPt->SetLineColor(3);
    hICHEP->SetMarkerStyle(20);
    hICHEP->SetMarkerColor(4);
    hICHEP->SetLineColor(4);
    //fS1YieldPt->Draw("p");
    //hICHEP->Draw("psame");
    legg = new TLegend(0.6,0.6,0.8,0.8);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("XSection Comparison For Ups(1S)");
    legge = legg->AddEntry(fS1YieldPt,  "Xsection ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hICHEP,  "Results with 3.1 pb^{-1}","p"); legge->SetTextColor(kBlack);
    //legg->Draw();
    
    //makeCanvas(1);
    //c1->Divide(3,1);
    //c1->cd(1);
    fAcceptanceProjPt_2S->GetXaxis()->SetTitle("P_{T}");
    fAcceptanceProjPt_2S->SetMinimum(0.2);
    fAcceptanceProjPt_2S->SetMaximum(1.05);
    fAcceptanceProjPt_2S->SetMarkerStyle(22);
    fAcceptanceProjPt_2S->SetMarkerColor(2);
    fAcceptanceProjPt_2S->SetLineColor(2);
    //fAcceptanceProjPt->Draw("p");
    legg = new TLegend(0.6,0.6,0.8,0.8);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Acceptance");
    legge = legg->AddEntry(fAcceptanceProjPt_2S,  "Y(2S) Acceptance ","p"); legge->SetTextColor(kBlack);
    //legg->Draw();
    //makeCanvas(1);
    fS2YieldPt->SetTitle("Differential XSection");
    fS2YieldPt->GetXaxis()->SetTitle("P_{T}");
    fS2YieldPt->SetMinimum(0.);
    fS2YieldPt->SetMaximum(0.4);
    fS2YieldPt->SetMarkerStyle(21);
    fS2YieldPt->SetMarkerColor(3);
    fS2YieldPt->SetLineColor(3);
    hICHEP->SetMarkerStyle(20);
    hICHEP->SetMarkerColor(4);
    hICHEP->SetLineColor(4);
    fS2YieldPt->Draw("p");
    hICHEP->Draw("psame");
    legg = new TLegend(0.4,0.6,0.7,0.8);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("XSection Comparison For Ups(2S)");
    legge = legg->AddEntry(fS2YieldPt,  "Xsection ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hICHEP,  "Results with 3.1 pb^{-1}","p"); legge->SetTextColor(kBlack);
    legg->Draw();
    

    //makeCanvas(1);
    //c1->Divide(3,1);
    //c1->cd(1);
    fAcceptanceProjPt_3S->GetXaxis()->SetTitle("P_{T}");
    fAcceptanceProjPt_3S->SetMinimum(0.2);
    fAcceptanceProjPt_3S->SetMaximum(1.05);
    fAcceptanceProjPt_3S->SetMarkerStyle(22);
    fAcceptanceProjPt_3S->SetMarkerColor(2);
    fAcceptanceProjPt_3S->SetLineColor(2);
    //fAcceptanceProjPt->Draw("p");
    legg = new TLegend(0.6,0.6,0.8,0.8);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Acceptance");
    legge = legg->AddEntry(fAcceptanceProjPt_3S,  "Y(3S) Acceptance ","p"); legge->SetTextColor(kBlack);
    //legg->Draw();
    //makeCanvas(1);
    fS3YieldPt->SetTitle("Differential XSection");
    fS3YieldPt->GetXaxis()->SetTitle("P_{T}");
    fS3YieldPt->SetMinimum(0.);
    fS3YieldPt->SetMaximum(0.15);
    fS3YieldPt->SetMarkerStyle(21);
    fS3YieldPt->SetMarkerColor(3);
    fS3YieldPt->SetLineColor(3);
    hICHEP->SetMarkerStyle(20);
    hICHEP->SetMarkerColor(4);
    hICHEP->SetLineColor(4);
    //fS3YieldPt->Draw("p");
    //hICHEP->Draw("psame");
    legg = new TLegend(0.4,0.6,0.7,0.8);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("XSection Comparison For Ups(3S)");
    legge = legg->AddEntry(fS3YieldPt,  "Xsection ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hICHEP,  "Results with 3.1 pb^{-1}","p"); legge->SetTextColor(kBlack);
    //legg->Draw();

    
    
    bin_contentAll=0; bin_contentYield=0;
    for (int i = 1; i <= fAcceptance->GetNbinsX(); ++i) {
      for (int j = 1; j <= fAcceptance->GetNbinsY(); ++j){
	
	bin_contentYield += fS1YieldCorrected->GetCellContent(i,j);
	
      }
      
      fS1YieldEta->SetBinContent(i,bin_contentYield/lumi);
      bin_contentAll=0;
      bin_contentYield=0;
      
    }
    
    for (int i = 1; i <= fAcceptance_2S->GetNbinsX(); ++i) {
      for (int j = 1; j <= fAcceptance_2S->GetNbinsY(); ++j){
	
	bin_contentYield += fS2YieldCorrected->GetCellContent(i,j);
	
      }
      
      fS2YieldEta->SetBinContent(i,bin_contentYield/lumi);
      bin_contentAll=0;
      bin_contentYield=0;
      
    }
    
    for (int i = 1; i <= fAcceptance_3S->GetNbinsX(); ++i) {
      for (int j = 1; j <= fAcceptance_3S->GetNbinsY(); ++j){
	
	bin_contentYield += fS3YieldCorrected->GetCellContent(i,j);
	
      }
      
      fS3YieldEta->SetBinContent(i,bin_contentYield/lumi);
      bin_contentAll=0;
      bin_contentYield=0;
      
    }    
    
    /*
    hICHEP_ETA->SetBinContent(1,1.79); hICHEP_ETA->SetBinError(1,0.3);
    hICHEP_ETA->SetBinContent(2,2.05); hICHEP_ETA->SetBinError(2,0.31);
    hICHEP_ETA->SetBinContent(3,2.05); hICHEP_ETA->SetBinError(3,0.31); 
    hICHEP_ETA->SetBinContent(4,1.79); hICHEP_ETA->SetBinError(4,0.3);    
    */    
    
    
    
    //makeCanvas(1);
    //c1->cd(1);
    fS1YieldEta->SetTitle("Yield Comparison");
    fS1YieldEta->GetXaxis()->SetTitle("Rapidity");
    fS1YieldEta->SetMinimum(1.2);
    fS1YieldEta->SetMaximum(2.4);
    fS1YieldEta->SetMarkerStyle(21);
    fS1YieldEta->SetMarkerColor(3);
    fS1YieldEta->SetLineColor(3);
    hICHEP_ETA->SetMarkerStyle(20);
    hICHEP_ETA->SetMarkerColor(4);
    hICHEP_ETA->SetLineColor(4);
    //fS1YieldEta->Draw("p");
    //hICHEP_ETA->Draw("psame");
    legg = new TLegend(0.6,0.7,0.8,0.9);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Yield Comparison For Ups(1S)");
    legge = legg->AddEntry(fS1YieldEta,  "Reconstructed Yield ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hICHEP_ETA,  "Results with 3.1 pb^{-1}","p"); legge->SetTextColor(kBlack);
    //legg->Draw();
    
    //makeCanvas(1);
    //c1->cd(1);
    fS2YieldEta->SetTitle("Yield Comparison");
    fS2YieldEta->GetXaxis()->SetTitle("Rapidity");
    fS2YieldEta->SetMinimum(0.0);
    fS2YieldEta->SetMaximum(1.0);
    fS2YieldEta->SetMarkerStyle(21);
    fS2YieldEta->SetMarkerColor(3);
    fS2YieldEta->SetLineColor(3);
    hICHEP_ETA->SetMarkerStyle(20);
    hICHEP_ETA->SetMarkerColor(4);
    hICHEP_ETA->SetLineColor(4);
    //fS2YieldEta->Draw("p");
    //hICHEP_ETA->Draw("psame");
    legg = new TLegend(0.6,0.7,0.8,0.9);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Yield Comparison For Ups(2S)");
    legge = legg->AddEntry(fS2YieldEta,  "Reconstructed Yield ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hICHEP_ETA,  "Results with 3.1 pb^{-1}","p"); legge->SetTextColor(kBlack);
    //legg->Draw();
    
    
    //makeCanvas(1);
    //c1->cd(1);
    fS3YieldEta->SetTitle("Yield Comparison");
    fS3YieldEta->GetXaxis()->SetTitle("Rapidity");
    fS3YieldEta->SetMinimum(0.0);
    fS3YieldEta->SetMaximum(0.4);
    fS3YieldEta->SetMarkerStyle(21);
    fS3YieldEta->SetMarkerColor(3);
    fS3YieldEta->SetLineColor(3);
    hICHEP_ETA->SetMarkerStyle(20);
    hICHEP_ETA->SetMarkerColor(4);
    hICHEP_ETA->SetLineColor(4);
    //fS3YieldEta->Draw("p");
    //hICHEP_ETA->Draw("psame");
    legg = new TLegend(0.6,0.7,0.8,0.9);
    legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
    legg->SetHeader("Yield Comparison For Ups(3S)");
    legge = legg->AddEntry(fS3YieldEta,  "Reconstructed Yield ","p"); legge->SetTextColor(kBlack);
    legge = legg->AddEntry(hICHEP_ETA,  "Results with 3.1 pb^{-1}","p"); legge->SetTextColor(kBlack);
    //legg->Draw();
     
    
    
  }
  
  
}

// ----------------------------------------------------------------------
void anaXS::projections() {

  zone();
  TH2D *hpt = fPtMctNeg->get2dHist("hpt", "hpt", 1); 

  //  hpt->Print("all"); 

  // -- momentum
  TH1D *hPmct = new TH1D("hPmct", "", hpt->GetNbinsY(), hpt->GetYaxis()->GetXbins()->GetArray()); 
  TH1D *hPmcp = new TH1D("hPmcp", "", hpt->GetNbinsY(), hpt->GetYaxis()->GetXbins()->GetArray()); 
  TH1D *hPtnp = new TH1D("hPtnp", "", hpt->GetNbinsY(), hpt->GetYaxis()->GetXbins()->GetArray()); 
  TH1D *hPmmb = new TH1D("hPmmb", "", hpt->GetNbinsY(), hpt->GetYaxis()->GetXbins()->GetArray()); 

  hPmct->SetMaximum(1.1); 
  setHist(hPmct, kRed, 25, 1.4); 
  setHist(hPmcp, kBlue, 24, 1.4); 
  setHist(hPtnp, kBlack, 20, 1.4); 
  setHist(hPmmb, kBlack, 26, 1.4); 

  double lo(0.), hi(0.); 
  int maxPad (0); 
  if (hpt->GetNbinsX() <= 4) {
    zone(2,2);
    maxPad = 4; 
  } else if (hpt->GetNbinsX() <= 6) {
    zone(2,3);
    maxPad = 6; 
  } else if (hpt->GetNbinsX() <= 9) {
    zone(3,3);
    maxPad = 9; 
  } else if (hpt->GetNbinsX() <= 12) {
    zone(3,4);
    maxPad = 12; 
  } else if (hpt->GetNbinsX() <= 16) {
    zone(4,4);
    maxPad = 16; 
  } else if (hpt->GetNbinsX() <= 20) {
    zone(4,5);
    maxPad = 20; 
  } else if (hpt->GetNbinsX() <= 24) {
    zone(4,6);
    maxPad = 24;
  }
  
  for (int i = 1; i <= hpt->GetNbinsX(); ++i) {
    c0->cd(i); shrinkPad(0.20, 0.15); 
    lo = hpt->GetXaxis()->GetBinLowEdge(i); 
    hi = hpt->GetXaxis()->GetBinLowEdge(i+1); 

    cout << "lo: " << lo << " hi: " << hi << endl;
    fPtMctNeg->projectP(hPmct, lo, hi);
    fPtMcpNeg->projectP(hPmcp, lo, hi);
    fPtTnpNeg->projectP(hPtnp, lo, hi);
    fPtMmbNeg->projectP(hPmmb, lo, hi);

    
    hPmct->SetMinimum(0.15); 
    setTitles(hPmct, "p_{T} [GeV]", "Efficiency", 0.09, 1.1, 0.8, 0.09);

    hPmct->DrawCopy("e");
    hPmcp->DrawCopy("samee");
    hPtnp->DrawCopy("samee");
    hPmmb->DrawCopy("samee");

    tl->DrawLatex(0.2, 0.92, Form("%3.1f < #eta < %3.1f", lo, hi)); 

    c0->Modified();
    c0->Update();
  }

  c0->cd(maxPad); 
  legg = new TLegend(0.3,0.3,0.8,0.8);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
  legg->SetHeader("negative muons");
  legge = legg->AddEntry(hPmct, " MC truth","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hPmcp,  " MC probe","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hPtnp,  " tag and probe 1","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hPmmb,  " tag and probe 2","p"); legge->SetTextColor(kBlack);
  legg->Draw();
  c0->Modified();
  c0->Update();

  c0->SaveAs(Form("%s/projections-%s-neg-pt.eps", fPtDirectory.c_str(), fSample.c_str())); 


  if (hpt->GetNbinsX() <= 4) {
    zone(2,2);
    maxPad = 4; 
  } else if (hpt->GetNbinsX() <= 6) {
    zone(2,3);
    maxPad = 6; 
  } else if (hpt->GetNbinsX() <= 9) {
    zone(3,3);
    maxPad = 9; 
  } else if (hpt->GetNbinsX() <= 12) {
    zone(3,4);
    maxPad = 12; 
  } else if (hpt->GetNbinsX() <= 16) {
    zone(4,4);
    maxPad = 16; 
  } else if (hpt->GetNbinsX() <= 20) {
    zone(4,5);
    maxPad = 20; 
  } else if (hpt->GetNbinsX() <= 24) {
    zone(4,6);
    maxPad = 24; 
  }
   

  for (int i = 1; i <= hpt->GetNbinsX(); ++i) {
    c0->cd(i); shrinkPad(0.20, 0.15); 
    lo = hpt->GetXaxis()->GetBinLowEdge(i); 
    hi = hpt->GetXaxis()->GetBinLowEdge(i+1); 

    cout << "lo: " << lo << " hi: " << hi << endl;
    fPtMctPos->projectP(hPmct, lo, hi);
    fPtMcpPos->projectP(hPmcp, lo, hi);
    fPtTnpPos->projectP(hPtnp, lo, hi);
    fPtMmbPos->projectP(hPmmb, lo, hi);

    hPmct->SetMinimum(0.15); 
    setTitles(hPmct, "p_{T} [GeV]", "Efficiency", 0.09, 1.1, 0.8, 0.09);

    hPmct->DrawCopy();
    hPmcp->DrawCopy("same");
    hPtnp->DrawCopy("samee");

    tl->DrawLatex(0.2, 0.92, Form("%3.1f < #eta < %3.1f", lo, hi)); 

    c0->Modified();
    c0->Update();

  }

  c0->cd(maxPad); 
  legg = new TLegend(0.3,0.3,0.8,0.8);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
  legg->SetHeader("positive muons");
  legge = legg->AddEntry(hPmct, " MC truth","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hPmcp,  " MC probe","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hPtnp,  " tag and probe 1","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hPmmb,  " tag and probe 2","p"); legge->SetTextColor(kBlack);
  legg->Draw();
  c0->Modified();
  c0->Update();

  c0->SaveAs(Form("%s/projections-%s-pos-pt.eps", fPtDirectory.c_str(), fSample.c_str())); 


  // -- eta
  TH1D *hTmct = new TH1D("hTmct", "", hpt->GetNbinsX(), hpt->GetXaxis()->GetXbins()->GetArray()); 
  TH1D *hTmcp = new TH1D("hTmcp", "", hpt->GetNbinsX(), hpt->GetXaxis()->GetXbins()->GetArray()); 
  TH1D *hTtnp = new TH1D("hTtnp", "", hpt->GetNbinsX(), hpt->GetXaxis()->GetXbins()->GetArray()); 
  TH1D *hTmmb = new TH1D("hTmmb", "", hpt->GetNbinsX(), hpt->GetXaxis()->GetXbins()->GetArray()); 

  setHist(hTmct, kRed, 25, 1.4); 
  setHist(hTmcp, kBlue, 24, 1.4); 
  setHist(hTtnp, kBlack, 20, 1.4); 
  setHist(hTmmb, kBlack, 26, 1.4); 

  if (hpt->GetNbinsY() <= 4) {
    zone(2,2);
    maxPad = 4; 
  } else if (hpt->GetNbinsY() <= 6) {
    zone(2,3);
    maxPad = 6; 
  } else if (hpt->GetNbinsY() <= 9) {
    zone(3,3);
    maxPad = 9; 
  } else if (hpt->GetNbinsY() <= 12) {
    zone(3,4);
    maxPad = 12; 
  } else if (hpt->GetNbinsY() <= 16) {
    zone(4,4);
    maxPad = 16; 
  } else if (hpt->GetNbinsY() <= 20) {
    zone(4,5);
    maxPad = 20; 
  }

  for (int i = 2; i <= hpt->GetNbinsY(); ++i) {
    c0->cd(i-1); shrinkPad(0.20, 0.15); 
    lo = hpt->GetYaxis()->GetBinLowEdge(i); 
    hi = hpt->GetYaxis()->GetBinLowEdge(i+1); 

    cout << "lo: " << lo << " hi: " << hi << endl;
    fPtMctNeg->projectT(hTmct, lo, hi);
    fPtMcpNeg->projectT(hTmcp, lo, hi);
    fPtTnpNeg->projectT(hTtnp, lo, hi);
    fPtMmbNeg->projectT(hTmmb, lo, hi);

    hTmct->SetMinimum(0.15); 
    hTmct->SetMaximum(1.1); 
    setTitles(hTmct, "#eta", "Efficiency", 0.09, 1.1, 0.8, 0.09);

    hTmct->DrawCopy("e");
    hTmcp->DrawCopy("samee");
    hTtnp->DrawCopy("samee");
    hTmmb->DrawCopy("samee");

    tl->DrawLatex(0.2, 0.92, Form("%3.1f < p_{T} < %3.1f GeV", lo, hi)); 

    c0->Modified();
    c0->Update();

  }

  c0->cd(maxPad); 
  legg = new TLegend(0.3,0.3,0.8,0.8);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
  legg->SetHeader("negative muons");
  legge = legg->AddEntry(hTmct, " MC truth","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hTmcp,  " MC probe","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hTtnp,  " tag and probe 1","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hTmmb,  " tag and probe 2","p"); legge->SetTextColor(kBlack);
  legg->Draw();
  c0->Modified();
  c0->Update();

  c0->SaveAs(Form("%s/projections-%s-neg-eta.eps", fPtDirectory.c_str(), fSample.c_str())); 

  if (hpt->GetNbinsY() <= 4) {
    zone(2,2);
    maxPad = 4; 
  } else if (hpt->GetNbinsY() <= 6) {
    zone(2,3);
    maxPad = 6; 
  } else if (hpt->GetNbinsY() <= 9) {
    zone(3,3);
    maxPad = 9; 
  } else if (hpt->GetNbinsY() <= 12) {
    zone(3,4);
    maxPad = 12; 
  } else if (hpt->GetNbinsY() <= 16) {
    zone(4,4);
    maxPad = 16; 
  } else if (hpt->GetNbinsY() <= 20) {
    zone(4,5);
    maxPad = 20; 
  }
 
  for (int i = 2; i <= hpt->GetNbinsY(); ++i) {
    c0->cd(i-1); shrinkPad(0.20, 0.15); 
    lo = hpt->GetYaxis()->GetBinLowEdge(i); 
    hi = hpt->GetYaxis()->GetBinLowEdge(i+1); 

    cout << "lo: " << lo << " hi: " << hi << endl;
    fPtMctPos->projectT(hTmct, lo, hi);
    fPtMcpPos->projectT(hTmcp, lo, hi);
    fPtTnpPos->projectT(hTtnp, lo, hi);
    fPtMmbPos->projectT(hTmmb, lo, hi);

    hTmct->SetMinimum(0.15); 
    hTmct->SetMaximum(1.1); 
    setTitles(hTmct, "#eta", "Efficiency", 0.09, 1.1, 0.8, 0.09);
    //    hTmct->SetTitle(Form("positive muons (%3.1f < p_{T} < %3.1f GeV, %s)", lo, hi, fSample.c_str())); 

    hTmct->DrawCopy("e");
    hTmcp->DrawCopy("samee");
    hTtnp->DrawCopy("samee");
    hTmmb->DrawCopy("samee");

    tl->DrawLatex(0.2, 0.92, Form("%3.1f < p_{T} < %3.1f GeV", lo, hi)); 

    c0->Modified();
    c0->Update();

  }
  c0->cd(maxPad); 
  legg = new TLegend(0.3,0.3,0.8,0.8);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.08); legg->SetTextFont(132); 
  legg->SetHeader("positive muons");
  legge = legg->AddEntry(hTmct, " MC truth","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hTmcp,  " MC probe","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hTtnp,  " tag and probe 1","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(hTmmb,  " tag and probe 2","p"); legge->SetTextColor(kBlack);
  legg->Draw();
  c0->Modified();
  c0->Update();
  c0->SaveAs(Form("%s/projections-%s-pos-eta.eps", fPtDirectory.c_str(), fSample.c_str())); 

}


void anaXS::ReadHistogramsDATA0(TFile *f, const char *s2, const char *s3, const char *s4, const char *binning) {  
  
  cout << "====> Reading histograms from " << f->GetName() << endl;
  
  TObject *obj;
  TKey    *key;
  TIter next(f->GetListOfKeys());

  TH1D  *h1, *h2; 
  TH2D  *h3; 
  float etamin, etamax, ptmin, ptmax; 
  int   n; 
  char searchString2[2000], searchString3[2000], searchString4[2000];
  char sp[] = "%"; 
  char sf[] = "f";
  sprintf(searchString2, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s2, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString2: " << searchString2 << endl;
  sprintf(searchString3, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s3, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString3: " << searchString3 << endl;   
  sprintf(searchString4, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s4, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString4: " << searchString4 << endl;    
  
  while ((key = (TKey*)next())) {
    n = -1; 
    obj = key->ReadObj();
    if (obj->InheritsFrom(TH2D::Class()) && !strcmp(obj->GetName(), binning)) {
      cout << "Extracting binning histogram from " << obj->GetName() << endl;
      fHbinning = (TH2D*)obj;
    }
    
    if (obj->InheritsFrom(TH2D::Class())) {
      h3 = (TH2D*)f->Get(obj->GetName());
      if ( !strcmp(h3->GetName(),"AllGenRes_1S")) {
	fAllGenRes = (TH2D*)h3->Clone("AllGenRes_1S");
	cout << "--> : " << fAllGenRes->GetName() << ", # entries = " << fAllGenRes->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"RecoGenRes_1S")) {
	fRecoGenRes = (TH2D*)h3->Clone("RecoGenRes_1S");
	cout << "--> : " << fRecoGenRes->GetName() << ", # entries = " << fRecoGenRes->GetSumOfWeights() << endl;
      }	
      
      if ( !strcmp(h3->GetName(),"PreSel_afterVtx_1S")) {
	fPreSelAfter = (TH2D*)h3->Clone("PreSel_afterVtx_1S");
	cout << "--> : " << fPreSelAfter->GetName() << ", # entries = " << fPreSelAfter->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"PreSel_beforeVtx_1S")) {
	fPreSelBefore = (TH2D*)h3->Clone("PreSel_beforeVtx_1S");
	cout << "--> : " << fPreSelBefore->GetName() << ", # entries = " << fPreSelBefore->GetSumOfWeights() << endl;
      }	
      
      if ( !strcmp(h3->GetName(),"TrueYield_1S")) {
	fTrueYield_1S = (TH2D*)h3->Clone("TrueYield_1S");
	cout << "--> : " << fTrueYield_1S->GetName() << ", # entries = " << fTrueYield_1S->GetSumOfWeights() << endl;
      }      
      
      
      if ( !strcmp(h3->GetName(),"AllGenRes_2S")) {
	fAllGenRes_2S = (TH2D*)h3->Clone("AllGenRes_2S");
	cout << "--> : " << fAllGenRes_2S->GetName() << ", # entries = " << fAllGenRes_2S->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"RecoGenRes_2S")) {
	fRecoGenRes_2S = (TH2D*)h3->Clone("RecoGenRes_2S");
	cout << "--> : " << fRecoGenRes_2S->GetName() << ", # entries = " << fRecoGenRes_2S->GetSumOfWeights() << endl;
      }	
      
      if ( !strcmp(h3->GetName(),"PreSel_afterVtx_2S")) {
	fPreSelAfter_2S = (TH2D*)h3->Clone("PreSel_afterVtx_2S");
	cout << "--> : " << fPreSelAfter_2S->GetName() << ", # entries = " << fPreSelAfter_2S->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"PreSel_beforeVtx_2S")) {
	fPreSelBefore_2S = (TH2D*)h3->Clone("PreSel_beforeVtx_2S");
	cout << "--> : " << fPreSelBefore_2S->GetName() << ", # entries = " << fPreSelBefore_2S->GetSumOfWeights() << endl;
      }	
      
      if ( !strcmp(h3->GetName(),"TrueYield_2S")) {
	fTrueYield_2S = (TH2D*)h3->Clone("TrueYield_2S");
	cout << "--> : " << fTrueYield_2S->GetName() << ", # entries = " << fTrueYield_2S->GetSumOfWeights() << endl;
      }      


      
      if ( !strcmp(h3->GetName(),"AllGenRes_3S")) {
	fAllGenRes_3S = (TH2D*)h3->Clone("AllGenRes_3S");
	cout << "--> : " << fAllGenRes_3S->GetName() << ", # entries = " << fAllGenRes_3S->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"RecoGenRes_3S")) {
	fRecoGenRes_3S = (TH2D*)h3->Clone("RecoGenRes_3S");
	cout << "--> : " << fRecoGenRes_3S->GetName() << ", # entries = " << fRecoGenRes_3S->GetSumOfWeights() << endl;
      }	
      
      if ( !strcmp(h3->GetName(),"PreSel_afterVtx_3S")) {
	fPreSelAfter_3S = (TH2D*)h3->Clone("PreSel_afterVtx_3S");
	cout << "--> : " << fPreSelAfter_3S->GetName() << ", # entries = " << fPreSelAfter_3S->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"PreSel_beforeVtx_3S")) {
	fPreSelBefore_3S = (TH2D*)h3->Clone("PreSel_beforeVtx_3S");
	cout << "--> : " << fPreSelBefore_3S->GetName() << ", # entries = " << fPreSelBefore_3S->GetSumOfWeights() << endl;
      }	
      
      if ( !strcmp(h3->GetName(),"TrueYield_3S")) {
	fTrueYield_3S = (TH2D*)h3->Clone("TrueYield_3S");
	cout << "--> : " << fTrueYield_3S->GetName() << ", # entries = " << fTrueYield_3S->GetSumOfWeights() << endl;
      }      

      
    }
    
    if (obj->InheritsFrom(TH1D::Class())) {
      h1 = (TH1D*)f->Get(obj->GetName());
	
      // -- sample AnaEff_1S
      n = sscanf(h1->GetName(), searchString2, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s2: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s2:%s", h1->GetName())); 
	fS2Vector.push_back(*h2); 
      }
      
      // -- sample AnaEff_2S
      n = sscanf(h1->GetName(), searchString3, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s3: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s3:%s", h1->GetName())); 
	fS3Vector.push_back(*h2); 
      }
      
      // -- sample AnaEff_3S
      n = sscanf(h1->GetName(), searchString4, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s4: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s4:%s", h1->GetName())); 
	fS4Vector.push_back(*h2); 
      }      
      
    }
    
  }
  
  fAnaEff = new TH2D("fAnaEff", "fAnaEff", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );
  
  fAnaEff_2S = new TH2D("fAnaEff_2S", "fAnaEff_2S", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );
  
  fAnaEff_3S = new TH2D("fAnaEff_3S", "fAnaEff_3S", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );
  
  fAcceptance = new TH2D("fAcceptance", "fAcceptance", 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 );
  
  fAcceptance_2S = new TH2D("fAcceptance_2S", "fAcceptance_2S", 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 );  

  fAcceptance_3S = new TH2D("fAcceptance_3S", "fAcceptance_3S", 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 );  
  
  fAcceptanceProjPt = new TH1D("Acceptance", "Acceptance", 
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );
  
  fAcceptanceProjPt_2S = new TH1D("Acceptance_2S", "Acceptance_2S", 
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );  
  
  fAcceptanceProjPt_3S = new TH1D("Acceptance_3S", "Acceptance_3S", 
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );
  
  fPreSelEff = new TH2D("fPreSelEff", "fPreSelEff", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);
  
  fPreSelEff_2S = new TH2D("fPreSelEff_2S", "fPreSelEff_2S", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);  
  
  fPreSelEff_3S = new TH2D("fPreSelEff_3S", "fPreSelEff_3S", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);    
  
  fAllGenResPt = new TH1D("AllGenResPt", "AllGenResPt", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);  
  
  fAllGenResPt_2S = new TH1D("AllGenResPt_2S", "AllGenResPt_2S", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);    
  
  fAllGenResPt_3S = new TH1D("AllGenResPt_3S", "AllGenResPt_3S", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);    
   
  fAllGenResEta = new TH1D("AllGenResEta", "AllGenResEta", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			);    
  
  fAllGenResEta_2S = new TH1D("AllGenResEta_2S", "AllGenResEta_2S", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			);    
  
  fAllGenResEta_3S = new TH1D("AllGenResEta_3S", "AllGenResEta_3S", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			);    
  
  
}


void anaXS::ReadHistogramsDATA1(TFile *f, const char *s1, const char *s5, const char *s8, const char *binning) {  
  
  cout << "====> Reading histograms from " << f->GetName() << endl;
  
  TObject *obj;
  TKey    *key;
  TIter next(f->GetListOfKeys());

  TH1D  *h1, *h2; 
  float etamin, etamax, ptmin, ptmax; 
  int   n; 
  char searchString1[2000], searchString5[2000] ,searchString8[2000];
  char sp[] = "%"; 
  char sf[] = "f";
  sprintf(searchString1, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s1, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString1: " << searchString1 << endl;
  
  sprintf(searchString5, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s5, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString5: " << searchString5 << endl;
  
  sprintf(searchString8, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s8, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString8: " << searchString8 << endl;

  while ((key = (TKey*)next())) {
    n = -1; 
    obj = key->ReadObj();
    if (obj->InheritsFrom(TH2D::Class()) && !strcmp(obj->GetName(), binning)) {
      cout << "Extracting binning histogram from " << obj->GetName() << endl;
      fHbinning = (TH2D*)obj;
    }
    
    if (obj->InheritsFrom(TH1D::Class())) {
      h1 = (TH1D*)f->Get(obj->GetName());
	
      // -- sample UpsilonMass
      n = sscanf(h1->GetName(), searchString1, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s1: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s1:%s", h1->GetName()));
	fS1Vector.push_back(*h2); 
	
      }
      
      // -- sample MuIDEff_1S
      n = sscanf(h1->GetName(), searchString5, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s5: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s5:%s", h1->GetName())); 
	fS5Vector.push_back(*h2); 
      }
      
      // -- sample TrigEff_1S
      n = sscanf(h1->GetName(), searchString8, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s8: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s8:%s", h1->GetName())); 
	fS8Vector.push_back(*h2); 
      }
      
      
    }
    
  }
  
  fS1Yield = new TH2D("fS1Yield", 
		      Form("%s fS1Yield", s1), 
		      fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		      fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		      );
  
  fS2Yield = new TH2D("fS2Yield", 
		      Form("%s fS2Yield", s1), 
		      fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		      fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		      );
  fS3Yield = new TH2D("fS3Yield", 
		      Form("%s fS3Yield", s1), 
		      fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		      fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		      );
  
  falpha = new TH2D("falpha", 
		      Form("%s falpha", s1), 
		      fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		      fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		      );    
  
  fn = new TH2D("fn", 
		      Form("%s fn", s1), 
		      fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		      fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		      );    

  fMuIDEff = new TH2D("fMuIDEff", "fMuIdEff", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );  
  
  fTrigEff = new TH2D("fTrigEff", "fTrigEff", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );   
  
  
  fS1YieldPt = new TH1D("S1YieldPt", "S1YieldPt", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);
  
  fS2YieldPt = new TH1D("S2YieldPt", "S2YieldPt", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);  
  
  fS3YieldPt = new TH1D("S3YieldPt", "S3YieldPt", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);
  
  fS1YieldEta = new TH1D("S1YieldEta", "S1YieldEta", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			);  
  
  fS2YieldEta = new TH1D("S2YieldEta", "S2YieldEta", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			); 
  
  fS3YieldEta = new TH1D("S3YieldEta", "S3YieldEta", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			);  
  
  fS1YieldCorrected = new TH2D("fS1YieldCorrected", 
			       Form("%s fS1YieldCorrected", s1), 
			       fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );
  
  fS2YieldCorrected = new TH2D("fS2YieldCorrected", 
			       Form("%s fS2YieldCorrected", s1), 
			       fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );  
  
  fS3YieldCorrected = new TH2D("fS3YieldCorrected", 
			       Form("%s fS3YieldCorrected", s1), 
			       fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );
  
  fS1YieldComparison = new TH2D("fS1YieldComparison", 
				Form("%s fS1YieldComparison", s1), 
				fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
				fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
				);
  
    
}



void anaXS::ReadHistograms(TFile *f, const char *s1, const char *s2, const char *s3, const char *s4, const char *s5, const char *s6, const char *s7, const char *s8, const char *s9, const char *s10, const char *binning) {  
  
  cout << "====> Reading histograms from " << f->GetName() << endl;
  
  TObject *obj;
  TKey    *key;
  TIter next(f->GetListOfKeys());

  TH1D  *h1, *h2; 
  TH2D  *h3; 
  float etamin, etamax, ptmin, ptmax; 
  int   n; 
  char searchString1[2000], searchString2[2000] ,searchString3[2000], searchString4[2000], searchString5[2000] ;
  char searchString6[2000], searchString7[2000] ,searchString8[2000], searchString9[2000], searchString10[2000];
  char sp[] = "%"; 
  char sf[] = "f";
  sprintf(searchString1, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s1, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString1: " << searchString1 << endl;
  
  sprintf(searchString2, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s2, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString2: " << searchString2 << endl;
  
  sprintf(searchString3, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s3, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString3: " << searchString3 << endl;

  sprintf(searchString4, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s4, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString4: " << searchString4 << endl;
  
  sprintf(searchString5, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s5, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString5: " << searchString5 << endl;
  
  sprintf(searchString6, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s6, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString6: " << searchString6 << endl;
  
  sprintf(searchString7, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s7, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString7: " << searchString7 << endl;

  sprintf(searchString8, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s8, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString8: " << searchString8 << endl;

  sprintf(searchString9, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s9, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString9: " << searchString9 << endl;

  sprintf(searchString10, "%s,rapidity%s%s_%s%s,pt%s%s_%s%s", s10, sp, sf, sp, sf, sp, sf, sp, sf); 
  cout << "searchString10: " << searchString10 << endl;
    
  while ((key = (TKey*)next())) {
    n = -1; 
    obj = key->ReadObj();
    if (obj->InheritsFrom(TH2D::Class()) && !strcmp(obj->GetName(), binning)) {
      cout << "Extracting binning histogram from " << obj->GetName() << endl;
      fHbinning = (TH2D*)obj;
    }
    
    if (obj->InheritsFrom(TH2D::Class())) {
      h3 = (TH2D*)f->Get(obj->GetName());
      if ( !strcmp(h3->GetName(),"AllGenRes_1S")) {
	fAllGenRes = (TH2D*)h3->Clone("AllGenRes_1S");
	cout << "--> : " << fAllGenRes->GetName() << ", # entries = " << fAllGenRes->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"RecoGenRes_1S")) {
	fRecoGenRes = (TH2D*)h3->Clone("RecoGenRes_1S");
	cout << "--> : " << fRecoGenRes->GetName() << ", # entries = " << fRecoGenRes->GetSumOfWeights() << endl;
      }	
      
      if ( !strcmp(h3->GetName(),"PreSel_afterVtx_1S")) {
	fPreSelAfter = (TH2D*)h3->Clone("PreSel_afterVtx_1S");
	cout << "--> : " << fPreSelAfter->GetName() << ", # entries = " << fPreSelAfter->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"PreSel_beforeVtx_1S")) {
	fPreSelBefore = (TH2D*)h3->Clone("PreSel_beforeVtx_1S");
	cout << "--> : " << fPreSelBefore->GetName() << ", # entries = " << fPreSelBefore->GetSumOfWeights() << endl;
      }	
      
      if ( !strcmp(h3->GetName(),"TrueYield_1S")) {
	fTrueYield_1S = (TH2D*)h3->Clone("TrueYield_1S");
	cout << "--> : " << fTrueYield_1S->GetName() << ", # entries = " << fTrueYield_1S->GetSumOfWeights() << endl;
      }       
      
      
      
      if ( !strcmp(h3->GetName(),"AllGenRes_2S")) {
	fAllGenRes_2S = (TH2D*)h3->Clone("AllGenRes_2S");
	cout << "--> : " << fAllGenRes_2S->GetName() << ", # entries = " << fAllGenRes_2S->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"RecoGenRes_2S")) {
	fRecoGenRes_2S = (TH2D*)h3->Clone("RecoGenRes_2S");
	cout << "--> : " << fRecoGenRes_2S->GetName() << ", # entries = " << fRecoGenRes_2S->GetSumOfWeights() << endl;
      }	
      
      if ( !strcmp(h3->GetName(),"PreSel_afterVtx_2S")) {
	fPreSelAfter_2S = (TH2D*)h3->Clone("PreSel_afterVtx_2S");
	cout << "--> : " << fPreSelAfter_2S->GetName() << ", # entries = " << fPreSelAfter_2S->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"PreSel_beforeVtx_2S")) {
	fPreSelBefore_2S = (TH2D*)h3->Clone("PreSel_beforeVtx_2S");
	cout << "--> : " << fPreSelBefore_2S->GetName() << ", # entries = " << fPreSelBefore_2S->GetSumOfWeights() << endl;
      }	

      if ( !strcmp(h3->GetName(),"TrueYield_2S")) {
	fTrueYield_2S = (TH2D*)h3->Clone("TrueYield_2S");
	cout << "--> : " << fTrueYield_2S->GetName() << ", # entries = " << fTrueYield_2S->GetSumOfWeights() << endl;
      }

      
      
      if ( !strcmp(h3->GetName(),"AllGenRes_3S")) {
	fAllGenRes_3S = (TH2D*)h3->Clone("AllGenRes_3S");
	cout << "--> : " << fAllGenRes_3S->GetName() << ", # entries = " << fAllGenRes_3S->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"RecoGenRes_3S")) {
	fRecoGenRes_3S = (TH2D*)h3->Clone("RecoGenRes_3S");
	cout << "--> : " << fRecoGenRes_3S->GetName() << ", # entries = " << fRecoGenRes_3S->GetSumOfWeights() << endl;
      }	
      
      if ( !strcmp(h3->GetName(),"PreSel_afterVtx_3S")) {
	fPreSelAfter_3S = (TH2D*)h3->Clone("PreSel_afterVtx_3S");
	cout << "--> : " << fPreSelAfter_3S->GetName() << ", # entries = " << fPreSelAfter_3S->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"PreSel_beforeVtx_3S")) {
	fPreSelBefore_3S = (TH2D*)h3->Clone("PreSel_beforeVtx_3S");
	cout << "--> : " << fPreSelBefore_3S->GetName() << ", # entries = " << fPreSelBefore_3S->GetSumOfWeights() << endl;
      }
      
      if ( !strcmp(h3->GetName(),"TrueYield_3S")) {
	fTrueYield_3S = (TH2D*)h3->Clone("TrueYield_3S");
	cout << "--> : " << fTrueYield_3S->GetName() << ", # entries = " << fTrueYield_3S->GetSumOfWeights() << endl;
      }      
      
      
      

      
    }
    
    if (obj->InheritsFrom(TH1D::Class())) {
      h1 = (TH1D*)f->Get(obj->GetName());
	
      // -- sample UpsilonMass
      n = sscanf(h1->GetName(), searchString1, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s1: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s1:%s", h1->GetName()));
	fS1Vector.push_back(*h2); 
	
      }
      
      // -- sample AnaEff_1S
      n = sscanf(h1->GetName(), searchString2, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s2: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s2:%s", h1->GetName())); 
	fS2Vector.push_back(*h2); 
      }
      
      // -- sample AnaEff_2S
      n = sscanf(h1->GetName(), searchString3, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s3: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s3:%s", h1->GetName())); 
	fS3Vector.push_back(*h2); 
      }
      
      // -- sample AnaEff_3S
      n = sscanf(h1->GetName(), searchString4, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s4: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s4:%s", h1->GetName())); 
	fS4Vector.push_back(*h2); 
      }	
      
      // -- sample MuIDEff_1S
      n = sscanf(h1->GetName(), searchString5, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s5: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s5:%s", h1->GetName())); 
	fS5Vector.push_back(*h2); 
      }
      
      // -- sample MuIDEff_2S
      n = sscanf(h1->GetName(), searchString6, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s6: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s6:%s", h1->GetName())); 
	fS6Vector.push_back(*h2); 
      }
      
      // -- sample MuIDEff_3S
      n = sscanf(h1->GetName(), searchString7, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s7: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s7:%s", h1->GetName())); 
	fS7Vector.push_back(*h2); 
      }	
      
      // -- sample TrigEff_1S
      n = sscanf(h1->GetName(), searchString8, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s8: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s8:%s", h1->GetName())); 
	fS8Vector.push_back(*h2); 
      }
      
      // -- sample TrigEff_2S
      n = sscanf(h1->GetName(), searchString9, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s9: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s9:%s", h1->GetName())); 
	fS9Vector.push_back(*h2); 
      }
      
      // -- sample TrigEff_3S
      n = sscanf(h1->GetName(), searchString10, &etamin, &etamax, &ptmin, &ptmax);
      if (n > 0) {
	cout << "s10: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s10:%s", h1->GetName())); 
	fS10Vector.push_back(*h2); 
      }	      
      
      
    }
    
  }
  
  fS1Yield = new TH2D("fS1Yield", 
		      Form("%s fS1Yield", s1), 
		      fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		      fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		      );
  
  fS2Yield = new TH2D("fS2Yield", 
		      Form("%s fS2Yield", s1), 
		      fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		      fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		      );
  
  fS3Yield = new TH2D("fS3Yield", 
		      Form("%s fS3Yield", s1), 
		      fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		      fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		      );  
  
  falpha = new TH2D("falpha", 
		      Form("%s falpha", s1), 
		      fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		      fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		      );    
  
  fn = new TH2D("fn", 
		      Form("%s fn", s1), 
		      fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		      fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		      );    
  
  fAnaEff = new TH2D("fAnaEff", "fAnaEff", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );
  
  fAnaEff_2S = new TH2D("fAnaEff_2S", "fAnaEff_2S", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );  
  
  fAnaEff_3S = new TH2D("fAnaEff_3S", "fAnaEff_3S", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );    
  
  fMuIDEff = new TH2D("fMuIDEff", "fMuIdEff", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     ); 
  
  fMuIDEff_2 = new TH2D("fMuIDEff_2", "fMuIdEff_2", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );   
  
  fMuIDEff_3 = new TH2D("fMuIDEff_3", "fMuIdEff_3", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     ); 
  
  fTrigEff = new TH2D("fTrigEff", "fTrigEff", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );  
  
  fTrigEff_2 = new TH2D("fTrigEff_2", "fTrigEff_2", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );    
  
  fTrigEff_3 = new TH2D("fTrigEff_3", "fTrigEff_3", 
		     fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		     fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		     );  
  
  fAcceptance = new TH2D("fAcceptance", "fAcceptance", 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 );
  
  fAcceptance_2S = new TH2D("fAcceptance_2S", "fAcceptance_2S", 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 ); 
  
  fAcceptance_3S = new TH2D("fAcceptance_3S", "fAcceptance_3S", 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 );  
  
  fPreSelEff = new TH2D("fPreSelEff", "fPreSelEff", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);
  
  fPreSelEff_2S = new TH2D("fPreSelEff_2S", "fPreSelEff_2S", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);  
  
  fPreSelEff_3S = new TH2D("fPreSelEff_3S", "fPreSelEff_3S", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);    

  
  fAcceptanceProjPt = new TH1D("Acceptance", "Acceptance", 
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );
  
  fAcceptanceProjPt_2S = new TH1D("Acceptance_2S", "Acceptance_2S", 
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );
  
  fAcceptanceProjPt_3S = new TH1D("Acceptance_3S", "Acceptance_3S", 
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );
  
  
  fS1YieldPt = new TH1D("S1YieldPt", "S1YieldPt", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);
  
  fS2YieldPt = new TH1D("S2YieldPt", "S2YieldPt", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);

  fS3YieldPt = new TH1D("S3YieldPt", "S3YieldPt", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);


  fS1YieldEta = new TH1D("S1YieldEta", "S1YieldEta", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			);  
  
  fS2YieldEta = new TH1D("S2YieldEta", "S2YieldEta", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			);  
  
  fS3YieldEta = new TH1D("S3YieldEta", "S3YieldEta", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			);  
  
  fAllGenResPt = new TH1D("AllGenResPt", "AllGenResPt", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);  
  
  fAllGenResPt_2S = new TH1D("AllGenResPt_2S", "AllGenResPt_2S", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);  
  
  fAllGenResPt_3S = new TH1D("AllGenResPt_3S", "AllGenResPt_3S", 
			fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			);  
  
  
  fAllGenResEta = new TH1D("AllGenResEta", "AllGenResEta", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			);    
  
  fAllGenResEta_2S = new TH1D("AllGenResEta_2S", "AllGenResEta_2S", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			);    
  
  fAllGenResEta_3S = new TH1D("AllGenResEta_3S", "AllGenResEta_3S", 
			fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray()
			);    
  
  
  fS1YieldCorrected = new TH2D("fS1YieldCorrected", 
			       Form("%s fS1YieldCorrected", s1), 
			       fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );
  
  fS2YieldCorrected = new TH2D("fS2YieldCorrected", 
			       Form("%s fS2YieldCorrected", s1), 
			       fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );
  
  fS3YieldCorrected = new TH2D("fS3YieldCorrected", 
			       Form("%s fS3YieldCorrected", s1), 
			       fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			       );
  
  fS1YieldComparison = new TH2D("fS1YieldComparison", 
				Form("%s fS1YieldComparison", s1), 
				fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
				fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
				);
  
  fS2YieldComparison = new TH2D("fS2YieldComparison", 
				Form("%s fS2YieldComparison", s1), 
				fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
				fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
				);  
  
  fS3YieldComparison = new TH2D("fS3YieldComparison", 
				Form("%s fS3YieldComparison", s1), 
				fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
				fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
				);
  
}

// ----------------------------------------------------------------------
void anaXS::readHistograms(TFile *f, 
			     const char *s1, const char *s2, const char *s3, 
			     const char *sm, const char *binning) {  

  cout << "====> Reading histograms from " << f->GetName() << endl;
  
  TObject *obj;
  TKey    *key;
  TIter next(f->GetListOfKeys());

  TH1D  *h1, *h2; 
  TH2D  *h3; 
  float etamin, etamax, ptmin, ptmax; 
  int   charge, n; 
  char searchString1[2000], searchString2[2000], searchString3[2000], 
    searchString1MC[2000], searchString2MC[2000]; 
  char sp[] = "%"; 
  char sf[] = "f"; 
  char sd[] = "d"; 
  sprintf(searchString1, "%s,eta%s%s_%s%s,pt%s%s_%s%s,Q%s%s", s1, sp, sf, sp, sf, sp, sf, sp, sf, sp, sd); 
  cout << "searchString1: " << searchString1 << endl;
  sprintf(searchString1MC, "%s%s,eta%s%s_%s%s,pt%s%s_%s%s,Q%s%s", s1, sm, sp, sf, sp, sf, sp, sf, sp, sf, sp, sd); 
  cout << "searchString1MC: " << searchString1MC << endl;

  sprintf(searchString2, "%s,eta%s%s_%s%s,pt%s%s_%s%s,Q%s%s", s2, sp, sf, sp, sf, sp, sf, sp, sf, sp, sd); 
  cout << "searchString2: " << searchString2 << endl;
  sprintf(searchString2MC, "%s%s,eta%s%s_%s%s,pt%s%s_%s%s,Q%s%s", s2, sm, sp, sf, sp, sf, sp, sf, sp, sf, sp, sd); 
  cout << "searchString2MC: " << searchString2MC << endl;

  sprintf(searchString3, "%s,eta%s%s_%s%s,pt%s%s_%s%s,Q%s%s", s3, sp, sf, sp, sf, sp, sf, sp, sf, sp, sd); 
  cout << "searchString3: " << searchString3 << endl;


  while ((key = (TKey*)next())) {
    n = -1; 
    obj = key->ReadObj();
    if (obj->InheritsFrom(TH2D::Class()) && !strcmp(obj->GetName(), binning)) {
      cout << "Extracting binning histogram from " << obj->GetName() << endl;
      fHbinning = (TH2D*)obj;
    }
    if (obj->InheritsFrom(TH1D::Class())) {
      h1 = (TH1D*)f->Get(obj->GetName());

      // -- sample 1
      n = sscanf(h1->GetName(), searchString1, &etamin, &etamax, &ptmin, &ptmax, &charge);
      if (n > 0) {
	cout << "s1: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s1:%s", h1->GetName())); 
	if (charge > 0) {
	  fS1VectorPos.push_back(*h2); 
	} else {
	  fS1VectorNeg.push_back(*h2); 
	}
      }

      // -- sample 1 MC probe
      n = sscanf(h1->GetName(), searchString1MC, &etamin, &etamax, &ptmin, &ptmax, &charge);
      if (n > 0) {
	cout << "mcp s1: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s1Mcp:%s", h1->GetName())); 
	if (charge > 0) {
	  fS1VectorMcpPos.push_back(*h2); 
	} else {
	  fS1VectorMcpNeg.push_back(*h2); 
	}
      }

      // -- sample 2
      n = sscanf(h1->GetName(), searchString2, &etamin, &etamax, &ptmin, &ptmax, &charge);
      if (n > 0) {
	cout << "s2: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s2:%s", h1->GetName())); 
	if (charge > 0) {
	  fS2VectorPos.push_back(*h2); 
	} else {
	  fS2VectorNeg.push_back(*h2); 
	}
      }

      // -- sample 2 MC probe
      n = sscanf(h1->GetName(), searchString2MC, &etamin, &etamax, &ptmin, &ptmax, &charge);
      if (n > 0) {
	cout << "mcp s2: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s1:%s", h1->GetName())); 
	if (charge > 0) {
	  fS2VectorMcpPos.push_back(*h2); 
	} else {
	  fS2VectorMcpNeg.push_back(*h2); 
	}
      }

      // -- sample 3 (mmbar)
      n = sscanf(h1->GetName(), searchString3, &etamin, &etamax, &ptmin, &ptmax, &charge);
      if (n > 0) {
	cout << "s3: " << h1->GetName() << ", # entries = " << h1->Integral(1, h1->GetNbinsX()) << endl;
	h2 = (TH1D*)h1->Clone(Form("s3:%s", h1->GetName())); 
	if (charge > 0) {
	  fS3VectorPos.push_back(*h2); 
	} else {
	  fS3VectorNeg.push_back(*h2); 
	}
      }

    }

    if (obj->InheritsFrom(TH2D::Class())) {
      h3 = (TH2D*)f->Get(obj->GetName());

      // -- MC truth
      if (!strcmp(h3->GetName(), "mEtaPt_neg")) {
	cout << "--> mct s1: " << h3->GetName() << ", # entries = " << h3->GetSumOfWeights() << endl;
	fS1MctNeg = (TH2D*)h3->Clone("fS1MctNeg"); 
	fS1MctNeg->SetTitle(Form("s1MctNeg:%s", h3->GetName())); 
      }
      if (!strcmp(h3->GetName(), "mEtaPt_pos")) {
	cout << "--> mct s1: " << h3->GetName() << ", # entries = " << h3->GetSumOfWeights() << endl;
	fS1MctPos = (TH2D*)h3->Clone("fS1MctPos"); 
	fS1MctPos->SetTitle(Form("s1MctPos:%s", h3->GetName())); 
      }

      if (!strcmp(h3->GetName(), "tEtaPt_neg")) {
	cout << "--> mct s2: " << h3->GetName() << ", # entries = " << h3->GetSumOfWeights() << endl;
	fS2MctNeg = (TH2D*)h3->Clone("fS2MctNeg"); 
	fS2MctNeg->SetTitle(Form("s2MctNeg:%s", h3->GetName())); 
      }
      if (!strcmp(h3->GetName(), "tEtaPt_pos")) {
	cout << "--> mct s2: " << h3->GetName() << ", # entries = " << h3->GetSumOfWeights() << endl;
	fS2MctPos = (TH2D*)h3->Clone("fS2MctPos"); 
	fS2MctPos->SetTitle(Form("s2MctPos:%s", h3->GetName())); 
      }
      
    }
  }


  fS1YieldPos = new TH2D("fS1YieldPos", 
			 Form("%s fS1YieldPos", s1), 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 );

  fS1YieldNeg = new TH2D("fS1YieldNeg", 
			 Form("%s fS1YieldNeg", s1), 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 );

  fS1McpPos = new TH2D("fS1McpPos", 
		       Form("%s fS1McpPos", s1), 
		       fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		       );
  
  fS1McpNeg = new TH2D("fS1McpNeg", 
		       Form("%s fS1McpNeg", s1), 
		       fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		       );



  fS2YieldPos = new TH2D("fS2YieldPos", 
			 Form("%s fS2YieldPos", s2), 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 );

  fS2YieldNeg = new TH2D("fS2YieldNeg", 
			 Form("%s fS2YieldNeg", s2), 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 );

  fS2McpPos = new TH2D("fS2McpPos", 
		       Form("%s fS2McpPos", s1), 
		       fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		       );
  
  fS2McpNeg = new TH2D("fS2McpNeg", 
		       Form("%s fS2McpNeg", s1), 
		       fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
		       fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
		       );



  fS3YieldPos = new TH2D("fS3YieldPos", 
			 Form("%s fS3YieldPos", s3), 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 );

  fS3YieldNeg = new TH2D("fS3YieldNeg", 
			 Form("%s fS3YieldNeg", s3), 
			 fHbinning->GetNbinsY(), fHbinning->GetYaxis()->GetXbins()->GetArray(),
			 fHbinning->GetNbinsX(), fHbinning->GetXaxis()->GetXbins()->GetArray()
			 );
  
}


// ----------------------------------------------------------------------
void anaXS::addBackground(vector<TH1D> &vec, double sb, double p0, double p1) {  
  TH1D *h;
  int s; 

  f0->SetParameters(-0.5, 1.); 
  f0->SetParameters(1., 0.); 
  for (unsigned int i = 0; i < fS1Vector.size(); ++i) {
    h = &(vec[i]);
    s = h->Integral(1, h->GetNbinsX());  
    h->FillRandom("f0", sb*s);
  }

}

void anaXS::integerEntries(TH1D  *h){

  int bincontent;
  double BinContent, remainder;
  
  for (Int_t i= 1 ; i <= h->GetNbinsX() ; i++){ 		
    
    bincontent = h->GetBinContent(i);  		
    BinContent = h->GetBinContent(i);
    remainder = BinContent - bincontent;
		
    if ( remainder >= 0.5 ){
      h->SetBinContent(i,bincontent+1);  
    } else if ( remainder < 0.5){
      h->SetBinContent(i,bincontent);
    }
  }

  
}

void anaXS::GetAnaEff(){
  
  TH1D *h;
  double pt, eta; 
  double yield, yieldE;
  double errN, errD, D;
  int    nbin;
  // AnaEff_1S
  for (unsigned int i = 0; i < fS2Vector.size(); ++i) {
    
    h = &(fS2Vector[i]);
    cout << h->GetName()  << endl;
    cout << h->GetBinContent(8)  << " "  << h->GetBinContent(2) << endl;
    yield = h->GetBinContent(8)/h->GetBinContent(2);
    D = h->GetBinContent(2);
    errN = TMath::Sqrt(h->GetBinContent(8));
    errD = TMath::Sqrt(h->GetBinContent(2));
    yieldE = TMath::Sqrt(((errN*errN)/(D*D)) + (errN*errN*errD*errD)/(D*D*D*D));
    cout << " AnaEff Ups(1S) =  "  << yield << "+/-" << yieldE << endl;
    GetBinCenters(h->GetName(), eta, pt);
    nbin = fAnaEff->FindBin(eta, pt); 
    cout << nbin  << endl;
    fAnaEff->SetBinContent(nbin, yield); 
    fAnaEff->SetBinError(nbin, yieldE); 
    
  }
  // AnaEff_2S
  for (unsigned int i = 0; i < fS3Vector.size(); ++i) {
    
    h = &(fS3Vector[i]);
    cout << h->GetName()  << endl;
    cout << h->GetBinContent(8)  << " "  << h->GetBinContent(2) << endl;
    yield = h->GetBinContent(8)/h->GetBinContent(2);
    D = h->GetBinContent(2);
    errN = TMath::Sqrt(h->GetBinContent(8));
    errD = TMath::Sqrt(h->GetBinContent(2));
    yieldE = TMath::Sqrt(((errN*errN)/(D*D)) + (errN*errN*errD*errD)/(D*D*D*D));
    cout << " AnaEff Ups(2S) =  "  << yield << "+/-" << yieldE << endl;
    GetBinCenters(h->GetName(), eta, pt);
    nbin = fAnaEff_2S->FindBin(eta, pt); 
    cout << nbin  << endl;
    fAnaEff_2S->SetBinContent(nbin, yield); 
    fAnaEff_2S->SetBinError(nbin, yieldE); 
    
  }
  // AnaEff_3S
  for (unsigned int i = 0; i < fS4Vector.size(); ++i) {
    
    h = &(fS4Vector[i]);
    cout << h->GetName()  << endl;
    cout << h->GetBinContent(8)  << " "  << h->GetBinContent(2) << endl;
    yield = h->GetBinContent(8)/h->GetBinContent(2);
    D = h->GetBinContent(2);
    errN = TMath::Sqrt(h->GetBinContent(8));
    errD = TMath::Sqrt(h->GetBinContent(2));
    yieldE = TMath::Sqrt(((errN*errN)/(D*D)) + (errN*errN*errD*errD)/(D*D*D*D));
    cout << " AnaEff  Ups(3S) =  "  << yield << "+/-" << yieldE << endl;
    GetBinCenters(h->GetName(), eta, pt);
    nbin = fAnaEff_3S->FindBin(eta, pt); 
    cout << nbin  << endl;
    fAnaEff_3S->SetBinContent(nbin, yield); 
    fAnaEff_3S->SetBinError(nbin, yieldE); 
    
  }  
    
  makeCanvas(1);
  c1->Divide(3,1);
  c1->cd(1);
  fAnaEff->Draw("colz");
  c1->cd(2);
  fAnaEff_2S->Draw("colz");  
  c1->cd(3);
  fAnaEff_3S->Draw("colz"); 
 
}

void anaXS::GetMuIDEff(int mode){
  
  /////////////////
  //TFile *f = new TFile("MuID.root", "RECREATE");
  /////////////////
  
  TH1D *h;
  double pt, eta; 
  double yield, yieldE;
  int    nbin;
  double corr(-1);
  
  if ( mode == 1 ){
    for (unsigned int i = 0; i < fS5Vector.size(); ++i) {
      
      h = &(fS5Vector[i]);
      cout << h->GetName()  << endl;
      yield = h->GetMean();
      yieldE = h->GetMeanError(); 
      GetBinCenters(h->GetName(), eta, pt);
      nbin = fMuIDEff->FindBin(eta, pt); 
      corr = fPtMuidCorr->effD(pt, eta, 0.);
      cout << nbin  << endl;
      yield*=corr;
      cout << " MuIDEff Ups(1S) " << yield << " +/- " << yieldE << endl;
      fMuIDEff->SetBinContent(nbin, yield); 
      fMuIDEff->SetBinError(nbin, yieldE); 
    }
    
    for (unsigned int i = 0; i < fS6Vector.size(); ++i) {
      
      h = &(fS6Vector[i]);
      cout << h->GetName()  << endl;
      yield = h->GetMean();
      yieldE = h->GetMeanError(); 
      GetBinCenters(h->GetName(), eta, pt);
      nbin = fMuIDEff_2->FindBin(eta, pt); 
      corr = fPtMuidCorr->effD(pt, eta, 0.);
      cout << nbin  << endl;
      yield*=corr;
      cout << " MuIDEff Ups(2S) " << yield << " +/- " << yieldE << endl;
      fMuIDEff_2->SetBinContent(nbin, yield); 
      fMuIDEff_2->SetBinError(nbin, yieldE); 
    }    
    
    for (unsigned int i = 0; i < fS7Vector.size(); ++i) {
      
      h = &(fS7Vector[i]);
      cout << h->GetName()  << endl;
      yield = h->GetMean();
      yieldE = h->GetMeanError(); 
      GetBinCenters(h->GetName(), eta, pt);
      nbin = fMuIDEff_3->FindBin(eta, pt); 
      corr = fPtMuidCorr->effD(pt, eta, 0.);
      cout << nbin  << endl;
      yield*=corr;
      cout << " MuIDEff Ups(3S) " << yield << " +/- " << yieldE << endl;
      fMuIDEff_3->SetBinContent(nbin, yield); 
      fMuIDEff_3->SetBinError(nbin, yieldE); 
    }     
    
    makeCanvas(1);
    c1->Divide(3,1);
    c1->cd(1);
    fMuIDEff->Draw("colz");
    c1->cd(2);
    fMuIDEff_2->Draw("colz");
    c1->cd(3);
    fMuIDEff_3->Draw("colz");
  }
  
  if ( mode == 2 ){
    for (unsigned int i = 0; i < fS5Vector.size(); ++i) {
      
      h = &(fS5Vector[i]);
      cout << h->GetName()  << endl;
      yield = h->GetMean();
      yieldE = h->GetMeanError(); 
      GetBinCenters(h->GetName(), eta, pt);
      nbin = fMuIDEff->FindBin(eta, pt); 
      corr = fPtMuidCorr->effD(pt, eta, 0.);
      cout << nbin  << endl;
      yield*=corr;
      cout << " MuIDEff Ups " << yield << " +/- " << yieldE << endl;
      fMuIDEff->SetBinContent(nbin, yield); 
      fMuIDEff->SetBinError(nbin, yieldE); 
    }
    makeCanvas(1);
    fMuIDEff->Draw("colz");
  }  
  
  
  ////////////
  //fMuIDEff_3->Write();
  ////////////
}

void anaXS::GetTrigEff(int mode){
  
  ////////////////
  //TFile *f = new TFile("Trig.root", "RECREATE");
  ////////////////
  TH1D *h;
  double pt, eta; 
  double yield, yieldE;
  int    nbin;
  double corr(-1);
  
  if ( mode == 1 ){
    for (unsigned int i = 0; i < fS8Vector.size(); ++i) {
      
      h = &(fS8Vector[i]);
      cout << h->GetName()  << endl;
      yield = h->GetMean(); 
      yieldE = h->GetMeanError(); 
      GetBinCenters(h->GetName(), eta, pt);
      nbin = fTrigEff->FindBin(eta, pt);
      corr = fPtTrigCorr->effD(pt, eta, 0.);
      cout << nbin  << endl;
      yield*=corr;
      cout << " TrigEff Ups(1S) " << yield << " +/- " << yieldE << endl;
      fTrigEff->SetBinContent(nbin, yield); 
      fTrigEff->SetBinError(nbin, yieldE); 
      
    }
    
    for (unsigned int i = 0; i < fS9Vector.size(); ++i) {
      
      h = &(fS9Vector[i]);
      cout << h->GetName()  << endl;
      yield = h->GetMean(); 
      yieldE = h->GetMeanError(); 
      GetBinCenters(h->GetName(), eta, pt);
      nbin = fTrigEff_2->FindBin(eta, pt);
      corr = fPtTrigCorr->effD(pt, eta, 0.);
      cout << nbin  << endl;
      yield*=corr;
      cout << " TrigEff Ups(2S) " << yield << " +/- " << yieldE << endl;
      fTrigEff_2->SetBinContent(nbin, yield); 
      fTrigEff_2->SetBinError(nbin, yieldE); 
      
    }
    
    for (unsigned int i = 0; i < fS10Vector.size(); ++i) {
      
      h = &(fS10Vector[i]);
      cout << h->GetName()  << endl;
      yield = h->GetMean(); 
      yieldE = h->GetMeanError(); 
      GetBinCenters(h->GetName(), eta, pt);
      nbin = fTrigEff_3->FindBin(eta, pt);
      corr = fPtTrigCorr->effD(pt, eta, 0.);
      cout << nbin  << endl;
      yield*=corr;
      cout << " TrigEff Ups(3S) " << yield << " +/- " << yieldE << endl;
      fTrigEff_3->SetBinContent(nbin, yield); 
      fTrigEff_3->SetBinError(nbin, yieldE); 
      
    }    
    
    makeCanvas(1);
    c1->Divide(3,1);
    c1->cd(1);
    fTrigEff->Draw("colz");
    c1->cd(2);
    fTrigEff_2->Draw("colz");  
    c1->cd(3);
    fTrigEff_3->Draw("colz");  
  }
  
  if ( mode == 2 ){
    for (unsigned int i = 0; i < fS8Vector.size(); ++i) {
      
      h = &(fS8Vector[i]);
      cout << h->GetName()  << endl;
      yield = h->GetMean(); 
      yieldE = h->GetMeanError(); 
      GetBinCenters(h->GetName(), eta, pt);
      nbin = fTrigEff->FindBin(eta, pt);
      corr = fPtTrigCorr->effD(pt, eta, 0.);
      cout << nbin  << endl;
      yield*=corr;
      cout << " TrigEff Ups " << yield << " +/- " << yieldE << endl;
      fTrigEff->SetBinContent(nbin, yield); 
      fTrigEff->SetBinError(nbin, yieldE); 
      
    }
    
    makeCanvas(1);
    fTrigEff->Draw("colz");  
    
  }
  
  ///////////
  //fTrigEff_3->Write();
  ////////////
}

void anaXS::GetPreSelEff(){
  double deno(-1.); double numa(-1.); double eff(-1);
  

  // Ups(1S) 
  for ( int iy = 1; iy <= fPreSelEff->GetNbinsX(); ++iy ){
    for ( int ipt = 1; ipt <= fPreSelEff->GetNbinsY(); ++ipt ){
      deno = fPreSelBefore->GetBinContent(iy,ipt);
      numa = fPreSelAfter->GetBinContent(iy,ipt);
      eff = numa/deno;
      if ( eff > 1. ) {
	cout << " !!! PreSelEff for UPS(1S) is higer than 1 !!! with " << eff  << endl;
	eff = 1;
      }
      fPreSelEff->SetBinContent(iy,ipt,eff);
    }
  }
  
  // Ups(2S) 
  for ( int iy = 1; iy <= fPreSelEff_2S->GetNbinsX(); ++iy ){
    for ( int ipt = 1; ipt <= fPreSelEff_2S->GetNbinsY(); ++ipt ){
      deno = fPreSelBefore_2S->GetBinContent(iy,ipt);
      numa = fPreSelAfter_2S->GetBinContent(iy,ipt);
      eff = numa/deno;
      if ( eff > 1. ) {
	cout << " !!! PreSelEff for UPS(2S) is higer than 1 !!! with " << eff  << endl;
	eff = 1;
      }
      fPreSelEff_2S->SetBinContent(iy,ipt,eff);
    }
  }  

  // Ups(3S) 
  for ( int iy = 1; iy <= fPreSelEff_3S->GetNbinsX(); ++iy ){
    for ( int ipt = 1; ipt <= fPreSelEff_3S->GetNbinsY(); ++ipt ){
      deno = fPreSelBefore_3S->GetBinContent(iy,ipt);
      numa = fPreSelAfter_3S->GetBinContent(iy,ipt);
      eff = numa/deno;
      if ( eff > 1. ) {
	cout << " !!! PreSelEff for UPS(3S) is higer than 1 !!! with " << eff  << endl;
	eff = 1;
      }
      fPreSelEff_3S->SetBinContent(iy,ipt,eff);
    }
  }  
  
  makeCanvas(1);
  c1->Divide(3,1);
  c1->cd(1);
  fPreSelEff->Draw("colz");
  c1->cd(2);
  fPreSelEff_2S->Draw("colz");
  c1->cd(3);
  fPreSelEff_3S->Draw("colz");  
}


void anaXS::Pull(int mode){
  
  
  if ( mode == 1 ){
    
    TH1D *h, *H, *h1, *h2, *h3, *h4;
    string fopt("LLE"); 
    
    double pt, eta; 
    int nbin(-99);
    double yield_1S(0.), yieldE_1S(0.), yield_1S_forE(0.), yield_2S_forE(0.); 
    double yield_2S(0.), yieldE_2S(0.);
    double yield_3S(0.), yieldE_3S(0.);
    double scale(1.033);
    
    double TrueYield_1S(0.), TrueYield_2S(0.);
    double pull_1s(0.), pull_2s(0.); double pull_err_1S(-99), pull_err_2S(-99);
    
    double sigma(-99);
    int status(0);
    const char* Status;
    gStyle->SetOptStat(0000000000000); 
    gStyle->SetOptFit(111111111);
    makeCanvas(1); 
    c1->Clear();
    
    TFile *f = new TFile("Pull_v1.root", "RECREATE");
    double content(-99);
    
    //for (unsigned int i = 0; i < fS1Vector.size(); ++i) {
    for (unsigned int i = 10; i < 12; ++i) {  
      
      c1->cd(1); shrinkPad(0.15, 0.26); 
      H = &(fS1Vector[i]);
      H->SetMinimum(0.); setTitles(H, "Mass_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
      
      TString frag(H->GetName()); 
      frag.ReplaceAll("s1:UpsilonMass,", "");
      
      h1 = new TH1D(Form("Pull_1S_%s",frag.Data() ), Form("Pull_1S_%s",frag.Data() ), 60, -6, 6);
      h2 = new TH1D(Form("Pull_2S_%s",frag.Data() ), Form("Pull_2S_%s",frag.Data() ), 60, -6, 6);
      
      h3 = new TH1D(Form("YieldCheck_1S_%s",frag.Data() ), Form("YieldCheck_1S_%s",frag.Data() ), 300, 0., 3000.);
      h4 = new TH1D(Form("PullError_1S_%s",frag.Data() ), Form("PullError_1S_%s",frag.Data() ), 100, 0., 100);
      
      for (unsigned int i = 0; i < 100; ++i) {
	setFunctionParameters(H, f13, 6);
	H->Fit(f13, fopt.c_str());
	status = 0;
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	cout << Status[0] << endl;
	if ( Status[0] == 'S' || Status[0] == 'P' ){
	  status = 1;
	} else if ( Status[0] == 'F' ){ 
	  status = -1;
	}
	cout << status << endl;
	///////////////
	if ( status == -1 || status == 0 ){
	  f13->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) , f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(7), f13->GetParameter(8), f13->GetParameter(9), f13->GetParameter(10));
	  H->Fit(f13, fopt.c_str());
	}
	
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	cout << Status[0] << endl;
	
	if ( Status[0] == 'S' || Status[0] == 'P' ){
	  status = 1;
	} else if ( Status[0] == 'F' ){ 
	  status = -1;
	}
	
	if ( Status[0] == 'S' || Status[0] == 'P' ) break;
      }
      
      // Ups 1S
      f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) );
      // Ups1S yield for Error Calculation
      yield_1S_forE  = f10->Integral(8.7,11.2)/H->GetBinWidth(1);
      cout << "yield_1S_forE = " << yield_1S_forE << endl;
      // Ups 2S
      f10->SetParameters( f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(7) );
      // Ups2S yield for Error Calculation
      yield_2S_forE  = f10->Integral(8.7,11.2)/H->GetBinWidth(1);
      cout << "yield_2S_forE = " << yield_2S_forE << endl;
      
      for (unsigned int i = 0; i < 10; ++i) {
	
	h = (TH1D*)H->Clone();
	
	for (int i=0;i<h->GetNbinsX();i++) {
	  content = h->GetBinContent(i);
	  sigma = h->GetBinError(i);
	  h->SetBinContent(i,content+gRandom->Gaus(0,sigma));
	}
	
	setFunctionParameters(h, f13, 6);
	h->Fit(f13, fopt.c_str());
	status = 0;
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	cout << Status[0] << endl;
	if ( Status[0] == 'S' || Status[0] == 'P' ){
	  status = 1;
	} else if ( Status[0] == 'F' ){ 
	  status = -1;
	}
	cout << status << endl;
	///////////////
	if ( status == -1 || status == 0 ){
	  f13->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) , f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(7), f13->GetParameter(8), f13->GetParameter(9), f13->GetParameter(10));
	  h->Fit(f13, fopt.c_str());
	}
      
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	cout << Status[0] << endl;
	
	if ( Status[0] == 'S' || Status[0] == 'P' ){
	  status = 1;
	} else if ( Status[0] == 'F' ){ 
	  status = -1;
	}
	///////////////
	cout << status << endl;
	// Ups 1S
	f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) );
	yield_1S  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_1S = TMath::Sqrt(yield_1S);
	cout << " Ups(1S) yield =  " << yield_1S << "+/-" << yieldE_1S  << endl;
	f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4)+f13->GetParError(4) );
	double yield_1S_e  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_1S = yield_1S - yield_1S_e;
	yieldE_1S = TMath::Abs(yieldE_1S);
	cout << " Ups(1S) yield =  " << yield_1S << "+/-" << yieldE_1S  << endl;
	if ( Status[0] == 'S' ) h3->Fill(yield_1S);
	
	// Ups 2S
	f10->SetParameters( f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(7) );
	yield_2S  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_2S = TMath::Sqrt(yield_2S);
	cout << " Ups(2S) yield =  " << yield_2S << "+/-" << yieldE_2S  << endl;
	f10->SetParameters( f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(7)+f13->GetParError(7) );
	double yield_2S_e  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_2S = yield_2S - yield_2S_e;
	yieldE_2S = TMath::Abs(yieldE_2S);
	cout << " Ups(2S) yield =  " << yield_2S << "+/-" << yieldE_2S  << endl;
	
	// Ups 3S
	f10->SetParameters( scale*f13->GetParameter(5), scale*f13->GetParameter(6), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(8) );
	yield_3S  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_3S = TMath::Sqrt(yield_3S);    
	
	
	GetBinCenters(h->GetName(), eta, pt);
	cout << " eta = " << eta << " pt = " << pt  << endl;
	nbin = fTrueYield_1S->FindBin(eta,pt);
	TrueYield_1S = fTrueYield_1S->GetBinContent(nbin);
	TrueYield_2S = fTrueYield_2S->GetBinContent(nbin);
	cout << "TrueYield_1S = " << TrueYield_1S << endl;
	cout << "TrueYield_2S = " << TrueYield_2S << endl;
	pull_err_1S = TMath::Sqrt((yieldE_1S*yieldE_1S)-yield_1S_forE);
	pull_err_2S = TMath::Sqrt((yieldE_2S*yieldE_2S)-yield_2S_forE);
	pull_1s = (yield_1S - TrueYield_1S)/pull_err_1S;
	pull_2s = (yield_2S*1.66 - TrueYield_2S)/pull_err_2S;
	cout << "pull_1s = " << pull_1s << endl;
	cout << "pull_2s = " << pull_2s << endl;
	if ( Status[0] == 'S' ) h4->Fill(pull_err_1S);
	if ( Status[0] == 'S' ) h1->Fill(pull_1s);
	if ( Status[0] == 'S' ) h2->Fill(pull_2s);
	
	
	c1->Modified();
	c1->Update();
	
	h->Reset();
	
      }
      
      h1->Write();
      h2->Write();
      h3->Write();
      h4->Write();
      
    }
    
  }
  
  if ( mode == 2 ){
    
    TH1D *h, *H, *h1, *h2, *h3, *h4;
    string fopt("LLE"); 
    
    double pt, eta; 
    int nbin(-99);
    double yield_1S(0.), yieldE_1S(0.);
    double yield_2S(0.), yieldE_2S(0.);
    double yield_3S(0.), yieldE_3S(0.);
    double scale(1.033);
    
    double TrueYield_1S(0.), TrueYield_2S(0.);
    double pull_1s(0.), pull_2s(0.);
    double integral_bg(-99), integral_1S(-99), integral_2S(-99), integral_3S(-99), integral_All(-99); 
    double sigma(-99); double pull_err_1S(-99), pull_err_2S(-99);
    int status(0);
    const char* Status;
    gStyle->SetOptStat(0000000000000); 
    gStyle->SetOptFit(111111111);
    makeCanvas(1); 
    c1->Clear();
    
    TFile *f = new TFile("Pull_v2.root", "RECREATE");
    double content(-99);
    
    
    //for (unsigned int i = 0; i < fS1Vector.size(); ++i) {
    for (unsigned int i = 10; i < 12; ++i) {  
      
      c1->cd(1); shrinkPad(0.15, 0.26); 
      H = &(fS1Vector[i]);
      H->SetMinimum(0.); setTitles(H, "Mass_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
      
      TString frag(H->GetName()); 
      frag.ReplaceAll("s1:UpsilonMass,", "");
      
      h1 = new TH1D(Form("Pull_1S_%s",frag.Data() ), Form("Pull_1S_%s",frag.Data() ), 60, -6, 6);
      h2 = new TH1D(Form("Pull_2S_%s",frag.Data() ), Form("Pull_2S_%s",frag.Data() ), 60, -6, 6);
      
      h3 = new TH1D(Form("YieldCheck_1S_%s",frag.Data() ), Form("YieldCheck_1S_%s",frag.Data() ), 300, 0., 3000.);
      h4 = new TH1D(Form("PullError_1S_%s",frag.Data() ), Form("PullError_1S_%s",frag.Data() ), 100, 0., 10.);
      
      for (unsigned int i = 0; i < 100; ++i) {
	setFunctionParameters(H, f13, 6);
	H->Fit(f13, fopt.c_str());
	status = 0;
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	if ( Status[0] == 'S' ){
	  status = 1;
	} else if ( Status[0] == 'F' || Status[0] == 'P' ){   
	  status = -1;
			     }
	if ( status == -1 || status == 0 ){
	  f13->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) , f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(7), f13->GetParameter(8), f13->GetParameter(9), f13->GetParameter(10));
	  H->Fit(f13, fopt.c_str());
	}
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	cout << Status[0] << endl;
	
	if ( Status[0] == 'S' ){
	  status = 1;
	} else if ( Status[0] == 'F' || Status[0] == 'P' ){
	  status = -1;
	}
	
	c1->Modified();
	c1->Update();
	
	if ( Status[0] == 'S' || Status[0] == 'P' ) break;
      }
      
      double p0(-99), p1(-99), p2(-99), p3(-99), p4(-99), p5(-99), p6(-99), p7(-99), p8(-99), p9(-99), p10(-99);
      p0 = f13->GetParameter(0); p1 = f13->GetParameter(1); p2 = f13->GetParameter(2); p3 = f13->GetParameter(3); 
      p4 = f13->GetParameter(4); p5 = f13->GetParameter(5); p6 = f13->GetParameter(6); p7 = f13->GetParameter(7);
      p8 = f13->GetParameter(8); p9 = f13->GetParameter(9); p10 = f13->GetParameter(10); 
      double e4(-99), e7(-99), e8(-99); 
      e4 = f13->GetParError(4); e7 = f13->GetParError(7); e8 = f13->GetParError(8); 
      
      for (unsigned int i = 0; i < 10; ++i) {
	
	cout << " i = " << i << endl;
	h = (TH1D*)H->Clone();
	
	h->Reset();
	
	/*
	f10->SetParameters( p0, p1, p2, p3, p4 );
	integral_1S = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	h->FillRandom("f10", integral_1S);
	cout << "Ups1S Yield = "  << integral_1S << endl;
	f10->SetParameters( p5, p6, p2, p3, p7 );
	integral_2S = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	h->FillRandom("f10", integral_2S);
	cout << "Ups2S Yield = "  << integral_2S << endl;
	f10->SetParameters( scale*p5, scale*p6, p2, p3, p8 );
	integral_3S = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	h->FillRandom("f10", integral_3S);	
	cout << "Ups3S Yield = "  << integral_3S << endl;
	f0->SetParameters( p9, p10);
	integral_bg = f0->Integral(8.7,11.2)/h->GetBinWidth(1);
	h->FillRandom("f0", integral_bg);		
	*/
	
	f13->SetParameters( p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10);
	integral_All = f13->Integral(8.7,11.2)/h->GetBinWidth(1);
	h->FillRandom("f13", integral_All);
	f10->SetParameters( p0, p1, p2, p3, p4 );
	integral_1S = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	f10->SetParameters( p5, p6, p2, p3, p7 );
	integral_2S = f10->Integral(8.7,11.2)/h->GetBinWidth(1);	
	
	
	setFunctionParameters(h, f13, 6);
	h->Fit(f13, fopt.c_str());
	status = 0;
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	cout << Status[0] << endl;
	if ( Status[0] == 'S' || Status[0] == 'P' ){
	  status = 1;
	} else if ( Status[0] == 'F' ){ 
	  status = -1;
	}
	cout << status << endl;
	///////////////
	if ( status == -1 || status == 0 ){
	  f13->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) , f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(7), f13->GetParameter(8), f13->GetParameter(9), f13->GetParameter(10));
	  h->Fit(f13, fopt.c_str());
	}
      
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	cout << Status[0] << endl;
	
	if ( Status[0] == 'S' || Status[0] == 'P' ){
	  status = 1;
	} else if ( Status[0] == 'F' ){ 
	  status = -1;
	}
	///////////////
	cout << status << endl;
	// Ups 1S
	f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) );
	yield_1S  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_1S = TMath::Sqrt(yield_1S);
	cout << " Ups(1S) yield =  " << yield_1S << "+/-" << yieldE_1S  << endl;
	f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4)+f13->GetParError(4) );
	double yield_1S_e  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_1S = yield_1S - yield_1S_e;
	yieldE_1S = TMath::Abs(yieldE_1S);
	cout << " Ups(1S) yield =  " << yield_1S << "+/-" << yieldE_1S  << endl;
	double check = (f13->GetParError(4)/f13->GetParameter(4))*yield_1S;
	cout << " check = "  << check << " yieldE_1S = "  << yieldE_1S  << endl;
	if ( Status[0] == 'S' ) h3->Fill(yield_1S);
	
	
	// Ups 2S
	f10->SetParameters( f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(7) );
	yield_2S  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_2S = TMath::Sqrt(yield_2S);
	cout << " Ups(2S) yield =  " << yield_2S << "+/-" << yieldE_2S  << endl;
	f10->SetParameters( f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(7)+f13->GetParError(7) );
	double yield_2S_e  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_2S = yield_2S - yield_2S_e;
	yieldE_2S = TMath::Abs(yieldE_2S);
	cout << " Ups(2S) yield =  " << yield_2S << "+/-" << yieldE_2S  << endl;
	
	// Ups 3S
	f10->SetParameters( scale*f13->GetParameter(5), scale*f13->GetParameter(6), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(8) );
	yield_3S  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_3S = TMath::Sqrt(yield_3S);    
	
	
	GetBinCenters(h->GetName(), eta, pt);
	cout << " eta = " << eta << " pt = " << pt  << endl;
	nbin = fTrueYield_1S->FindBin(eta,pt);
	TrueYield_1S = fTrueYield_1S->GetBinContent(nbin);
	TrueYield_2S = fTrueYield_2S->GetBinContent(nbin);
	cout << "TrueYield_1S = " << TrueYield_1S << endl;
	cout << "TrueYield_2S = " << TrueYield_2S << endl;
	pull_err_1S = TMath::Sqrt((yieldE_1S*yieldE_1S)-integral_1S);
	pull_err_2S = TMath::Sqrt((yieldE_2S*yieldE_2S)-integral_2S);
	pull_1s = (yield_1S - TrueYield_1S)/pull_err_1S;
	pull_2s = (yield_2S*1.66 - TrueYield_2S)/pull_err_2S;
	cout << "pull_1s = " << pull_1s << endl;
	cout << "pull_2s = " << pull_2s << endl;
	if ( Status[0] == 'S' ) h4->Fill(pull_err_1S);
	if ( Status[0] == 'S' ) h1->Fill(pull_1s);
	if ( Status[0] == 'S' ) h2->Fill(pull_2s);
	
	
	c1->Modified();
	c1->Update();
	
	h->Reset();
	
      }
      
      h1->Write();
      h2->Write();
      h3->Write();
      h4->Write();
      
    }
    
  }
  
}

void anaXS::FITUpsilon(int mode){
  
  if ( mode == 1 ){
    
    int PRINT(1); 
    double PRINTX(0.5);
    
    
    
    TH1D *h; 
    
    //  string fopt("LLIEMQ"); 
    string fopt("LLE"); 
    
    double pt, eta; 
    
    double yield_1S(0.), yieldE_1S(0.);
    double yield_2S(0.), yieldE_2S(0.);
    double yield_3S(0.), yieldE_3S(0.);
    double alpha(0.), n(0.);
    double scale(1.033);
    
    int    nbin;
    int fitted(0);
    
    int status(0);
    const char* Status;
    
    gStyle->SetOptStat(0000000000000); 
    gStyle->SetOptFit(00000000000000);
    //gStyle->SetOptStat(111111111);
    gStyle->SetOptFit(111111111);
    makeCanvas(1); 
    c1->Clear();
    
    TFile *f = new TFile("Yield.root", "RECREATE");
    
    
    for (unsigned int i = 0; i < fS1Vector.size(); ++i) {
      
      // -- positive charge
      c1->cd(1); shrinkPad(0.15, 0.26); 
      h = &(fS1Vector[i]);
      h->SetMinimum(0.); setTitles(h, "Mass_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
      if (h->GetSumOfWeights() > 10.) {
	setFunctionParameters(h, f13, 6);
	h->Fit(f13, fopt.c_str());
	status = 0;
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	cout << Status[0] << endl;
	if ( Status[0] == 'S' || Status[0] == 'P' ){
	  status = 1;
	} else if ( Status[0] == 'F' ){ 
	  status = -1;
	}
	cout << status << endl;
	///////////////
	if ( status == -1 || status == 0 ){
	  f13->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) , f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(7), f13->GetParameter(8), f13->GetParameter(9), f13->GetParameter(10));
	  h->Fit(f13, fopt.c_str());
	}
	
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	cout << Status[0] << endl;
	
	if ( Status[0] == 'S' || Status[0] == 'P' ){
	  status = 1;
	} else if ( Status[0] == 'F' ){ 
	  status = -1;
	}
	///////////////
	cout << status << endl;
	// Ups 1S
	f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) );
	yield_1S  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4)+f13->GetParError(4) );
	double yield_1S_e  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_1S = yield_1S - yield_1S_e;
	yieldE_1S = TMath::Abs(yieldE_1S);
	cout << " Ups(1S) yield =  " << yield_1S << "+/-" << yieldE_1S  << endl;
	alpha = f13->GetParameter(2);
	n = f13->GetParameter(3);
	
	// Ups 2S
	f10->SetParameters( f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(7) );
	yield_2S  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	f10->SetParameters( f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(7)+f13->GetParError(7) );
	double yield_2S_e  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_2S = yield_2S - yield_2S_e;
	yieldE_2S = TMath::Abs(yieldE_2S);
	cout << " Ups(2S) yield =  " << yield_2S << "+/-" << yieldE_2S  << endl;
	// Ups 3S
	
	f10->SetParameters( scale*f13->GetParameter(5), scale*f13->GetParameter(6), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(8) );
	yield_3S  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	f10->SetParameters( f13->GetParameter(5), f13->GetParameter(6), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(8)+f13->GetParError(8) );
	double yield_3S_e  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_3S = yield_3S - yield_3S_e;
	yieldE_3S = TMath::Abs(yieldE_3S);
	cout << " Ups(3S) yield =  " << yield_3S << "+/-" << yieldE_3S  << endl;
	
		
	///
	fitted = 1;
      } else {
	h->Draw();
	yield_1S  = h->GetSumOfWeights();
	yieldE_1S = TMath::Sqrt(h->GetSumOfWeights());
	fitted = 0;
      }
      if (PRINT) tl->DrawLatex(PRINTX, 0.40, Form("Ups(1S) Yield = %4.1f +/- %4.1f", yield_1S,yieldE_1S)); 
      if (PRINT) tl->DrawLatex(PRINTX, 0.33, Form("Ups(2S) Yield = %4.1f +/- %4.1f", yield_2S,yieldE_2S));
      if (PRINT) tl->DrawLatex(PRINTX, 0.26, Form("Ups(3S) Yield = %4.1f +/- %4.1f", yield_3S,yieldE_3S));
      cout << " --> " << h->GetName() << ", Ups(1S) Yield = " << yield_1S << "+/-" << yieldE_1S << endl;
      cout << " --> " << h->GetName() << ", Ups(2S) Yield = " << yield_2S << "+/-" << yieldE_2S << endl;
      cout << " --> " << h->GetName() << ", Ups(3S) Yield = " << yield_3S << "+/-" << yieldE_3S << endl;
      GetBinCenters(h->GetName(), eta, pt); 
      cout << "/////" << endl;
      cout << eta << "  " << pt << endl;
      nbin = fS1Yield->FindBin(eta, pt); 
      cout << nbin  << endl;
      fS1Yield->SetBinContent(nbin, yield_1S); 
      fS1Yield->SetBinError(nbin, yieldE_1S);
      fS2Yield->SetBinContent(nbin, yield_2S); 
      fS2Yield->SetBinError(nbin, yieldE_2S);
      fS3Yield->SetBinContent(nbin, yield_3S); 
      fS3Yield->SetBinError(nbin, yieldE_3S);
      
      falpha->SetBinContent(nbin, alpha);
      fn->SetBinContent(nbin, n);
      cout << "/////" << endl;
      c1->Modified();
      c1->Update();
      
      
      TString frag(h->GetName()); 
      frag.ReplaceAll("s3:mmbar,", ""); 
      frag.ReplaceAll(",Q1", ""); 
      frag.ReplaceAll(".", "_");
      frag.ReplaceAll(",", "_");
      frag.ReplaceAll(":", "_");
      
      c1->SaveAs(Form("%s/massfits_%s_%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data())); 
    }
    
    fS1Yield->Write();
    fS2Yield->Write();
    fS3Yield->Write();
    falpha->Write();
    fn->Write();
    c1->Clear();
  }
  
  if ( mode == 2  ){
    int PRINT(1); 
    double PRINTX(0.5);
    
    TH1D *h; 
    
    //  string fopt("LLIEMQ"); 
    string fopt("LLE"); 
    
    double pt, eta; 
    
    double yield_1S(0.), yieldE_1S(0.);
    double alpha(0.), n(0.);
    
    int    nbin;
    int fitted(0);
    
    int status(0);
    const char* Status;
    
    gStyle->SetOptStat(0000000000000); 
    gStyle->SetOptFit(00000000000000);	
    //gStyle->SetOptStat(111111111);
    gStyle->SetOptFit(111111111);
    makeCanvas(1); 
    c1->Clear();
    
    TFile *f = new TFile("FitParameters.root", "RECREATE");
    
    
    for (unsigned int i = 0; i < fS1Vector.size(); ++i) {
      
      // -- positive charge
      c1->cd(1); shrinkPad(0.15, 0.08); 
      h = &(fS1Vector[i]);
      h->SetMinimum(0.);  setTitles(h, "Mass_{#mu #mu} [GeV]", "Entries/Bin", 0.07, 0.75, 0.58, 0.07);
      if (h->GetSumOfWeights() > 10.) {
	
	for (unsigned int i = 0; i < h->GetNbinsX(); ++i) {
	  if ( h->GetBinContent(i) == 0  ){
	    h->GetBinError(i,1);
	  }
	}
	
	setFunctionParameters(h, f6, 5); 
	h->Fit(f6, fopt.c_str());
	status = 0;
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	cout << Status[0] << endl;
	if ( Status[0] == 'S' || Status[0] == 'P' ){
	  status = 1;
	} else if ( Status[0] == 'F' ){ 
	  status = -1;
	}
	cout << status << endl;
	///////////////
	if ( status == -1 || status == 0 ){
	  f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4), f6->GetParameter(5), f6->GetParameter(6));
	  h->Fit(f6, fopt.c_str());
	}
	cout << gMinuit->fCstatu.Data() << endl;
	Status = gMinuit->fCstatu.Data();
	cout << Status[0] << endl;
	
	if ( Status[0] == 'S' || Status[0] == 'P' ){
	  status = 1;
	} else if ( Status[0] == 'F' ){ 
	  status = -1;
	}
	///////////////
	cout << status << endl;
	yield_1S  = f6->Integral(8.7,11.2)/h->GetBinWidth(1);
	yieldE_1S = TMath::Sqrt(yield_1S);
	alpha = f6->GetParameter(2);
	n = f6->GetParameter(3);
	fitted = 1;
      } else {
	h->Draw();
	yield_1S  = h->GetSumOfWeights();
	yieldE_1S = TMath::Sqrt(h->GetSumOfWeights());
	fitted = 0;
      }
      if (PRINT) tl->DrawLatex(PRINTX, 0.40, Form("Ups Yield = %4.1f +/- %4.1f", yield_1S,yieldE_1S)); 
      cout << " --> " << h->GetName() << ", Ups Yield = " << yield_1S << "+/-" << yieldE_1S << endl;
      GetBinCenters(h->GetName(), eta, pt); 
      cout << "/////" << endl;
      cout << eta << "  " << pt << endl;
      nbin = fS1Yield->FindBin(eta, pt); 
      cout << nbin  << endl;
      fS1Yield->SetBinContent(nbin, yield_1S); 
      fS1Yield->SetBinError(nbin, yieldE_1S);
      falpha->SetBinContent(nbin, alpha);
      fn->SetBinContent(nbin, n);
      cout << "/////" << endl;
      c1->Modified();
      c1->Update();
      
      
      TString frag(h->GetName()); 
      frag.ReplaceAll("s3:mmbar,", ""); 
      frag.ReplaceAll(",Q1", ""); 
      frag.ReplaceAll(".", "_");
      frag.ReplaceAll(",", "_");
      frag.ReplaceAll(":", "_");
      
      c1->SaveAs(Form("%s/massfits_%s_%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data())); 
    }
    
    fn->Write();
    falpha->Write();
    fS1Yield->Write();
    c1->Clear();
    
    
  }
}
// ----------------------------------------------------------------------
void anaXS::fitUpsilon(int mode) {  


  if ( mode == 0 ){   // fit All three Upsilons with pol + Gaussians

  //    o function for fitting
  //    o do in there something about low statistics cases

	int PRINT(1); 
	double PRINTX(0.3);

	TH1D *h; 

	  //  string fopt("LLIEMQ"); 
	string fopt("LLE"); 

  	double pt, eta; 
  	int    Q; 
	
  	double yield, yieldE; 
  	int    nbin; 

  	gStyle->SetOptStat(PRINT); 
  	gStyle->SetOptFit(PRINT); 
  	makeCanvas(1); 
  	c1->Clear();
  	c1->Divide(3,1); 
  	//  for (int i = 1; i < 4; ++i) {
  	for (unsigned int i = 0; i < fS1VectorPos.size(); ++i) {

  		  // -- positive charge
    		c1->cd(1); shrinkPad(0.15, 0.26); 
    		h = &(fS1VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f2, 2); 
    		  h->Fit(f2, fopt.c_str());
    		  yield  = f2->GetParameter(0)/h->GetBinWidth(1);
    		  yieldE = f2->GetParError(0)/h->GetBinWidth(1);
    		} else {
    		  h->Draw();
    		  yield  = h->GetSumOfWeights();
    		  yieldE = TMath::Sqrt(h->GetSumOfWeights());
    		}
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS1YieldPos->FindBin(eta, pt); 
    		fS1YieldPos->SetBinContent(nbin, yield); 
    		fS1YieldPos->SetBinError(nbin, yieldE); 
		
    		c1->Modified();
    		c1->Update();
		
    		c1->cd(2); shrinkPad(0.15, 0.26); 
    		h = &(fS2VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f2, 2); 
    		  h->Fit(f2, fopt.c_str());
    		  yield  = f2->GetParameter(0)/h->GetBinWidth(1);
    		  yieldE = f2->GetParError(0)/h->GetBinWidth(1);
    		} else {
    		  h->Draw();
    		  yield  = h->GetSumOfWeights();
    		  yieldE = TMath::Sqrt(h->GetSumOfWeights());
    		}
		
 		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS2YieldPos->FindBin(eta, pt); 
    		fS2YieldPos->SetBinContent(nbin, yield); 
    		fS2YieldPos->SetBinError(nbin, yieldE); 
		
    		c1->Modified();
    		c1->Update();
		
    		c1->cd(3);  shrinkPad(0.15, 0.26); 
    		h = &(fS3VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f2, 2); 
     	  	  h->Fit(f2, fopt.c_str());
      	  	  yield  = f2->GetParameter(0)/h->GetBinWidth(1);
      	  	  yieldE = f2->GetParError(0)/h->GetBinWidth(1);
    		} else {
      		  h->Draw();
      		  yield  = h->GetSumOfWeights();
      		  yieldE = TMath::Sqrt(h->GetSumOfWeights());
    		}
	
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS3YieldPos->FindBin(eta, pt); 
    		double mm = fS1YieldPos->GetBinContent(nbin); 
    		double mmE= fS1YieldPos->GetBinError(nbin); 
    		fS3YieldPos->SetBinContent(nbin, yield+mm); 
    		fS3YieldPos->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 
		
   		TString frag(h->GetName()); 
    		frag.ReplaceAll("s3:mmbar,", ""); 
    		frag.ReplaceAll(",Q1", ""); 
    		frag.ReplaceAll(".", ":"); 
	
    		c1->SaveAs(Form("%s/massfits-pos-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data())); 
  	}
  	
  
  		c1->Clear();
  		c1->Divide(3,1); 
  	for (unsigned int i = 0; i < fS1VectorNeg.size(); ++i) {

  		  // -- negative charge
    		c1->cd(1); shrinkPad(0.15, 0.26);  
    		h = &(fS1VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f2, 2); 
    		  h->Fit(f2, fopt.c_str());
    		  yield  = f2->GetParameter(0)/h->GetBinWidth(1);
    		  yieldE = f2->GetParError(0)/h->GetBinWidth(1);
    		} else {
    	  	  h->Draw();
    	  	  yield  = h->GetSumOfWeights();
    	  	  yieldE = TMath::Sqrt(h->GetSumOfWeights());
    		}
	
		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS1YieldNeg->FindBin(eta, pt); 
    		fS1YieldNeg->SetBinContent(nbin, yield); 
    		fS1YieldNeg->SetBinError(nbin, yieldE); 
	
    		c1->Modified();
    		c1->Update();
	
    		c1->cd(2); shrinkPad(0.15, 0.26); 
    		h = &(fS2VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    	  	  setFunctionParameters(h, f2, 2); 
    	  	  h->Fit(f2, fopt.c_str());
    	    	  yield  = f2->GetParameter(0)/h->GetBinWidth(1);
    	   	  yieldE = f2->GetParError(0)/h->GetBinWidth(1);
    	 	} else {
    	  	  h->Draw();
    	   	  yield  = h->GetSumOfWeights();
    	   	  yieldE = TMath::Sqrt(h->GetSumOfWeights());
    		}
	
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS2YieldNeg->FindBin(eta, pt); 
    		fS2YieldNeg->SetBinContent(nbin, yield); 
    		fS2YieldNeg->SetBinError(nbin, yieldE); 
	
    		c1->cd(3);  shrinkPad(0.15, 0.26); 
    		h = &(fS3VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f2, 2); 
    	  	  h->Fit(f2, fopt.c_str());
    	  	  yield  = f2->GetParameter(0)/h->GetBinWidth(1);
    	  	  yieldE = f2->GetParError(0)/h->GetBinWidth(1);
    		} else {
    	  	  h->Draw();
    	  	  yield  = h->GetSumOfWeights();
    	  	  yieldE = TMath::Sqrt(h->GetSumOfWeights());
    		}

    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS3YieldNeg->FindBin(eta, pt); 
    		double mm = fS1YieldNeg->GetBinContent(nbin); 
    		double mmE= fS1YieldNeg->GetBinError(nbin); 
    		fS3YieldNeg->SetBinContent(nbin, yield+mm); 
    		fS3YieldNeg->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 
		
		
    		c1->Modified();
    		c1->Update();
	
    		TString frag(h->GetName()); 
    		frag.ReplaceAll("s3:mmbar,", ""); 
    		frag.ReplaceAll(",Q-1", ""); 
    		frag.ReplaceAll(".", ":"); 
		
    		c1->SaveAs(Form("%s/massfits-neg-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data())); 

  	}

  }
  
  if ( mode == 1 ){   // fit Upsilon(1S) with pol + CB

  //    o function for fitting
  //    o do in there something about low statistics cases

	int PRINT(1); 
	double PRINTX(0.3);

	TH1D *h; 

	  //  string fopt("LLIEMQ"); 
	string fopt("LLE"); 

  	double pt, eta; 
  	int    Q; 
	
  	double yield, yieldE, width, widthE; 
  	int    nbin;
	int fitted(0);
	 
  	int status(0);
	const char* Status;
		
	TH1D *hChisq_mm    = new TH1D("hChisq/ndof_mm"   , "hChisq/ndof_mm"   , 40, -20.,20.);
	TH1D *hChisq_mt    = new TH1D("hChisq/ndof_mt"   , "hChisq/ndof_mt"   , 40, -20.,20.); 
	TH1D *hChisq_mmbar = new TH1D("hChisq/ndof_mmbar", "hChisq/ndof_mmbar", 40, -20.,20.);
	
 	TH1D *hSigma_mm    = new TH1D("hSigma_mm"   , "hSigma_mm"   , 28, -0.14,0.14);
	TH1D *hSigma_mt    = new TH1D("hSigma_mt"   , "hSigma_mt"   , 28, -0.14,0.14); 
	TH1D *hSigma_mmbar = new TH1D("hSigma_mmbar", "hSigma_mmbar", 28, -0.14,0.14);
	
 	TH1D *hAlpha_mm    = new TH1D("hAlpha_mm"   , "hAlpha_mm"   , 30, -3.,3.);
	TH1D *hAlpha_mt    = new TH1D("hAlpha_mt"   , "hAlpha_mt"   , 30, -3.,3.); 
	TH1D *hAlpha_mmbar = new TH1D("hAlpha_mmbar", "hAlpha_mmbar", 30, -3.,3.);		
	
 	TH1D *hn_mm        = new TH1D("hn_mm"   , "hn_mm"   , 40, -4.,4.);
	TH1D *hn_mt        = new TH1D("hn_mt"   , "hn_mt"   , 40, -4.,4.); 
	TH1D *hn_mmbar     = new TH1D("hn_mmbar", "hn_mmbar", 40, -4.,4.);
	
  	gStyle->SetOptStat(PRINT); 
  	gStyle->SetOptFit(PRINT); 
  	makeCanvas(1); 
  	c1->Clear();
  	c1->Divide(3,1);
	 
  	for (unsigned int i = 0; i < fS1VectorPos.size(); ++i) {

  		  // -- positive charge
    		c1->cd(1); shrinkPad(0.15, 0.26); 
    		h = &(fS1VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f6, 5); 
    		  h->Fit(f6, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  ///////////////
		  if ( status == -1 || status == 0 ){
			f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
			f6->GetParameter(5), f6->GetParameter(6));
			h->Fit(f6, fopt.c_str());
		  }
			
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
			
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	    	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  /////////////////////////	
		  if ( status == 1 ){
			hChisq_mm->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
			hSigma_mm->Fill(status*f6->GetParameter(1));
			hAlpha_mm->Fill(status*f6->GetParameter(2));
			hn_mm->Fill(status*f6->GetParameter(3));
		  }
		  f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
    		  yield  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
    		  yieldE = TMath::Sqrt(yield);
		  fitted = 1;
    		} else {
    		  h->Draw();
    		  yield  = h->GetSumOfWeights();
    		  yieldE = TMath::Sqrt(h->GetSumOfWeights());
		  fitted = 0;
    		}
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS1YieldPos->FindBin(eta, pt); 
    		fS1YieldPos->SetBinContent(nbin, yield); 
    		fS1YieldPos->SetBinError(nbin, yieldE); 
		
    		c1->Modified();
    		c1->Update();
		
    		c1->cd(2); shrinkPad(0.15, 0.26); 
    		h = &(fS2VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f6, 5); 
    		  h->Fit(f6, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	     	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  ///////////////
		  if ( status == -1 || status == 0 ){
			f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
			f6->GetParameter(5), f6->GetParameter(6));
			h->Fit(f6, fopt.c_str());
		  }
		
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  /////////////////////////
		  if ( status == 1 ){
			hChisq_mt->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
			hSigma_mt->Fill(status*f6->GetParameter(1));
			hAlpha_mt->Fill(status*f6->GetParameter(2));
			hn_mt->Fill(status*f6->GetParameter(3));
		  }
		  f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
		  yield  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
      		  yieldE = TMath::Sqrt(yield);
      		  width  = f6->GetParameter(1); 
      		  widthE = f6->GetParError(1); 
      		  fitted = 1; 
    		} else {
    		  h->Draw();
    		  yield  = h->GetSumOfWeights();
    		  yieldE = TMath::Sqrt(h->GetSumOfWeights());
		  fitted = 0;
    		}
		
 		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS2YieldPos->FindBin(eta, pt); 
    		fS2YieldPos->SetBinContent(nbin, yield); 
    		fS2YieldPos->SetBinError(nbin, yieldE); 
		
    		c1->Modified();
    		c1->Update();
		
    		c1->cd(3);  shrinkPad(0.15, 0.26); 
    		h = &(fS3VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f6, 5); 
     	  	  h->Fit(f6, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  ///////////////
		  if ( status == -1 || status == 0 ){
			f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
			f6->GetParameter(5), f6->GetParameter(6));
			h->Fit(f6, fopt.c_str());
		  }
		
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  /////////////////////////
		  if ( status == 1 ){
			hChisq_mmbar->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
			hSigma_mmbar->Fill(status*f6->GetParameter(1));
			hAlpha_mmbar->Fill(status*f6->GetParameter(2));
			hn_mmbar->Fill(status*f6->GetParameter(3));
		  }		
		  f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
		  yield  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
      		  yieldE = TMath::Sqrt(yield);		
      		  fitted = 1; 
    		} else {
      		  h->Draw();
      		  yield  = h->GetSumOfWeights();
      		  yieldE = TMath::Sqrt(h->GetSumOfWeights());
		  fitted = 0;
    		}
	
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS3YieldPos->FindBin(eta, pt); 
    		double mm = fS1YieldPos->GetBinContent(nbin); 
    		double mmE= fS1YieldPos->GetBinError(nbin); 
    		fS3YieldPos->SetBinContent(nbin, yield+mm); 
    		fS3YieldPos->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 
		
   		TString frag(h->GetName()); 
    		frag.ReplaceAll("s3:mmbar,", ""); 
    		frag.ReplaceAll(",Q1", ""); 
    		frag.ReplaceAll(".", ":"); 
	
    		c1->SaveAs(Form("%s/massfits-pos-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data())); 
  	}
  	
  
  		c1->Clear();
  		c1->Divide(3,1); 
  	for (unsigned int i = 0; i < fS1VectorNeg.size(); ++i) {

  		  // -- negative charge
    		c1->cd(1); shrinkPad(0.15, 0.26);  
    		h = &(fS1VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f6, 5); 
    		  h->Fit(f6, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  ///////////////
		  if ( status == -1 || status == 0 ){
			f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
			f6->GetParameter(5), f6->GetParameter(6));
			h->Fit(f6, fopt.c_str());
		  }
		
	  	  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  /////////////////////////
		  if ( status == 1 ){
			hChisq_mm->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
			hSigma_mm->Fill(status*f6->GetParameter(1));
			hAlpha_mm->Fill(status*f6->GetParameter(2));
			hn_mm->Fill(status*f6->GetParameter(3));
		  }		
		  f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
		  yield  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
      		  yieldE = TMath::Sqrt(yield);		
    		} else {
    	  	  h->Draw();
    	  	  yield  = h->GetSumOfWeights();
    	  	  yieldE = TMath::Sqrt(h->GetSumOfWeights());
    		}
	
		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS1YieldNeg->FindBin(eta, pt); 
    		fS1YieldNeg->SetBinContent(nbin, yield); 
    		fS1YieldNeg->SetBinError(nbin, yieldE); 
	
    		c1->Modified();
    		c1->Update();
	
    		c1->cd(2); shrinkPad(0.15, 0.26); 
    		h = &(fS2VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    	  	  setFunctionParameters(h, f6, 5); 
    	  	  h->Fit(f6, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
		}
		cout << status << endl;
		///////////////
		if ( status == -1 || status == 0 ){
			f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
			f6->GetParameter(5), f6->GetParameter(6));
			h->Fit(f6, fopt.c_str());
		}
		
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		/////////////////////////
		if ( status == 1 ){
			hChisq_mt->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
			hSigma_mt->Fill(status*f6->GetParameter(1));
			hAlpha_mt->Fill(status*f6->GetParameter(2));
			hn_mt->Fill(status*f6->GetParameter(3));
		}		
		f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
		yield  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);		
      		width  = f6->GetParameter(1); 
      		widthE = f6->GetParError(1); 
		fitted = 1;  
    	 	} else {
    	  	  h->Draw();
    	   	  yield  = h->GetSumOfWeights();
    	   	  yieldE = TMath::Sqrt(h->GetSumOfWeights());
		  fitted = 0;
    		}
	
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS2YieldNeg->FindBin(eta, pt); 
    		fS2YieldNeg->SetBinContent(nbin, yield); 
    		fS2YieldNeg->SetBinError(nbin, yieldE); 
	
    		c1->cd(3);  shrinkPad(0.15, 0.26); 
    		h = &(fS3VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f6, 5); 
    	  	  h->Fit(f6, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  ///////////////
		  if ( status == -1 || status == 0 ){
			f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
			f6->GetParameter(5), f6->GetParameter(6));
			h->Fit(f6, fopt.c_str());
		  }
		
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  /////////////////////////
		  if ( status == 1 ){
			hChisq_mmbar->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
			hSigma_mmbar->Fill(status*f6->GetParameter(1));
			hAlpha_mmbar->Fill(status*f6->GetParameter(2));
			hn_mmbar->Fill(status*f6->GetParameter(3));
		  }		
		  f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
		  yield  = f10->Integral(8.7,11.2)/h->GetBinWidth(1);
      		  yieldE = TMath::Sqrt(yield);		
      		  width  = f6->GetParameter(1); 
      		  widthE = f6->GetParError(1); 
		  fitted = 1;
    		} else {
    	  	  h->Draw();
    	  	  yield  = h->GetSumOfWeights();
    	  	  yieldE = TMath::Sqrt(h->GetSumOfWeights());
		  fitted = 0;
    		}

    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS3YieldNeg->FindBin(eta, pt); 
    		double mm = fS1YieldNeg->GetBinContent(nbin); 
    		double mmE= fS1YieldNeg->GetBinError(nbin); 
    		fS3YieldNeg->SetBinContent(nbin, yield+mm); 
    		fS3YieldNeg->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 
		
		
    		c1->Modified();
    		c1->Update();
	
    		TString frag(h->GetName()); 
    		frag.ReplaceAll("s3:mmbar,", ""); 
    		frag.ReplaceAll(",Q-1", ""); 
    		frag.ReplaceAll(".", ":"); 
		
    		c1->SaveAs(Form("%s/massfits-neg-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data())); 

  	}

	TCanvas *c2 = new TCanvas("c2", "c2");
  	c2->Clear();
        c2->Divide(3,1);
        gStyle->SetOptStat(1111111111);
        c2->cd(1);  shrinkPad(0.1, 0.2);
        hChisq_mm->Draw();
        c2->cd(2);  shrinkPad(0.1, 0.2);	
        hChisq_mt->Draw();
        c2->cd(3);  shrinkPad(0.1, 0.2);	    
        hChisq_mmbar->Draw();
        c2->SaveAs("chi2_polCB.eps");
	    
        TCanvas *c3 = new TCanvas("c3", "c3");
        c3->Clear();
        c3->Divide(3,1);
        gStyle->SetOptStat(1111111111);
        c3->cd(1);  shrinkPad(0.1, 0.2);
        hSigma_mm->Draw();
        c3->cd(2);  shrinkPad(0.1, 0.2);
        hSigma_mt->Draw();
        c3->cd(3);  shrinkPad(0.1, 0.2);	    
        hSigma_mmbar->Draw();
        c3->SaveAs("Sigma_polCB.eps");	    	
	    	
        TCanvas *c4 = new TCanvas("c4", "c4");
        c4->Clear();
        c4->Divide(3,1);
        gStyle->SetOptStat(1111111111);
        c4->cd(1);  shrinkPad(0.1, 0.2);
        hAlpha_mm->Draw();
        c4->cd(2);  shrinkPad(0.1, 0.2);
        hAlpha_mt->Draw();
        c4->cd(3);  shrinkPad(0.1, 0.2);	    
        hAlpha_mmbar->Draw();
        c4->SaveAs("Alpha_polCB.eps");
		
        TCanvas *c5 = new TCanvas("c5", "c5");
        c5->Clear();
        c5->Divide(3,1);
        gStyle->SetOptStat(1111111111);
        c5->cd(1);  shrinkPad(0.1, 0.2);
        hn_mm->Draw();
        c5->cd(2);  shrinkPad(0.1, 0.2);
        hn_mt->Draw();
        c5->cd(3);  shrinkPad(0.1, 0.2);	    
        hn_mmbar->Draw();
        c5->SaveAs("n_polCB.eps");
	
  }
  
  if ( mode == 2 ){   // fit ALL Upsilons with pol + 3CB

  //    o function for fitting
  //    o do in there something about low statistics cases

	int PRINT(1); 
	double PRINTX(0.3);

	TH1D *h; 

	  //  string fopt("LLIEMQ"); 
	string fopt("LLE"); 

  	double pt, eta; 
  	int    Q; 
	
  	double yield, yieldE, width, widthE; 
  	int    nbin;
	int fitted(0);
	 
  	int status(0);
	const char* Status;
	
	int Npt = 6;
	double PTbin[] = {0., 2., 3., 4., 5., 6., 15.};				
								
	int Neta = 5;					
	double Etabin[] = {-2.4, -1.2, -0.4, 0.4, 1.2, 2.4};
	
	int bin(-1);
	
	TH2D *hmean1_pos_mm = new TH2D("hmean_pos_mm", "hmean1_pos_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hmean1_pos_mt = new TH2D("hmean1_pos_mt", "hmean1_pos_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hmean1_pos_mmbar = new TH2D("hmean1_pos_mmbar", "hmean1_pos_mmbar", Neta, Etabin, Npt, PTbin);
	
	TH2D *hmean2_pos_mm = new TH2D("hmean2_pos_mm", "hmean2_pos_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hmean2_pos_mt = new TH2D("hmean2_pos_mt", "hmean2_pos_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hmean2_pos_mmbar = new TH2D("hmean2_pos_mmbar", "hmean2_pos_mmbar", Neta, Etabin, Npt, PTbin);
	
	TH2D *hsigma1_pos_mm = new TH2D("hsigma1_pos_mm", "hsigma1_pos_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hsigma1_pos_mt = new TH2D("hsigma1_pos_mt", "hsigma1_pos_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hsigma1_pos_mmbar = new TH2D("hsigma_pos1_mmbar", "hsigma1_pos_mmbar", Neta, Etabin, Npt, PTbin);
	
	TH2D *hsigma2_pos_mm = new TH2D("hsigma2_pos_mm", "hsigma2_pos_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hsigma2_pos_mt = new TH2D("hsigma2_pos_mt", "hsigma2_pos_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hsigma2_pos_mmbar = new TH2D("hsigma2_pos_mmbar", "hsigma2_pos_mmbar", Neta, Etabin, Npt, PTbin);		
	
	TH2D *hratio12_pos_mm = new TH2D("hratio12_pos_mm", "hratio12_pos_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hratio12_pos_mt = new TH2D("hratio12_pos_mt", "hratio12_pos_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hratio12_pos_mmbar = new TH2D("hratio12_pos_mmbar", "hratio12_pos_mmbar", Neta, Etabin, Npt, PTbin);
	
	TH2D *hratio13_pos_mm = new TH2D("hratio13_pos_mm", "hratio13_pos_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hratio13_pos_mt = new TH2D("hratio13_pos_mt", "hratio13_pos_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hratio13_pos_mmbar = new TH2D("hratio13_pos_mmbar", "hratio13_pos_mmbar", Neta, Etabin, Npt, PTbin);



	TH2D *hmean1_neg_mm = new TH2D("hmean_neg_mm", "hmean1_neg_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hmean1_neg_mt = new TH2D("hmean1_neg_mt", "hmean1_neg_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hmean1_neg_mmbar = new TH2D("hmean1_neg_mmbar", "hmean1_neg_mmbar", Neta, Etabin, Npt, PTbin);
	
	TH2D *hmean2_neg_mm = new TH2D("hmean2_neg_mm", "hmean2_neg_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hmean2_neg_mt = new TH2D("hmean2_neg_mt", "hmean2_neg_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hmean2_neg_mmbar = new TH2D("hmean2_neg_mmbar", "hmean2_neg_mmbar", Neta, Etabin, Npt, PTbin);
	
	TH2D *hsigma1_neg_mm = new TH2D("hsigma1_neg_mm", "hsigma1_neg_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hsigma1_neg_mt = new TH2D("hsigma1_neg_mt", "hsigma1_neg_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hsigma1_neg_mmbar = new TH2D("hsigma_neg1_mmbar", "hsigma1_neg_mmbar", Neta, Etabin, Npt, PTbin);
	
	TH2D *hsigma2_neg_mm = new TH2D("hsigma2_neg_mm", "hsigma2_neg_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hsigma2_neg_mt = new TH2D("hsigma2_neg_mt", "hsigma2_neg_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hsigma2_neg_mmbar = new TH2D("hsigma2_neg_mmbar", "hsigma2_neg_mmbar", Neta, Etabin, Npt, PTbin);		
	
	TH2D *hratio12_neg_mm = new TH2D("hratio12_neg_mm", "hratio12_neg_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hratio12_neg_mt = new TH2D("hratio12_neg_mt", "hratio12_neg_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hratio12_neg_mmbar = new TH2D("hratio12_neg_mmbar", "hratio12_neg_mmbar", Neta, Etabin, Npt, PTbin);
	
	TH2D *hratio13_neg_mm = new TH2D("hratio13_neg_mm", "hratio13_neg_mm", Neta, Etabin, Npt, PTbin);
	TH2D *hratio13_neg_mt = new TH2D("hratio13_neg_mt", "hratio13_neg_mt", Neta, Etabin, Npt, PTbin);
	TH2D *hratio13_neg_mmbar = new TH2D("hratio13_neg_mmbar", "hratio13_neg_mmbar", Neta, Etabin, Npt, PTbin);
	
	
				
	TH1D *hChisq_mm    = new TH1D("hChisq/ndof_mm"   , "hChisq/ndof_mm"   , 40, -20.,20.);
	TH1D *hChisq_mt    = new TH1D("hChisq/ndof_mt"   , "hChisq/ndof_mt"   , 40, -20.,20.); 
	TH1D *hChisq_mmbar = new TH1D("hChisq/ndof_mmbar", "hChisq/ndof_mmbar", 40, -20.,20.);
	
	
 	TH1D *hSigma_mm    = new TH1D("hSigma_mm"   , "hSigma_mm"   , 28, -0.07,0.07);
	TH1D *hSigma_mt    = new TH1D("hSigma_mt"   , "hSigma_mt"   , 28, -0.07,0.07); 
	TH1D *hSigma_mmbar = new TH1D("hSigma_mmbar", "hSigma_mmbar", 28, -0.07,0.07);
	
 	TH1D *hAlpha_mm    = new TH1D("hAlpha_mm"   , "hAlpha_mm"   , 30, -3.,3.);
	TH1D *hAlpha_mt    = new TH1D("hAlpha_mt"   , "hAlpha_mt"   , 30, -3.,3.); 
	TH1D *hAlpha_mmbar = new TH1D("hAlpha_mmbar", "hAlpha_mmbar", 30, -3.,3.);		
	
 	TH1D *hn_mm        = new TH1D("hn_mm"   , "hn_mm"   , 40, -4.,4.);
	TH1D *hn_mt        = new TH1D("hn_mt"   , "hn_mt"   , 40, -4.,4.); 
	TH1D *hn_mmbar     = new TH1D("hn_mmbar", "hn_mmbar", 40, -4.,4.);
	
//  	gStyle->SetOptStat(PRINT); 
//  	gStyle->SetOptFit(PRINT);
  	gStyle->SetOptStat(0000000000000); 
  	gStyle->SetOptFit(00000000000000);		 
  	makeCanvas(1); 
  	c1->Clear();
  	c1->Divide(3,1);
	 
  	for (unsigned int i = 0; i < fS1VectorPos.size(); ++i) {

  		  // -- positive charge
    		c1->cd(1); shrinkPad(0.15, 0.26); 
    		h = &(fS1VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f13, 6); 
    		  h->Fit(f13, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  ///////////////
		  if ( status == -1 || status == 0 ){
			f13->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) ,
			f13->GetParameter(5), f13->GetParameter(6));
			h->Fit(f13, fopt.c_str());
		  }
			
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
			
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	    	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  /////////////////////////	
		  if ( status == 1 ){
			hChisq_mm->Fill(status*(f13->GetChisquare()/f13->GetNDF()));
			hSigma_mm->Fill(status*f13->GetParameter(1));
			hAlpha_mm->Fill(status*f13->GetParameter(2));
			hn_mm->Fill(status*f13->GetParameter(3));
			
			///
			getBinCenters(h->GetName(), eta, pt, Q);
			bin = hmean1_pos_mm->FindBin(eta, pt);
			hmean1_pos_mm->SetBinContent(bin, status*f13->GetParameter(0));
			hmean2_pos_mm->SetBinContent(bin, status*f13->GetParameter(5));
			hsigma1_pos_mm->SetBinContent(bin, status*f13->GetParameter(1));
			hsigma2_pos_mm->SetBinContent(bin, status*f13->GetParameter(6));
			hratio12_pos_mm->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(7)));
			hratio13_pos_mm->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(8)));
			///
		  }
		  f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) );
    		  yield  = f10->Integral(8.7,9.7)/h->GetBinWidth(1);
    		  yieldE = TMath::Sqrt(yield);
		  fitted = 1;
    		} else {
    		  h->Draw();
    		  yield  = h->GetSumOfWeights();
    		  yieldE = TMath::Sqrt(h->GetSumOfWeights());
		  fitted = 0;
    		}
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS1YieldPos->FindBin(eta, pt); 
    		fS1YieldPos->SetBinContent(nbin, yield); 
    		fS1YieldPos->SetBinError(nbin, yieldE); 
		
    		c1->Modified();
    		c1->Update();
		
    		c1->cd(2); shrinkPad(0.15, 0.26); 
    		h = &(fS2VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f13, 6); 
    		  h->Fit(f13, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	     	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  ///////////////
		  if ( status == -1 || status == 0 ){
			f13->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) ,
			f13->GetParameter(5), f13->GetParameter(6));
			h->Fit(f13, fopt.c_str());
		  }
		
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  /////////////////////////
		  if ( status == 1 ){
			hChisq_mt->Fill(status*(f13->GetChisquare()/f13->GetNDF()));
			hSigma_mt->Fill(status*f13->GetParameter(1));
			hAlpha_mt->Fill(status*f13->GetParameter(2));
			hn_mt->Fill(status*f13->GetParameter(3));
			
			///
			getBinCenters(h->GetName(), eta, pt, Q);
			bin = hmean1_pos_mt->FindBin(eta, pt);
			hmean1_pos_mt->SetBinContent(bin, status*f13->GetParameter(0));
			hmean2_pos_mt->SetBinContent(bin, status*f13->GetParameter(5));
			hsigma1_pos_mt->SetBinContent(bin, status*f13->GetParameter(1));
			hsigma2_pos_mt->SetBinContent(bin, status*f13->GetParameter(6));
			hratio12_pos_mt->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(7)));
			hratio13_pos_mt->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(8)));
			///
			
		  }
		  f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) );
		  yield  = f10->Integral(8.7,9.7)/h->GetBinWidth(1);
      		  yieldE = TMath::Sqrt(yield);
      		  width  = f13->GetParameter(1); 
      		  widthE = f13->GetParError(1); 
      		  fitted = 1; 
    		} else {
    		  h->Draw();
    		  yield  = h->GetSumOfWeights();
    		  yieldE = TMath::Sqrt(h->GetSumOfWeights());
		  fitted = 0;
    		}
		
 		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS2YieldPos->FindBin(eta, pt); 
    		fS2YieldPos->SetBinContent(nbin, yield); 
    		fS2YieldPos->SetBinError(nbin, yieldE); 
		
    		c1->Modified();
    		c1->Update();
		
    		c1->cd(3);  shrinkPad(0.15, 0.26); 
    		h = &(fS3VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f13, 6); 
     	  	  h->Fit(f13, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  ///////////////
		  if ( status == -1 || status == 0 ){
			f13->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) ,
			f13->GetParameter(5), f13->GetParameter(6));
			h->Fit(f13, fopt.c_str());
		  }
		
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  /////////////////////////
		  if ( status == 1 ){
			hChisq_mmbar->Fill(status*(f13->GetChisquare()/f13->GetNDF()));
			hSigma_mmbar->Fill(status*f13->GetParameter(1));
			hAlpha_mmbar->Fill(status*f13->GetParameter(2));
			hn_mmbar->Fill(status*f13->GetParameter(3));
			
			///
			getBinCenters(h->GetName(), eta, pt, Q);
			bin = hmean1_pos_mmbar->FindBin(eta, pt);
			hmean1_pos_mmbar->SetBinContent(bin, status*f13->GetParameter(0));
			hmean2_pos_mmbar->SetBinContent(bin, status*f13->GetParameter(5));
			hsigma1_pos_mmbar->SetBinContent(bin, status*f13->GetParameter(1));
			hsigma2_pos_mmbar->SetBinContent(bin, status*f13->GetParameter(6));
			hratio12_pos_mmbar->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(7)));
			hratio13_pos_mmbar->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(8)));
			///
			
		  }		
		  f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) );
		  yield  = f10->Integral(8.7,9.7)/h->GetBinWidth(1);
      		  yieldE = TMath::Sqrt(yield);		
      		  fitted = 1; 
    		} else {
      		  h->Draw();
      		  yield  = h->GetSumOfWeights();
      		  yieldE = TMath::Sqrt(h->GetSumOfWeights());
		  fitted = 0;
    		}
	
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS3YieldPos->FindBin(eta, pt); 
    		double mm = fS1YieldPos->GetBinContent(nbin); 
    		double mmE= fS1YieldPos->GetBinError(nbin); 
    		fS3YieldPos->SetBinContent(nbin, yield+mm); 
    		fS3YieldPos->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 
		
   		TString frag(h->GetName()); 
    		frag.ReplaceAll("s3:mmbar,", ""); 
    		frag.ReplaceAll(",Q1", ""); 
    		frag.ReplaceAll(".", ":"); 
	
    		c1->SaveAs(Form("%s/massfits-pos-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data())); 
  	}
  	
  
  		c1->Clear();
  		c1->Divide(3,1); 
  	for (unsigned int i = 0; i < fS1VectorNeg.size(); ++i) {

  		  // -- negative charge
    		c1->cd(1); shrinkPad(0.15, 0.26);  
    		h = &(fS1VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f13, 6); 
    		  h->Fit(f13, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  ///////////////
		  if ( status == -1 || status == 0 ){
			f13->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) ,
			f13->GetParameter(5), f13->GetParameter(6));
			h->Fit(f13, fopt.c_str());
		  }
		
	  	  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  /////////////////////////
		  if ( status == 1 ){
			hChisq_mm->Fill(status*(f13->GetChisquare()/f13->GetNDF()));
			hSigma_mm->Fill(status*f13->GetParameter(1));
			hAlpha_mm->Fill(status*f13->GetParameter(2));
			hn_mm->Fill(status*f13->GetParameter(3));
			
			///
			getBinCenters(h->GetName(), eta, pt, Q);
			bin = hmean1_neg_mm->FindBin(eta, pt);
			hmean1_neg_mm->SetBinContent(bin, status*f13->GetParameter(0));
			hmean2_neg_mm->SetBinContent(bin, status*f13->GetParameter(5));
			hsigma1_neg_mm->SetBinContent(bin, status*f13->GetParameter(1));
			hsigma2_neg_mm->SetBinContent(bin, status*f13->GetParameter(6));
			hratio12_neg_mm->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(7)));
			hratio13_neg_mm->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(8)));
			///
			
		  }		
		  f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) );
		  yield  = f10->Integral(8.7,9.7)/h->GetBinWidth(1);
      		  yieldE = TMath::Sqrt(yield);		
    		} else {
    	  	  h->Draw();
    	  	  yield  = h->GetSumOfWeights();
    	  	  yieldE = TMath::Sqrt(h->GetSumOfWeights());
    		}
	
		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS1YieldNeg->FindBin(eta, pt); 
    		fS1YieldNeg->SetBinContent(nbin, yield); 
    		fS1YieldNeg->SetBinError(nbin, yieldE); 
	
    		c1->Modified();
    		c1->Update();
	
    		c1->cd(2); shrinkPad(0.15, 0.26); 
    		h = &(fS2VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    	  	  setFunctionParameters(h, f13, 6); 
    	  	  h->Fit(f13, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
		}
		cout << status << endl;
		///////////////
		if ( status == -1 || status == 0 ){
			f13->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) ,
			f13->GetParameter(5), f13->GetParameter(6));
			h->Fit(f13, fopt.c_str());
		}
		
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		/////////////////////////
		if ( status == 1 ){
			hChisq_mt->Fill(status*(f13->GetChisquare()/f13->GetNDF()));
			hSigma_mt->Fill(status*f13->GetParameter(1));
			hAlpha_mt->Fill(status*f13->GetParameter(2));
			hn_mt->Fill(status*f13->GetParameter(3));
			
			///
			getBinCenters(h->GetName(), eta, pt, Q);
			bin = hmean1_neg_mt->FindBin(eta, pt);
			hmean1_neg_mt->SetBinContent(bin, status*f13->GetParameter(0));
			hmean2_neg_mt->SetBinContent(bin, status*f13->GetParameter(5));
			hsigma1_neg_mt->SetBinContent(bin, status*f13->GetParameter(1));
			hsigma2_neg_mt->SetBinContent(bin, status*f13->GetParameter(6));
			hratio12_neg_mt->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(7)));
			hratio13_neg_mt->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(8)));
			///
			
		}		
		f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) );
		yield  = f10->Integral(8.7,9.7)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);		
      		width  = f13->GetParameter(1); 
      		widthE = f13->GetParError(1); 
		fitted = 1;  
    	 	} else {
    	  	  h->Draw();
    	   	  yield  = h->GetSumOfWeights();
    	   	  yieldE = TMath::Sqrt(h->GetSumOfWeights());
		  fitted = 0;
    		}
	
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS2YieldNeg->FindBin(eta, pt); 
    		fS2YieldNeg->SetBinContent(nbin, yield); 
    		fS2YieldNeg->SetBinError(nbin, yieldE); 
	
    		c1->cd(3);  shrinkPad(0.15, 0.26); 
    		h = &(fS3VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		if (h->GetSumOfWeights() > 10.) {
    		  setFunctionParameters(h, f13, 6); 
    	  	  h->Fit(f13, fopt.c_str());
		  status = 0;
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  ///////////////
		  if ( status == -1 || status == 0 ){
			f13->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) ,
			f13->GetParameter(5), f13->GetParameter(6));
			h->Fit(f13, fopt.c_str());
		  }
		
		  cout << gMinuit->fCstatu.Data() << endl;
		  Status = gMinuit->fCstatu.Data();
		  cout << Status[0] << endl;
		
		  if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	  } else if ( Status[0] == 'F' ){ 
	      		status = -1;
		  }
		  cout << status << endl;
		  /////////////////////////
		  if ( status == 1 ){
			hChisq_mmbar->Fill(status*(f13->GetChisquare()/f13->GetNDF()));
			hSigma_mmbar->Fill(status*f13->GetParameter(1));
			hAlpha_mmbar->Fill(status*f13->GetParameter(2));
			hn_mmbar->Fill(status*f13->GetParameter(3));
			
			///
			getBinCenters(h->GetName(), eta, pt, Q);
			bin = hmean1_neg_mmbar->FindBin(eta, pt);
			hmean1_neg_mmbar->SetBinContent(bin, status*f13->GetParameter(0));
			hmean2_neg_mmbar->SetBinContent(bin, status*f13->GetParameter(5));
			hsigma1_neg_mmbar->SetBinContent(bin, status*f13->GetParameter(1));
			hsigma2_neg_mmbar->SetBinContent(bin, status*f13->GetParameter(6));
			hratio12_neg_mmbar->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(7)));
			hratio13_neg_mmbar->SetBinContent(bin, status*(f13->GetParameter(4)/f13->GetParameter(8)));
			///
			
		  }		
		  f10->SetParameters( f13->GetParameter(0), f13->GetParameter(1), f13->GetParameter(2), f13->GetParameter(3), f13->GetParameter(4) );
		  yield  = f10->Integral(8.7,9.7)/h->GetBinWidth(1);
      		  yieldE = TMath::Sqrt(yield);		
      		  width  = f13->GetParameter(1); 
      		  widthE = f13->GetParError(1); 
		  fitted = 1;
    		} else {
    	  	  h->Draw();
    	  	  yield  = h->GetSumOfWeights();
    	  	  yieldE = TMath::Sqrt(h->GetSumOfWeights());
		  fitted = 0;
    		}

    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
    		cout << " --> " << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS3YieldNeg->FindBin(eta, pt); 
    		double mm = fS1YieldNeg->GetBinContent(nbin); 
    		double mmE= fS1YieldNeg->GetBinError(nbin); 
    		fS3YieldNeg->SetBinContent(nbin, yield+mm); 
    		fS3YieldNeg->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 
		
		
    		c1->Modified();
    		c1->Update();
	
    		TString frag(h->GetName()); 
    		frag.ReplaceAll("s3:mmbar,", ""); 
    		frag.ReplaceAll(",Q-1", ""); 
    		frag.ReplaceAll(".", ":"); 
		
    		c1->SaveAs(Form("%s/massfits-neg-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data())); 

  	}

	TCanvas *c2 = new TCanvas("c2", "c2");
  	c2->Clear();
        c2->Divide(3,1);
        gStyle->SetOptStat(1111111111);
        c2->cd(1);  shrinkPad(0.1, 0.2);
        hChisq_mm->Draw();
        c2->cd(2);  shrinkPad(0.1, 0.2);	
        hChisq_mt->Draw();
        c2->cd(3);  shrinkPad(0.1, 0.2);	    
        hChisq_mmbar->Draw();
        c2->SaveAs("chi2_polCB.eps");
	    
        TCanvas *c3 = new TCanvas("c3", "c3");
        c3->Clear();
        c3->Divide(3,1);
        gStyle->SetOptStat(1111111111);
        c3->cd(1);  shrinkPad(0.1, 0.2);
        hSigma_mm->Draw();
        c3->cd(2);  shrinkPad(0.1, 0.2);
        hSigma_mt->Draw();
        c3->cd(3);  shrinkPad(0.1, 0.2);	    
        hSigma_mmbar->Draw();
        c3->SaveAs("Sigma_polCB.eps");	    	
	    	
        TCanvas *c4 = new TCanvas("c4", "c4");
        c4->Clear();
        c4->Divide(3,1);
        gStyle->SetOptStat(1111111111);
        c4->cd(1);  shrinkPad(0.1, 0.2);
        hAlpha_mm->Draw();
        c4->cd(2);  shrinkPad(0.1, 0.2);
        hAlpha_mt->Draw();
        c4->cd(3);  shrinkPad(0.1, 0.2);	    
        hAlpha_mmbar->Draw();
        c4->SaveAs("Alpha_polCB.eps");
		
        TCanvas *c5 = new TCanvas("c5", "c5");
        c5->Clear();
        c5->Divide(3,1);
        gStyle->SetOptStat(1111111111);
        c5->cd(1);  shrinkPad(0.1, 0.2);
        hn_mm->Draw();
        c5->cd(2);  shrinkPad(0.1, 0.2);
        hn_mt->Draw();
        c5->cd(3);  shrinkPad(0.1, 0.2);	    
        hn_mmbar->Draw();
        c5->SaveAs("n_polCB.eps");
	
	TCanvas *c6 = new TCanvas("c6", "c6");
        c6->Clear();
        c6->Divide(3,1);
	hmean1_pos_mm->SetMinimum(9.41); hmean1_pos_mm->SetMaximum(9.51);
	hmean1_pos_mt->SetMinimum(9.41); hmean1_pos_mt->SetMaximum(9.51);
	hmean1_pos_mmbar->SetMinimum(9.41); hmean1_pos_mmbar->SetMaximum(9.51);
        gStyle->SetOptStat(1111111111);
        c6->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean1_pos_mm->SetStats(false);
        hmean1_pos_mm->Draw("colz");
        c6->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean1_pos_mt->SetStats(false);
        hmean1_pos_mt->Draw("colz");
        c6->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean1_pos_mmbar->SetStats(false);
        hmean1_pos_mmbar->Draw("colz");
        c6->SaveAs("mean1_pos.eps");	
	
	TCanvas *c7 = new TCanvas("c7", "c7");
        c7->Clear();
        c7->Divide(3,1);
	hmean2_pos_mm->SetMinimum(9.9); hmean2_pos_mm->SetMaximum(10.1);
	hmean2_pos_mt->SetMinimum(9.9); hmean2_pos_mt->SetMaximum(10.1);
	hmean2_pos_mmbar->SetMinimum(9.9); hmean2_pos_mmbar->SetMaximum(10.1);
        gStyle->SetOptStat(1111111111);
        c7->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean2_pos_mm->SetStats(false);
        hmean2_pos_mm->Draw("colz");
        c7->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean2_pos_mt->SetStats(false);
        hmean2_pos_mt->Draw("colz");
        c7->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean2_pos_mmbar->SetStats(false);
        hmean2_pos_mmbar->Draw("colz");
        c7->SaveAs("mean2_pos.eps");	
	
	TCanvas *c8 = new TCanvas("c8", "c8");
        c8->Clear();
        c8->Divide(3,1);
	hsigma1_pos_mm->SetMinimum(0.04); hsigma1_pos_mm->SetMaximum(0.14);
	hsigma1_pos_mt->SetMinimum(0.04); hsigma1_pos_mt->SetMaximum(0.14);
	hsigma1_pos_mmbar->SetMinimum(0.04); hsigma1_pos_mmbar->SetMaximum(0.14);
        gStyle->SetOptStat(1111111111);
        c8->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma1_pos_mm->SetStats(false);
        hsigma1_pos_mm->Draw("colz");
        c8->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma1_pos_mt->SetStats(false);
        hsigma1_pos_mt->Draw("colz");
        c8->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma1_pos_mmbar->SetStats(false);
        hsigma1_pos_mmbar->Draw("colz");
        c8->SaveAs("sigma1_pos.eps");	
	
	TCanvas *c9 = new TCanvas("c9", "c9");
        c9->Clear();
        c9->Divide(3,1);
	hsigma2_pos_mm->SetMinimum(0.06); hsigma2_pos_mm->SetMaximum(0.20);
	hsigma2_pos_mt->SetMinimum(0.06); hsigma2_pos_mt->SetMaximum(0.20);
	hsigma2_pos_mmbar->SetMinimum(0.06); hsigma2_pos_mmbar->SetMaximum(0.20);
        gStyle->SetOptStat(1111111111);
        c9->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma2_pos_mm->SetStats(false);
        hsigma2_pos_mm->Draw("colz");
        c9->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma2_pos_mt->SetStats(false);
        hsigma2_pos_mt->Draw("colz");
        c9->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma2_pos_mmbar->SetStats(false);
        hsigma2_pos_mmbar->Draw("colz");
        c9->SaveAs("sigma2_pos.eps");
	
	TCanvas *c10 = new TCanvas("c10", "c10");
        c10->Clear();
        c10->Divide(3,1);
	hratio12_pos_mm->SetMinimum(1.); hratio12_pos_mm->SetMaximum(6.);
	hratio12_pos_mt->SetMinimum(1.); hratio12_pos_mt->SetMaximum(6.);
	hratio12_pos_mmbar->SetMinimum(1.); hratio12_pos_mmbar->SetMaximum(6.);
        gStyle->SetOptStat(1111111111);
        c10->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio12_pos_mm->SetStats(false);
        hratio12_pos_mm->Draw("colz");
        c10->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio12_pos_mt->SetStats(false);
        hratio12_pos_mt->Draw("colz");
        c10->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio12_pos_mmbar->SetStats(false);
        hratio12_pos_mmbar->Draw("colz");
        c10->SaveAs("ratio12_pos.eps");
	
	TCanvas *c11 = new TCanvas("c11", "c11");
        c11->Clear();
        c11->Divide(3,1);
	hratio13_pos_mm->SetMinimum(1.); hratio13_pos_mm->SetMaximum(15.);
	hratio13_pos_mt->SetMinimum(1.); hratio13_pos_mt->SetMaximum(15.);
	hratio13_pos_mmbar->SetMinimum(1.); hratio13_pos_mmbar->SetMaximum(15.);
        gStyle->SetOptStat(1111111111);
        c11->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio13_pos_mm->SetStats(false);
        hratio13_pos_mm->Draw("colz");
        c11->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio13_pos_mt->SetStats(false);
        hratio13_pos_mt->Draw("colz");
        c11->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio13_pos_mmbar->SetStats(false);
        hratio13_pos_mmbar->Draw("colz");
        c11->SaveAs("ratio13_pos.eps");			
		
	TCanvas *c12 = new TCanvas("c12", "c12");
        c12->Clear();
        c12->Divide(3,1);
	hmean1_neg_mm->SetMinimum(9.41); hmean1_neg_mm->SetMaximum(9.51);
	hmean1_neg_mt->SetMinimum(9.41); hmean1_neg_mt->SetMaximum(9.51);
	hmean1_neg_mmbar->SetMinimum(9.41); hmean1_neg_mmbar->SetMaximum(9.51);
        gStyle->SetOptStat(1111111111);
        c12->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean1_neg_mm->SetStats(false);
        hmean1_neg_mm->Draw("colz");
        c12->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean1_neg_mt->SetStats(false);
        hmean1_neg_mt->Draw("colz");
        c12->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean1_neg_mmbar->SetStats(false);
        hmean1_neg_mmbar->Draw("colz");
        c12->SaveAs("mean1_neg.eps");	
	
	TCanvas *c13 = new TCanvas("c13", "c13");
        c13->Clear();
        c13->Divide(3,1);
	hmean2_neg_mm->SetMinimum(9.9); hmean2_neg_mm->SetMaximum(10.1);
	hmean2_neg_mt->SetMinimum(9.9); hmean2_neg_mt->SetMaximum(10.1);
	hmean2_neg_mmbar->SetMinimum(9.9); hmean2_neg_mmbar->SetMaximum(10.1);
        gStyle->SetOptStat(1111111111);
        c13->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean2_neg_mm->SetStats(false);
        hmean2_neg_mm->Draw("colz");
        c13->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean2_neg_mt->SetStats(false);
        hmean2_neg_mt->Draw("colz");
        c13->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hmean2_neg_mmbar->SetStats(false);
        hmean2_neg_mmbar->Draw("colz");
        c13->SaveAs("mean2_neg.eps");	
	
	TCanvas *c14 = new TCanvas("c14", "c14");
        c14->Clear();
        c14->Divide(3,1);
	hsigma1_neg_mm->SetMinimum(0.04); hsigma1_neg_mm->SetMaximum(0.14);
	hsigma1_neg_mt->SetMinimum(0.04); hsigma1_neg_mt->SetMaximum(0.14);
	hsigma1_neg_mmbar->SetMinimum(0.04); hsigma1_neg_mmbar->SetMaximum(0.14);
        gStyle->SetOptStat(1111111111);
        c14->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma1_neg_mm->SetStats(false);
        hsigma1_neg_mm->Draw("colz");
        c14->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma1_neg_mt->SetStats(false);
        hsigma1_neg_mt->Draw("colz");
        c14->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma1_neg_mmbar->SetStats(false);
        hsigma1_neg_mmbar->Draw("colz");
        c14->SaveAs("sigma1_neg.eps");	
	
	TCanvas *c15 = new TCanvas("c15", "c15");
        c15->Clear();
        c15->Divide(3,1);
	hsigma2_neg_mm->SetMinimum(0.06); hsigma2_neg_mm->SetMaximum(0.20);
	hsigma2_neg_mt->SetMinimum(0.06); hsigma2_neg_mt->SetMaximum(0.20);
	hsigma2_neg_mmbar->SetMinimum(0.06); hsigma2_neg_mmbar->SetMaximum(0.20);
        gStyle->SetOptStat(1111111111);
        c15->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma2_neg_mm->SetStats(false);
        hsigma2_neg_mm->Draw("colz");
        c15->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma2_neg_mt->SetStats(false);
        hsigma2_neg_mt->Draw("colz");
        c15->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hsigma2_neg_mmbar->SetStats(false);
        hsigma2_neg_mmbar->Draw("colz");
        c15->SaveAs("sigma2_neg.eps");
	
	TCanvas *c16 = new TCanvas("c16", "c16");
        c16->Clear();
        c16->Divide(3,1);
	hratio12_neg_mm->SetMinimum(1.); hratio12_neg_mm->SetMaximum(6.);
	hratio12_neg_mt->SetMinimum(1.); hratio12_neg_mt->SetMaximum(6.);
	hratio12_neg_mmbar->SetMinimum(1.); hratio12_neg_mmbar->SetMaximum(6.);
        gStyle->SetOptStat(1111111111);
        c16->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio12_neg_mm->SetStats(false);
        hratio12_neg_mm->Draw("colz");
        c16->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio12_neg_mt->SetStats(false);
        hratio12_neg_mt->Draw("colz");
        c16->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio12_neg_mmbar->SetStats(false);
        hratio12_neg_mmbar->Draw("colz");
        c16->SaveAs("ratio12_neg.eps");
	
	TCanvas *c17 = new TCanvas("c17", "c17");
        c17->Clear();
        c17->Divide(3,1);
	hratio13_neg_mm->SetMinimum(1.); hratio13_neg_mm->SetMaximum(15.);
	hratio13_neg_mt->SetMinimum(1.); hratio13_neg_mt->SetMaximum(15.);
	hratio13_neg_mmbar->SetMinimum(1.); hratio13_neg_mmbar->SetMaximum(15.);
        gStyle->SetOptStat(1111111111);
        c17->cd(1);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio13_neg_mm->SetStats(false);
        hratio13_neg_mm->Draw("colz");
        c17->cd(2);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio13_neg_mt->SetStats(false);
        hratio13_neg_mt->Draw("colz");
        c17->cd(3);  shrinkPad(0.1, 0.2); gPad->SetRightMargin(0.15);
	hratio13_neg_mmbar->SetStats(false);
        hratio13_neg_mmbar->Draw("colz");
        c17->SaveAs("ratio13_neg.eps");			

  }

}

// ----------------------------------------------------------------------
void anaXS::fitJpsi(int mode) {  

  if ( mode == 0 ){   // pol + Gauss	
  
    	int PRINT(1); 
  	double PRINTX(0.3); 
  	// -- FIXME: 
  	//    o function for fitting

  	TH1D *h; 

  	//  string fopt("LLIEMQ"); 
  	string fopt("LLE"); 

  	double pt, eta; 
  	int    Q; 

  	double yield, yieldE, width, widthE; 
  	int    nbin; 
  	int fitted(0); 
  	
	int status(0);
	const char* Status;
	
  	int MINFIT(20);
	 
 	TH1D *hChisq_mm    = new TH1D("hChisq/ndof_mm"   , "hChisq/ndof_mm"   , 40, -20.,20.);
	TH1D *hChisq_mt    = new TH1D("hChisq/ndof_mt"   , "hChisq/ndof_mt"   , 40, -20.,20.); 
	TH1D *hChisq_mmbar = new TH1D("hChisq/ndof_mmbar", "hChisq/ndof_mmbar", 40, -20.,20.);

  	gStyle->SetOptStat(PRINT); 
  	gStyle->SetOptFit(PRINT); 
  	makeCanvas(1); 
  	c1->Clear();
  	c1->Divide(3,1); 
  	// -- positive charge
  	for (unsigned int i = 0; i < fS1VectorPos.size(); ++i) {
  	  c1->cd(1);  shrinkPad(0.15, 0.26); 
  	  h = &(fS1VectorPos[i]);
  	  h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
  	  setFunctionParameters(h, f1, 1); 
  	  if (h->GetSumOfWeights() > MINFIT) {
  	    h->Fit(f1, fopt.c_str());
	    cout << gMinuit->fCstatu.Data() << endl;
	    Status = gMinuit->fCstatu.Data();
	    cout << Status[0] << endl;
	    if ( Status[0] == 'S' || Status[0] == 'P' ){
	    	status = 1;
	    } else if  ( Status[0] == 'F' ) { 
	    	status = -1;
	    	}
	    cout << status << endl;
	    hChisq_mm->Fill(status*(f1->GetChisquare()/f1->GetNDF()));
  	    yield  = f1->GetParameter(0)/h->GetBinWidth(1);
  	    yieldE = f1->GetParError(0)/h->GetBinWidth(1);
  	    fitted = 1; 
  	  } else {
  	    h->Draw();
  	    int ib = h->FindBin(3.1); 
  	    yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
  	    yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
  	    yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
  	    if (yield < 0) yield = 1.; 
	
      	yieldE =TMath::Sqrt(yield); 
      	fitted = 0; 
    	}

    	cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
	    if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    	getBinCenters(h->GetName(), eta, pt, Q); 
    	nbin = fS1YieldPos->FindBin(eta, pt); 
    	fS1YieldPos->SetBinContent(nbin, yield); 
    	fS1YieldPos->SetBinError(nbin, yieldE); 

    	c1->cd(2);  shrinkPad(0.15, 0.26); 
    	h = &(fS2VectorPos[i]);
    	h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    	setFunctionParameters(h, f1, 1); 
    	if (h->GetSumOfWeights() > MINFIT) {
    	  h->Fit(f1, fopt.c_str());
	  cout << gMinuit->fCstatu.Data() << endl;
	  Status = gMinuit->fCstatu.Data();
	  cout << Status[0] << endl;
	  if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	status = 1;
	  } else if ( Status[0] == 'F' ) { 
	    	status = -1;
	    	}
	  cout << status << endl;
	  hChisq_mt->Fill(status*(f1->GetChisquare()/f1->GetNDF()));
    	  yield  = f1->GetParameter(0)/h->GetBinWidth(1);
    	  yieldE = f1->GetParError(0)/h->GetBinWidth(1);
    	  width  = f1->GetParameter(2); 
    	  widthE = f1->GetParError(2); 
    	  fitted = 1; 
    	} else {
    	  int ib = h->FindBin(3.1); 
    	  h->Draw();
    	  yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
    	  yieldE =TMath::Sqrt(yield); 
    	  yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
    	  yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
    	  if (yield < 0) yield = 1.; 
    	  width  = 0.030; 
    	  width  = 0.010; 
    	  fitted = 0; 
    	}
    
    	cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    	if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    	getBinCenters(h->GetName(), eta, pt, Q); 
    	nbin = fS2YieldPos->FindBin(eta, pt); 
    	fS2YieldPos->SetBinContent(nbin, yield); 
    	fS2YieldPos->SetBinError(nbin, yieldE); 

    	c1->cd(3);  shrinkPad(0.15, 0.26); 
    	h = &(fS3VectorPos[i]);
    	h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    	setFunctionParameters(h, f1, 1); 
    	if (h->GetSumOfWeights() > MINFIT) {
    	  h->Fit(f1, fopt.c_str());
	  cout << gMinuit->fCstatu.Data() << endl;
	  Status = gMinuit->fCstatu.Data();
	  cout << Status[0] << endl;
	  if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	status = 1;
	  } else if ( Status[0] == 'F' ){ 
	  	status = -1;
		}
	  cout << status << endl;
	  hChisq_mmbar->Fill(status*(f1->GetChisquare()/f1->GetNDF()));
    	  yield  = f1->GetParameter(0)/h->GetBinWidth(1);
    	  yieldE = f1->GetParError(0)/h->GetBinWidth(1);
    	  fitted = 1; 
    	} else {
    	  h->Draw();
    	  int ib = h->FindBin(3.1); 
    	  yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
    	  yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
    	  yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
    	  if (yield < 0) yield = 1.; 
    	  yieldE =TMath::Sqrt(yield); 
    	  fitted = 0; 
    	}

    	cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    	if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

   	 getBinCenters(h->GetName(), eta, pt, Q); 
   	 nbin = fS3YieldPos->FindBin(eta, pt); 
   	 double mm = fS1YieldPos->GetBinContent(nbin); 
   	 double mmE= fS1YieldPos->GetBinError(nbin); 
   	 fS3YieldPos->SetBinContent(nbin, yield+mm); 
   	 fS3YieldPos->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 
	
    	c1->Update();
   	c1->Draw();

    	TString frag(h->GetName()); 
   	frag.ReplaceAll("s3:mmbar,", ""); 
  	frag.ReplaceAll(",Q1", ""); 
  	frag.ReplaceAll(".", ":"); 

  	  c1->SaveAs(Form("%s/massfits-pos-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data()));
	  
	   
  	}


  	c1->Clear();
  	c1->Divide(3,1); 
  	// -- negative charge
  	for (unsigned int i = 0; i < fS1VectorNeg.size(); ++i) {
  	  c1->cd(1);  shrinkPad(0.15, 0.26); 
  	  h = &(fS1VectorNeg[i]);
  	  h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
  	  setFunctionParameters(h, f1, 1); 
  	  if (h->GetSumOfWeights() > MINFIT) {
  	    h->Fit(f1, fopt.c_str());
	    cout << gMinuit->fCstatu.Data() << endl;
	    Status = gMinuit->fCstatu.Data();
	    cout << Status[0] << endl;
	    if ( Status[0] == 'S' || Status[0] == 'P' ){
	    	status = 1;
	    } else if ( Status[0] == 'F' ){ 
	    	status = -1;
		}
	    cout << status << endl;
	    hChisq_mm->Fill(status*(f1->GetChisquare()/f1->GetNDF()));
  	    yield  = f1->GetParameter(0)/h->GetBinWidth(1);
  	    yieldE = f1->GetParError(0)/h->GetBinWidth(1);
  	  } else {
  	    h->Draw();
  	    int ib = h->FindBin(3.1); 
  	    yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
  	    yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
  	    yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
  	    if (yield < 0) yield = 1.; 
  	    yieldE =TMath::Sqrt(yield); 
  	  }
  	  cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
	    if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 
	
	    getBinCenters(h->GetName(), eta, pt, Q); 
	    nbin = fS1YieldNeg->FindBin(eta, pt); 
	    fS1YieldNeg->SetBinContent(nbin, yield); 
	    fS1YieldNeg->SetBinError(nbin, yieldE); 
	
	    c1->cd(2);  shrinkPad(0.15, 0.26); 
	    h = &(fS2VectorNeg[i]);
	    h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
	    setFunctionParameters(h, f1, 1); 
	    if (h->GetSumOfWeights() > MINFIT) {
	      h->Fit(f1, fopt.c_str());
	      cout << gMinuit->fCstatu.Data() << endl;
	      Status = gMinuit->fCstatu.Data();
	      cout << Status[0] << endl;
	      if ( Status[0] == 'S' || Status[0] == 'P' ){
	    	status = 1;
	      } else if ( Status[0] == 'F' ){ 
	      	status = -1;
		}
	      cout << status << endl;
	      hChisq_mt->Fill(status*(f1->GetChisquare()/f1->GetNDF()));
	      yield  = f1->GetParameter(0)/h->GetBinWidth(1);
	      yieldE = f1->GetParError(0)/h->GetBinWidth(1);
	      width  = f1->GetParameter(2); 
	      widthE = f1->GetParError(2); 
	    } else {
	      h->Draw();
	      int ib = h->FindBin(3.1); 
	      yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
	      yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
	      yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
	      if (yield < 0) yield = 1.; 
	      yieldE =TMath::Sqrt(yield); 
	      width  = 0.030; 
	      widthE = 0.010; 
    	    }

  	  cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
	    if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

	    getBinCenters(h->GetName(), eta, pt, Q); 
	    nbin = fS2YieldNeg->FindBin(eta, pt); 
	    fS2YieldNeg->SetBinContent(nbin, yield); 
	    fS2YieldNeg->SetBinError(nbin, yieldE); 
	
	    c1->cd(3);  shrinkPad(0.15, 0.2); 
	    h = &(fS3VectorNeg[i]);
	    h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
	    setFunctionParameters(h, f1, 1); 
	    if (h->GetSumOfWeights() > MINFIT) {
	      h->Fit(f1, fopt.c_str());
	      cout << gMinuit->fCstatu.Data() << endl;
	      Status = gMinuit->fCstatu.Data();
	      cout << Status[0] << endl;
	      if ( Status[0] == 'S' || Status[0] == 'P' ){
	    	status = 1;
	      } else  if ( Status[0] == 'F' ){ 
	      	status = -1;
		}
	      cout << status << endl;
	      hChisq_mmbar->Fill(status*(f1->GetChisquare()/f1->GetNDF()));
	      yield  = f1->GetParameter(0)/h->GetBinWidth(1);
	      yieldE = f1->GetParError(0)/h->GetBinWidth(1);
	      fitted = 1; 
	    } else {
	      h->Draw();
	      int ib = h->FindBin(3.1); 
	      yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
	      yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
	      yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
	      if (yield < 0) yield = 1.; 
	      yieldE =TMath::Sqrt(yield); 
	      fitted = 0; 
	    }

	    cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
	    if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

	    getBinCenters(h->GetName(), eta, pt, Q); 
	    nbin = fS3YieldNeg->FindBin(eta, pt); 
	    double mm = fS1YieldNeg->GetBinContent(nbin); 
	    double mmE= fS1YieldNeg->GetBinError(nbin); 
	    fS3YieldNeg->SetBinContent(nbin, yield+mm); 
	    fS3YieldNeg->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 
	
	    c1->Update();
	    c1->Draw();
	
	    TString frag(h->GetName()); 
	    frag.ReplaceAll("s3:mmbar,", ""); 
	    frag.ReplaceAll(",Q-1", ""); 
	    frag.ReplaceAll(".", ":"); 
	
	    c1->SaveAs(Form("%s/massfits-neg-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data()));
	    
	    TCanvas *c2 = new TCanvas("c2", "c2");
  	    c2->Clear();
  	    c2->Divide(3,1);
	    gStyle->SetOptStat(1111111111);
	    c2->cd(1);  shrinkPad(0.1, 0.2);
	    hChisq_mm->Draw();
	    c2->cd(2);  shrinkPad(0.1, 0.2);
	    hChisq_mt->Draw();
	    c2->cd(3);  shrinkPad(0.1, 0.2);	    
	    hChisq_mmbar->Draw();
	    c2->SaveAs(Form("chi2_polGaus.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data()));	     
  	}

 	 c0->Clear(); 

 	 c0->Divide(2,2);
 	 c0->cd(1);
 	 fS1YieldPos->Draw("colz");
 	 c0->cd(2);
 	 fS2YieldPos->Draw("colz");

 	 c0->cd(3);
 	 fS1YieldNeg->Draw("colz");
 	 c0->cd(4);
 	 fS2YieldNeg->Draw("colz");	
	  
  }
  
  if ( mode == 1 ){   // pol + CB
  
  	int PRINT(1); 
  	double PRINTX(0.3); 
  	// -- FIXME: 
  	//    o function for fitting

  	TH1D *h; 

  	//  string fopt("LLIEMQ"); 
  	string fopt("LLE"); 

  	double pt, eta; 
  	int    Q; 
	
  	double yield, yieldE, width, widthE; 
  	int    nbin; 
  	int fitted(0);
	 
  	int status(0);
	const char* Status;

  	int MINFIT(20);
	 
 	TH1D *hChisq_mm    = new TH1D("hChisq/ndof_mm"   , "hChisq/ndof_mm"   , 40, -20.,20.);
	TH1D *hChisq_mt    = new TH1D("hChisq/ndof_mt"   , "hChisq/ndof_mt"   , 40, -20.,20.); 
	TH1D *hChisq_mmbar = new TH1D("hChisq/ndof_mmbar", "hChisq/ndof_mmbar", 40, -20.,20.);
	
 	TH1D *hSigma_mm    = new TH1D("hSigma_mm"   , "hSigma_mm"   , 20, -0.07,0.07);
	TH1D *hSigma_mt    = new TH1D("hSigma_mt"   , "hSigma_mt"   , 20, -0.07,0.07); 
	TH1D *hSigma_mmbar = new TH1D("hSigma_mmbar", "hSigma_mmbar", 20, -0.07,0.07);
	
 	TH1D *hAlpha_mm    = new TH1D("hAlpha_mm"   , "hAlpha_mm"   , 30, -3.,3.);
	TH1D *hAlpha_mt    = new TH1D("hAlpha_mt"   , "hAlpha_mt"   , 30, -3.,3.); 
	TH1D *hAlpha_mmbar = new TH1D("hAlpha_mmbar", "hAlpha_mmbar", 30, -3.,3.);		
	
 	TH1D *hn_mm        = new TH1D("hn_mm"   , "hn_mm"   , 40, -4.,4.);
	TH1D *hn_mt        = new TH1D("hn_mt"   , "hn_mt"   , 40, -4.,4.); 
	TH1D *hn_mmbar     = new TH1D("hn_mmbar", "hn_mmbar", 40, -4.,4.);
		
  	//gStyle->SetOptStat(PRINT); 
  	//gStyle->SetOptFit(PRINT);
  	gStyle->SetOptStat(0000000000000); 
  	gStyle->SetOptFit(00000000000000);		 
  	makeCanvas(1); 
  	c1->Clear();
  	c1->Divide(3,1); 
  	// -- positive charge  
  	for (unsigned int i = 0; i < fS1VectorPos.size(); ++i) {
    		c1->cd(1);  shrinkPad(0.15, 0.26); 
    		h = &(fS1VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f6, 3); 
    		if (h->GetSumOfWeights() > MINFIT) {
    		  	h->Fit(f6, fopt.c_str());
			status = 0;
			cout << gMinuit->fCstatu.Data() << endl;
			Status = gMinuit->fCstatu.Data();
			cout << Status[0] << endl;
		        if ( Status[0] == 'S' || Status[0] == 'P' ){
	   		    	status = 1;
	      		} else if ( Status[0] == 'F' ){ 
	      			status = -1;
				}
			cout << status << endl;
			///////////////
			if ( status == -1 || status == 0 ){
				f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
				f6->GetParameter(5), f6->GetParameter(6));
				h->Fit(f6, fopt.c_str());
			}
			
			cout << gMinuit->fCstatu.Data() << endl;
			Status = gMinuit->fCstatu.Data();
			cout << Status[0] << endl;
			
			if ( Status[0] == 'S' || Status[0] == 'P' ){
	   		    	status = 1;
	      		} else if ( Status[0] == 'F' ){ 
	      			status = -1;
				}
			cout << status << endl;
			/////////////////////////	
			if ( status == 1 ){
				hChisq_mm->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
				hSigma_mm->Fill(status*f6->GetParameter(1));
				hAlpha_mm->Fill(status*f6->GetParameter(2));
				hn_mm->Fill(status*f6->GetParameter(3));
			}
			f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
			yield  = f10->Integral(2.8,3.4)/h->GetBinWidth(1);
      			yieldE = TMath::Sqrt(yield);
      			fitted = 1; 
    		} else {
    			h->Draw();
   			int ib = h->FindBin(3.1); 
      			yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      			yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      			yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
     			if (yield < 0) yield = 1.; 

      			yieldE =TMath::Sqrt(yield); 
      			fitted = 0; 
    		}
		
    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

	        getBinCenters(h->GetName(), eta, pt, Q); 
		nbin = fS1YieldPos->FindBin(eta, pt); 
    	        fS1YieldPos->SetBinContent(nbin, yield); 
   		fS1YieldPos->SetBinError(nbin, yieldE); 

    		c1->cd(2);  shrinkPad(0.15, 0.26); 
    		h = &(fS2VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f6, 3); 
   		if (h->GetSumOfWeights() > MINFIT) {
      		h->Fit(f6, fopt.c_str());
		status = 0;
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		///////////////
		if ( status == -1 || status == 0 ){
			f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
			f6->GetParameter(5), f6->GetParameter(6));
			h->Fit(f6, fopt.c_str());
		}
		
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		/////////////////////////
		if ( status == 1 ){
			hChisq_mt->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
			hSigma_mt->Fill(status*f6->GetParameter(1));
			hAlpha_mt->Fill(status*f6->GetParameter(2));
			hn_mt->Fill(status*f6->GetParameter(3));
		}
		f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
		yield  = f10->Integral(2.8,3.4)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);
      		width  = f6->GetParameter(1); 
      		widthE = f6->GetParError(1); 
      		fitted = 1; 
    		} else {
      		int ib = h->FindBin(3.1); 
      		h->Draw();
      		yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      		yieldE =TMath::Sqrt(yield); 
      		yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      		yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
      		if (yield < 0) yield = 1.; 
      		width  = 0.030; 
      		width  = 0.010; 
      		fitted = 0; 
    		}	
    		
    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
	  	if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS2YieldPos->FindBin(eta, pt); 
    		fS2YieldPos->SetBinContent(nbin, yield); 
    		fS2YieldPos->SetBinError(nbin, yieldE); 
	
    		c1->cd(3);  shrinkPad(0.15, 0.26); 
    		h = &(fS3VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f6, 3); 
    		if (h->GetSumOfWeights() > MINFIT) {
      		h->Fit(f6, fopt.c_str());
		status = 0;
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		///////////////
		if ( status == -1 || status == 0 ){
			f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
			f6->GetParameter(5), f6->GetParameter(6));
			h->Fit(f6, fopt.c_str());
		}
		
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		/////////////////////////
		if ( status == 1 ){
			hChisq_mmbar->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
			hSigma_mmbar->Fill(status*f6->GetParameter(1));
			hAlpha_mmbar->Fill(status*f6->GetParameter(2));
			hn_mmbar->Fill(status*f6->GetParameter(3));
		}		
		f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
		yield  = f10->Integral(2.8,3.4)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);		
      		fitted = 1; 
    		} else {
    	 	h->Draw();
      		int ib = h->FindBin(3.1); 
      		yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      		yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      		yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
      		if (yield < 0) yield = 1.; 
      		yieldE =TMath::Sqrt(yield); 
      		fitted = 0; 
    		}	

    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
	 	<< h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS3YieldPos->FindBin(eta, pt); 
    		double mm = fS1YieldPos->GetBinContent(nbin); 
    		double mmE= fS1YieldPos->GetBinError(nbin); 
    		fS3YieldPos->SetBinContent(nbin, yield+mm); 
    		fS3YieldPos->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 

    		c1->Update();
    		c1->Draw();

    		TString frag(h->GetName()); 
    		frag.ReplaceAll("s3:mmbar,", ""); 
    		frag.ReplaceAll(",Q1", ""); 
    		frag.ReplaceAll(".", ":"); 
	
    		c1->SaveAs(Form("%s/massfits-pos-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data())); 
  		}


  		c1->Clear();
  		c1->Divide(3,1); 
  		// -- negative charge
  		for (unsigned int i = 0; i < fS1VectorNeg.size(); ++i) {
    		c1->cd(1);  shrinkPad(0.15, 0.26); 
    		h = &(fS1VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f6, 3); 
    		if (h->GetSumOfWeights() > MINFIT) {
    		h->Fit(f6, fopt.c_str());
		status = 0;
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		///////////////
		if ( status == -1 || status == 0 ){
			f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
			f6->GetParameter(5), f6->GetParameter(6));
			h->Fit(f6, fopt.c_str());
		}
		
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		/////////////////////////
		if ( status == 1 ){
			hChisq_mm->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
			hSigma_mm->Fill(status*f6->GetParameter(1));
			hAlpha_mm->Fill(status*f6->GetParameter(2));
			hn_mm->Fill(status*f6->GetParameter(3));
		}		
		f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
		yield  = f10->Integral(2.8,3.4)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);		
    		} else {
      		h->Draw();
      		int ib = h->FindBin(3.1); 
      		yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      		yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      		yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
      		if (yield < 0) yield = 1.; 
      		yieldE =TMath::Sqrt(yield); 
    		}
    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
	 	<< h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS1YieldNeg->FindBin(eta, pt); 
    		fS1YieldNeg->SetBinContent(nbin, yield); 
    		fS1YieldNeg->SetBinError(nbin, yieldE); 

    		c1->cd(2);  shrinkPad(0.15, 0.26); 
    		h = &(fS2VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f6, 3); 
    		if (h->GetSumOfWeights() > MINFIT) {
    		h->Fit(f6, fopt.c_str());
		status = 0;
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		///////////////
		if ( status == -1 || status == 0 ){
			f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
			f6->GetParameter(5), f6->GetParameter(6));
			h->Fit(f6, fopt.c_str());
		}
		
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		/////////////////////////
		if ( status == 1 ){
			hChisq_mt->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
			hSigma_mt->Fill(status*f6->GetParameter(1));
			hAlpha_mt->Fill(status*f6->GetParameter(2));
			hn_mt->Fill(status*f6->GetParameter(3));
		}		
		f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
		yield  = f10->Integral(2.8,3.4)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);		
      		width  = f6->GetParameter(1); 
      		widthE = f6->GetParError(1); 
    		} else {
      		h->Draw();
      		int ib = h->FindBin(3.1); 
      		yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      		yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      		yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
      		if (yield < 0) yield = 1.; 
      		yieldE =TMath::Sqrt(yield); 
      		width  = 0.030; 
      		widthE = 0.010; 
    		}

    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS2YieldNeg->FindBin(eta, pt); 
    		fS2YieldNeg->SetBinContent(nbin, yield); 
    		fS2YieldNeg->SetBinError(nbin, yieldE); 

    		c1->cd(3);  shrinkPad(0.15, 0.2); 
    		h = &(fS3VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f6, 3); 
    		if (h->GetSumOfWeights() > MINFIT) {
      		h->Fit(f6, fopt.c_str());
		status = 0;
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		///////////////
		if ( status == -1 || status == 0 ){
			f6->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) ,
			f6->GetParameter(5), f6->GetParameter(6));
			h->Fit(f6, fopt.c_str());
		}
		
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	  	    	status = 1;
	      	} else if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		/////////////////////////
		if ( status == 1 ){
			hChisq_mmbar->Fill(status*(f6->GetChisquare()/f6->GetNDF()));
			hSigma_mmbar->Fill(status*f6->GetParameter(1));
			hAlpha_mmbar->Fill(status*f6->GetParameter(2));
			hn_mmbar->Fill(status*f6->GetParameter(3));
		}		
		f10->SetParameters( f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), f6->GetParameter(3), f6->GetParameter(4) );
		yield  = f10->Integral(2.8,3.4)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);		
      		fitted = 1; 
    		} else {
      		h->Draw();
      		int ib = h->FindBin(3.1); 
      		yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      		yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      		yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
     		if (yield < 0) yield = 1.; 
     		yieldE =TMath::Sqrt(yield); 
     		fitted = 0; 
    		}

    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
	 	<< h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS3YieldNeg->FindBin(eta, pt); 
    		double mm = fS1YieldNeg->GetBinContent(nbin); 
    		double mmE= fS1YieldNeg->GetBinError(nbin); 
    		fS3YieldNeg->SetBinContent(nbin, yield+mm); 
    		fS3YieldNeg->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 

    		c1->Update();
    		c1->Draw();
	
    		TString frag(h->GetName()); 
    		frag.ReplaceAll("s3:mmbar,", ""); 
    		frag.ReplaceAll(",Q-1", ""); 
    		frag.ReplaceAll(".", ":"); 

    		c1->SaveAs(Form("%s/massfits-neg-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data()));
		
 		}
				 
	        TCanvas *c2 = new TCanvas("c2", "c2");
  	        c2->Clear();
  	        c2->Divide(3,1);
	        gStyle->SetOptStat(1111111111);
	        c2->cd(1);  shrinkPad(0.1, 0.2);
	        hChisq_mm->Draw();
	        c2->cd(2);  shrinkPad(0.1, 0.2);
	        hChisq_mt->Draw();
	        c2->cd(3);  shrinkPad(0.1, 0.2);	    
	        hChisq_mmbar->Draw();
	        c2->SaveAs("chi2_polCB.eps");
	    
	        TCanvas *c3 = new TCanvas("c3", "c3");
  	        c3->Clear();
  	        c3->Divide(3,1);
	        gStyle->SetOptStat(1111111111);
	        c3->cd(1);  shrinkPad(0.1, 0.2);
	        hSigma_mm->Draw();
	        c3->cd(2);  shrinkPad(0.1, 0.2);
	        hSigma_mt->Draw();
	        c3->cd(3);  shrinkPad(0.1, 0.2);	    
	        hSigma_mmbar->Draw();
	        c3->SaveAs("Sigma_polCB.eps");	    	
	    	
	        TCanvas *c4 = new TCanvas("c4", "c4");
  	        c4->Clear();
  	        c4->Divide(3,1);
	        gStyle->SetOptStat(1111111111);
	        c4->cd(1);  shrinkPad(0.1, 0.2);
	        hAlpha_mm->Draw();
	        c4->cd(2);  shrinkPad(0.1, 0.2);
	        hAlpha_mt->Draw();
	        c4->cd(3);  shrinkPad(0.1, 0.2);	    
	        hAlpha_mmbar->Draw();
	        c4->SaveAs("Alpha_polCB.eps");
		
	        TCanvas *c5 = new TCanvas("c5", "c5");
  	        c5->Clear();
  	        c5->Divide(3,1);
	        gStyle->SetOptStat(1111111111);
	        c5->cd(1);  shrinkPad(0.1, 0.2);
	        hn_mm->Draw();
	        c5->cd(2);  shrinkPad(0.1, 0.2);
	        hn_mt->Draw();
	        c5->cd(3);  shrinkPad(0.1, 0.2);	    
	        hn_mmbar->Draw();
	        c5->SaveAs("n_polCB.eps");		

  		c0->Clear(); 

  		c0->Divide(2,2);
  		c0->cd(1);
  		fS1YieldPos->Draw("colz");
  		c0->cd(2);
  		fS2YieldPos->Draw("colz");
	
  		c0->cd(3);
  		fS1YieldNeg->Draw("colz");
  		c0->cd(4);
  		fS2YieldNeg->Draw("colz");
	
  
  }
  
  if ( mode == 2 ){   // pol + CB + Gaussian
  
  	int PRINT(1); 
  	double PRINTX(0.3); 
  	// -- FIXME: 
  	//    o function for fitting

  	TH1D *h; 

  	//  string fopt("LLIEMQ"); 
  	string fopt("LLE"); 

  	double pt, eta; 
  	int    Q; 
	
  	double yield, yieldE, width, widthE; 
  	int    nbin; 
  	int fitted(0); 
  	
	int status(0);
	const char* Status;
	
  	int MINFIT(50);
	 
 	TH1D *hChisq_mm    = new TH1D("hChisq/ndof_mm"   , "hChisq/ndof_mm"   , 40, -20.,20.);
	TH1D *hChisq_mt    = new TH1D("hChisq/ndof_mt"   , "hChisq/ndof_mt"   , 40, -20.,20.); 
	TH1D *hChisq_mmbar = new TH1D("hChisq/ndof_mmbar", "hChisq/ndof_mmbar", 40, -20.,20.);
	
  	gStyle->SetOptStat(PRINT); 
  	gStyle->SetOptFit(PRINT); 
  	makeCanvas(1); 
  	c1->Clear();
  	c1->Divide(3,1); 
  	// -- positive charge  
  	for (unsigned int i = 0; i < fS1VectorPos.size(); ++i) {
    		c1->cd(1);  shrinkPad(0.15, 0.26); 
    		h = &(fS1VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f12, 4); 
    		if (h->GetSumOfWeights() > MINFIT) {
    		  	h->Fit(f12, fopt.c_str());
			cout << gMinuit->fCstatu.Data() << endl;
			Status = gMinuit->fCstatu.Data();
			cout << Status[0] << endl;
			if ( Status[0] == 'S' || Status[0] == 'P' ){
	   		   	status = 1;
	      		} else  if ( Status[0] == 'F' ){ 
	      			status = -1;
				}
			cout << status << endl;			
			hChisq_mm->Fill(status*(f12->GetChisquare()/f12->GetNDF()));
			f11->SetParameters( f12->GetParameter(0), f12->GetParameter(1), f12->GetParameter(2), f12->GetParameter(3), f12->GetParameter(4), f12->GetParameter(5), f12->GetParameter(6),
			f12->GetParameter(7) );
			yield  = f11->Integral(2.8,3.4)/h->GetBinWidth(1);
      			yieldE = TMath::Sqrt(yield);
      			fitted = 1; 
    		} else {
    			h->Draw();
   			int ib = h->FindBin(3.1); 
      			yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      			yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      			yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
     			if (yield < 0) yield = 1.; 

      			yieldE =TMath::Sqrt(yield); 
      			fitted = 0; 
    		}

    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

	        getBinCenters(h->GetName(), eta, pt, Q); 
		nbin = fS1YieldPos->FindBin(eta, pt); 
    	        fS1YieldPos->SetBinContent(nbin, yield); 
   		fS1YieldPos->SetBinError(nbin, yieldE); 

    		c1->cd(2);  shrinkPad(0.15, 0.26); 
    		h = &(fS2VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f12, 4); 
   		if (h->GetSumOfWeights() > MINFIT) {
      		h->Fit(f12, fopt.c_str());
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else  if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		hChisq_mt->Fill(status*(f12->GetChisquare()/f12->GetNDF()));
		f11->SetParameters( f12->GetParameter(0), f12->GetParameter(1), f12->GetParameter(2), f12->GetParameter(3), f12->GetParameter(4), f12->GetParameter(5), f12->GetParameter(6),
		f12->GetParameter(7) );
		yield  = f11->Integral(2.8,3.4)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);
      		width  = f12->GetParameter(1); 
      		widthE = f12->GetParError(1); 
      		fitted = 1; 
    		} else {
      		int ib = h->FindBin(3.1); 
      		h->Draw();
      		yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      		yieldE =TMath::Sqrt(yield); 
      		yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      		yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
      		if (yield < 0) yield = 1.; 
      		width  = 0.030; 
      		width  = 0.010; 
      		fitted = 0; 
    		}	
    
    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
	  	if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS2YieldPos->FindBin(eta, pt); 
    		fS2YieldPos->SetBinContent(nbin, yield); 
    		fS2YieldPos->SetBinError(nbin, yieldE); 
	
    		c1->cd(3);  shrinkPad(0.15, 0.26); 
    		h = &(fS3VectorPos[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f12, 4); 
    		if (h->GetSumOfWeights() > MINFIT) {
      		h->Fit(f12, fopt.c_str());
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else  if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		hChisq_mmbar->Fill(status*(f12->GetChisquare()/f12->GetNDF()));
		f11->SetParameters( f12->GetParameter(0), f12->GetParameter(1), f12->GetParameter(2), f12->GetParameter(3), f12->GetParameter(4), f12->GetParameter(5), f12->GetParameter(6),
		f12->GetParameter(7) );
		yield  = f11->Integral(2.8,3.4)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);		
      		fitted = 1; 
    		} else {
    	 	h->Draw();
      		int ib = h->FindBin(3.1); 
      		yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      		yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      		yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
      		if (yield < 0) yield = 1.; 
      		yieldE =TMath::Sqrt(yield); 
      		fitted = 0; 
    		}	

    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
	 	<< h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS3YieldPos->FindBin(eta, pt); 
    		double mm = fS1YieldPos->GetBinContent(nbin); 
    		double mmE= fS1YieldPos->GetBinError(nbin); 
    		fS3YieldPos->SetBinContent(nbin, yield+mm); 
    		fS3YieldPos->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 

    		c1->Update();
    		c1->Draw();

    		TString frag(h->GetName()); 
    		frag.ReplaceAll("s3:mmbar,", ""); 
    		frag.ReplaceAll(",Q1", ""); 
    		frag.ReplaceAll(".", ":"); 
	
    		c1->SaveAs(Form("%s/massfits-pos-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data())); 
  		}


  		c1->Clear();
  		c1->Divide(3,1); 
  		// -- negative charge
  		for (unsigned int i = 0; i < fS1VectorNeg.size(); ++i) {
    		c1->cd(1);  shrinkPad(0.15, 0.26); 
    		h = &(fS1VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f12, 4); 
    		if (h->GetSumOfWeights() > MINFIT) {
    		h->Fit(f12, fopt.c_str());		
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else  if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		hChisq_mm->Fill(status*(f12->GetChisquare()/f12->GetNDF()));
		f11->SetParameters( f12->GetParameter(0), f12->GetParameter(1), f12->GetParameter(2), f12->GetParameter(3), f12->GetParameter(4), f12->GetParameter(5), f12->GetParameter(6),
		f12->GetParameter(7) );
		yield  = f11->Integral(2.8,3.4)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);		
    		} else {
      		h->Draw();
      		int ib = h->FindBin(3.1); 
      		yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      		yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      		yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
      		if (yield < 0) yield = 1.; 
      		yieldE =TMath::Sqrt(yield); 
    		}
    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
	 	<< h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS1YieldNeg->FindBin(eta, pt); 
    		fS1YieldNeg->SetBinContent(nbin, yield); 
    		fS1YieldNeg->SetBinError(nbin, yieldE); 

    		c1->cd(2);  shrinkPad(0.15, 0.26); 
    		h = &(fS2VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f12, 4); 
    		if (h->GetSumOfWeights() > MINFIT) {
    		h->Fit(f12, fopt.c_str());
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else  if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		hChisq_mt->Fill(status*(f12->GetChisquare()/f12->GetNDF()));
		f11->SetParameters( f12->GetParameter(0), f12->GetParameter(1), f12->GetParameter(2), f12->GetParameter(3), f12->GetParameter(4) , f12->GetParameter(5), f12->GetParameter(6),
		f12->GetParameter(7) );
		yield  = f11->Integral(2.8,3.4)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);		
      		width  = f12->GetParameter(1); 
      		widthE = f12->GetParError(1); 
    		} else {
      		h->Draw();
      		int ib = h->FindBin(3.1); 
      		yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      		yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      		yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
      		if (yield < 0) yield = 1.; 
      		yieldE =TMath::Sqrt(yield); 
      		width  = 0.030; 
      		widthE = 0.010; 
    		}

    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
		 << h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS2YieldNeg->FindBin(eta, pt); 
    		fS2YieldNeg->SetBinContent(nbin, yield); 
    		fS2YieldNeg->SetBinError(nbin, yieldE); 

    		c1->cd(3);  shrinkPad(0.15, 0.2); 
    		h = &(fS3VectorNeg[i]);
    		h->SetMinimum(0.); setTitles(h, "m_{#mu #bar{#mu} } [GeV]", "Entries/Bin", 0.08, 0.9, 1.8, 0.07);
    		setFunctionParameters(h, f12, 4); 
    		if (h->GetSumOfWeights() > MINFIT) {
      		h->Fit(f12, fopt.c_str());
		cout << gMinuit->fCstatu.Data() << endl;
		Status = gMinuit->fCstatu.Data();
		cout << Status[0] << endl;
		if ( Status[0] == 'S' || Status[0] == 'P' ){
	   	   	status = 1;
	      	} else  if ( Status[0] == 'F' ){ 
	      		status = -1;
			}
		cout << status << endl;
		hChisq_mmbar->Fill(status*(f12->GetChisquare()/f12->GetNDF()));
		f11->SetParameters( f12->GetParameter(0), f12->GetParameter(1), f12->GetParameter(2), f12->GetParameter(3), f12->GetParameter(4), f12->GetParameter(5), f12->GetParameter(6),
		f12->GetParameter(7) );
		yield  = f11->Integral(2.8,3.4)/h->GetBinWidth(1);
      		yieldE = TMath::Sqrt(yield);		
      		fitted = 1; 
    		} else {
      		h->Draw();
      		int ib = h->FindBin(3.1); 
      		yield  = h->Integral(ib - 0.2*h->GetNbinsX(),  ib + 0.2*h->GetNbinsX()); 
      		yield  -=h->Integral(h->GetNbinsX() - 0.2*h->GetNbinsX(), h->GetNbinsX()); 
      		yield  -=h->Integral(1, 0.2*h->GetNbinsX()); 
     		if (yield < 0) yield = 1.; 
     		yieldE =TMath::Sqrt(yield); 
     		fitted = 0; 
    		}

    		cout << " --> " << fitted << " " << h->GetSumOfWeights() << " .. " << yield << "  " 
	 	<< h->GetName() << "  " << yield << "+/-" << yieldE << endl;
    		if (PRINT) tl->DrawLatex(PRINTX, 0.85, Form("%4.1f", yield)); 

    		getBinCenters(h->GetName(), eta, pt, Q); 
    		nbin = fS3YieldNeg->FindBin(eta, pt); 
    		double mm = fS1YieldNeg->GetBinContent(nbin); 
    		double mmE= fS1YieldNeg->GetBinError(nbin); 
    		fS3YieldNeg->SetBinContent(nbin, yield+mm); 
    		fS3YieldNeg->SetBinError(nbin, TMath::Sqrt(yieldE*yieldE + mmE*mmE)); 

    		c1->Update();
    		c1->Draw();
	
    		TString frag(h->GetName()); 
    		frag.ReplaceAll("s3:mmbar,", ""); 
    		frag.ReplaceAll(",Q-1", ""); 
    		frag.ReplaceAll(".", ":"); 

    		c1->SaveAs(Form("%s/massfits-neg-%s-%s.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data()));
		 
	        TCanvas *c2 = new TCanvas("c2", "c2");
  	        c2->Clear();
  	        c2->Divide(3,1);
	        gStyle->SetOptStat(1111111111);
	        c2->cd(1);  shrinkPad(0.1, 0.2);
	        hChisq_mm->Draw();
	        c2->cd(2);  shrinkPad(0.1, 0.2);
	        hChisq_mt->Draw();
	        c2->cd(3);  shrinkPad(0.1, 0.2);	    
	        hChisq_mmbar->Draw();
	        c2->SaveAs(Form("chi2_polCBGauss.eps", fPtDirectory.c_str(), fSample.c_str(), frag.Data()));
		
  		}

  		c0->Clear(); 

  		c0->Divide(2,2);
  		c0->cd(1);
  		fS1YieldPos->Draw("colz");
  		c0->cd(2);
  		fS2YieldPos->Draw("colz");
	
  		c0->cd(3);
  		fS1YieldNeg->Draw("colz");
  		c0->cd(4);
  		fS2YieldNeg->Draw("colz");
	
  
  }
  
  	
}


// ----------------------------------------------------------------------
void anaXS::McpYields() {

  int nbin(0); 
  double yield(0.), yieldE(0.); 
  double pt(0.), eta(0.); 
  int Q(0); 

  TH1D *h; 

  cout << "fS1VectorMcpPos.size() = " << fS1VectorMcpPos.size() << endl;
  for (unsigned int i = 0; i < fS1VectorMcpPos.size(); ++i) {
    h = &(fS1VectorMcpPos[i]);
    yield  = h->GetSumOfWeights(); 
    yieldE = TMath::Sqrt(yield); 
    getBinCenters(h->GetName(), eta, pt, Q); 
    nbin = fS1McpPos->FindBin(eta, pt); 
    fS1McpPos->SetBinContent(nbin, yield); 
    fS1McpPos->SetBinError(nbin, yieldE); 

    h = &(fS1VectorMcpNeg[i]);
    yield  = h->GetSumOfWeights(); 
    yieldE = TMath::Sqrt(yield); 
    getBinCenters(h->GetName(), eta, pt, Q); 
    nbin = fS1McpNeg->FindBin(eta, pt); 
    fS1McpNeg->SetBinContent(nbin, yield); 
    fS1McpNeg->SetBinError(nbin, yieldE); 

    h = &(fS2VectorMcpPos[i]);
    yield  = h->GetSumOfWeights(); 
    yieldE = TMath::Sqrt(yield); 
    getBinCenters(h->GetName(), eta, pt, Q); 
    nbin = fS2McpPos->FindBin(eta, pt); 
    fS2McpPos->SetBinContent(nbin, yield); 
    fS2McpPos->SetBinError(nbin, yieldE); 

    h = &(fS2VectorMcpNeg[i]);
    yield  = h->GetSumOfWeights(); 
    yieldE = TMath::Sqrt(yield); 
    getBinCenters(h->GetName(), eta, pt, Q); 
    nbin = fS2McpNeg->FindBin(eta, pt); 
    fS2McpNeg->SetBinContent(nbin, yield); 
    fS2McpNeg->SetBinError(nbin, yieldE); 

  }

}


// ----------------------------------------------------------------------
void anaXS::totalMass() {

  TH1D *h1 = &(fS1VectorPos[0]); 
  fS1TotalPos = new TH1D("S1TotalPos", "", h1->GetNbinsX(), h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()+1)); 
  fS1TotalNeg = new TH1D("S1TotalNeg", "", h1->GetNbinsX(), h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()+1)); 
  fS2TotalPos = new TH1D("S2TotalPos", "", h1->GetNbinsX(), h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()+1)); 
  fS2TotalNeg = new TH1D("S2TotalNeg", "", h1->GetNbinsX(), h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX()+1)); 

  
  for (unsigned int i = 0; i < fS1VectorPos.size(); ++i) {
    fS1TotalPos->Add(&(fS1VectorPos[i]));
    fS1TotalNeg->Add(&(fS1VectorNeg[i]));
  }

  for (unsigned int i = 0; i < fS2VectorPos.size(); ++i) {
    fS2TotalPos->Add(&(fS2VectorPos[i]));
    fS2TotalNeg->Add(&(fS2VectorNeg[i]));
  }
    
  zone(1); 
  shrinkPad(0.20, 0.25);
  fS1TotalPos->SetMinimum(0.); 
  setHist(fS1TotalPos, kBlue, 24); 
  setHist(fS1TotalNeg, kRed, 25); 

  setTitles(fS1TotalPos, "m_{#mu #mu} [GeV]", "Entries/Bin", 0.08, 1.2, 1.7, 0.07); 

  fS1TotalPos->Draw("e");
  fS1TotalNeg->Draw("esame");

  legg = new TLegend(0.7,0.6,0.85,0.85);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05); legg->SetTextFont(132); 
  legge = legg->AddEntry(fS1TotalPos, " positive tag muons","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(fS2TotalNeg,  " negative tag muons","p"); legge->SetTextColor(kBlack);
  legg->Draw();

  c0->SaveAs(Form("%s/totalMass-s1-%s.eps", fPtDirectory.c_str(), fSample.c_str())); 


  zone(1); 
  shrinkPad(0.20, 0.25);
  fS2TotalPos->SetMinimum(0.); 
  setHist(fS2TotalPos, kBlue, 24); 
  setHist(fS2TotalNeg, kRed, 25); 

  setTitles(fS2TotalPos, "m_{#mu t} [GeV]", "Entries/Bin", 0.08, 1.0, 1.7, 0.07); 

  fS2TotalPos->Draw("e");
  fS2TotalNeg->Draw("esame");

  legg = new TLegend(0.4,0.2,0.8,0.4);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05); legg->SetTextFont(132); 
  legge = legg->AddEntry(fS2TotalPos, " positive tag muons","p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(fS2TotalNeg,  " negative tag muons","p"); legge->SetTextColor(kBlack);
  legg->Draw();

  c0->SaveAs(Form("%s/totalMass-s2-%s.eps", fPtDirectory.c_str(), fSample.c_str())); 

  
}

// ----------------------------------------------------------------------
void anaXS::validation() {

  zone(2,2);
  fPtMctNeg->drawEff();
  c0->cd(2);
  fPtMctPos->drawEff();
  c0->cd(3);
  fPtMcpNeg->drawEff();
  c0->cd(4);
  fPtMcpPos->drawEff();

  // -- 
  makeCanvas(4);
  c3->cd();
  c3->Clear();
  c3->Divide(2,2);
  c3->cd(1);
  fS1McpPos->Draw("colz");

  c3->cd(2);
  fS2McpPos->Draw("colz");

  c3->cd(3);
  fS1McpNeg->Draw("colz");

  c3->cd(4);
  fS2McpNeg->Draw("colz");

  // -- 
  makeCanvas(16);
  c5->cd();
  c5->Clear();
  c5->Divide(2,2);

  c5->cd(1);
  fPtTnpPos->drawEff();

  c5->cd(2);
  fPtTnpNeg->drawEff();

  c5->cd(3);
  fPtMcpPos->drawEff();

  c5->cd(4);
  fPtMcpNeg->drawEff();


}

// ----------------------------------------------------------------------
void anaXS::biasPlots(const char *fname, const char *psname,  int mode) {

  TFile *f = TFile::Open(fname);

  //   TFile *fPsi = TFile::Open("/data/ursl/muonID-rootfiles/biased/upsilon/ups_biasstudy.root");
  //   TFile *fPsi = TFile::Open("/data/ursl/muonID-rootfiles/biased/jpsi/jpsi_biasstudy.root");

  f->cd();

  TH1D *ht = new TH1D("ht", "", 48, -2.4, 2.4); 

  TH1D *h0 = new TH1D("h1", "", 48, -2.4, 2.4); setFilledHist(h0, kBlack, kYellow, 1000, 2);
  TH1D *h1 = new TH1D("h1", "", 48, -2.4, 2.4); setFilledHist(h1, kRed, kRed, 3005, 2);
  TH1D *h2 = new TH1D("h2", "", 48, -2.4, 2.4); setFilledHist(h2, kCyan, kCyan, 1000, 2);
  TH1D *h3 = new TH1D("h3", "", 48, -2.4, 2.4); setFilledHist(h3, kBlue, kBlue, 3004, 2);
  TH1D *h4 = new TH1D("h4", "", 48, -2.4, 2.4); setFilledHist(h4, kBlack, 1, 3006, 2);
  TH1D *h5 = new TH1D("h5", "", 48, -2.4, 2.4); setFilledHist(h5, kMagenta, 1, 3006, 2);

  setTitles(h0, "#eta", "muons/bin", 0.06, 1.2, 1.7); 
  h0->SetMinimum(0);   

  TTree *t = (TTree*)f->Get("events"); 

  float pe, pp, te, tp; 
  t->SetBranchAddress("pe", &pe);
  t->SetBranchAddress("te", &te);
  t->SetBranchAddress("pp", &pp);
  t->SetBranchAddress("tp", &tp);

  Int_t nentries = Int_t(t->GetEntries());

  Int_t nb(0), nbytes(0);
  //  for (Int_t jentry = 0; jentry < 20000; jentry++) {
  for (Int_t jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);   nbytes += nb;
    
    //    if (pp < 3 || pp > 6) continue; 
    //    if (pp > 3 && pp < 6) continue; 

    ht->Fill(te); 
    h0->Fill(pe); 
    if (te*pe > 0.) { h1->Fill(pe); }
    if (TMath::Abs(te - pe) < 0.3) { h2->Fill(pe); }
    if (tp > 4) { h3->Fill(pe); }
    if (tp > 7) { h4->Fill(pe); }
    if (tp > 6) { h5->Fill(pe); }
    //    if (tp < 4 && TMath::Abs(te - pe)<0.2) { h4->Fill(pe); }
    //    if (TMath::Abs(pp-tp) < 3) { h4->Fill(pe); }
    //    if (TMath::Abs(te) < 2.2 && tp<5 && te*pe > 0.) { h4->Fill(pe); }
  }  


  //  if (TMath::Abs(te) - TMath::Abs(pe) < 0.2 && tp>4.) { h4->Fill(pe); }
  //  if (TMath::Abs(te - pe) < 0.1) { h2->Fill(pe); }
  //  if (TMath::Abs(te - pe) < 0.2) { h2->Fill(pe); }
  //  if (te > 5.0) { h2->Fill(pe); }
  //  if (te > 6.0) { h2->Fill(pe); }
  //  if (te > 7.0) { h2->Fill(pe); }

  h0->SetMaximum(1.4*h0->GetMaximum()); 

  gStyle->SetOptStat(0); 
  zone();
  shrinkPad(0.15, 0.25); 
  h0->Draw("");
  ht->Draw("esame");
  if (mode & 1) h1->Draw("same");
  if (mode & 2) h2->Draw("same");
  if (mode & 4) h3->Draw("same");
  if (mode & 8) h4->Draw("same");



  legg = new TLegend(0.47, 0.7, 0.6, 0.88);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.025); legg->SetTextFont(132); 
  legge = legg->AddEntry(ht, "all tags", "p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(h0, "all probes", "f"); legge->SetTextColor(kBlack);
  if (mode & 1) legge = legg->AddEntry(h1, " #eta_{tag} #times #eta_{probe} > 0", "f"); legge->SetTextColor(kBlack);
  if (mode & 2) legge = legg->AddEntry(h2, " |#eta_{tag} - #eta_{probe}| < 0.3", "f"); legge->SetTextColor(kBlack);
  if (mode & 4) legge = legg->AddEntry(h3, " p_{T, tag} > 4 GeV", "f"); legge->SetTextColor(kBlack);
  if (mode & 8) legge = legg->AddEntry(h4, " p_{T, tag} > 7 GeV", "f"); legge->SetTextColor(kBlack);

  cout << "=============> " << f->GetName() << endl;
  cout << "ht: " << ht->GetSumOfWeights() << endl;
  cout << "h0: " << h0->GetSumOfWeights() << endl;
  cout << "h1: " << h1->GetSumOfWeights() << endl;
  cout << "h2: " << h2->GetSumOfWeights() << endl;
  cout << "h3: " << h3->GetSumOfWeights() << endl;
  cout << "h4: " << h4->GetSumOfWeights() << endl;
  cout << "h5: " << h5->GetSumOfWeights() << endl;
  cout << "==========================" << endl;

  legg->Draw();


  c0->SaveAs(Form("%s/biasPlots-%s", fPtDirectory.c_str(), psname)); 

}

// ----------------------------------------------------------------------
bool anaXS::GetBinCenters(string hname, double &eta, double &pT) {
  bool ok(false); 

  float etamin, etamax, ptmin, ptmax; 
  int   n; 
  char searchString[2000], hchar[2000]; 

  size_t ipos = hname.find(string("rapidity")); 
  sprintf(hchar, "%s", hname.substr(ipos).c_str()); 

  char sp[] = "%"; 
  char sf[] = "f"; 
  sprintf(searchString, "rapidity%s%s_%s%s,pt%s%s_%s%s", sp, sf, sp, sf, sp, sf, sp, sf, sp); 
 
  n = sscanf(hchar, searchString, &etamin, &etamax, &ptmin, &ptmax);

  if (n > 0) {
    ok = true; 
    eta = 0.5*(etamin + etamax);
    pT  = 0.5*(ptmin + ptmax);
  } else {
    cout << "could not parse " << hname << endl;
  }
	 

  return ok;
}




// ----------------------------------------------------------------------
bool anaXS::getBinCenters(string hname, double &eta, double &pT, int &Q) {
  bool ok(false); 

  float etamin, etamax, ptmin, ptmax; 
  int   charge, n; 
  char searchString[2000], hchar[2000]; 

  size_t ipos = hname.find(string("eta")); 
  sprintf(hchar, "%s", hname.substr(ipos).c_str()); 

  char sp[] = "%"; 
  char sf[] = "f"; 
  char sd[] = "d"; 
  sprintf(searchString, "eta%s%s_%s%s,pt%s%s_%s%s,Q%s%s", sp, sf, sp, sf, sp, sf, sp, sf, sp, sd); 
 
  n = sscanf(hchar, searchString, &etamin, &etamax, &ptmin, &ptmax, &charge);

  if (n > 0) {
    ok = true; 
    eta = 0.5*(etamin + etamax);
    pT  = 0.5*(ptmin + ptmax);
    Q   = charge;
  } else {
    cout << "could not parse " << hname << endl;
  }
	 

  return ok;
}



// ----------------------------------------------------------------------
void anaXS::setFunctionParameters(TH1D *h, TF1 *f, int mode) {
  const int EDG(4), NB(EDG+1); 
  double p0(0.), p1(0.), g0(0.), g1(0.), g2(0.), g3(0.), g4(0.), g5(0.), g6(0.);
  double c0(0.), c1(0.), c2(0.), c3(0.), c4(0.), c5(0.), c6(0.), c7(0.), c8(0.), c9(0.), c10(0.); 
  double lo(0.), hi(0.), dx(0.);
  double lo_2(0.), hi_2(0.);
  double p0E(0.), p1E(0.); 

  // -- (background) pol1 parameters from end points
  lo = h->Integral(1, 1+EDG); 
  hi = h->Integral(h->GetNbinsX()-EDG, h->GetNbinsX());
  dx = h->GetBinLowEdge(h->GetNbinsX()+1) - h->GetBinLowEdge(1);
  p1 = (hi-lo)/NB/dx;
  p0 = lo/NB - p1*h->GetBinLowEdge(1);
  f0->SetParameters(p0, p1); 
  p0E = TMath::Sqrt(lo)/lo; 
  p1E = TMath::Sqrt(hi)/hi; 
  p1E = TMath::Sqrt(p1E*p1E + p0E*p0E); 

  // -- pol1
  if (0 == mode) {
    f->SetParameters(p0, p1); 
  }
  
  // -- pol1 + Gauss for J/Psi
  if (1 == mode) {
    g1 = 3.1;
    g2 = 0.04;
    lo = h->GetMaximumBin()-2.*g2/h->GetBinWidth(1); 
    hi = h->GetMaximumBin()+2.*g2/h->GetBinWidth(1);
    g0 = (h->Integral(lo, hi) - f0->Integral(lo, hi))*h->GetBinWidth(1);  
    if (g0 < 0) g0 = h->Integral(lo, hi)*h->GetBinWidth(1);
    f->ReleaseParameter(0); 
    f->ReleaseParameter(1); 
    f->ReleaseParameter(2); 
    f->ReleaseParameter(3); 
    f->ReleaseParameter(4); 

    f->SetParameters(g0, g1, g2, p0, p1); 
    f->SetParLimits(1, 3.050, 3.150); 
    f->SetParLimits(2, 0.015, 0.040); 
  }


  // -- pol1 + Gauss for All three Upsilons
  if (2 == mode) {
    g1 = 9.45;
    g2 = 0.06;
    lo = h->GetMaximumBin()-2.*g2/h->GetBinWidth(1); 
    hi = h->GetMaximumBin()+2.*g2/h->GetBinWidth(1);
    g0 = (h->Integral(lo, hi) - f0->Integral(lo, hi))*h->GetBinWidth(1);  
    if (g0 < 0) g0 = h->Integral(lo, hi)*h->GetBinWidth(1);

    g3 = 0.6*g0; 
    g4 = 10.0;
    g5 = 0.08;

    g6 = 0.1*g0; 

    f->ReleaseParameter(0); 
    f->ReleaseParameter(1); 
    f->ReleaseParameter(2); 
    f->ReleaseParameter(3); 
    f->ReleaseParameter(4); 
    f->ReleaseParameter(5); 
    f->ReleaseParameter(6); 
    f->ReleaseParameter(7); 
    f->ReleaseParameter(8); 

    f->SetParameters(g0, g1, g2, g3, g4, g5, g6, p0, p1); 
    f->SetParLimits(0, 0., 10000000.); 
    f->SetParLimits(2, 0.04, 0.14); 
    f->SetParLimits(3, 0., 10000000.); 
    f->SetParLimits(4, 9.9, 10.1); 
    f->SetParLimits(5, 0.06, 0.20); 
    f->SetParLimits(6, 0., 10000000.); 

  }
  
  // -- pol1 + CB for J/Psi
  if (3 == mode) {
    c0 = 3.1;
    c1 = 0.03;
    c2 = 1.5;
    lo = h->GetMaximumBin()-2.*c1/h->GetBinWidth(1); 
    hi = h->GetMaximumBin()+2.*c1/h->GetBinWidth(1);
    c4 = (h->Integral(lo, hi) - f0->Integral(lo, hi))*h->GetBinWidth(1);
    if (c4 < 0) c4 = h->Integral(lo, hi)*h->GetBinWidth(1);
    f->ReleaseParameter(0); 
    f->ReleaseParameter(1); 
    f->ReleaseParameter(2); 
    f->ReleaseParameter(3); 
    f->ReleaseParameter(4);
    f->ReleaseParameter(5); 
    f->ReleaseParameter(6);     

    f->SetParameters(c0, c1, c2, c3, c4, p0, p1);     
    f->SetParLimits(0, 3.050, 3.150); 
    f->SetParLimits(1, 0.015, 0.07);
    f->FixParameter(2, 1.5);
    f->FixParameter(3, 3.);
     
  }
    
  // -- pol1 + CB + Gauss for J/Psi
  if (4 == mode) {
    g1 = 3.1;
    g2 = 0.04;
    lo = h->GetMaximumBin()-2.*g2/h->GetBinWidth(1); 
    hi = h->GetMaximumBin()+2.*g2/h->GetBinWidth(1);
    g0 = (h->Integral(lo, hi) - f0->Integral(lo, hi))*h->GetBinWidth(1);  
    if (g0 < 0) g0 = h->Integral(lo, hi)*h->GetBinWidth(1);
    c0 = g1;
    c1 = g2;
    c2 = 10;
    c3 = 0.03;
    c4 = g0;
    f->ReleaseParameter(0); 
    f->ReleaseParameter(1); 
    f->ReleaseParameter(2); 
    f->ReleaseParameter(3); 
    f->ReleaseParameter(4);
    f->ReleaseParameter(5); 
    f->ReleaseParameter(6); 
    f->ReleaseParameter(7); 
    f->ReleaseParameter(8); 
    f->ReleaseParameter(9);     

    f->SetParameters(c0, c1, c2, c3, c4, g0, g1, g2, p0, p1);
    f->SetParLimits(0, 3.000, 3.200); 
    f->SetParLimits(1, 0.015, 0.075);
    f->SetParLimits(2, 1., 100);
    f->SetParLimits(3, 0.1, 0.9);
    f->SetParLimits(4, 0., g0*h->GetBinWidth(1));     
    f->SetParLimits(6, 3.050, 3.150); 
    f->SetParLimits(7, 0.015, 0.040); 
  }
  
  // -- pol1 + CB for Ups(1S)
  if (5 == mode) {
    c0 = 9.46;				    
    c1 = 0.06;
    //c0 = 10.02;				    
    //c1 = 0.08;    
    //c0 = 10.35;				    
    //c1 = 0.10;        
    c2 = 1.9;
    c3 = 1.6;
    lo = h->GetMaximumBin()-2.*c1/h->GetBinWidth(1); 
    hi = h->GetMaximumBin()+2.*c1/h->GetBinWidth(1);
    c4 = (h->Integral(lo, hi) - f0->Integral(lo, hi))/h->GetBinWidth(1);
    if (c4 < 0) c4 = h->Integral(lo, hi)*h->GetBinWidth(1);
    f->ReleaseParameter(0); 
    f->ReleaseParameter(1); 
    f->ReleaseParameter(2); 
    f->ReleaseParameter(3); 
    f->ReleaseParameter(4);
    f->ReleaseParameter(5); 
    f->ReleaseParameter(6);     
    
    f->SetParameters(c0, c1, c2, c3, c4, p0, p1);     
    f->SetParLimits(0, 9.410, 9.510); 
    f->SetParLimits(1, 0.04, 0.14);
    //f->SetParLimits(0, 9.9, 10.1); 
    //f->SetParLimits(1, 0.06, 0.20);      
    //f->SetParLimits(0, 10.22, 10.43); 
    //f->SetParLimits(1, 0.06, 0.20);    
    //f->SetParLimits(2, 1.6, 2.2);
    //f->SetParLimits(3, 1., 2.2);
    f->FixParameter(2, 1.9);
    f->FixParameter(3, 1.6);
    f->SetParLimits(4, 0, 10000000);
    f->FixParameter(5, 0.);
    f->FixParameter(6, 0.);
  }
  
    // -- pol1 + CB for All Upsilons
  if (6 == mode) {
    c0 = 9.46;
    c1 = 0.08; // 0.06 -> 0.08
    c2 = 1.9;
    c3 = 1.6;
    lo = h->GetMaximumBin()-2.*c1/h->GetBinWidth(1); 
    hi = h->GetMaximumBin()+2.*c1/h->GetBinWidth(1);
    c4 = (h->Integral(lo, hi) - f0->Integral(lo, hi))/h->GetBinWidth(1);    //  * -> /
    if (c4 < 0) c4 = h->Integral(lo, hi)*h->GetBinWidth(1);
    c5 = 10.02;
    c6 = 0.10;
    h->GetXaxis()->SetRangeUser(9.8,10.2);
    lo_2 = h->GetMaximumBin()-2.*c6/h->GetBinWidth(1); 
    hi_2 = h->GetMaximumBin()+2.*c6/h->GetBinWidth(1);
    c7 = (h->Integral(lo_2, hi_2) - f0->Integral(lo_2, hi_2))/h->GetBinWidth(1);
    h->GetXaxis()->SetRangeUser(8.7,11.2);
    c8 = c7;  
    f->ReleaseParameter(0); 
    f->ReleaseParameter(1); 
    f->ReleaseParameter(2); 
    f->ReleaseParameter(3); 
    f->ReleaseParameter(4);
    f->ReleaseParameter(5); 
    f->ReleaseParameter(6);
    f->ReleaseParameter(7); 
    f->ReleaseParameter(8);
    f->ReleaseParameter(9); 
    f->ReleaseParameter(10);         

    f->SetParameters(c0, c1, c2, c3, c4, c5, c6, c7, c8, p0, p1);     
    f->SetParLimits(0, 9.410, 9.510); 
    f->SetParLimits(1, 0.04, 0.14);
    f->FixParameter(2, 1.9);
    //f->SetParLimits(2, 1.6, 2.2);
    f->FixParameter(3, 1.6);
    //f->SetParLimits(3, 1., 2.2);
    f->SetParLimits(4, 0, 10000000);
    f->SetParLimits(5, 9.9, 10.1);
    f->SetParLimits(6, 0.06, 0.20);
    f->SetParLimits(7, 0, 10000000);
    f->SetParLimits(8, 0, 10000000);
    //f->FixParameter(9, 0.);
    //f->FixParameter(10, 0.);
     
  }
  
  
  
}


// ----------------------------------------------------------------------
void anaXS::makeCanvas(int i) {
  if (i & 16) { 
    c5 = new TCanvas("c5", "c5", 210,   0, 800, 1000);
    c5->ToggleEventStatus();
  }
  if (i & 8) { 
    c4 = new TCanvas("c4", "c4", 210,   0, 800, 600);
    c4->ToggleEventStatus();
  }
  if (i & 4) {
    c3 = new TCanvas("c3", "c3", 200,  20, 800, 800);
    c3->ToggleEventStatus();
  }
  if (i & 1) {
    //    c1 = new TCanvas("c1", "c1", 20,  60, 1200, 400);
    c1 = new TCanvas("c1", "c1", 20,  60, 1000, 400);
    c1->ToggleEventStatus();
  }
  if (i & 2) { 
    c2 = new TCanvas("c2", "c2", 300, 200, 800, 400);
    c2->ToggleEventStatus();
  }
}

// ----------------------------------------------------------------------
int  anaXS::wait() {
  cout << " Continue [<RET>|q]?  "; 
  char x;
  x = getchar();
  if ((x == 'q') || (x == 'Q')) return 1;
  return 0;
}
