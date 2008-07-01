#include "MyMakePlots.h"

#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TObject.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TKey.h"
#include "TImageDump.h"
#include "TLegend.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>

std::vector<std::string> getAllKeys (const TDirectory* fDir, const std::string& fClassName) {
  std::cout << "getAllKeys-> " << fDir->GetName() << ", " <<  fClassName << std::endl;
  //  fDir->ls();
  std::vector<std::string> result;
  TIter next (fDir->GetListOfKeys ());
  for (TKey* key = 0; (key = (TKey *) next());) {
    std::cout << "key from list: " << key->GetName()  << '/' << key->GetClassName () << std::endl;
    if (fClassName == key->GetClassName ()) {
      result.push_back (std::string (key->GetName ()));
      //std::cout << key->GetName() << std::endl;
    } 
  }
  return result;
}

TObject* getObject (TDirectory* fDir, const std::vector <std::string>& fObjectName) {
  TObject* result = 0; // nothing so far
  TDirectory* dir = fDir;
  for (unsigned i = 0; i < fObjectName.size (); ++i) {
    dir->GetObject (fObjectName[i].c_str(), result);
    if (result) {
      if (i < fObjectName.size () - 1) {
	dir = (TDirectory*) result;
	result = 0;
      }
    }
    else {
      std::cerr << "[make_plots] getObject-> Can not find (sub)dir/object " << fObjectName[i] << " in directory " << dir->GetName () << std::endl;
      return 0;
    }
  }
  return result;
}

TH1* HistoCompare(TH1 *h, TH1 *refHisto) {

  std::cout << "Comparing " << h->GetName() << " and " << refHisto->GetName() << std::endl;
  std::cout << " entries " << h->GetEntries() << " and " << refHisto->GetEntries() << std::endl;

  if (h->GetEntries() == 0) {
    TH1 *resHisto = (TH1*) refHisto->Clone(TString("Reference"));
    resHisto->SetTitle("Reference");
    return resHisto;
  }
  //create a residual histogram
  //if (resHisto_) delete resHisto_;
  TH1 *resHisto = (TH1*) h->Clone(TString(h->GetName())+"_residuals");
  resHisto->Reset();

  // clone input histogram
  TH1 *htemp = (TH1*) h->Clone(TString(h->GetName()));

  
  // normalize histograms
  //htemp->Sumw2();
  //refHisto->Sumw2();

  //htemp->GetBinContent(25);
  //  ******Temp for Higgs plots ***** ///
  //refHisto->Scale(htemp->GetBinContent(25)/refHisto->GetBinContent(25));


  if (refHisto->Integral() != 0) refHisto->Scale(1./refHisto->Integral());
  if (htemp->Integral() != 0 ) htemp->Scale(1./htemp->Integral());
  refHisto ->SetMarkerStyle(1);
        
  resHisto->Add( htemp, refHisto, 1., -1.);

  resHisto->SetTitle("Residuals");
  resHisto->SetYTitle("");
 
  return resHisto;
  
}

MyMakePlots::MyMakePlots() {

	extension = "png";
	compare = false;
	compare_filename = "";
	logaxis = false;
	
}

void MyMakePlots::Draw() {
  std::cout << "MyMakePlots" << std::endl;
  // get a file
  TFile *rel_file = new TFile(root_filename);
  TFile *ref_file = 0;

  ofstream HistFail(fileprefix+"HistFail.txt");


  //HistFail.open(fileprefix+"HistFail.txt");
	
  std::string HistoType[6] =  {"TH1D", "TH1F", "TH2F", "TH3F", "TProfile", "TProfile2D"} ;
  std::vector<std::string> dirList1 = getAllKeys(rel_file, "TDirectory");
  std::cout << int(dirList1.size()) << std::endl;
 for(unsigned int q = 0 ; q < 6; ++q)
    {
  if(dirList1.size() == 0){
    std::vector<std::string> histKeys = getAllKeys(gDirectory, HistoType[q]);
    std::cout << histKeys.size() << std::endl;
    for(unsigned ihist = 0; ihist < histKeys.size(); ++ihist)
      {
	TH1* hist = 0;
	gDirectory->GetObject(histKeys[ihist].c_str(), hist);
	if(hist)
	  {
	    std::vector<std::string>histPathName;
	    histPathName.push_back(histKeys[ihist]);
	    TH1* hist1 = (TH1*)getObject(rel_file, histPathName);
	    TString thename = hist1->GetTitle();
	    
	    TLatex *label = new TLatex(0.01,0.01,thename);
	    label->SetNDC();
	    label->SetTextColor(15);
	    label->SetTextSize(0.02);
	    TCanvas *canvas1 = new TCanvas(TString(hist->GetName()), TString(hist->GetName()), 800,800);
	    //if (compare){
	    canvas1->Divide(1,2);
	    canvas1->cd(1);
	    TH1 *htemp1 = (TH1*) hist1->Clone(TString(hist1->GetName()));
	    htemp1->Sumw2();
	    if (htemp1->Integral() != 0 ) htemp1->Scale(1./htemp1->Integral());
	    htemp1->SetLineColor(kRed);
	    htemp1->SetLineWidth(4);
	    htemp1->Draw();
	    canvas1->Update();
	    if (compare){
	    ref_file = new TFile(compare_filename);
	    TH1* ref_hist = (TH1*)getObject(ref_file, histPathName);
	    if(ref_hist)
	      {
		ref_hist->Sumw2();
		if (ref_hist->Integral() != 0) ref_hist->Scale(1./ref_hist->Integral());
		//canvas1->cd(1);
		ref_hist->SetLineColor(kBlue);
		ref_hist->SetLineWidth(2);
		ref_hist->SetLineWidth(2);
		ref_hist->Draw("SAME");
		canvas1->Update();
		double Chi2Test = 0.; double KSTest = 0.;
		Chi2Test = ref_hist->Chi2Test(hist1,"UU");
		KSTest = ref_hist->KolmogorovTest(hist1,"UO");
		TH1* res_hist = HistoCompare(htemp1, ref_hist);
		TLegend *leg = new TLegend(0.1,0.8,0.25,0.9,"","NDC");
		leg->AddEntry(htemp1, "Rel", "L");
		leg->AddEntry(ref_hist, "Ref", "L");
		leg->Draw();
		canvas1->cd(2);
		res_hist->Draw("");
		
		//TPad *npad = new TPad("npad", "", 0.1, 0.5, 0.4, 0.8);
		//npad->Draw();
		//npad->cd();
		//res_hist->Draw();
		
		char buf[1024];
		sprintf (buf, "#chi^{2}= %1.2f", Chi2Test);
		TLatex *labelchi2 = new TLatex(0.4,0.91,buf);
		labelchi2->SetNDC();
		labelchi2->SetTextSize(0.035);
		labelchi2->Draw();
		canvas1->Update();
		sprintf (buf, "KS= %2.2f", KSTest);
		TLatex *labelkg = new TLatex(0.6,0.91,buf);
		labelkg->SetNDC();
		if (Chi2Test < 0.8 && KSTest < .8) {gPad->SetFillColor(2); HistFail << histKeys[ihist] << ":" << root_filename << ":" << compare_filename << std::endl;}
		labelkg->SetTextSize(0.035);
		labelkg->Draw();
	
	      }else{
		std::cout << " no reference plot " << histKeys[ihist] << std::endl;
		TLatex *label = new TLatex(0.2,0.5,"no reference plot available");
		label->SetNDC();
		label->SetTextSize(0.07);
		label->Draw();
	      }
	    }else{
	       hist->Draw();
	       //label->Draw();
	    }
	    canvas1->cd();
	    canvas1->Print(webpath+"/"+TString(hist->GetName())+"."+extension);
	  }else{
	  std::cout << "Can not get histograms" << std::endl;
	  }
      }
  }else
    {
    
    for(unsigned idir = 0; idir < dirList1.size(); ++idir)
      {
	TDirectory* dir1 = 0;
	rel_file->GetObject(dirList1[idir].c_str(), dir1);
	if(dir1)
	  {
	    std::vector<std::string> histKeys = getAllKeys(dir1, HistoType[q]);
	    if(histKeys.size() != 0)
	      {
		for(unsigned ihist = 0; ihist < histKeys.size(); ++ihist)
		{
		  TH1* hist = 0;
		  dir1->GetObject(histKeys[ihist].c_str(), hist);
		  if(hist)
		    {
		      TCanvas *canvas1 = new TCanvas(TString(hist->GetName()), TString(hist->GetName()), 800,800);
		      if (compare){
			canvas1->Divide(1,2);
			canvas1->cd();
			//hist->SetLineColor(kBlue);
			//hist->SetLineSize(2);
			hist->Draw();
			
			std::vector<std::string>histPathName;
			histPathName.push_back(histKeys[ihist]);
			TH1* ref_hist = (TH1*)getObject(ref_file, histPathName);
			if(ref_hist)
			  {
			    //ref_hist->SetLineColor(kRed);
			    //ref_hist->SetLineWidth(2);
			    ref_hist->Draw("SAME");
			    double Chi2Test = 0.; double KSTest = 0.;
			    Chi2Test = ref_hist->Chi2Test(hist, "UFOF");
			    KSTest = ref_hist->KolmogorovTest(hist, "UO");
			    TH1* res_hist = HistoCompare(hist, ref_hist);
			    canvas1->cd(2);
			    res_hist->Draw();
			    char buf[1024];
			    sprintf (buf, "#chi^{2}= %1.2f", Chi2Test);
			    TLatex *labelchi2 = new TLatex(0.4,0.91,buf);
			    labelchi2->SetNDC();
			    labelchi2->SetTextSize(0.035);
			    labelchi2->Draw();
			    canvas1->Update();
			    sprintf (buf, "KG= %2.2f", KSTest);
			    TLatex *labelkg = new TLatex(0.6,0.91,buf);
			    labelkg->SetNDC();
			    if (KSTest<0.8) {gPad->SetFillColor(2); HistFail << histKeys[ihist] << ":" << root_filename << ":" << compare_filename << std::endl;}
			    labelkg->SetTextSize(0.035);
			    labelkg->Draw();
			  }else{
			    std::cout << " no reference plot " << histKeys[ihist] << std::endl;
			    TLatex *label = new TLatex(0.2,0.5,"no reference plot available");
			    label->SetNDC();
			    label->SetTextSize(0.07);
			    label->Draw();
			  }
		      }else{
			hist->Draw();
			//label->Draw();
		      }
		      canvas1->cd();
		      canvas1->Print(webpath+"/"+TString(hist->GetName())+"."+extension);
		    }
		}
	      }
	  }
      
      }
    }
    
    
    }
  HistFail.close();
      
}

		
	   
			     
  
