#ifndef STACK_COMMON_H
#define STACK_COMMON_H

#include <iostream>

const string zmmName = "Analysis_zmm_2_4.root";     
const string wmnName = "Analysis_wmn_2_4.root";     
const string ttbarName = "Analysis_ttbar_2_4.root"; 
const string qcdName = "Analysis_qcd_2_4_all.root";  
const string dataName = "data13XXXX.root";  

const int canvasSizeX = 500;
const int canvasSizeY = 500;

const Color_t zLineColor = kAzure+2;
const Color_t zFillColor = kAzure+6;
const Color_t wLineColor = kAzure+4;
const Color_t wFillColor = kAzure+2;
const Color_t qcdLineColor = kGreen+3;
const Color_t qcdFillColor = kGreen+1;
const Color_t ttFillColor =  kRed+1;
const Color_t ttLineColor = kRed+3;

const double lumi =.017 ;
const double lumiZ = 10. ;
const double lumiW = 10.;
const double lumiQ = 10.;
const double lumiT =10.;

const double mMin = 60;
const double mMax = 120;

TFile z(zmmName.c_str()) ; 
TFile w( wmnName.c_str()) ; 
TFile tt(ttbarName.c_str()) ; 
TFile qcd(qcdName.c_str()) ; 
TFile data(dataName.c_str()) ; 

TCanvas *c1 = new TCanvas("c1","Stack plot",10,10,canvasSizeX, canvasSizeY);

void setHisto(TH1 * h, Color_t fill, Color_t line, double scale, int rebin) {
  h->SetFillColor(fill);
  h->SetLineColor(line);
  h->Scale(scale);
  h->Rebin(rebin);  
}

void stat(TH1 * h1, TH1 * h2, TH1 * h3, TH1 * h4, TH1 * hdata, int rebin) {
  double a = mMin/rebin, b = mMax/rebin;
  double i1 = h1->Integral(a, b);
  double i2 = h2->Integral(a, b);
  double i3 = h3->Integral(a, b);
  double i4 = h4->Integral(a, b);
  double idata = hdata != 0 ? hdata->Integral(mMin, mMax) : 0;
  std::cout.setf(0,ios::floatfield);
  std::cout.setf(ios::fixed,ios::floatfield);
  std::cout.precision(1);
  std::cout <<"Zmm (" << mMin << ", " << mMax << ") = ";
  std::cout.precision(8);
  std::cout << i1 <<std::endl;
  std::cout.precision(1);
  std::cout <<"Wmn (" << mMin << ", " << mMax << ") = ";
  std::cout.precision(8);
  std::cout << i2 <<std::endl;
  std::cout.precision(1);
  std::cout <<"QCD (" << mMin << ", " << mMax << ") = ";
  std::cout.precision(8);
  std::cout << i4 <<std::endl;
  std::cout.precision(1);
  std::cout <<"tt~ (" << mMin << ", " << mMax << ") = ";
  std::cout.precision(8);
  std::cout << i3 <<std::endl; 
  std::cout.precision(1);
  std::cout <<"dat (" << mMin << ", " << mMax << ") = ";
  std::cout.precision(8);
  std::cout << idata <<std::endl; 
}

void makeStack(TH1 * h1, TH1 * h2, TH1 * h3, TH1 * h4, TH1 * hdata,
	       double min, int rebin) {
  setHisto(h1, zFillColor, zLineColor, lumi/lumiZ, rebin);
  setHisto(h2, wFillColor, wLineColor, lumi/lumiW, rebin);
  setHisto(h3, ttFillColor, ttLineColor, lumi/lumiT, rebin);
  setHisto(h4, qcdFillColor, qcdLineColor, lumi/lumiQ, rebin);

  THStack * hs = new THStack("hs","");
  hs->Add(h4);
  hs->Add(h3);
  hs->Add(h2);
  hs->Add(h1);

  hs->Draw("HIST");
  if(hdata != 0) {
    hdata->SetMarkerStyle(20);
    hdata->SetMarkerSize(1.4);
    hdata->SetMarkerColor(kBlack);
    hdata->SetLineWidth(2);
    hdata->SetLineColor(kBlack);
    hdata->Draw("epsame"); 
    hdata->GetXaxis()->SetLabelSize(0);
    hdata->GetYaxis()->SetLabelSize(0);
    hs->SetMaximum( hdata->GetMaximum() + 1);
  }
  hs->SetMinimum(min);
  hs->GetXaxis()->SetTitle("m_{#mu #mu} (GeV/c^{2})");

  std::string yTag = "";
  switch(rebin) {
  case 1: yTag = "events/(GeV/c^{2})"; break;
  case 2: yTag = "events/(2 GeV/c^{2})"; break;
  case 3: yTag = "events/(3 GeV/c^{2})"; break;
  case 4: yTag = "events/(4 GeV/c^{2})"; break;
  case 5: yTag = "events/(5 GeV/c^{2})"; break;
  case 10: yTag = "events/(10 GeV/c^{2})"; break;
  default:
    std::cerr << ">>> ERROR: set y tag for rebin = " << rebin << std::endl;
  };

  hs->GetYaxis()->SetTitle(yTag.c_str());
  hs->GetXaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetTitleSize(0.05);
  hs->GetXaxis()->SetTitleOffset(1.0);
  hs->GetYaxis()->SetTitleOffset(1.3);
  hs->GetYaxis()->SetLabelOffset(0);
  hs->GetXaxis()->SetLabelSize(.05);
  hs->GetYaxis()->SetLabelSize(.05);
  leg = new TLegend(0.8,0.8,0.95,0.95);
  if(hdata != 0)
    leg->AddEntry(hdata,"data");
  leg->AddEntry(h1,"Z#rightarrow#mu #mu","f");
  leg->AddEntry(h2,"W#rightarrow#mu #nu","f");
  leg->AddEntry(h4,"QCD","f");
  leg->AddEntry(h3,"t#bar{t}","f"); 
  leg->SetFillColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->SetShadowColor(kWhite);
  leg->Draw();
  c1->SetLogy();

  stat(h1, h2, h3, h4, hdata, rebin);
}

void makePlots(const char * name, int rebin, const char * plot,
	       double min = 0.001,
	       bool doData = false) {
  TH1F *h1 = (TH1F*)z.Get(name);
  TH1F *h2 = (TH1F*)w.Get(name);
  TH1F *h3 = (TH1F*)tt.Get(name);
  TH1F *h4 = (TH1F*)qcd.Get(name);
  //  TH1F *hdata = doData? (TH1F*)data.Get(name) : 0;
  if (doData) { 
 TH1F *hdat = (TH1F*)data.Get(name) ;
TH1F *hdata = new TH1F ("hdata", "hdata", 200, 0, 200);
  for ( size_t i =60 ; i<120 ; i++  ){
 int val = hdat->GetBinContent(i) ;
    hdata->SetBinContent(i , val );
  }
  } else {TH1F * hdata=0;}
  hdata->Rebin(rebin);
  makeStack(h1, h2, h3, h4, hdata, min, rebin);
  c1->SaveAs((std::string(plot)+".eps").c_str());
  c1->SaveAs((std::string(plot)+".gif").c_str());
}

#endif
