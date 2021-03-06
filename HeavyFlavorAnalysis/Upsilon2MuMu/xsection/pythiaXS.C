#include "TMath.h"
#include <math.h>
#include "TF1.h"
#include "TH1.h"

void XS(){
  
  TFile *f = new TFile("~/root/xs/upsilon/101201.fl10.mm.COMBINED.xsReader_10ptbins_ClosureTest.root");
  TH1D *S1; S1 = (TH1D*)gFile->Get("RecoGenRes_1S");
  TH1D *S11; S11 = (TH1D*)gFile->Get("RecoGenRes_2S");  
  TH1D *S21; S21 = (TH1D*)gFile->Get("RecoGenRes_3S");  
  
  TFile *f1 = new TFile("~/root/xs/upsilon/Acc_All_0_50.xsReader_3Sbin.default.root");
  TH1D *S2; S2 = (TH1D*)gFile->Get("UG_AllGenRes_1S");  
  TH1D *S3; S3 = (TH1D*)gFile->Get("UG_RecoGenRes_1S");   
  TH1D *S12; S12 = (TH1D*)gFile->Get("UG_AllGenRes_2S");  
  TH1D *S13; S13 = (TH1D*)gFile->Get("UG_RecoGenRes_2S");   
  TH1D *S22; S22 = (TH1D*)gFile->Get("UG_AllGenRes_2S");  
  TH1D *S23; S23 = (TH1D*)gFile->Get("UG_RecoGenRes_2S");
  
  TH1D *S4 = (TH1D*)S1->Clone();
  TH1D *S5 = (TH1D*)S11->Clone();
  TH1D *S6 = (TH1D*)S21->Clone();
  
  S4Pt = new TH1D("Pt_1S", "Pt_1S", S4->GetNbinsY(), S4->GetYaxis()->GetXbins()->GetArray());  
  S5Pt = new TH1D("Pt_2S", "Pt_2S", S5->GetNbinsY(), S5->GetYaxis()->GetXbins()->GetArray());  
  S6Pt = new TH1D("Pt_3S", "Pt_3S", S6->GetNbinsY(), S6->GetYaxis()->GetXbins()->GetArray());  
  S4Y = new TH1D("Y_1S", "Y_1S", S4->GetNbinsX(), S4->GetXaxis()->GetXbins()->GetArray());    
  S5Y = new TH1D("Y_2S", "Y_2S", S5->GetNbinsX(), S5->GetXaxis()->GetXbins()->GetArray());    
  S6Y = new TH1D("Y_3S", "Y_3S", S6->GetNbinsX(), S6->GetXaxis()->GetXbins()->GetArray());  
  
  S4Pt_Nor = new TH1D("Pt_1S_Nor", "Pt_1S_Nor", S4->GetNbinsY(), S4->GetYaxis()->GetXbins()->GetArray());  
  S5Pt_Nor = new TH1D("Pt_2S_Nor", "Pt_2S_Nor", S5->GetNbinsY(), S5->GetYaxis()->GetXbins()->GetArray());  
  S6Pt_Nor = new TH1D("Pt_3S_Nor", "Pt_3S_Nor", S6->GetNbinsY(), S6->GetYaxis()->GetXbins()->GetArray());  
  S4Y_Nor = new TH1D("Y_1S_Nor", "Y_1S_Nor", S4->GetNbinsX(), S4->GetXaxis()->GetXbins()->GetArray());    
  S5Y_Nor = new TH1D("Y_2S_Nor", "Y_2S_Nor", S5->GetNbinsX(), S5->GetXaxis()->GetXbins()->GetArray());    
  S6Y_Nor = new TH1D("Y_3S_Nor", "Y_3S_Nor", S6->GetNbinsX(), S6->GetXaxis()->GetXbins()->GetArray());    
    
  double reco(-1.); double reco2(-1.); double reco3(-1.); 
  double deno(-1.); double deno2(-1.); double deno3(-1.);
  double numa(-1.); double numa2(-1.); double numa3(-1.);
  double xs(-1.); double xs2(-1.); double xs3(-1.);  
  double lumi(162500.); double lumi2(190200.); double lumi3(566500.);
  double n1(2.); double n2(2.7); double n3(1.4);  
  double totxs(0.); double totxs2(0.); double totxs3(0.); 
  
  for ( int iy = 1; iy <= S1->GetNbinsX(); ++iy ){
    for ( int ipt = 1; ipt <= S1->GetNbinsY(); ++ipt ){
      reco = S1->GetBinContent(iy,ipt);
      deno = S2->GetBinContent(iy,ipt);
      numa = S3->GetBinContent(iy,ipt);
      xs = (reco/lumi)*(deno/numa);
      S4->SetBinContent(iy,ipt, xs);
      totxs += xs;
      
      reco2 = S11->GetBinContent(iy,ipt);
      deno2 = S12->GetBinContent(iy,ipt);
      numa2 = S13->GetBinContent(iy,ipt);
      xs2 = (reco2/lumi2)*(deno2/numa2);
      S5->SetBinContent(iy,ipt, xs2);      
      totxs2 += xs2;
      
      reco3 = S21->GetBinContent(iy,ipt);
      deno3 = S22->GetBinContent(iy,ipt);
      numa3 = S23->GetBinContent(iy,ipt);
      xs3 = (reco3/lumi3)*(deno3/numa3);
      S6->SetBinContent(iy,ipt, xs3);
      totxs3 += xs3;
      
    }
  }
  
  cout << "totxs = " << totxs << " totxs2 = " << totxs2 << " totxs3 = " << totxs3 << endl;
  
  double s4_pt(0.); double s5_pt(0.); double s6_pt(0.); 
  
  for (int j = 1; j <= S4->GetNbinsY(); ++j){
    for (int i = 1; i <= S4->GetNbinsX(); ++i) {
      s4_pt += S4->GetCellContent(i,j);
      s5_pt += S5->GetCellContent(i,j);
      s6_pt += S6->GetCellContent(i,j);
    }
    
    S4Pt->SetBinContent(j,s4_pt/S4Pt->GetBinWidth(j));
    S5Pt->SetBinContent(j,s5_pt/S5Pt->GetBinWidth(j));
    S6Pt->SetBinContent(j,s6_pt/S6Pt->GetBinWidth(j));
    
    S4Pt_Nor->SetBinContent(j,(s4_pt/S4Pt->GetBinWidth(j))/n1);
    S5Pt_Nor->SetBinContent(j,(s5_pt/S5Pt->GetBinWidth(j))/n2);
    S6Pt_Nor->SetBinContent(j,(s6_pt/S6Pt->GetBinWidth(j))/n3);    
    
    s4_pt=0.; s5_pt=0.; s6_pt=0.;
    
  }
  
  
  double s4_y(0.); double s5_y(0.); double s6_y(0.); 
  
  for (int i = 1; i <= S4->GetNbinsX(); ++i){
    for (int j = 1; j <= S4->GetNbinsY(); ++j) {
      s4_y += S4->GetCellContent(i,j);
      s5_y += S5->GetCellContent(i,j);
      s6_y += S6->GetCellContent(i,j);
    }
    
    S4Y->SetBinContent(i,s4_y/S4Y->GetBinWidth(i));
    S5Y->SetBinContent(i,s5_y/S5Y->GetBinWidth(i));
    S6Y->SetBinContent(i,s6_y/S6Y->GetBinWidth(i));
    
    S4Y_Nor->SetBinContent(i,(s4_y/S4Y->GetBinWidth(i))/n1);
    S5Y_Nor->SetBinContent(i,(s5_y/S5Y->GetBinWidth(i))/n2);
    S6Y_Nor->SetBinContent(i,(s6_y/S6Y->GetBinWidth(i))/n3);   
    
    s4_y=0.; s5_y=0.; s6_y=0.;
    
  }
  
  TFile *f = new TFile("Pythia_XS.root", "RECREATE");
  S4Pt->Write(); S4Pt_Nor->Write();
  S5Pt->Write(); S5Pt_Nor->Write();
  S6Pt->Write(); S6Pt_Nor->Write();
  S4Y->Write(); S4Y_Nor->Write();
  S5Y->Write(); S5Y_Nor->Write();
  S6Y->Write(); S6Y_Nor->Write();
}

void Overlay(){
  
  gStyle->SetOptStat(00000000000);
  TFile *f = new TFile("Pythia_XS.root");
  TH1D *h1;
  h1 = (TH1D*)gFile->Get("Pt_1S");
  TH1D *h2;
  h2 = (TH1D*)gFile->Get("Pt_2S");
  TH1D *h3;
  h3 = (TH1D*)gFile->Get("Pt_3S"); 
  
  TFile *f = new TFile("Final1S.root");
  TGraphAsymmErrors *S1;
  S1 = (TGraphAsymmErrors*)gFile->Get("Ups1S");
  TFile *f = new TFile("Final2S.root");
  TGraphAsymmErrors *S2;
  S2 = (TGraphAsymmErrors*)gFile->Get("Ups2S");
  TFile *f = new TFile("Final3S.root");
  TGraphAsymmErrors *S3;
  S3 = (TGraphAsymmErrors*)gFile->Get("Ups3S");  
  
  TCanvas *c1 = new TCanvas("c1", "c1", 800,600);
  c1->SetLogy();
  S1->SetLineColor(1); S1->SetMarkerColor(1); S2->SetLineColor(2); S2->SetMarkerColor(2); S3->SetLineColor(4); S3->SetMarkerColor(4);
  S1->SetMarkerStyle(20); S2->SetMarkerStyle(21); S3->SetMarkerStyle(22);
  S1->GetXaxis()->SetTitle("p_{T}^{#Upsilon}(GeV/c)");
  S1->GetYaxis()->SetTitle("d#sigma/dp_{T}#times Br(#mu#mu)");
  S1->SetTitle("");
  S1->SetMinimum(0.0001);
  S1->SetMaximum(2.5);
  S1->Draw("AP");
  h1->Draw("same");
  S2->Draw("P");
  h2->SetLineColor(2);
  h2->Draw("same");
  S3->Draw("P");
  h3->SetLineColor(4);
  h3->Draw("same");  
  legg = new TLegend(0.5,0.5,0.7,0.9);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05); legg->SetTextFont(132); 
  legg->SetHeader("");
  legge = legg->AddEntry(S1, "Data #Upsilon(1S)" ,"p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(h1, "Pythia #Upsilon(1S)" ,"l"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(S2, "Data #Upsilon(2S)","p"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(h2, "Pythia #Upsilon(2S)","l"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(S3, "Data #Upsilon(3S)","p"); legge->SetTextColor(kBlue);
  legge = legg->AddEntry(h3, "Pythia #Upsilon(3S)","l"); legge->SetTextColor(kBlue);
  legg->Draw();
  c1->SaveAs("Overlay_pythia.pdf");
}

void Overlay_Nor(){
  
  gStyle->SetOptStat(00000000000);
  TFile *f = new TFile("Pythia_XS.root");
  TH1D *h1;
  h1 = (TH1D*)gFile->Get("Pt_1S_Nor");
  TH1D *h2;
  h2 = (TH1D*)gFile->Get("Pt_2S_Nor");
  TH1D *h3;
  h3 = (TH1D*)gFile->Get("Pt_3S_Nor"); 
  
  TFile *f = new TFile("Final1S.root");
  TGraphAsymmErrors *S1;
  S1 = (TGraphAsymmErrors*)gFile->Get("Ups1S");
  TFile *f = new TFile("Final2S.root");
  TGraphAsymmErrors *S2;
  S2 = (TGraphAsymmErrors*)gFile->Get("Ups2S");
  TFile *f = new TFile("Final3S.root");
  TGraphAsymmErrors *S3;
  S3 = (TGraphAsymmErrors*)gFile->Get("Ups3S");  
  
  TCanvas *c1 = new TCanvas("c1", "c1", 800,600);
  c1->SetLogy();
  S1->SetLineColor(1); S1->SetMarkerColor(1); S2->SetLineColor(2); S2->SetMarkerColor(2); S3->SetLineColor(4); S3->SetMarkerColor(4);
  S1->SetMarkerStyle(20); S2->SetMarkerStyle(21); S3->SetMarkerStyle(22);
  S1->GetXaxis()->SetTitle("p_{T}^{#Upsilon}(GeV/c)");
  S1->GetYaxis()->SetTitle("d#sigma/dp_{T}#times Br(#mu#mu)");
  S1->SetTitle("");
  S1->SetMinimum(0.0001);
  S1->SetMaximum(2.5);
  S1->Draw("AP");
  h1->Draw("same");
  S2->Draw("P");
  h2->SetLineColor(2);
  h2->Draw("same");
  S3->Draw("P");
  h3->SetLineColor(4);
  h3->Draw("same");  
  legg = new TLegend(0.45,0.5,0.7,0.9);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05); legg->SetTextFont(132); 
  legg->SetHeader("");
  legge = legg->AddEntry(S1, "Data #Upsilon(1S)" ,"p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(h1, "Pythia(Normalized) #Upsilon(1S)" ,"l"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(S2, "Data #Upsilon(2S)","p"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(h2, "Pythia(Normalized #Upsilon(2S)","l"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(S3, "Data #Upsilon(3S)","p"); legge->SetTextColor(kBlue);
  legge = legg->AddEntry(h3, "Pythia(Normalized #Upsilon(3S)","l"); legge->SetTextColor(kBlue);
  legg->Draw();
  c1->SaveAs("Overlay_pythia_Nor.pdf");
}

void Overlay_rap(){
  
  gStyle->SetOptStat(00000000000);
  TFile *f = new TFile("Pythia_XS.root");
  TH1D *h1;
  h1 = (TH1D*)gFile->Get("Y_1S");
  TH1D *h2;
  h2 = (TH1D*)gFile->Get("Y_2S");
  TH1D *h3;
  h3 = (TH1D*)gFile->Get("Y_3S"); 
  
  TFile *f = new TFile("Final1S_rap.root");
  TGraphAsymmErrors *S1;
  S1 = (TGraphAsymmErrors*)gFile->Get("Ups1S");
  TFile *f = new TFile("Final2S_rap.root");
  TGraphAsymmErrors *S2;
  S2 = (TGraphAsymmErrors*)gFile->Get("Ups2S");
  TFile *f = new TFile("Final3S_rap.root");
  TGraphAsymmErrors *S3;
  S3 = (TGraphAsymmErrors*)gFile->Get("Ups3S");  
  
  TCanvas *c2 = new TCanvas("c2", "c2", 800,600);
  S1->SetLineColor(1); S1->SetMarkerColor(1); S2->SetLineColor(2); S2->SetMarkerColor(2); S3->SetLineColor(4); S3->SetMarkerColor(4);
  S1->SetMarkerStyle(20); S2->SetMarkerStyle(21); S3->SetMarkerStyle(22);
  S1->GetXaxis()->SetTitle("|y^{#Upsilon}|");
  S1->GetYaxis()->SetTitle("d#sigma/dp_{T}#times Br(#mu#mu)");
  S1->SetTitle("");
  S1->SetMinimum(0.0001);
  S1->SetMaximum(9.5);
  S1->Draw("AP");
  h1->Draw("same");
  S2->Draw("P");
  h2->SetLineColor(2);
  h2->Draw("same");
  S3->Draw("P");
  h3->SetLineColor(4);
  h3->Draw("same");  
  legg = new TLegend(0.71,0.62,0.85,0.9);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04); legg->SetTextFont(132); 
  legg->SetHeader("");
  legge = legg->AddEntry(S1, "Data #Upsilon(1S)" ,"p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(h1, "Pythia #Upsilon(1S)" ,"l"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(S2, "Data #Upsilon(2S)","p"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(h2, "Pythia #Upsilon(2S)","l"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(S3, "Data #Upsilon(3S)","p"); legge->SetTextColor(kBlue);
  legge = legg->AddEntry(h3, "Pythia #Upsilon(3S)","l"); legge->SetTextColor(kBlue);
  legg->Draw();
  c2->SaveAs("Overlay_pythia_rap.pdf");
}

void Overlay_rap_Nor(){
  
  gStyle->SetOptStat(00000000000);
  TFile *f = new TFile("Pythia_XS.root");
  TH1D *h1;
  h1 = (TH1D*)gFile->Get("Y_1S_Nor");
  TH1D *h2;
  h2 = (TH1D*)gFile->Get("Y_2S_Nor");
  TH1D *h3;
  h3 = (TH1D*)gFile->Get("Y_3S_Nor"); 
  
  TFile *f = new TFile("Final1S_rap.root");
  TGraphAsymmErrors *S1;
  S1 = (TGraphAsymmErrors*)gFile->Get("Ups1S");
  TFile *f = new TFile("Final2S_rap.root");
  TGraphAsymmErrors *S2;
  S2 = (TGraphAsymmErrors*)gFile->Get("Ups2S");
  TFile *f = new TFile("Final3S_rap.root");
  TGraphAsymmErrors *S3;
  S3 = (TGraphAsymmErrors*)gFile->Get("Ups3S");  
  
  TCanvas *c2 = new TCanvas("c2", "c2", 800,600);
  S1->SetLineColor(1); S1->SetMarkerColor(1); S2->SetLineColor(2); S2->SetMarkerColor(2); S3->SetLineColor(4); S3->SetMarkerColor(4);
  S1->SetMarkerStyle(20); S2->SetMarkerStyle(21); S3->SetMarkerStyle(22);
  S1->GetXaxis()->SetTitle("|y^{#Upsilon}|");
  S1->GetYaxis()->SetTitle("d#sigma/dp_{T}#times Br(#mu#mu)");
  S1->SetTitle("");
  S1->SetMinimum(0.0001);
  S1->SetMaximum(9.5);
  S1->Draw("AP");
  h1->Draw("same");
  S2->Draw("P");
  h2->SetLineColor(2);
  h2->Draw("same");
  S3->Draw("P");
  h3->SetLineColor(4);
  h3->Draw("same");  
  legg = new TLegend(0.41,0.62,0.85,0.9);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.04); legg->SetTextFont(132); 
  legg->SetHeader("");
  legge = legg->AddEntry(S1, "Data #Upsilon(1S)" ,"p"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(h1, "Pythia(Normalized) #Upsilon(1S)" ,"l"); legge->SetTextColor(kBlack);
  legge = legg->AddEntry(S2, "Data #Upsilon(2S)","p"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(h2, "Pythia(Normalized) #Upsilon(2S)","l"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(S3, "Data #Upsilon(3S)","p"); legge->SetTextColor(kBlue);
  legge = legg->AddEntry(h3, "Pythia(Normalized) #Upsilon(3S)","l"); legge->SetTextColor(kBlue);
  legg->Draw();
  c2->SaveAs("Overlay_pythia_rap_Nor.pdf");
}

void Ratio(){
  
  gStyle->SetOptStat(00000000000);
  TFile *f = new TFile("Final1S.root");
  TH1D *S1h; TH1D *S1l;
  S1h = (TH1D*)gFile->Get("Ups1S_ratioh");
  S1l = (TH1D*)gFile->Get("Ups1S_ratiol");
  TFile *f = new TFile("Final2S.root");
  TH1D *S2h; TH1D *S2l;
  S2h = (TH1D*)gFile->Get("Ups2S_ratioh");
  S2l = (TH1D*)gFile->Get("Ups2S_ratiol");  
  TFile *f = new TFile("Final3S.root");
  TH1D *S3h; TH1D *S3l;
  S3h = (TH1D*)gFile->Get("Ups3S_ratioh");
  S3l = (TH1D*)gFile->Get("Ups3S_ratiol");   
  
  TFile *f = new TFile("Pythia_XS.root");
  TH1D *h1;
  h1 = (TH1D*)gFile->Get("Pt_1S");
  TH1D *h2;
  h2 = (TH1D*)gFile->Get("Pt_2S");
  TH1D *h3;
  h3 = (TH1D*)gFile->Get("Pt_3S"); 
  TH1D *h21 = (TH1D*)h1->Clone();
  TH1D *h31 = (TH1D*)h1->Clone();
  
  for (int i = 1; i <= h1->GetNbinsX(); ++i) {
    
    double bincontent21 = h2->GetBinContent(i)/h1->GetBinContent(i); 
    h21->SetBinContent(i,bincontent21);
    double bincontent31 = h3->GetBinContent(i)/h1->GetBinContent(i); 
    h31->SetBinContent(i,bincontent31);   
    
  }
  
  double s21[10]; double s31[10]; double e21_h[10]; double e31_h[10]; double e21_l[10]; double e31_l[10];
  double term1_h(0.), term2_h(0.),term3_h(0.), term21_h(0.), term31_h(0.);
  double term1_l(0.), term2_l(0.),term3_l(0.), term21_l(0.), term31_l(0.);
  for (int i = 1; i <= S1h->GetNbinsX(); ++i) {
    
   s21[i-1] = S2h->GetBinContent(i)/S1h->GetBinContent(i);
   s31[i-1] = S3h->GetBinContent(i)/S1h->GetBinContent(i);
      
   term1_h = (S1h->GetBinError(i)/S1h->GetBinContent(i))*(S1h->GetBinError(i)/S1h->GetBinContent(i));
   term2_h = (S2h->GetBinError(i)/S2h->GetBinContent(i))*(S2h->GetBinError(i)/S2h->GetBinContent(i));
   term3_h = (S3h->GetBinError(i)/S3h->GetBinContent(i))*(S3h->GetBinError(i)/S3h->GetBinContent(i));
   term21_h = term1_h+term2_h;
   term31_h = term1_h+term3_h;
   e21_h[i-1] = s21[i-1]*TMath::Sqrt(term21_h);
   e31_h[i-1] = s31[i-1]*TMath::Sqrt(term31_h);
      
   term1_l = (S1l->GetBinError(i)/S1l->GetBinContent(i))*(S1l->GetBinError(i)/S1l->GetBinContent(i));
   term2_l = (S2l->GetBinError(i)/S2l->GetBinContent(i))*(S2l->GetBinError(i)/S2l->GetBinContent(i));
   term3_l = (S3l->GetBinError(i)/S3l->GetBinContent(i))*(S3l->GetBinError(i)/S3l->GetBinContent(i));
   term21_l = term1_l+term2_l;
   term31_l = term1_l+term3_l;
   e21_l[i-1] = s21[i-1]*TMath::Sqrt(term21_l);
   e31_l[i-1] = s31[i-1]*TMath::Sqrt(term31_l);   
   
   cout << e21_h[i-1] << endl; cout << e21_l[i-1] << endl;
   cout << e31_h[i-1] << endl; cout << e31_l[i-1] << endl;
   cout << endl;
  }  
  
  double xbin[10] = {1.27, 2.95, 5.72, 8.92, 11.30, 14.32, 17.73, 22.06, 27.16, 36.25};
  double xl[10] = {1.27, 0.95, 1.72, 0.92, 1.30, 1.32, 1.73, 2.06, 2.16, 6.25};
  double xh[10] = {0.73, 1.05, 2.28, 1.08, 1.70, 1.68, 2.27, 2.94, 2.84, 13.75};

  S21 = new TGraphAsymmErrors(10,xbin,s21,xl,xh,e21_l,e21_h);
  S21->SetName("Ups1S_2S");
  S21->SetMarkerColor(2);
  S21->SetMarkerStyle(21);
    
  S31 = new TGraphAsymmErrors(10,xbin,s31,xl,xh,e31_l,e31_h);
  S31->SetName("Ups1S_3S");
  S31->SetMarkerColor(4);
  S31->SetMarkerStyle(22);
    
  S21->GetXaxis()->SetTitle("p_{T}^{#Upsilon}(GeV/c)");
  S21->GetYaxis()->SetTitle("Ratio");
  S21->SetTitle("");
  S21->SetMinimum(0.);
  S21->SetMaximum(1.);
  h21->SetLineColor(2);
  h31->SetLineColor(4);
  S21->Draw("AP");
  S31->Draw("P");
  h21->Draw("same");
  h31->Draw("same");  
  legg = new TLegend(0.1,0.6,0.5,0.9);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05); legg->SetTextFont(132); 
  legg->SetHeader("");
  legge = legg->AddEntry(S21, "Data #Upsilon(2S)/#Upsilon(1S)" ,"p"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(S31, "Data #Upsilon(3S)/#Upsilon(1S)","p"); legge->SetTextColor(kBlue);
  legge = legg->AddEntry(h21, "Pythia #Upsilon(2S)/#Upsilon(1S)" ,"l"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(h31, "Pythia #Upsilon(3S)/#Upsilon(1S)","l"); legge->SetTextColor(kBlue);  
  legg->Draw();
  c1->SaveAs("Ratio_pythia.pdf");
  
}


void Ratio_rap(){
  
  gStyle->SetOptStat(00000000000);
  TFile *f = new TFile("Final1S_rap.root");
  TH1D *S1h; TH1D *S1l;
  S1h = (TH1D*)gFile->Get("Ups1S_ratioh");
  S1l = (TH1D*)gFile->Get("Ups1S_ratiol");
  TFile *f = new TFile("Final2S_rap.root");
  TH1D *S2h; TH1D *S2l;
  S2h = (TH1D*)gFile->Get("Ups2S_ratioh");
  S2l = (TH1D*)gFile->Get("Ups2S_ratiol");  
  TFile *f = new TFile("Final3S_rap.root");
  TH1D *S3h; TH1D *S3l;
  S3h = (TH1D*)gFile->Get("Ups3S_ratioh");
  S3l = (TH1D*)gFile->Get("Ups3S_ratiol");   
    
  TFile *f = new TFile("Pythia_XS.root");
  TH1D *h1;
  h1 = (TH1D*)gFile->Get("Y_1S");
  TH1D *h2;
  h2 = (TH1D*)gFile->Get("Y_2S");
  TH1D *h3;
  h3 = (TH1D*)gFile->Get("Y_3S"); 
  TH1D *h21 = (TH1D*)h1->Clone();
  TH1D *h31 = (TH1D*)h1->Clone();
  
  for (int i = 1; i <= h1->GetNbinsX(); ++i) {
    
    double bincontent21 = h2->GetBinContent(i)/h1->GetBinContent(i); 
    h21->SetBinContent(i,bincontent21);
    double bincontent31 = h3->GetBinContent(i)/h1->GetBinContent(i); 
    h31->SetBinContent(i,bincontent31);   
    
  }
  
  double s21[6]; double s31[6]; double e21_h[6]; double e31_h[6]; double e21_l[6]; double e31_l[6];
  double term1_h(0.), term2_h(0.),term3_h(0.), term21_h(0.), term31_h(0.);
  double term1_l(0.), term2_l(0.),term3_l(0.), term21_l(0.), term31_l(0.);
  for (int i = 1; i <= S1h->GetNbinsX(); ++i) {
    
   s21[i-1] = S2h->GetBinContent(i)/S1h->GetBinContent(i);
   s31[i-1] = S3h->GetBinContent(i)/S1h->GetBinContent(i);
      
   term1_h = (S1h->GetBinError(i)/S1h->GetBinContent(i))*(S1h->GetBinError(i)/S1h->GetBinContent(i));
   term2_h = (S2h->GetBinError(i)/S2h->GetBinContent(i))*(S2h->GetBinError(i)/S2h->GetBinContent(i));
   term3_h = (S3h->GetBinError(i)/S3h->GetBinContent(i))*(S3h->GetBinError(i)/S3h->GetBinContent(i));
   term21_h = term1_h+term2_h;
   term31_h = term1_h+term3_h;
   e21_h[i-1] = s21[i-1]*TMath::Sqrt(term21_h);
   e31_h[i-1] = s31[i-1]*TMath::Sqrt(term31_h);
      
   term1_l = (S1l->GetBinError(i)/S1l->GetBinContent(i))*(S1l->GetBinError(i)/S1l->GetBinContent(i));
   term2_l = (S2l->GetBinError(i)/S2l->GetBinContent(i))*(S2l->GetBinError(i)/S2l->GetBinContent(i));
   term3_l = (S3l->GetBinError(i)/S3l->GetBinContent(i))*(S3l->GetBinError(i)/S3l->GetBinContent(i));
   term21_l = term1_l+term2_l;
   term31_l = term1_l+term3_l;
   e21_l[i-1] = s21[i-1]*TMath::Sqrt(term21_l);
   e31_l[i-1] = s31[i-1]*TMath::Sqrt(term31_l);   
   
   cout << e21_h[i-1] << endl; cout << e21_l[i-1] << endl;
   cout << e31_h[i-1] << endl; cout << e31_l[i-1] << endl;
   cout << endl;
  }  
  
  double xbin[6] = {0.208, 0.611, 1.004, 1.395, 1.786, 2.126};
  double xh[6] = {0.192, 0.189, 0.196, 0.205, 0.214, 0.274};
  double xl[6] = {0.208, 0.211, 0.204, 0.195, 0.186, 0.126};

  S21 = new TGraphAsymmErrors(6,xbin,s21,xl,xh,e21_l,e21_h);
  S21->SetName("Ups1S_2S");
  S21->SetMarkerColor(2);
  S21->SetMarkerStyle(21);
    
  S31 = new TGraphAsymmErrors(6,xbin,s31,xl,xh,e31_l,e31_h);
  S31->SetName("Ups1S_3S");
  S31->SetMarkerColor(4);
  S31->SetMarkerStyle(22);
    
  S21->GetXaxis()->SetTitle("|y^{#Upsilon}|");
  S21->GetYaxis()->SetTitle("Ratio");
  S21->SetTitle("");
  S21->SetMinimum(0.);
  S21->SetMaximum(0.6);
  h21->SetLineColor(2);
  h31->SetLineColor(4);
  S21->Draw("AP");
  S31->Draw("P");
  h21->Draw("same");
  h31->Draw("same"); 
  legg = new TLegend(0.1,0.6,0.5,0.9);
  legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05); legg->SetTextFont(132); 
  legg->SetHeader("");
  legge = legg->AddEntry(S21, "#Upsilon(2S)/#Upsilon(1S)" ,"p"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(S31, "#Upsilon(3S)/#Upsilon(1S)","p"); legge->SetTextColor(kBlue);
  legge = legg->AddEntry(h21, "Pythia #Upsilon(2S)/#Upsilon(1S)" ,"l"); legge->SetTextColor(kRed);
  legge = legg->AddEntry(h31, "Pythia #Upsilon(3S)/#Upsilon(1S)","l"); legge->SetTextColor(kBlue);
  legg->Draw();
  c1->SaveAs("Ratio_pythia_rap.pdf");

}
