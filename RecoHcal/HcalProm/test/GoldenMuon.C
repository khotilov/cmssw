{
#include "TCanvas.h"
#include <TStyle.h>
#include "TFile.h"



      printf("z1: gROOT Reset \n");
        gROOT->Reset();
        gROOT->SetStyle("Plain");
        gStyle->SetTitleFont(42);           // title font
        gStyle->SetTitleFontSize(0.07);     // title font size
        gStyle->SetLabelSize(0.05,"xy");     //
        gStyle->SetStatFontSize(0.05);     //for box size
        gStyle->SetStatX(0.97);     // for box x position
        gStyle->SetStatW(0.25);     // for box width
        gStyle->SetOptFit(1111); // pcev

   TFile *hfile = new TFile("/tmp/efe/hf_almall.root", "READ");     //open file 

   prom->cd();
   Muon->cd();

   TCanvas *c1 = new TCanvas("c1","Check that Muon passes though DT in 2 places");

string title;
//below is your output directory
title = "/tmp/efe/";

   c1->Divide(2,2);
   c1->cd(1);
   gPad -> SetLogy();
   //     TAxis *axis = histo->GetXaxis();
   //     axis->SetLabelSize(0.06);
   hAngleMuonHB2DT ->Draw();

gPad->SetGridx();
gPad->SetGridy();
   c1->cd(2);
   hDeltaImpXYHB2DT->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c1->cd(3);
   hDeltaZImpXYHB2DT->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c1->cd(4);
   hPhiDeltaTowerHB2DT->Draw();
gPad->SetGridx();
gPad->SetGridy();

string title1 = title + "pic1.gif";
char ctitle1[1024];
strcpy(ctitle1,title1.c_str());
   c1->Print(ctitle1);

// wait
//   c1->WaitPrimitive();
/////////////////////////////////////////////////
   TCanvas *c2 = new TCanvas("c2","E of Muon passes though DT in 2 places, NTowerEta<5");
   c2->Divide(2,2);
   c2->cd(1);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTTopMinus.Fit("gaus","E","",0.8.,2.4);
   c2->cd(2);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTTopPlus.Fit("gaus","E","",0.8.,2.4);
   c2->cd(3);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotMinus.Fit("gaus","E","",0.8.,2.4);
   c2->cd(4);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotPlus.Fit("gaus","E","",0.8.,2.4);
string title2 = title + "pic2.gif";
char ctitle2[1024];
strcpy(ctitle2,title2.c_str());
   c2->Print(ctitle2);
// wait
//   c2.WaitPrimitive();
                                                                                                                                      
/////////////////////////////////////////////////
   TCanvas *c21 = new TCanvas("c21","Time of Muon passes though DT in 2 places, NTowerEta<5");
   c21->Divide(2,2);
   c21->cd(1);
   hTimeMuonHB2DTTopMinus->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c21->cd(2);
   hTimeMuonHB2DTTopPlus->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c21->cd(3);
   hTimeMuonHB2DTBotMinus->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c21->cd(4);
   hTimeMuonHB2DTBotPlus->Draw();
gPad->SetGridx();
gPad->SetGridy();
string title21 = title + "pic21.gif";
char ctitle21[1024];
strcpy(ctitle21,title21.c_str());
   c21->Print(ctitle21);
//  c21.WaitPrimitive();
/////////////////////////////////////////////////
   TCanvas *c22 = new TCanvas("c22","Number Eta Towers of Muon passes though DT in 2 places, NTowerEta<5");
   c22->Divide(2,2);
   c22->cd(1);
   hNumTowerMuonHB2DTTopMinus ->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c22->cd(2);
   hNumTowerMuonHB2DTTopPlus->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c22->cd(3);
   hNumTowerMuonHB2DTBotMinus->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c22->cd(4);
   hNumTowerMuonHB2DTBotPlus->Draw();
gPad->SetGridx();
gPad->SetGridy();
string title22 = title + "pic22.gif";
char ctitle22[1024];
strcpy(ctitle22,title22.c_str());
   c22->Print(ctitle22);
// wait
//  c22.WaitPrimitive();
/////////////////////////////////////////////////
   TCanvas *c4 = new TCanvas("c4","Number Eta Towers of Muon passes though DT in 2 places, NTowerEta<5");
   c4->Divide(2,1);
   c4->cd(1);
   hIdPhiTowerHB2DTMinus->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c4->cd(2);
   hIdPhiTowerHB2DTPlus->Draw();
gPad->SetGridx();
gPad->SetGridy();                                                                                                                      
//   c4->cd(3);
//   hIdPhiTowerHB2DTsevMinus->Draw();
//   c4->cd(4);
//   hIdPhiTowerHB2DTsevPlus->Draw();
string title4 = title + "pic4.gif";
char ctitle4[1024];
strcpy(ctitle4,title4.c_str());
   c4->Print(ctitle4);
// wait
//  c4.WaitPrimitive();
                                                                                                                                      
/////////////////////////////////////////////////
   TCanvas *c51 = new TCanvas("c51","Number Eta Towers of Muon passes though DT in 2 places, NTowerEta<5");
   c51->Divide(2,2);
   c51->cd(1);
   hNumHitsHB2DT->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c51->cd(2);
   hLengthMuonDTHB2DT->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c51->cd(3);
   hNumHitsHB2DT2->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c51->cd(4);
   hLengthMuonDTHB2DT2->Draw();
gPad->SetGridx();
gPad->SetGridy();
string title51 = title + "pic51.gif";
char ctitle51[1024];
strcpy(ctitle51,title51.c_str());
   c51->Print(ctitle51);
// wait
//  c51.WaitPrimitive();
/////////////////////////////////////////////////
   TCanvas *c61 = new TCanvas("c61","Number Eta Towers of Muon passes though DT in 2 places, NTowerEta<5");
   c61->Divide(2,2);
   c61->cd(1);
   //hEmuonHB2DTTopMinusTimeMinus->Draw();
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTTopMinusTimeMinus.Fit("gaus","E","",0.8.,2.4);
   c61->cd(2);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTTopPlusTimeMinus.Fit("gaus","E","",0.8.,2.4);
   c61->cd(3);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotMinusTimeMinus.Fit("gaus","E","",0.8.,2.4);
   c61->cd(4);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotPlusTimeMinus.Fit("gaus","E","",0.8.,2.4);
string title61 = title + "pic61.gif";
char ctitle61[1024];
strcpy(ctitle61,title61.c_str());
   c61->Print(ctitle61);
// wait
//  c61.WaitPrimitive();
/////////////////////////////////////////////////
   TCanvas *c62 = new TCanvas("c62","Number Eta Towers of Muon passes though DT in 2 places, NTowerEta<5");
   c62->Divide(2,2);
   c62->cd(1);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTTopMinusTimePlus->Fit("gaus","E","",0.8.,2.4);
   c62->cd(2);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTTopPlusTimePlus->Fit("gaus","E","",0.8.,2.4);
   c62->cd(3);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotMinusTimePlus->Fit("gaus","E","",0.8.,2.4);
   c62->cd(4);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotPlusTimePlus->Fit("gaus","E","",0.8.,2.4);
string title62 = title + "pic62.gif";
char ctitle62[1024];
strcpy(ctitle62,title62.c_str());
   c62->Print(ctitle62);
// wait
//  c62.WaitPrimitive();
/////////////////////////////////////////////////
   TCanvas *c63 = new TCanvas("c63","Number Eta Towers of Muon passes though DT in 2 places, NTowerEta<5");
   c63->Divide(2,2);
   c63->cd(1);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotPlusTimeMinus1->Fit("gaus","E","",0.8.,2.4);
   c63->cd(2);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotPlusTimeMinus2->Fit("gaus","E","",0.8.,2.4);
   c63->cd(3);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotPlusTimeMinus3->Fit("gaus","E","",0.8.,2.4);
   c63->cd(4);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotPlusTimeMinus4->Fit("gaus","E","",0.8.,2.4);
string title63 = title + "pic63.gif";
char ctitle63[1024];
strcpy(ctitle63,title63.c_str());
   c63->Print(ctitle63);
 
// wait
//  c63.WaitPrimitive();
/////////////////////////////////////////////////
   TCanvas *c64 = new TCanvas("c64","Number Eta Towers of Muon passes though DT in 2 places, NTowerEta<5");
   c64->Divide(2,2);
   c64->cd(1);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotPlusTimePlus1->Fit("gaus","E","",0.8.,2.4);
   c64->cd(2);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotPlusTimePlus2->Fit("gaus","E","",0.8.,2.4);
   c64->cd(3);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotPlusTimePlus3->Fit("gaus","E","",0.8.,2.4);
   c64->cd(4);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonHB2DTBotPlusTimePlus4->Fit("gaus","E","",0.8.,2.4);
string title64 = title + "pic64.gif";
char ctitle64[1024];
strcpy(ctitle64,title64.c_str());
   c64->Print(ctitle64);
// wait
//  c64.WaitPrimitive();
/////////////////////////////////////////////////
   TCanvas *c65 = new TCanvas("c65","Number Eta Towers of Muon passes though DT in 2 places, NTowerEta<5");
   c65->Divide(3,2);
   c65->cd(1);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonPhiDetaTower1->Fit("gaus","E","",0.8.,2.4);
   c65->cd(2);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonPhiDetaTower2->Fit("gaus","E","",0.8.,2.4);
   c65->cd(3);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonPhiDetaTower3->Fit("gaus","E","",0.8.,2.4);
   c65->cd(4);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonPhiDetaTower4->Fit("gaus","E","",0.8.,2.4);
   c65->cd(5);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonPhiDetaTower5->Fit("gaus","E","",0.8.,2.4);
   c65->cd(6);
gPad->SetGridx();
gPad->SetGridy();
   hEmuonPhiDetaTower6->Fit("gaus","E","",0.8.,2.4);
string title65 = title + "pic65.gif";
char ctitle65[1024];
strcpy(ctitle65,title65.c_str());
   c65->Print(ctitle65);
// wait
//  c65.WaitPrimitive();
                                                                                                                                      
////////////////////////////////////////////////
   TCanvas *c7 = new TCanvas("c7","Number Eta Towers of Muon passes though DT in 2 places, NTowerEta<5");
   c7->Divide(2,2);
   c7->cd(1);
   hImpXYHB2DT->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c7->cd(2);
   hZImpXYHB2DT->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c7->cd(3);
   hLmuonDTImpXY->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c7->cd(4);
   hImpXYvsZ->Draw();
gPad->SetGridx();
gPad->SetGridy();
string title7 = title + "pic7.gif";
char ctitle7[1024];
strcpy(ctitle7,title7.c_str());
   c7->Print(ctitle7);
// wait
//  c7.WaitPrimitive();
////////////////////////////////////////////////
        gStyle->SetOptStat(0);   //  no statistics _or_
        //gStyle->SetGrid();   //  no statistics _or_
////////////////////////////////////////////////
   TCanvas *c8 = new TCanvas("c8","Time vs Idphi");
   c8->Divide(2,2);
   c8->cd(1);
   hProfTimeAsIdPhiMinus->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c8->cd(2);
   hProfTimeAsIdPhiPlus->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c8->cd(3);
//   hIdPhiPlusVsE->Draw();
   hProfTimeAsIdEtaTop->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c8->cd(4);
//   hIdPhiMinusVsE->Draw();
   hProfTimeAsIdEtaBot->Draw();
gPad->SetGridx();
gPad->SetGridy();
string title8 = title + "pic8.gif";
char ctitle8[1024];
strcpy(ctitle8,title8.c_str());
   c8->Print(ctitle8);
// wait
//  c8.WaitPrimitive();
////////////////////////////////////////////////
   TCanvas *c81 = new TCanvas("c81","IdEtaTower vs E");
   c81->Divide(1,2);
   c81->cd(1);
gPad->SetGridx();
gPad->SetGridy();
   hIdEtaTopVsE -> Draw();
   c81->cd(2);
   hIdEtaBotVsE -> Draw();
gPad->SetGridx();
gPad->SetGridy();
string title81 = title + "pic81.gif";
char ctitle81[1024];
strcpy(ctitle81,title81.c_str());
   c81->Print(ctitle81);
// wait
  //c81.WaitPrimitive();

////////////////////////////////////////////////
// make histo projection for hIdPhiPlusVsE

   Int_t nbiny = 74; // number or IdPhiTowers + 2
   Int_t nbinx = 60; // number of bins for Energy
   TH1F *hEmuon = new TH1F("hEmuon","Emuom for one phi slice", 60, -2,10);
   TH1F *hEmuonGaussVsIdPhiTowersPlus = 
       new TH1F("hEmuonGaussVsIdPhiTowersPlus","Emuom mean Gauss as IdPhiTowers, HB+", 72, 0.5,72.5);
   TH1F *hEmuonGaussVsIdPhiTowersMinus = 
       new TH1F("hEmuonGaussVsIdPhiTowersMinus","Emuom mean Gauss as IdPhiTowers, HB-", 72, 0.5,72.5);
   TH1F *hEmuonGaussSigmaVsIdPhiTowersPlus = 
       new TH1F("hEmuonGaussSigmaVsIdPhiTowersPlus","Emuom sigma Gauss as IdPhiTowers, HB+", 72, 0.5,72.5);
   TH1F *hEmuonGaussSigmaVsIdPhiTowersMinus = 
       new TH1F("hEmuonGaussSigmaVsIdPhiTowersMinus","Emuom sigma Gauss as IdPhiTowers, HB-", 72, 0.5,72.5);
   TH1F *hEmeanGaussTopPlus = new TH1F("hEmeanGaussTopPlus","Emean Gauss Top HB+, IdPhi: 11-27",40,1.5,2.3); 
   TH1F *hEmeanGaussBotPlus = new TH1F("hEmeanGaussBotPlus","Emean Gauss Bottom HB+, IdPhi: 47-62",40,1.5,2.3); 
   TH1F *hEmeanGaussTopMinus = new TH1F("hEmeanGaussTopMinus","Emean Gauss Top HB-, IdPhi: 11-27",40,1.5,2.3); 
   TH1F *hEmeanGaussBotMinus = new TH1F("hEmeanGaussBotMinus","Emean Gauss Bottom HB-, IdPhi: 47-62",40,1.5,2.3); 
   // HB+
   for(Int_t ik=2; ik<=nbiny-1; ik++){ // make cicle by IdPhiTowers
      // creat histo
      hEmuon->Clear(); //Clean historgram contents = 0;
      for(Int_t jk=1; jk<=nbinx; jk++){ // make cicle by Energy bins 
          Int_t bin2 = hIdPhiPlusVsE->GetBin(jk,ik);// get index for 2d histo
          Double_t binContent = hIdPhiPlusVsE->GetBinContent(bin2); // get content    
          Double_t binError = hIdPhiPlusVsE->GetBinError(bin2); // get error
          hEmuon->SetBinContent(jk,binContent);    
          hEmuon->SetBinError(jk,binError);    
 
      } //end cicle by Energy bins   
      hEmuon->Fit("gaus","E","",0.8.,2.6);
      Double_t Emean = gaus->GetParameter(1);// value of the 2nd parameter (mean) 
      Double_t errEmean = gaus->GetParError(1);// value of the 2nd parameter (mean) 
      Double_t Esigma = gaus->GetParameter(2);// value of the 3nd parameter (mean) 
      Double_t errEsigma = gaus->GetParError(2);// value of the 3nd parameter (mean) 
      if(errEmean>fabs(Emean)/5){errEmean=0;Emean=0;Esigma=0;errEsigma=0;} 
      if(errEsigma>fabs(Esigma)/2){errEmean=0;Emean=0;Esigma=0;errEsigma=0;} 
       cout << " idPhi = "<< ik-1<<"par mean = "<< Emean <<"+/-"<<errEmean << std::endl;
      hEmuonGaussVsIdPhiTowersPlus->SetBinContent(ik-1,Emean); // ik-1 - IdPhiTower 
      hEmuonGaussVsIdPhiTowersPlus->SetBinError(ik-1,errEmean); 
      hEmuonGaussSigmaVsIdPhiTowersPlus->SetBinContent(ik-1,Esigma); 
      hEmuonGaussSigmaVsIdPhiTowersPlus->SetBinError(ik-1,errEsigma);
      if((ik-1)>=11&&(ik-1)<=27)hEmeanGaussTopPlus->Fill(Emean); 
      if((ik-1)>=47&&(ik-1)<=62)hEmeanGaussBotPlus->Fill(Emean); 
      //c9.WaitPrimitive();
   } // end  cicle by IdPhiTowers 
   // HB-
   for(Int_t ik=2; ik<=nbiny-1; ik++){ // make cicle by IdPhiTowers
      // creat histo
      hEmuon->Clear(); //Clean historgram contents = 0;
      for(Int_t jk=1; jk<=nbinx; jk++){ // make cicle by Energy bins
          Int_t bin2 = hIdPhiMinusVsE->GetBin(jk,ik);// get index for 2d histo
          Double_t binContent = hIdPhiMinusVsE->GetBinContent(bin2); // get content
          Double_t binError = hIdPhiMinusVsE->GetBinError(bin2); // get error
          hEmuon->SetBinContent(jk,binContent);
          hEmuon->SetBinError(jk,binError);
                                                                                                       
      } //end cicle by Energy bins
      hEmuon->Fit("gaus","E","",0.8.,2.6);
      Double_t Emean = gaus->GetParameter(1);// value of the 2nd parameter (mean)
      Double_t errEmean = gaus->GetParError(1);// value of the 2nd parameter (mean)
      Double_t Esigma = gaus->GetParameter(2);// value of the 3nd parameter (mean) 
      Double_t errEsigma = gaus->GetParError(2);// value of the 3nd parameter (mean) 
      if(errEmean>fabs(Emean)/5){errEmean=0;Emean=0;Esigma=0;errEsigma=0;}
      if(errEsigma>fabs(Esigma)/2){errEmean=0;Emean=0;Esigma=0;errEsigma=0;} 
       cout << " idPhi = "<< ik-1<<"par mean = "<< Emean <<"+/-"<<errEmean << std::endl;
      hEmuonGaussVsIdPhiTowersMinus->SetBinContent(ik-1,Emean);
      hEmuonGaussVsIdPhiTowersMinus->SetBinError(ik-1,errEmean);
      hEmuonGaussSigmaVsIdPhiTowersMinus->SetBinContent(ik-1,Esigma); 
      hEmuonGaussSigmaVsIdPhiTowersMinus->SetBinError(ik-1,errEsigma); 
      if((ik-1)>=11&&(ik-1)<=27&&(ik-1)!=24)hEmeanGaussTopMinus->Fill(Emean); // without problematic point 24 
      if((ik-1)>=47&&(ik-1)<=62)hEmeanGaussBotMinus->Fill(Emean); 
     // c8.WaitPrimitive();
   } // end  cicle by IdPhiTowers
   TH1F *hEmuonGaussVsIdEtaTowersTop = 
       new TH1F("hEmuonGaussVsIdEtaTowersTop","Emuom mean Gauss as IdEtaTowers, Top HB", 29, -14.5, 14.5);
   TH1F *hEmuonGaussVsIdEtaTowersBot = 
       new TH1F("hEmuonGaussVsIdEtaTowersBot","Emuom mean Gauss as IdEtaTowers, Bottom HB", 29, -14.5, 14.5);
//   Int_t nbinyEta = 59;
//   Int_t nbinxEta = 60;
   Int_t nbinyEta = 59;
   Int_t nbinxEta = 30;
   for(Int_t ik=1; ik<=nbinyEta; ik++){ // make cicle by IdEtaTowers
      // creat histo
      hEmuon->Clear(); //Clean historgram contents = 0;
      for(Int_t jk=1; jk<=nbinxEta; jk++){ // make cicle by Energy bins
          Int_t bin2 = hIdEtaTopVsE->GetBin(jk,ik);// get index for 2d histo
          Double_t binContent = hIdEtaTopVsE->GetBinContent(bin2); // get content
          Double_t binError = hIdEtaTopVsE->GetBinError(bin2); // get error
          hEmuon->SetBinContent(jk,binContent);
          hEmuon->SetBinError(jk,binError);
                                                                                                                          
      } //end cicle by Energy bins
      hEmuon->Fit("gaus","E","",0.8.,2.6);
      Double_t Emean = gaus->GetParameter(1);// value of the 2nd parameter (mean)
      Double_t errEmean = gaus->GetParError(1);// value of the 2nd parameter (mean)
      Double_t Esigma = gaus->GetParameter(2);// value of the 3nd parameter (mean)
      Double_t errEsigma = gaus->GetParError(2);// value of the 3nd parameter (mean)
      if(errEmean>fabs(Emean)/5||errEmean>0.5){errEmean=0;Emean=0;Esigma=0;errEsigma=0;}
      if(errEsigma>fabs(Esigma)/2){errEmean=0;Emean=0;Esigma=0;errEsigma=0;}
      cout << " idEta = "<< (float(ik)-29)/2 <<"par mean = "<< Emean <<"+/-"<<errEmean << std::endl;
      hEmuonGaussVsIdEtaTowersTop->SetBinContent(ik,Emean);
      hEmuonGaussVsIdEtaTowersTop->SetBinError(ik,errEmean);
      //c81.WaitPrimitive();
   } // end cicle by IdEtaTowers
   for(Int_t ik=1; ik<=nbinyEta; ik++){ // make cicle by IdEtaTowers
      // creat histo
      hEmuon->Clear(); //Clean historgram contents = 0;
      for(Int_t jk=1; jk<=nbinxEta; jk++){ // make cicle by Energy bins
          Int_t bin2 = hIdEtaBotVsE->GetBin(jk,ik);// get index for 2d histo
          Double_t binContent = hIdEtaBotVsE->GetBinContent(bin2); // get content
          Double_t binError = hIdEtaBotVsE->GetBinError(bin2); // get error
          hEmuon->SetBinContent(jk,binContent);
          hEmuon->SetBinError(jk,binError);
                                                                                                                          
      } //end cicle by Energy bins
      hEmuon->Fit("gaus","E","",0.8.,2.6);
      Double_t Emean = gaus->GetParameter(1);// value of the 2nd parameter (mean)
      Double_t errEmean = gaus->GetParError(1);// value of the 2nd parameter (mean)
      Double_t Esigma = gaus->GetParameter(2);// value of the 3nd parameter (mean)
      Double_t errEsigma = gaus->GetParError(2);// value of the 3nd parameter (mean)
      if(errEmean>fabs(Emean)/5||errEmean>0.5){errEmean=0;Emean=0;Esigma=0;errEsigma=0;}
      if(errEsigma>fabs(Esigma)/2){errEmean=0;Emean=0;Esigma=0;errEsigma=0;}
      cout << " idEta = "<< (float(ik)-29)/2 <<"par mean = "<< Emean <<"+/-"<<errEmean << std::endl;
      hEmuonGaussVsIdEtaTowersBot->SetBinContent(ik,Emean);
      hEmuonGaussVsIdEtaTowersBot->SetBinError(ik,errEmean);
//      c81.WaitPrimitive();
   } // end cicle by IdEtaTowers

cout<<"next histos"<<endl;

   TCanvas *c9 = new TCanvas("c9","Emuon Gauss mean for different IdPhiTower, HB+");
   c9->Divide(1,2);
   c9->cd(1); 
   hEmuonGaussVsIdPhiTowersPlus->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c9->cd(2); 
   hEmuonGaussVsIdPhiTowersMinus->Draw();
gPad->SetGridx();
gPad->SetGridy();
string title9 = title + "EmuonGaussMeanForHB.gif";
char ctitle9[1024];
strcpy(ctitle9,title9.c_str());
   c9->Print(ctitle9);
//      c9.WaitPrimitive();
/////////////////////////////////////
   TCanvas *c10 = new TCanvas("c10","Emuon Gauss sigma for different IdPhiTower, HB+");
   c10->Divide(1,2);
   c10->cd(1); 
   hEmuonGaussSigmaVsIdPhiTowersPlus->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c10->cd(2); 
   hEmuonGaussSigmaVsIdPhiTowersMinus->Draw();
gPad->SetGridx();
gPad->SetGridy();
string title10 = title + "EmuonGaussSigmaForDiffPhiTower.gif";
char ctitle10[1024];
strcpy(ctitle10,title10.c_str());
c10->Print(ctitle10);
   
//      c10.WaitPrimitive();
/////////////////////////////////////
   TCanvas *c11 = new TCanvas("c11","Emuon Gauss sigma for HB");
   c11->Divide(2,2);
   c11->cd(1); 
gPad->SetGridx();
gPad->SetGridy();
   hEmeanGaussTopMinus->Draw();
   c11->cd(2);
 gPad->SetGridx();
gPad->SetGridy();
   hEmeanGaussTopPlus->Draw();
   c11->cd(3);
 gPad->SetGridx();
gPad->SetGridy();
   hEmeanGaussBotMinus->Draw();
   c11->cd(4); 
gPad->SetGridx();
gPad->SetGridy();
   hEmeanGaussBotPlus->Draw();
string title11 = title + "EmuonGaussSigmaForHB.gif";
char ctitle11[1024];
strcpy(ctitle11,title11.c_str());
   c11->Print(ctitle11);
//      c11.WaitPrimitive();
/////////////////////////////////////
   TCanvas *c12 = new TCanvas("c12","Emuon Gauss mean for different IdEtaTower,HB");
   c12->SetGrid();   
   c12->Divide(1,2);
   c12->cd(1);

   hEmuonGaussVsIdEtaTowersTop->Draw();
gPad->SetGridx();
gPad->SetGridy();
   c12->cd(2);
   hEmuonGaussVsIdEtaTowersBot->Draw();
gPad->SetGridx();
gPad->SetGridy();
string title12 = title + "final.gif";
char ctitle12[1024];
strcpy(ctitle12,title12.c_str());
   c12->Print(ctitle12);



//      c12.WaitPrimitive();
/////////////////////////////////////

}

