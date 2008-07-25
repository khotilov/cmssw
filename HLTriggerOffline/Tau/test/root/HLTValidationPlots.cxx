run(string fname1,string fname2,string v1,string v2,string path)
{
  L1Val(fname1,fname2,v1,v2,"DQMData/HLT/HLTTAU/L1");
  L2Val(fname1,fname2,v1,v2,"DQMData/HLT/HLTTAU/"+path+"/L2");
  L25Val(fname1,fname2,v1,v2,"DQMData/HLT/HLTTAU/"+path+"/L25");
  if(path=="SingleTau"||path=="SingleTauMET")
    L3Val(fname1,fname2,v1,v2,"DQMData/HLT/HLTTAU/"+path+"/L3");

  if(path=="ElectronTau")
    ElectronVal(fname1,fname2,v1,v2,"DQMData/HLT/HLTTAU/"+path+"/Electron");

}

draw(string fname1,string fname2,string v1,string v2,string path)
{
  L1Draw(fname1,fname2,v1,v2,"DQMData/HLT/HLTTAU/L1");
  L2Draw(fname1,fname2,v1,v2,"DQMData/HLT/HLTTAU/"+path+"/L2");
  L25Draw(fname1,fname2,v1,v2,"DQMData/HLT/HLTTAU/"+path+"/L25");

  if(path=="SingleTau"||path=="SingleTauMET")
    L3Draw(fname1,fname2,v1,v2,"DQMData/HLT/HLTTAU/"+path+"/L3");

  if(path=="ElectronTau")
    ElectronDraw(fname1,fname2,v1,v2,"DQMData/HLT/HLTTAU/"+path+"/Electron");

}




L1Draw(string fname1,string fname2,string v1,string v2,string folder)
{
  TCanvas *l1 = new TCanvas;
  l1->Divide(4,5);
  l1->cd(1);
  DrawHistos(fname1,fname2,v1,v2,folder,"L1TauEt","L1 E_{t}","n #tau cands","L1 Tau inclusive E_{t}",0,false);
  l1->cd(2);
  DrawHistos(fname1,fname2,v1,v2,folder,"L1TauEta","L1 #eta","n #tau cands","L1 Tau inclusive #eta",0,false);
  l1->cd(3);
  DrawHistos(fname1,fname2,v1, v2,folder,"L1TauPhi","L1 #phi","n #tau cands","L1 Tau inclusive #phi",0,false);
  l1->cd(4);
  DrawHistos(fname1,fname2, v1, v2,folder,"L1Tau1Et","L1 E_{t}","n #tau cands","LeadingL1 Tau  E_{t}",0,false);
  l1->cd(5);
  DrawHistos(fname1,fname2, v1, v2,folder,"L1Tau1Eta","L1 #eta","n #tau cands","Leading L1 Tau  #eta",0,false);
  l1->cd(6);
  DrawHistos(fname1,fname2, v1, v2,folder,"L1Tau1Phi","L1 #phi","n #tau cands","Leading L1Tau  #phi",0,false);
  l1->cd(7);
  DrawHistos(fname1,fname2, v1, v2,folder,"L1Tau2Et","L1 E_{t}","n #tau cands","2nd L1 Tau  E_{t}",0,false);
  l1->cd(8);
  DrawHistos(fname1,fname2, v1, v2,folder,"L1Tau2Eta","L1 #eta","n #tau cands","2nd  L1 Tau  #eta",0,false);
  l1->cd(9);
  DrawHistos(fname1,fname2, v1, v2,folder,"L1Tau2Phi","L1 #phi","n #tau cands","2nd  L1Tau  #phi",0,false);
  l1->cd(10);
  DrawHistos(fname1,fname2, v1, v2,folder,"L1minusRefTauEt","#Delta E_{t}","n #tau cands","L1 Tau E_{t} - Ref E_{t}",0,false);
  l1->cd(11);
  DrawHistos(fname1,fname2, v1, v2,folder,"L1minusMCoverRefTauEt","#Delta E_{t}/E_{t}","n #tau cands","(L1 Tau E_{t} - Ref E_{t})/RefE_{t}",0,false);
  l1->cd(12);
  //Single Object Efficiency

  DrawEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauEt","RefTauHadEt","Ref E_{t}","Efficiency","Regional Tau Efficiency vs E_{t}");
  l1->cd(13);
  DrawEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauEta","RefTauHadEta","Ref #eta","Efficiency","Regional Tau Efficiency vs #eta");
  l1->cd(14);
  DrawEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauPhi","RefTauHadPhi","Ref #phi","Efficiency","Regional Tau Efficiency vs #phi");
  l1->cd(15);
  DrawEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauMuonEt","RefTauMuonEt","Ref E_{t}","Efficiency","Regional Muon Efficiency vs E_{t}");
  l1->cd(16);
  DrawEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauMuonEta","RefTauMuonEta","Ref #eta","Efficiency","Regional Muon Efficiency vs #eta");
  l1->cd(17);
  DrawEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauMuonPhi","RefTauMuonPhi","Ref #phi","Efficiency","Regional Muon Efficiency vs #phi");
  l1->cd(18);
  DrawEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauElecEt","RefTauElecEt","Ref E_{t}","Efficiency","Regional Eg Efficiency vs E_{t}");
  l1->cd(19);
  DrawEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauElecEta","RefTauElecEta","Ref #eta","Efficiency","Regional Eg Efficiency vs #eta");
  l1->cd(20);
  DrawEffHistos(fname1,fname2, v1, v2,folder,"EffRefMatchedTauElecPhi","RefTauElecPhi","Ref #phi","Efficiency","Regional Eg Efficiency vs #phi");

  //Path Threshold Efficiency



}


L1Val(string fname1,string fname2,string v1,string v2,string folder)
{
  TCanvas *l1 = new TCanvas;

  GetHistos(fname1,fname2,v1,v2,folder,"L1TauEt","L1 E_{t}","n #tau cands","L1 Tau inclusive E_{t}",0,false);
  GetHistos(fname1,fname2,v1,v2,folder,"L1TauEta","L1 #eta","n #tau cands","L1 Tau inclusive #eta",0,false);
  GetHistos(fname1,fname2,v1, v2,folder,"L1TauPhi","L1 #phi","n #tau cands","L1 Tau inclusive #phi",0,false);
  GetHistos(fname1,fname2, v1, v2,folder,"L1Tau1Et","L1 E_{t}","n #tau cands","LeadingL1 Tau  E_{t}",0,false);
  GetHistos(fname1,fname2, v1, v2,folder,"L1Tau1Eta","L1 #eta","n #tau cands","Leading L1 Tau  #eta",0,false);
  GetHistos(fname1,fname2, v1, v2,folder,"L1Tau1Phi","L1 #phi","n #tau cands","Leading L1Tau  #phi",0,false);
  GetHistos(fname1,fname2, v1, v2,folder,"L1Tau2Et","L1 E_{t}","n #tau cands","2nd L1 Tau  E_{t}",0,false);
  GetHistos(fname1,fname2, v1, v2,folder,"L1Tau2Eta","L1 #eta","n #tau cands","2nd  L1 Tau  #eta",0,false);
  GetHistos(fname1,fname2, v1, v2,folder,"L1Tau2Phi","L1 #phi","n #tau cands","2nd  L1Tau  #phi",0,false);
  GetHistos(fname1,fname2, v1, v2,folder,"L1minusRefTauEt","#Delta E_{t}","n #tau cands","L1 Tau E_{t} - Ref E_{t}",0,false);
  GetHistos(fname1,fname2, v1, v2,folder,"L1minusMCoverRefTauEt","#Delta E_{t}/E_{t}","n #tau cands","(L1 Tau E_{t} - Ref E_{t})/RefE_{t}",0,false);
  GetEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauEt","RefTauHadEt","Ref E_{t}","Efficiency","Regional Tau Efficiency vs E_{t}");
  GetEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauEta","RefTauHadEta","Ref #eta","Efficiency","Regional Tau Efficiency vs #eta");
  GetEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauPhi","RefTauHadPhi","Ref #phi","Efficiency","Regional Tau Efficiency vs #phi");
  GetEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauMuonEt","RefTauMuonEt","Ref E_{t}","Efficiency","Regional Muon Efficiency vs E_{t}");
  GetEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauMuonEta","RefTauMuonEta","Ref #eta","Efficiency","Regional Muon Efficiency vs #eta");
  GetEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauMuonPhi","RefTauMuonPhi","Ref #phi","Efficiency","Regional Muon Efficiency vs #phi");
  GetEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauElecEt","RefTauElecEt","Ref E_{t}","Efficiency","Regional Eg Efficiency vs E_{t}");
  GetEffHistos(fname1,fname2, v1, v2,folder,"L1RefMatchedTauElecEta","RefTauElecEta","Ref #eta","Efficiency","Regional Eg Efficiency vs #eta");
  GetEffHistos(fname1,fname2, v1, v2,folder,"EffRefMatchedTauElecPhi","RefTauElecPhi","Ref #phi","Efficiency","Regional Eg Efficiency vs #phi");

  //Path Threshold Efficiency



}


L2Val(string fname1,string fname2,string v1,string v2,string folder)
{

  GetHistos(fname1,fname2,v1,v2,folder,"L2tauCandEt","L2 jet E_{t}","n. #tau cands."," L2 jet E_{t}",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L2tauCandEta","L2 jet #eta","n. #tau cands."," L2 jet #eta",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L2ecalIsolEt","L2 ECAL Isol. E_{t}","n. #tau cands.","ECAL #Sigma E_{t} in annulus",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L2seedTowerEt","L2 Seed Tower E_{t}","n. #tau cands.","Seed Tower E_{t} in annulus",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L2towerIsolEt","L2 Tower isol E_{t}","n. #tau cands.","Tower #Sigma E_{t} in annulus",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L2nClusters","L2 Number of  Clusters","n. #tau cands.","Number of clusters",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L2clusterEtaRMS","L2 Cluster #eta RMS","n. #tau cands.","Cluster #eta RMS",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L2clusterPhiRMS","L2 Cluster #phi RMS","n. #tau cands.","Cluster #phi RMS",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L2clusterDeltaRRMS","L2 Cluster #DeltaR RMS","n. #tau cands.","Cluster #DeltaR RMS",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"MET","MET [GeV]","n. events","Missing E_{t}",0,false );
  GetEffHistos(fname1,fname2,v1,v2,folder,"L2EtEffNum","L2EtEffDenom","Ref. E_{t}[GeV]","Efficiency","L2 Efficiency vs E_{t}");
 


}

L25Val(string fname1,string fname2,string v1,string v2,string folder)
{

  GetHistos(fname1,fname2,v1,v2,folder,"L25jetEt","L25 Jet E_{t} [GeV]","n. #tau cands.","L25 Jet E_{t} [GeV]",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L25jetEta","L25 Jet #eta","n. #tau cands.","L25 Jet #eta",0,false );
  // GetHistos(fname1,fname2,v1,v2,folder,"L25EtEff","L25 Reco Jet E_{t} [GeV]","Efficiency","L25 Efficiency vs E_{t}",0,true );
  //GetHistos(fname1,fname2,v1,v2,folder,"L25EtaEff","L25 Efficiency vs #eta","n. #tau cands.","Reco Jet #eta",0,true );
  GetHistos(fname1,fname2,v1,v2,folder,"L25nTrksInJet","L25 Number of Pixel Tracks","n. #tau cands.","# Pixel Tracks in L25 Jet",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L25nQTrksInJet","L25 Number of Q Pixel Tracks","n. #tau cands.","# Quality Pixel Tracks in L25 Jet",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L25signalLeadTrkPt","L25 Lead Track P_{t}[GeV]","n. #tau cands.","Lead Track P_{t} [GeV]",0,false );
  GetEffHistos(fname1,fname2,v1,v2,folder,"L25IsoJetEt","L25jetEt","Ref. E_{t}[GeV]","Efficiency","L25 Efficiency vs E_{t}");
  GetEffHistos(fname1,fname2,v1,v2,folder,"L25IsoJetEta","L25jetEta","#eta","Efficiency","L25 Efficiency vs #eta");
  

}

L3Val(string fname1,string fname2,string v1,string v2,string folder)
{
  GetHistos(fname1,fname2,v1,v2,folder,"L3jetEt","L3 Jet E_{t} [GeV]","n. #tau cands.","L3 Jet E_{t} [GeV]",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L3jetEta","L3 Jet #eta","n. #tau cands.","L3 Jet #eta",0,false );
  // GetHistos(fname1,fname2,v1,v2,folder,"L3EtEff","L3 Jet E_{t} [GeV]","Efficiency","L3 Efficiency vs E_{t}",0,true );
  // GetHistos(fname1,fname2,v1,v2,folder,"L3EtaEff","L3 Efficiency vs #eta","n. #tau cands.","Reco Jet #eta",0,true );

  GetHistos(fname1,fname2,v1,v2,folder,"L3nTrksInJet","L3 Number of  Silicon Tracks","n. #tau cands.","# Si Tracks in L3 Jet",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L3nQTrksInJet","L3 Number of Quality Si Tracks","n. #tau cands.","# Quality Si. Tracks in L3 Jet",0,false );
  GetHistos(fname1,fname2,v1,v2,folder,"L3signalLeadTrkPt","L3 Lead Track P_{t}[GeV]","n. #tau cands.","Lead Track P_{t} [GeV]",0,false );
  GetEffHistos(fname1,fname2,v1,v2,folder,"L3IsoJetEt","L3jetEt","Ref. E_{t}[GeV]","Efficiency","L25 Efficiency vs E_{t}");
  GetEffHistos(fname1,fname2,v1,v2,folder,"L3IsoJetEta","L3jetEta","#eta","Efficiency","L25 Efficiency vs #eta");
 
}




ElectronVal(string fname1,string fname2,string v1,string v2,string f)
{
  string folder=f+"/Electron";
  gStyle->SetOptTitle(1);
  TFile* file1 =   new TFile(fname1.c_str());
  string Histo_path  = folder+"/total eff";
  TH1F*  hf1 = (TH1F*)file1->Get(Histo_path.c_str());

  std::vector<string> vars;
  vars.push_back("et");
  vars.push_back("eta");


  for(int i=0;i<vars.size();++i)
    for(int filter=1; filter < hf1->GetNbinsX() -2 ; filter++)
      {
	string name1(hf1->GetXaxis()->GetBinLabel(filter));
	string name2(hf1->GetXaxis()->GetBinLabel(filter+1));
	hname1=name1+vars[i];
	hname2=name2+vars[i];
	GetEffHistos(fname1,fname2,v1,v2,folder,hname2,hname1,vars[i].c_str(),"Efficiency",name1.c_str());
      }


}


ElectronDraw(string fname1,string fname2,string v1,string v2,string f)
{
  string folder=f+"/Electron";
  gStyle->SetOptTitle(1);

  TFile* file1 =   new TFile(fname1.c_str());
  string Histo_path  = folder+"/total eff";
  TH1F*  hf1 = (TH1F*)file1->Get(Histo_path.c_str());

  std::vector<string> vars;
  vars.push_back("et");
  vars.push_back("eta");
  TCanvas *el =new TCanvas("el","el");
  el->Divide(3,3);
  int pad=1;
  for(int i=0;i<vars.size();++i)
    for(int filter=1; filter < hf1->GetNbinsX() -2 ; filter++)
      {
	el->cd(pad);
	string name1(hf1->GetXaxis()->GetBinLabel(filter));
	string name2(hf1->GetXaxis()->GetBinLabel(filter+1));
	string hname1=name1+vars[i];
	string hname2=name2+vars[i];
	DrawEffHistos(fname1,fname2,v1,v2,folder,hname2,hname1,vars[i].c_str(),"Efficiency",name1.c_str());
	pad+=1;
      }


}






L2Draw(string fname1,string fname2,string v1,string v2,string folder)
{
  TCanvas *l2 = new TCanvas;
  l2->Divide(3,4);
  l2->cd(1);
  DrawHistos(fname1,fname2,v1,v2,folder,"L2tauCandEt","L2 jet E_{t}","n. #tau cands."," L2 jet E_{t}",0,false );
 l2->cd(2); 
 DrawHistos(fname1,fname2,v1,v2,folder,"L2tauCandEta","L2 jet #eta","n. #tau cands."," L2 jet #eta",0,false );
  l2->cd(3); 
 DrawHistos(fname1,fname2,v1,v2,folder,"L2ecalIsolEt","L2 ECAL Isol. E_{t}","n. #tau cands.","ECAL #Sigma E_{t} in annulus",0,false );
 l2->cd(4); 
 DrawHistos(fname1,fname2,v1,v2,folder,"L2seedTowerEt","L2 Seed Tower E_{t}","n. #tau cands.","Seed Tower E_{t} in annulus",0,false );
 l2->cd(5);  
 DrawHistos(fname1,fname2,v1,v2,folder,"L2towerIsolEt","L2 Tower isol E_{t}","n. #tau cands.","Tower #Sigma E_{t} in annulus",0,false );
 l2->cd(6);
 DrawHistos(fname1,fname2,v1,v2,folder,"L2nClusters","L2 Number of  Clusters","n. #tau cands.","Number of clusters",0,false );
 l2->cd(7);  
 DrawHistos(fname1,fname2,v1,v2,folder,"L2clusterEtaRMS","L2 Cluster #eta RMS","n. #tau cands.","Cluster #eta RMS",0,false );
 l2->cd(8);  
 DrawHistos(fname1,fname2,v1,v2,folder,"L2clusterPhiRMS","L2 Cluster #phi RMS","n. #tau cands.","Cluster #phi RMS",0,false );
 l2->cd(9);  
 DrawHistos(fname1,fname2,v1,v2,folder,"L2clusterDeltaRRMS","L2 Cluster #DeltaR RMS","n. #tau cands.","Cluster #DeltaR RMS",0,false );
 l2->cd(10);  
 DrawHistos(fname1,fname2,v1,v2,folder,"MET","MET [GeV]","n. events","Missing E_{t}",0,false );
 l2->cd(11);
 DrawEffHistos(fname1,fname2,v1,v2,folder,"L2EtEffNum","L2EtEffDenom","Ref. E_{t}[GeV]","Efficiency","L2 Efficiency vs E_{t}");


}



L25Draw(string fname1,string fname2,string v1,string v2,string folder)
{
  TCanvas *l25 = new TCanvas;
  l25->Divide(3,3);
  l25->cd(1);
  DrawHistos(fname1,fname2,v1,v2,folder,"L25jetEt","L25 Jet E_{t} [GeV]","n. #tau cands.","L25 Jet E_{t} [GeV]",0,false );
  l25->cd(2);
  DrawHistos(fname1,fname2,v1,v2,folder,"L25jetEta","L25 Jet #eta","n. #tau cands.","L25 Jet #eta",0,false );
  //  l25->cd(3);
  // DrawHistos(fname1,fname2,v1,v2,folder,"L25EtEff","L25 Reco Jet E_{t} [GeV]","Efficiency","L25 Efficiency vs E_{t}",0,true );
  //l25->cd(4);
  // DrawHistos(fname1,fname2,v1,v2,folder,"L25EtaEff","L25 Efficiency vs #eta","n. #tau cands.","Reco Jet #eta",0,true );
  l25->cd(5);
  DrawHistos(fname1,fname2,v1,v2,folder,"L25nTrksInJet","L25 Number of Pixel Tracks","n. #tau cands.","# Pixel Tracks in L25 Jet",0,false );
  l25->cd(6);
  DrawHistos(fname1,fname2,v1,v2,folder,"L25nQTrksInJet","L25 Number of Q Pixel Tracks","n. #tau cands.","# Quality Pixel Tracks in L25 Jet",0,false );
  l25->cd(7);
  DrawHistos(fname1,fname2,v1,v2,folder,"L25signalLeadTrkPt","L25 Lead Track P_{t}[GeV]","n. #tau cands.","Lead Track P_{t} [GeV]",0,false );
  l25->cd(3);
  DrawEffHistos(fname1,fname2,v1,v2,folder,"L25IsoJetEt","L25jetEt","Ref. E_{t}[GeV]","Efficiency","L25 Efficiency vs E_{t}");
  l25->cd(4);
  DrawEffHistos(fname1,fname2,v1,v2,folder,"L25IsoJetEta","L25jetEta","#eta","Efficiency","L25 Efficiency vs #eta");
  

}

L3Draw(string fname1,string fname2,string v1,string v2,string folder)
{
  TCanvas *l33 = new TCanvas;
  l33->Divide(3,3);
  l33->cd(1);
  DrawHistos(fname1,fname2,v1,v2,folder,"L3jetEt","L3 Jet E_{t} [GeV]","n. #tau cands.","L3 Jet E_{t} [GeV]",0,false );
  l33->cd(2);
  DrawHistos(fname1,fname2,v1,v2,folder,"L3jetEta","L3 Jet #eta","n. #tau cands.","L3 Jet #eta",0,false );
  l33->cd(3);
  DrawHistos(fname1,fname2,v1,v2,folder,"L3EtEff","L3 Jet E_{t} [GeV]","Efficiency","L3 Efficiency vs E_{t}",0,true );
  l33->cd(4);
  DrawHistos(fname1,fname2,v1,v2,folder,"L3EtaEff","L3 Efficiency vs #eta","n. #tau cands.","Reco Jet #eta",0,true );
  l33->cd(5);
  DrawHistos(fname1,fname2,v1,v2,folder,"L3nTrksInJet","L3 Number of  Silicon Tracks","n. #tau cands.","# Si Tracks in L3 Jet",0,false );
  l33->cd(6);
  DrawHistos(fname1,fname2,v1,v2,folder,"L3nQTrksInJet","L3 Number of Quality Si Tracks","n. #tau cands.","# Quality Si. Tracks in L3 Jet",0,false );
  l33->cd(7);
  DrawHistos(fname1,fname2,v1,v2,folder,"L3signalLeadTrkPt","L3 Lead Track P_{t}[GeV]","n. #tau cands.","Lead Track P_{t} [GeV]",0,false );
  l25->cd(3);
  DrawEffHistos(fname1,fname2,v1,v2,folder,"L25IsoJetEt","L25jetEt","Ref. E_{t}[GeV]","Efficiency","L25 Efficiency vs E_{t}");
  l25->cd(4);
  DrawEffHistos(fname1,fname2,v1,v2,folder,"L25IsoJetEta","L25jetEta","#eta","Efficiency","L25 Efficiency vs #eta");
 
}

GetHistos(string f1,string f2,string v1,string v2,string modName,string histo,char* xlabel = "",char *ylabel = "",char *title, double scale = 1,bool eff)
{
  TCanvas *c = new TCanvas;
  c->cd();
  TFile *f = new TFile(f1.c_str());
  TFile *ff = new TFile(f2.c_str());
  
  TH1F *h = ((TH1F*)f->Get((modName+"/"+histo).c_str()))->Clone();
  TH1F *hh = ((TH1F*)ff->Get((modName+"/"+histo).c_str()))->Clone();

  if(!eff)
    h->Sumw2();

  h->GetXaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetNdivisions(509);
  h->GetYaxis()->SetNdivisions(509);
  h->GetYaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(1.25);

   h->SetMarkerColor(kYellow);
  h->SetMarkerStyle(20);
  h->SetFillColor(kYellow);
  h->GetXaxis()->SetTitle(xlabel);
  h->GetYaxis()->SetTitle(ylabel);
  h->SetTitle(title);
  if(!eff)
    {
      if(h->Integral()>0)
	if(scale==0) h->Scale(1/h->Integral());
	else
	  h->Scale(scale);
     
      hh->Sumw2();
    }
  hh->GetXaxis()->SetLabelSize(0.06);
  hh->GetYaxis()->SetLabelSize(0.06);
  hh->GetXaxis()->SetTitleSize(0.06);
  hh->GetYaxis()->SetTitleSize(0.06);
  hh->SetLineColor(kRed);
  hh->SetMarkerColor(kRed);
  hh->SetMarkerStyle(20);
  hh->GetXaxis()->SetTitle(xlabel);
  hh->GetYaxis()->SetTitle(ylabel);
  hh->SetTitle(title);
  if(!eff)
    {
      if(hh->Integral()>0)
	if(scale==0) hh->Scale(1/hh->Integral());
	else
	  hh->Scale(scale);
      
	
    }

  h->Draw("HIST");
  hh->Draw("SAME");
  TLegend *l = new TLegend(0.7,0.5,0.9,0.7);
  l->AddEntry(h,v1.c_str());
  l->AddEntry(hh,v2.c_str());
  l->Draw();

  c->SaveAs((histo+".gif").c_str());
  delete c;
   f->Close();
   ff->Close();

}


GetEffHistos(string f1,string f2,string v1,string v2,string modName,string histo1,string histo2,char* xlabel = "",char *ylabel = "",char *title)
{
  TCanvas *c = new TCanvas;
  c->cd();
  TFile *f = new TFile(f1.c_str());
  TFile *ff = new TFile(f2.c_str());
  
  TH1F *h = ((TH1F*)f->Get((modName+"/"+histo1).c_str()))->Clone();
  TH1F *hd = ((TH1F*)f->Get((modName+"/"+histo2).c_str()))->Clone();

  TH1F *hh = ((TH1F*)ff->Get((modName+"/"+histo1).c_str()))->Clone();
  TH1F *hhd = ((TH1F*)ff->Get((modName+"/"+histo2).c_str()))->Clone();

  h->Sumw2();
  hd->Sumw2();
  hh->Sumw2();
  hhd->Sumw2();


  h->Divide(h,hd,1.,1.,"b");
  hh->Divide(hh,hhd,1.,1.,"b");


  h->GetXaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetNdivisions(509);
  h->GetYaxis()->SetNdivisions(509);
  h->GetYaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(1.25);

   h->SetMarkerColor(kYellow);
  h->SetMarkerStyle(20);
  h->SetFillColor(kYellow);
  h->GetXaxis()->SetTitle(xlabel);
  h->GetYaxis()->SetTitle(ylabel);
  h->SetTitle(title);
  hh->GetXaxis()->SetLabelSize(0.06);
  hh->GetYaxis()->SetLabelSize(0.06);
  hh->GetXaxis()->SetTitleSize(0.06);
  hh->GetYaxis()->SetTitleSize(0.06);
  hh->SetLineColor(kRed);
  hh->SetMarkerColor(kRed);
  hh->SetMarkerStyle(20);
  hh->GetXaxis()->SetTitle(xlabel);
  hh->GetYaxis()->SetTitle(ylabel);
  hh->SetTitle(title);
  h->Draw("HIST");
  hh->Draw("SAME");
  TLegend *l = new TLegend(0.7,0.5,0.9,0.7);
  l->AddEntry(h,v1.c_str());
  l->AddEntry(hh,v2.c_str());
  l->Draw();

  c->SaveAs((histo1+"Eff.gif").c_str());
  delete c;
   f->Close();
   ff->Close();

}


DrawHistos(string f1,string f2,string v1,string v2,string modName,string histo,char* xlabel = "",char *ylabel = "",char *title, double scale = 1,bool eff)
{

  TFile *f = new TFile(f1.c_str());
  TFile *ff = new TFile(f2.c_str());
  
  TH1F *h = ((TH1F*)f->Get((modName+"/"+histo).c_str()))->Clone();
  TH1F *hh = ((TH1F*)ff->Get((modName+"/"+histo).c_str()))->Clone();

  if(!eff)
    h->Sumw2();

  h->GetXaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetNdivisions(509);
  h->GetYaxis()->SetNdivisions(509);
  h->GetYaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(1.25);
  h->SetMarkerColor(kYellow);
  h->SetMarkerStyle(20);
  h->SetFillColor(kYellow);
  h->GetXaxis()->SetTitle(xlabel);
  h->GetYaxis()->SetTitle(ylabel);
  h->SetTitle(title);

  if(!eff)
    {
      if(h->Integral()>0)
	if(scale==0) h->Scale(1/h->Integral());
	else
	  h->Scale(scale);
     
      hh->Sumw2();
    }
  hh->GetXaxis()->SetLabelSize(0.06);
  hh->GetYaxis()->SetLabelSize(0.06);
  hh->GetXaxis()->SetTitleSize(0.06);
  hh->GetYaxis()->SetTitleSize(0.06);
  hh->SetLineColor(kRed);
  hh->SetMarkerColor(kRed);
  hh->SetMarkerStyle(20);
  hh->GetXaxis()->SetTitle(xlabel);
  hh->GetYaxis()->SetTitle(ylabel);
  hh->SetTitle(title);
  if(!eff)
    {
      if(hh->Integral()>0)
	if(scale==0) hh->Scale(1/hh->Integral());
	else
	  hh->Scale(scale);
      
	
    }

  h->Draw("HIST");
  hh->Draw("SAME");
  TLegend *l = new TLegend(0.7,0.5,0.9,0.7);
  l->AddEntry(h,v1.c_str());
  l->AddEntry(hh,v2.c_str());
  l->Draw();

}



DrawEffHistos(string f1,string f2,string v1,string v2,string modName,string histo1,string histo2,char* xlabel = "",char *ylabel = "",char *title)
{

  TFile *f = new TFile(f1.c_str());
  TFile *ff = new TFile(f2.c_str());
  
  TH1F *h = ((TH1F*)f->Get((modName+"/"+histo1).c_str()))->Clone();
  TH1F *hd = ((TH1F*)f->Get((modName+"/"+histo2).c_str()))->Clone();

  TH1F *hh = ((TH1F*)ff->Get((modName+"/"+histo1).c_str()))->Clone();
  TH1F *hhd = ((TH1F*)ff->Get((modName+"/"+histo2).c_str()))->Clone();

 h->Sumw2();
  hd->Sumw2();
  hh->Sumw2();
  hhd->Sumw2();

  TGraphAsymmErrors* g1 = new TGraphAsymmErrors;
  g1->SetName(h->GetName());
  g1->BayesDivide(h,hd);

  g1->GetXaxis()->SetLabelSize(0.06);
  g1->GetXaxis()->SetNdivisions(509);
  g1->GetYaxis()->SetNdivisions(509);
  g1->GetYaxis()->SetLabelSize(0.06);
  g1->GetXaxis()->SetTitleSize(0.06);
  g1->GetYaxis()->SetTitleSize(0.06);
  g1->GetYaxis()->SetTitleOffset(1.25);

  TGraphAsymmErrors* g2 = new TGraphAsymmErrors;
  g2->SetName(h->GetName());
  g2->BayesDivide(hh,hhd);

  g2->GetXaxis()->SetLabelSize(0.06);
  g2->GetXaxis()->SetNdivisions(509);
  g2->GetYaxis()->SetNdivisions(509);
  g2->GetYaxis()->SetLabelSize(0.06);
  g2->GetXaxis()->SetTitleSize(0.06);
  g2->GetYaxis()->SetTitleSize(0.06);
  g2->GetYaxis()->SetTitleOffset(1.25);




  g1->SetMarkerColor(kYellow);
  g1->SetMarkerStyle(20);
  g1->SetFillColor(kYellow);
  g1->GetXaxis()->SetTitle(xlabel);
  g1->GetYaxis()->SetTitle(ylabel);
  g1->SetTitle(title);
 

  g2->SetLineColor(kRed);
  g2->SetMarkerColor(kRed);
  g2->SetMarkerStyle(20);
  g2->GetXaxis()->SetTitle(xlabel);
  g2->GetYaxis()->SetTitle(ylabel);
  g2->SetTitle(title);


  g1->Draw("AP");
  g2->Draw("Psame");
  TLegend *l = new TLegend(0.7,0.5,0.9,0.7);
  l->AddEntry(g1,v1.c_str());
  l->AddEntry(g2,v2.c_str());
  l->Draw();
}

