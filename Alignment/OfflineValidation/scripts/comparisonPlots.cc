#include "comparisonPlots.h"
#include <string>
#include <sstream>

#include "TProfile.h"
#include "TList.h"
#include "TNtuple.h"
#include "TString.h"
#include <iostream>
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"

comparisonPlots::comparisonPlots(std::string filename, std::string outputDir, std::string outputFilename)
{
	
	_outputDir = outputDir;
	
	fin = new TFile(filename.c_str());
	fin->cd();
	
	output = new TFile((outputDir+outputFilename).c_str(),"recreate");
	output->cd();
	
	readTree();
	
}

void comparisonPlots::readTree(){
	
	data = (TTree*)fin->Get("alignTree");
	data->SetBranchAddress("id",&id_);
	data->SetBranchAddress("mid",&mid_);
	data->SetBranchAddress("level",&level_);
	data->SetBranchAddress("mlevel",&mlevel_);
	data->SetBranchAddress("sublevel",&sublevel_);
	data->SetBranchAddress("x",&x_);
	data->SetBranchAddress("y",&y_);
	data->SetBranchAddress("z",&z_);
	data->SetBranchAddress("alpha",&alpha_);
	data->SetBranchAddress("beta",&beta_);
	data->SetBranchAddress("gamma",&gamma_);
	data->SetBranchAddress("phi",&phi_);
	data->SetBranchAddress("eta",&eta_);
	data->SetBranchAddress("r",&r_);
	data->SetBranchAddress("dx",&dx_);
	data->SetBranchAddress("dy",&dy_);
	data->SetBranchAddress("dz",&dz_);
	data->SetBranchAddress("dphi",&dphi_);
	data->SetBranchAddress("dr",&dr_);	
	data->SetBranchAddress("dalpha",&dalpha_);
	data->SetBranchAddress("dbeta",&dbeta_);
	data->SetBranchAddress("dgamma",&dgamma_);
	data->SetBranchAddress("useDetId",&useDetId_);
	data->SetBranchAddress("detDim",&detDim_);
}

void comparisonPlots::plot3x5(TCut Cut, char* dirName, bool savePlot, std::string plotName, bool autolimits){
	
	// ---------  create directory for histograms ---------
	//const char* dirName = Cut;
	TDirectory* plotDir = output->mkdir( dirName );
	
	// ---------  get right limits for histogram ---------
	double minimumR, maximumR;
	double minimumZ, maximumZ;
	double minimumPhi, maximumPhi;
	double minimumDR, maximumDR;
	double minimumDZ, maximumDZ;
	double minimumRDPhi, maximumRDPhi;
	double minimumDX, maximumDX;
	double minimumDY, maximumDY;
	if (autolimits){
		// ---------  get right limits for histogram ---------
		TH1F* phr = new TH1F("phr", "phr", 200, 0, 200);
		TH1F* phz = new TH1F("phz", "phz", 400, -300, 300);
		TH1F* phphi = new TH1F("phphi", "phphi", 200, -3.15, 3.15);
		TH1F* phdr = new TH1F("phdr", "phdr", 2000, -10, 10);
		TH1F* phdz = new TH1F("phdz", "phdz", 2000, -10, 10);
		TH1F* phrdphi = new TH1F("phrdphi", "phrdphi", 200, -10, 10);
		TH1F* phdx = new TH1F("phdx", "phy", 2000, -10, 10);
		TH1F* phdy = new TH1F("phdy", "phy", 2000, -10, 10);
		data->Project("phr","r",Cut);
		data->Project("phz","z",Cut);
		data->Project("phphi","phi",Cut);
		data->Project("phdr","dr",Cut);
		data->Project("phdz","dz",Cut);
		data->Project("phrdphi","r*dphi",Cut);
		data->Project("phdx","dx",Cut);
		data->Project("phdy","dy",Cut);
		getHistMaxMin(phr, maximumR, minimumR, 0);
		getHistMaxMin(phz, maximumZ, minimumZ, 0);
		getHistMaxMin(phphi, maximumPhi, minimumPhi, 0);
		getHistMaxMin(phdr, maximumDR, minimumDR, 1);
		getHistMaxMin(phdz, maximumDZ, minimumDZ, 1);
		getHistMaxMin(phrdphi, maximumRDPhi, minimumRDPhi, 1);
		getHistMaxMin(phdx, maximumDX, minimumDX, 1);
		getHistMaxMin(phdy, maximumDY, minimumDY, 1);
	}
	else{
		minimumR = 0., maximumR = 200.; 
		minimumZ = -300., maximumZ = 300.; 
		minimumPhi = -3.15, maximumPhi = 3.15;
		minimumDR = -1, maximumDR = 1;
		minimumDZ = -1, maximumDZ = 1;
		minimumRDPhi = -1, maximumRDPhi = 1;
		minimumDX = -1, maximumDX = 1;
		minimumDY = -1, maximumDY = 1;
	}
	
	
	// ---------  declare histograms ---------
	TH1F* h_dr = new TH1F("h_dr", "#Delta r", 2000, minimumDR, maximumDR);
	TH1F* h_dz = new TH1F("h_dz", "#Delta z", 2000, minimumDZ, maximumDZ);
	TH1F* h_rdphi = new TH1F("h_rdphi", "r* #Delta #phi", 2000, minimumRDPhi, maximumRDPhi);
	TH2F* h_drVr = new TH2F("h_drVr","#Delta r vs. r",200,minimumR,maximumR,200,minimumDR,maximumDR);
	TH2F* h_dzVr = new TH2F("h_dzVr","#Delta z vs. r",200,minimumR,maximumR,200,minimumDZ,maximumDZ);
	TH2F* h_rdphiVr = new TH2F("h_rdphiVr","r#Delta #phi vs. r",200,minimumR,maximumR,200,minimumRDPhi,maximumRDPhi);
	TH2F* h_dxVr = new TH2F("h_dxVr","#Delta x vs. r", 200,minimumR,maximumR, 200,minimumDX,maximumDX);
	TH2F* h_dyVr = new TH2F("h_dyVr","#Delta y vs. r", 200,minimumR,maximumR, 200,minimumDY,maximumDY);
	TH2F* h_drVz = new TH2F("h_drVz","#Delta r vs. z", 200,minimumZ,maximumZ, 200,minimumDR,maximumDR);
	TH2F* h_dzVz = new TH2F("h_dzVz","#Delta z vs. z", 200,minimumZ,maximumZ, 200,minimumDZ,maximumDZ);
	TH2F* h_rdphiVz = new TH2F("h_rdphiVz","r#Delta #phi vs. z", 200,minimumZ,maximumZ, 200,minimumRDPhi,maximumRDPhi);
	TH2F* h_dxVz = new TH2F("h_dxVz","#Delta x vs. z", 200,minimumZ,maximumZ, 200,minimumDX,maximumDX);
	TH2F* h_dyVz = new TH2F("h_dyVz","#Delta y vs. z", 200,minimumZ,maximumZ, 200,minimumDY,maximumDY);
	TH2F* h_drVphi = new TH2F("h_drVphi","#Delta r vs. #phi", 200,minimumPhi,maximumPhi,200,minimumDR,maximumDR);
	TH2F* h_dzVphi = new TH2F("h_dzVphi","#Delta z vs. #phi", 200,minimumPhi,maximumPhi, 200,minimumDZ,maximumDZ);
	TH2F* h_rdphiVphi = new TH2F("h_rdphiVphi","r#Delta #phi vs. #phi", 200,minimumPhi,maximumPhi,200,minimumRDPhi,maximumRDPhi);
	TH2F* h_dxVphi = new TH2F("h_dxVphi","#Delta x vs. #phi", 200,minimumPhi,maximumPhi, 200,minimumDX,maximumDX);
	TH2F* h_dyVphi = new TH2F("h_dyVphi","#Delta y vs. #phi", 200,minimumPhi,maximumPhi, 200,minimumDY,maximumDY);
	
	h_drVr->SetMarkerStyle(6);
	h_dzVr->SetMarkerStyle(6);
	h_rdphiVr->SetMarkerStyle(6);
	h_dxVr->SetMarkerStyle(6);
	h_dyVr->SetMarkerStyle(6);
	h_drVz->SetMarkerStyle(6);
	h_dzVz->SetMarkerStyle(6);
	h_rdphiVz->SetMarkerStyle(6);
	h_dxVz->SetMarkerStyle(6);
	h_dyVz->SetMarkerStyle(6);
	h_drVphi->SetMarkerStyle(6);
	h_dzVphi->SetMarkerStyle(6);
	h_rdphiVphi->SetMarkerStyle(6);
	h_dxVphi->SetMarkerStyle(6);
	h_dyVphi->SetMarkerStyle(6);
	
	TCut zCut[2];
	zCut[0] = "z < 0";
	zCut[1] = "z >= 0";
	
	// ---------  project tree onto histograms ---------
	data->Project("h_dr","dr",Cut);
	data->Project("h_dz","dz",Cut);
	data->Project("h_rdphi","r*dphi",Cut);
	data->Project("h_drVr", "dr:r",Cut);
	TGraph* gr_drVr_Array[2];
	TMultiGraph* mgr_drVr=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dr:r",Cut);	  
		gr_drVr_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_drVr_Array[i]->SetMarkerColor(4);	  
		gr_drVr_Array[i]->SetMarkerStyle(6);	  
		mgr_drVr->Add(gr_drVr_Array[i],"p");
	}	
	
	data->Project("h_dzVr", "dz:r",Cut);
	TGraph* gr_dzVr_Array[2];
	TMultiGraph* mgr_dzVr=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dz:r",Cut+zCut[i]);	  
		gr_dzVr_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_dzVr_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_dzVr_Array[i]->SetMarkerStyle(6);	  
		mgr_dzVr->Add(gr_dzVr_Array[i],"p");
	}
	
	data->Project("h_rdphiVr", "r*dphi:r",Cut);
	TGraph* gr_rdphiVr_Array[2];
	TMultiGraph* mgr_rdphiVr=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("r*dphi:r",Cut+zCut[i]);	  
		gr_rdphiVr_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_rdphiVr_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_rdphiVr_Array[i]->SetMarkerStyle(6);	  
		mgr_rdphiVr->Add(gr_rdphiVr_Array[i],"p");
	}
	
	data->Project("h_dxVr", "dx:r",Cut);
	TGraph* gr_dxVr_Array[2];
	TMultiGraph* mgr_dxVr=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dx:r",Cut+zCut[i]);	  
		gr_dxVr_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_dxVr_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_dxVr_Array[i]->SetMarkerStyle(6);	  
		mgr_dxVr->Add(gr_dxVr_Array[i],"p");
	}
	
	data->Project("h_dyVr", "dy:r",Cut);
	TGraph* gr_dyVr_Array[2];
	TMultiGraph* mgr_dyVr=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dy:r",Cut+zCut[i]);	  
		gr_dyVr_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_dyVr_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_dyVr_Array[i]->SetMarkerStyle(6);	  
		mgr_dyVr->Add(gr_dyVr_Array[i],"p");
	}
	
	data->Project("h_drVz", "dr:z",Cut);
	TGraph* gr_drVz_Array[2];
	TMultiGraph* mgr_drVz=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dr:z",Cut+zCut[i]);	  
		gr_drVz_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_drVz_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_drVz_Array[i]->SetMarkerStyle(6);	  
		mgr_drVz->Add(gr_drVz_Array[i],"p");
	}
	
	data->Project("h_dzVz", "dz:z",Cut);
	TGraph* gr_dzVz_Array[2];
	TMultiGraph* mgr_dzVz=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dz:z",Cut+zCut[i]);	  
		gr_dzVz_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_dzVz_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_dzVz_Array[i]->SetMarkerStyle(6);	  
		mgr_dzVz->Add(gr_dzVz_Array[i],"p");
	}
	
	data->Project("h_rdphiVz", "r*dphi:z",Cut);
	TGraph* gr_rdphiVz_Array[2];
	TMultiGraph* mgr_rdphiVz=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("r*dphi:z",Cut+zCut[i]);	  
		gr_rdphiVz_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_rdphiVz_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_rdphiVz_Array[i]->SetMarkerStyle(6);	  
		mgr_rdphiVz->Add(gr_rdphiVz_Array[i],"p");
	}
	
	data->Project("h_dxVz", "dx:z",Cut);
	TGraph* gr_dxVz_Array[2];
	TMultiGraph* mgr_dxVz=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dx:z",Cut+zCut[i]);	  
		gr_dxVz_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_dxVz_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_dxVz_Array[i]->SetMarkerStyle(6);	  
		mgr_dxVz->Add(gr_dxVz_Array[i],"p");
	}
	
	data->Project("h_dyVz", "dy:z",Cut);
	TGraph* gr_dyVz_Array[2];
	TMultiGraph* mgr_dyVz=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dy:z",Cut+zCut[i]);	  
		gr_dyVz_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_dyVz_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_dyVz_Array[i]->SetMarkerStyle(6);	  
		mgr_dyVz->Add(gr_dyVz_Array[i],"p");
	}
	
	data->Project("h_drVphi", "dr:phi",Cut);
	TGraph* gr_drVphi_Array[2];
	TMultiGraph* mgr_drVphi=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dr:phi",Cut+zCut[i]);	  
		gr_drVphi_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_drVphi_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_drVphi_Array[i]->SetMarkerStyle(6);	  
		mgr_drVphi->Add(gr_drVphi_Array[i],"p");
	}
	
	data->Project("h_dzVphi", "dz:phi",Cut);
	TGraph* gr_dzVphi_Array[2];
	TMultiGraph* mgr_dzVphi=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dz:phi",Cut+zCut[i]);	  
		gr_dzVphi_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_dzVphi_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_dzVphi_Array[i]->SetMarkerStyle(6);	  
		mgr_dzVphi->Add(gr_dzVphi_Array[i],"p");
	}
	
	data->Project("h_rdphiVphi", "r*dphi:phi",Cut);
	TGraph* gr_rdphiVphi_Array[2];
	TMultiGraph* mgr_rdphiVphi=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("r*dphi:phi",Cut+zCut[i]);	  
		gr_rdphiVphi_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_rdphiVphi_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_rdphiVphi_Array[i]->SetMarkerStyle(6);	  
		mgr_rdphiVphi->Add(gr_rdphiVphi_Array[i],"p");
	}
	
	data->Project("h_dxVphi", "dx:phi",Cut);
	TGraph* gr_dxVphi_Array[2];
	TMultiGraph* mgr_dxVphi=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dx:phi",Cut+zCut[i]);	  
		gr_dxVphi_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_dxVphi_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_dxVphi_Array[i]->SetMarkerStyle(6);	  
		mgr_dxVphi->Add(gr_dxVphi_Array[i],"p");
	}
	
	data->Project("h_dyVphi", "dy:phi",Cut);
	TGraph* gr_dyVphi_Array[2];
	TMultiGraph* mgr_dyVphi=new TMultiGraph();
	for ( int i = 0; i < 2; i++) {
		data->Draw("dy:phi",Cut+zCut[i]);	  
		gr_dyVphi_Array[i] = new TGraph(data->GetSelectedRows(),data->GetV2(),data->GetV1());
		gr_dyVphi_Array[i]->SetMarkerColor((i+1)*2);	  
		gr_dyVphi_Array[i]->SetMarkerStyle(6);	  
		mgr_dyVphi->Add(gr_dyVphi_Array[i],"p");
	}
	
	
	
	// ---------  draw histograms ---------
	TCanvas* c0 = new TCanvas("c0", "c0", 200, 10, 900, 300);
	c0->SetFillColor(0);
	c0->Divide(3,1);
	c0->cd(1);
	h_dr->Draw();
	c0->cd(2);
	h_dz->Draw();
	c0->cd(3);
	h_rdphi->Draw();
	if (savePlot) c0->Print((_outputDir+"plot3x1_"+plotName).c_str());
	
	mgr_drVr->SetTitle("#Delta r vs. r");           
	mgr_dzVr->SetTitle("#Delta z vs. r");           
	mgr_rdphiVr->SetTitle("r#Delta #phi vs. r");    
	mgr_dxVr->SetTitle("#Delta x vs. r");           
	mgr_dyVr->SetTitle("#Delta y vs. r");           
	mgr_drVz->SetTitle("#Delta r vs. z");           
	mgr_dzVz->SetTitle("#Delta z vs. z");           
	mgr_rdphiVz->SetTitle("r#Delta #phi vs. z");    
	mgr_dxVz->SetTitle("#Delta x vs. z");           
	mgr_dyVz->SetTitle("#Delta y vs. z");           
	mgr_drVphi->SetTitle("#Delta r vs. #phi");      
	mgr_dzVphi->SetTitle("#Delta z vs. #phi");      
	mgr_rdphiVphi->SetTitle("r#Delta #phi vs. #phi");
	mgr_dxVphi->SetTitle("#Delta x vs. #phi");      
	mgr_dyVphi->SetTitle("#Delta y vs. #phi");   
	
	// ---------  draw histograms ---------
	TCanvas* c = new TCanvas("c", "c", 200, 10, 1200, 700);
	c->SetFillColor(0);
	data->SetMarkerSize(0.5);
	data->SetMarkerStyle(6);
	c->Divide(5,3);
	c->cd(1);
	if (!autolimits) h_drVr->Draw();
	if (autolimits) mgr_drVr->Draw("a");
	c->Update();
	c->cd(2);
	if (!autolimits) h_dzVr->Draw();
	if (autolimits) mgr_dzVr->Draw("a");
	c->Update();
	c->cd(3);
	if (!autolimits) h_rdphiVr->Draw();
	if (autolimits) mgr_rdphiVr->Draw("a");
	c->cd(4);
	if (!autolimits) h_dxVr->Draw();
	if (autolimits) mgr_dxVr->Draw("a");
	c->cd(5);
	if (!autolimits) h_dyVr->Draw();
	if (autolimits) mgr_dyVr->Draw("a");
	c->cd(6);
	if (!autolimits) h_drVz->Draw();
	if (autolimits) mgr_drVz->Draw("a");
	c->cd(7);
	if (!autolimits) h_dzVz->Draw();
	if (autolimits) mgr_dzVz->Draw("a");
	c->cd(8);
	if (!autolimits) h_rdphiVz->Draw();
	if (autolimits) mgr_rdphiVz->Draw("a");
	c->cd(9);
	if (!autolimits) h_dxVz->Draw();
	if (autolimits) mgr_dxVz->Draw("a");
	c->cd(10);
	if (!autolimits) h_dyVz->Draw();
	if (autolimits) mgr_dyVz->Draw("a");
	c->cd(11);
	if (!autolimits) h_drVphi->Draw();
	if (autolimits) mgr_drVphi->Draw("a");
	c->cd(12);
	if (!autolimits) h_dzVphi->Draw();
	if (autolimits) mgr_dzVphi->Draw("a");
	c->cd(13);
	if (!autolimits) h_rdphiVphi->Draw();
	if (autolimits) mgr_rdphiVphi->Draw("a");
	c->cd(14);
	if (!autolimits) h_dxVphi->Draw();
	if (autolimits) mgr_dxVphi->Draw("a");
	c->cd(15);
	if (!autolimits) h_dyVphi->Draw();
	if (autolimits) mgr_dyVphi->Draw("a");
	
	c->Update();
	
	// ---------  set output directory for histograms ---------
	plotDir->cd();
	h_dr->Write(); h_dz->Write(); h_rdphi->Write();
	h_drVr->Write(); h_dzVr->Write(); h_rdphiVr->Write(); h_dxVr->Write(); h_dyVr->Write();
	h_drVz->Write(); h_dzVz->Write(); h_rdphiVz->Write(); h_dxVz->Write(); h_dyVz->Write();
	h_drVphi->Write(); h_dzVphi->Write(); h_rdphiVphi->Write(); h_dxVphi->Write(); h_dyVphi->Write();
	
	if (savePlot) c->Print((_outputDir+"plot3x5_"+plotName).c_str());
	
}

void comparisonPlots::plot3x5Profile(TCut Cut, char* dirName, int nBins, bool savePlot, std::string plotName, bool autolimits){
	
	// ---------  create directory for histograms ---------
	//const char* dirName = Cut;
	string s;// = "profile";
	s = s + dirName;
	s.append("_profile");
	TDirectory* plotDir = output->mkdir( s.data() );
	
	
	double minimumR, maximumR;
	double minimumZ, maximumZ;
	double minimumPhi, maximumPhi;
	double minimumDR, maximumDR;
	double minimumDZ, maximumDZ;
	double minimumRDPhi, maximumRDPhi;
	double minimumDX, maximumDX;
	double minimumDY, maximumDY;
	if (autolimits){
		// ---------  get right limits for histogram ---------
		TH1F* phr = new TH1F("phr", "phr", 200, 0, 200);
		TH1F* phz = new TH1F("phz", "phz", 400, -300, 300);
		TH1F* phphi = new TH1F("phphi", "phphi", 200, -3.15, 3.15);
		TH1F* phdr = new TH1F("phdr", "phdr", 2000, -10, 10);
		TH1F* phdz = new TH1F("phdz", "phdz", 2000, -10, 10);
		TH1F* phrdphi = new TH1F("phrdphi", "phrdphi", 200, -10, 10);
		TH1F* phdx = new TH1F("phdx", "phy", 2000, -10, 10);
		TH1F* phdy = new TH1F("phdy", "phy", 2000, -10, 10);
		data->Project("phr","r",Cut);
		data->Project("phz","z",Cut);
		data->Project("phphi","phi",Cut);
		data->Project("phdr","dr",Cut);
		data->Project("phdz","dz",Cut);
		data->Project("phrdphi","r*dphi",Cut);
		data->Project("phdx","dx",Cut);
		data->Project("phdy","dy",Cut);
		getHistMaxMin(phr, maximumR, minimumR, 0);
		getHistMaxMin(phz, maximumZ, minimumZ, 0);
		getHistMaxMin(phphi, maximumPhi, minimumPhi, 0);
		getHistMaxMin(phdr, maximumDR, minimumDR, 1);
		getHistMaxMin(phdz, maximumDZ, minimumDZ, 1);
		getHistMaxMin(phrdphi, maximumRDPhi, minimumRDPhi, 1);
		getHistMaxMin(phdx, maximumDX, minimumDX, 1);
		getHistMaxMin(phdy, maximumDY, minimumDY, 1);
	}
	else{
		minimumR = 0., maximumR = 200.; 
		minimumZ = -300., maximumZ = 300.; 
		minimumPhi = -3.15, maximumPhi = 3.15;
		minimumDR = -1, maximumDR = 1;
		minimumDZ = -1, maximumDZ = 1;
		minimumRDPhi = -1, maximumRDPhi = 1;
		minimumDX = -1, maximumDX = 1;
		minimumDY = -1, maximumDY = 1;
	}
	
	// ---------  declare histograms ---------
	TProfile* hprof_drVr = new TProfile("hprof_drVr","#Delta r vs. r",nBins,minimumR,maximumR,minimumDR,maximumDR);
	TProfile* hprof_dzVr = new TProfile("hprof_dzVr","#Delta z vs. r",nBins,minimumR,maximumR,minimumDZ,maximumDZ);
	TProfile* hprof_rdphiVr = new TProfile("hprof_rdphiVr","r#Delta #phi vs. r",nBins,minimumR,maximumR,minimumRDPhi,maximumRDPhi);
	TProfile* hprof_dxVr = new TProfile("hprof_dxVr","#Delta x vs. r", nBins,minimumR,maximumR,minimumDX,maximumDX);
	TProfile* hprof_dyVr = new TProfile("hprof_dyVr","#Delta y vs. r", nBins,minimumR,maximumR,minimumDY,maximumDY);
	TProfile* hprof_drVz = new TProfile("hprof_drVz","#Delta r vs. z", nBins,minimumZ,maximumZ,minimumDR,maximumDR);
	TProfile* hprof_dzVz = new TProfile("hprof_dzVz","#Delta z vs. z", nBins,minimumZ,maximumZ,minimumDZ,maximumDZ);
	TProfile* hprof_rdphiVz = new TProfile("hprof_rdphiVz","r#Delta #phi vs. z", nBins,minimumZ,maximumZ,minimumRDPhi,maximumRDPhi);
	TProfile* hprof_dxVz = new TProfile("hprof_dxVz","#Delta x vs. z", nBins,minimumZ,maximumZ,minimumDX,maximumDX);
	TProfile* hprof_dyVz = new TProfile("hprof_dyVz","#Delta y vs. z", nBins,minimumZ,maximumZ,minimumDY,maximumDY);
	TProfile* hprof_drVphi = new TProfile("hprof_drVphi","#Delta r vs. #phi", nBins,minimumPhi,maximumPhi,minimumDR,maximumDR);
	TProfile* hprof_dzVphi = new TProfile("hprof_dzVphi","#Delta z vs. #phi", nBins,minimumPhi,maximumPhi,minimumDZ,maximumDZ);
	TProfile* hprof_rdphiVphi = new TProfile("hprof_rdphiVphi","r#Delta #phi vs. #phi", nBins,minimumPhi,maximumPhi,minimumRDPhi,maximumRDPhi);
	TProfile* hprof_dxVphi = new TProfile("hprof_dxVphi","#Delta x vs. #phi", nBins,minimumPhi,maximumPhi,minimumDX,maximumDX);
	TProfile* hprof_dyVphi = new TProfile("hprof_dyVphi","#Delta y vs. #phi", nBins,minimumPhi,maximumPhi,minimumDY,maximumDY);
	
	// ---------  project tree onto histograms ---------
	data->Project("hprof_drVr", "dr:r",Cut,"prof");
	data->Project("hprof_dzVr", "dz:r",Cut,"prof");
	data->Project("hprof_rdphiVr", "r*dphi:r",Cut,"prof");
	data->Project("hprof_dxVr", "dx:r",Cut,"prof");
	data->Project("hprof_dyVr", "dy:r",Cut,"prof");
	data->Project("hprof_drVz", "dr:z",Cut,"prof");
	data->Project("hprof_dzVz", "dz:z",Cut,"prof");
	data->Project("hprof_rdphiVz", "r*dphi:z",Cut,"prof");
	data->Project("hprof_dxVz", "dx:z",Cut,"prof");
	data->Project("hprof_dyVz", "dy:z",Cut,"prof");
	data->Project("hprof_drVphi", "dr:phi",Cut,"prof");
	data->Project("hprof_dzVphi", "dz:phi",Cut,"prof");
	data->Project("hprof_rdphiVphi", "r*dphi:phi",Cut,"prof");
	data->Project("hprof_dxVphi", "dx:phi",Cut,"prof");
	data->Project("hprof_dyVphi", "dy:phi",Cut,"prof");
	
	// ---------  draw histograms ---------
	TCanvas* cp = new TCanvas("cp", "cp", 200, 10, 1200, 700);
	cp->SetFillColor(0);
	data->SetMarkerSize(0.5);
	data->SetMarkerStyle(6);
	cp->Divide(5,3);
	cp->cd(1);
	hprof_drVr->Draw();
	cp->cd(2);
	hprof_dzVr->Draw();
	cp->cd(3);
	hprof_rdphiVr->Draw();
	cp->cd(4);
	hprof_dxVr->Draw();
	cp->cd(5);
	hprof_dyVr->Draw();
	cp->cd(6);
	hprof_drVz->Draw();
	cp->cd(7);
	hprof_dzVz->Draw();
	cp->cd(8);
	hprof_rdphiVz->Draw();
	cp->cd(9);
	hprof_dxVz->Draw();
	cp->cd(10);
	hprof_dyVz->Draw();
	cp->cd(11);
	hprof_drVphi->Draw();
	cp->cd(12);
	hprof_dzVphi->Draw();
	cp->cd(13);
	hprof_rdphiVphi->Draw();
	cp->cd(14);
	hprof_dxVphi->Draw();
	cp->cd(15);
	hprof_dyVphi->Draw();
	
	// ---------  set output directory for histograms ---------
	plotDir->cd();
	hprof_drVr->Write(); hprof_dzVr->Write(); hprof_rdphiVr->Write(); hprof_dxVr->Write(); hprof_dyVr->Write();
	hprof_drVz->Write(); hprof_dzVz->Write(); hprof_rdphiVz->Write(); hprof_dxVz->Write(); hprof_dyVz->Write();
	hprof_drVphi->Write(); hprof_dzVphi->Write(); hprof_rdphiVphi->Write(); hprof_dxVphi->Write(); hprof_dyVphi->Write();
	
	if (savePlot) cp->Print((_outputDir+"plot3x5Prof_"+plotName).c_str());
}


void comparisonPlots::getMaxMin(){
	
	data->GetEntry(0);
	minR = r_; maxR = r_;
	minZ = z_; maxZ = z_;
	minPhi = phi_; maxPhi = phi_;
	minDR = dr_; maxDR = dr_;
	minDZ = dz_; maxDZ = dz_;
	minRDPhi = r_*dphi_; maxRDPhi = r_*dphi_;
	minDX = dx_; maxDX = dx_;
	minDY = dy_; maxDY = dy_;
	
	int nEntries = data->GetEntries();
	for (int i = 1; i < nEntries; ++i){
		data->GetEntry(i);
		
		if (r_ < minR) minR = r_;
		if (r_ > maxR) maxR = r_;
		if (z_ < minZ) minZ = z_;
		if (z_ > maxZ) maxZ = z_;
		if (phi_ < minPhi) minPhi = phi_;
		if (phi_ > maxPhi) maxPhi = phi_;
		if (dr_ < minDR) minDR = dr_;
		if (dr_ > maxDR) maxDR = dr_;
		if (dz_ < minDZ) minDZ = dz_;
		if (dz_ > maxDZ) maxDZ = dz_;
		if (r_*dphi_ < minRDPhi) minRDPhi = r_*dphi_;
		if (r_*dphi_ > maxRDPhi) maxRDPhi = r_*dphi_;
		if (dx_ < minDX) minDX = dx_;
		if (dx_ > maxDX) maxDX = dx_;
		if (dy_ < minDY) minDY = dy_;
		if (dy_ > maxDY) maxDY = dy_;
	}
}


void comparisonPlots::getHistMaxMin(TH1* hist, double &max, double &min, int flag){
	
	int nBins = hist->GetNbinsX();
	for (int i = 0; i < nBins; ++i){
		double binContent = hist->GetBinContent(i);
		if (binContent > 0){
			//double binWidth = hist->GetBinLowEdge(i) - hist->GetBinLowEdge(i-1);
			//std::cout << "bin width1: " << hist->GetBinWidth(i) << ", bin width2: " << binWidth << std::endl;
			if (flag == 0) max = hist->GetBinLowEdge(i) + 2.*hist->GetBinWidth(i) + 0.1*hist->GetBinLowEdge(i);
			if (flag == 1) max = hist->GetBinLowEdge(i) + hist->GetBinWidth(i) + 0.1*hist->GetBinLowEdge(i);
		}
	}
	for (int i = (nBins-1); i >= 0; i--){
		double binContent = hist->GetBinContent(i);
		if (binContent > 0) min = hist->GetBinLowEdge(i) - 0.1*hist->GetBinLowEdge(i);
	}
	//std::cout << "max: " << max << ", min: " << min << std::endl;
}

void comparisonPlots::Write()
{
	output->Write();
}


