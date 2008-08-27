/*
    Author:  Jeremy Werner (jwerner@princeton.edu)
    Institution: Princeton University, USA
    Date: May, 2008
*/
#include "RecoLuminosity/CCCWBM/interface/CCCWBM.h"

using std::cout;
using std::endl;

CCCWBM::CCCWBM(const edm::ParameterSet& iConfig){

  listenPort  = iConfig.getUntrackedParameter<unsigned int>("SourcePort", 51001);
  DistribIP   = iConfig.getUntrackedParameter<std::string>("HLXIP","localhost");
  reconnTime  = iConfig.getUntrackedParameter<unsigned int>("ReconnectionTime",5);
  AquireMode  = iConfig.getUntrackedParameter<unsigned int>("AquireMode",0);
  outputDirectory = iConfig.getUntrackedParameter<std::string>("outputDir","/cms/mon/data/htdocs/lumi/LHCLumi/test/");

}

CCCWBM::~CCCWBM(){
  HT.Disconnect();
}


void CCCWBM::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  int errorCode;

  do{
    errorCode = HT.ReceiveLumiSection(lumiSection);

    while(errorCode !=1){
      HT.Disconnect();
      cout << "Connecting to TCPDistributor" << endl;
      errorCode = HT.Connect();
      if(errorCode != 1){
	cout << "*** Connection Failed: " << errorCode << " Attempting reconnect in " << reconnTime << " seconds." << endl;
	sleep(reconnTime);
      }
    }    
  }while(errorCode != 1);
  

  if(currentRunNumber > -1 && currentRunNumber != lumiSection.hdr.runNumber){
    this->resetHistos();
  }else{
    currentRunNumber= lumiSection.hdr.runNumber;
  }



  time_t now_time = time(NULL);
 
  int time_diff_min = int((now_time - start_time)/60); //for real

  instLumi[time_diff_min]= lumiSection.lumiSummary.InstantLumi; //for real

  std::ostringstream tmpLSIdentifier;
  tmpLSIdentifier <<", Current Run: "<<lumiSection.hdr.runNumber<<", LS: "<<lumiSection.hdr.sectionNumber<<", Orbit Range: "<<lumiSection.hdr.startOrbit<<"-"<<lumiSection.hdr.startOrbit+lumiSection.hdr.numOrbits - 1<<std::endl;
  IntLumiTitle= "Integrated Luminosity";
  IntLumiTitle+= tmpLSIdentifier.str();
  InstLumiTitle= "Instantaneous Luminosity";
  InstLumiTitle+= tmpLSIdentifier.str();
  LumiVsBXTitle= "Lumi Vs BX";
  LumiVsBXTitle+= tmpLSIdentifier.str();
  LumiSpreadTitle= "Lumi Spread";
  LumiSpreadTitle+= tmpLSIdentifier.str();

  outputfile->cd();
  int numExtraBins = 60;
  if(time_diff_min > 1380){ numExtraBins = 1440 - time_diff_min;}

  if(time_diff_min > 1439){
    //24 hours has elapsed... time to reset the histos
    std::cout<<"WARNING!!! Duration of current run has exceeded 24 hours:  Histograms will now be reset!"<<std::endl;
    this->resetHistos();

  }

  //gROOT->Reset();
  TStyle* gStyle= (TStyle*) new TStyle("Plain","Plain");
  gStyle->SetOptStat("");
  gROOT->SetBatch(kTRUE);
  gStyle->cd();
  //std::cout<<"time offset= "<<3600*timeinfo->tm_hour + 60*timeinfo->tm_min + timeinfo->tm_sec<<std::endl;
  Int_t wrongBy = 2;
  //if(timeinfo->tm_hour >= 12){ wrongBy= 2;}
  gStyle->SetTimeOffset(offset_time);
  //gStyle->SetTimeOffset(3600*timeinfo->tm_hour + 60*timeinfo->tm_min + timeinfo->tm_sec - 3*3600); //I don't know why the -3*3600 should be there... 


  //std::cout<<"really bout to init histos"<<std::endl;

  IntLumi= (TH1F*) new TH1F("IntLumi",IntLumiTitle.c_str(),time_diff_min + numExtraBins,0, 60*time_diff_min + 60*numExtraBins);
  InstLumi= (TH1F*) new TH1F("InstLumi",InstLumiTitle.c_str(),time_diff_min + numExtraBins,0, 60*time_diff_min + 60*numExtraBins); 
  LumiVsBX= (TH1F*) new TH1F("LumiVsBX",LumiVsBXTitle.c_str(),3564,0,3564);
  LumiSpread= (TH1F*) new TH1F("LumiSpread",LumiSpreadTitle.c_str(),100,0,0.2);
  
  //std::cout<<"done initing histos"<<std::endl;

  LumiVsBX->GetXaxis()->SetTitle("BX");
  LumiVsBX->GetYaxis()->SetTitle("% Inst Luminosity/BX");
  LumiVsBX->GetYaxis()->SetLabelSize(0.030);
  LumiVsBX->GetYaxis()->SetTitleOffset(1.3);
  LumiSpread->GetXaxis()->SetTitle("% of Inst Lumi");
  LumiSpread->GetYaxis()->SetTitle("#BX/(% of Inst Lumi)");
  LumiSpread->GetYaxis()->SetLabelSize(0.030);
  LumiSpread->GetYaxis()->SetTitleOffset(1.3);
  gStyle->cd();
  IntLumi->UseCurrentStyle();
  IntLumi->GetXaxis()->SetTitle(start_date.str().c_str());
  IntLumi->GetYaxis()->SetTitle("Integrated Lumi [10^31 cm-2]");
  IntLumi->GetYaxis()->SetLabelSize(0.030);
  IntLumi->GetYaxis()->SetTitleOffset(1.3);
  IntLumi->GetXaxis()->SetTimeDisplay(true);
  //IntLumi->GetXaxis()->SetTimeOffset(3600*timeinfo->tm_hour + 60*timeinfo->tm_min + timeinfo->tm_sec);
  IntLumi->GetXaxis()->SetTimeFormat("%H:%M");
  
  InstLumi->UseCurrentStyle();
  InstLumi->GetXaxis()->SetTitle(start_date.str().c_str());
  InstLumi->GetYaxis()->SetTitle("Inst. Lumi [10^31 cm-2 s-1]");
  InstLumi->GetYaxis()->SetLabelSize(0.030);
  InstLumi->GetYaxis()->SetTitleOffset(1.3);
  InstLumi->GetXaxis()->SetTimeDisplay(true);
  InstLumi->GetXaxis()->SetTimeFormat("%H:%M");

  //std::cout<<"initialized histos"<<std::endl;

  //fill the histos
  Double_t tmplumi=0;
  for(int j=0; j< lumiSection.hdr.numBunches; j++){
    tmplumi= tmplumi + lumiSection.lumiDetail.LHCLumi[j];
  }

  for(int j=0; j< lumiSection.hdr.numBunches; j++){
    LumiVsBX->SetBinContent(j+1, lumiSection.lumiDetail.LHCLumi[j]*100.0/tmplumi);
    LumiSpread->Fill(lumiSection.lumiDetail.LHCLumi[j]*100.0/tmplumi);
  }

  for(int j=0; j< time_diff_min; j++){
    InstLumi->SetBinContent(j+1, instLumi[j]); 
    intLumi[j]=InstLumi->Integral()*60.0; //inst lumis are per sec, so must multiply by 60
    IntLumi->SetBinContent(j+1, intLumi[j]); 
  }

  timeWasFilled[time_diff_min]= true;

  //std::cout<<"filled histos"<<std::endl;
  //draw histos
  TCanvas *c = new TCanvas;
  if(!c){
    std::cout<<"ERROR!!! Couldn't create canvas"<<std::endl;
    //return;
  }

  IntLumi->Draw();

  gSystem->ProcessEvents();
  TImage *img = TImage::Create();
  if (!img) {
    std::cout<<"ERROR!!!  Could not create an image"<<std::endl;
    //return;
  }
 
  //std::cout<<"created image"<<std::endl;                                                                      
  img->FromPad(c);
  //std::cout<<"got image from pad"<<std::endl;                                                                 
  

  std::string tmpFileName = outputDirectory + "IntLumi.png";
  img->WriteImage(tmpFileName.c_str());
  InstLumi->Draw();
  img->FromPad(c);
  tmpFileName = outputDirectory + "InstLumi.png";
  img->WriteImage(tmpFileName.c_str());
  LumiVsBX->Draw();
  img->FromPad(c);
  tmpFileName = outputDirectory + "LumiVsBX.png";
  img->WriteImage(tmpFileName.c_str());
  LumiSpread->Draw();
  img->FromPad(c);
  tmpFileName = outputDirectory + "LumiSpread.png";
  img->WriteImage(tmpFileName.c_str());

  outputfile->Write();

  //std::cout<<"wrote output file"<<std::endl;
  delete IntLumi;
  delete InstLumi;
  delete LumiVsBX;
  delete LumiSpread;
  delete gStyle;
  delete c;
  img->Delete();

}  //analyze done 

void CCCWBM::beginJob(const edm::EventSetup&){
 
  outputfile = (TFile*) new TFile("output.root","RECREATE");
  this->resetHistos();

  int errorCode;

  HT.SetIP(DistribIP);

  errorCode = HT.SetPort(listenPort);
  errorCode = HT.SetMode(AquireMode);
  
   do{
    errorCode = HT.Connect();
    if(errorCode != 1){
      cout << "Error code= "<<errorCode<<" Attempting to reconnect in " << reconnTime << " seconds." << endl;
      sleep(reconnTime);
    }
   }while(errorCode != 1);

}

void CCCWBM::endJob() {

  HT.Disconnect();
}

void CCCWBM::resetHistos(){

  start_time = time(NULL);
  timeinfo = localtime(&start_time);
  offset_time = 3600*timeinfo->tm_hour + 60*timeinfo->tm_min + timeinfo->tm_sec - 2*3600;
  //std::cout<<"reset time= "<<start_time<<std::endl;
  //std::cout<<"...Or local time and date: "<<asctime(timeinfo)<<std::endl;
  start_date.clear();
  start_date.str().clear();
  start_date<<"time (starting from "<<asctime(timeinfo)<<")"<<std::endl;
  //std::cout<<"start_date= "<<start_date.str().c_str()<<std::endl;
  

  //init the arrays
  for(int j=0; j<1440; j++){
    intLumi[j]= 0;
    instLumi[j]= 0;
    timeWasFilled[j]= false;
  }
  
  currentRunNumber= -1;

}

DEFINE_FWK_MODULE(CCCWBM);
