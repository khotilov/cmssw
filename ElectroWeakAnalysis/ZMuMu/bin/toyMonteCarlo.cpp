/*************************/
/*                       */
/* author: Pasquale Noli */
/* INFN Naples           */
/* Toy Montecarlo        */
/*                       */
/*************************/

//root include
#include "TRandom3.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TDirectory.h"
//std include
#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <sstream>
#include <cmath>
using namespace std;


void fillRandom(int N, TH1F *pdf, TH1F * histo){
  double m =0;
  int i=0; 
  do{
    m=pdf->GetRandom();
    if(m>=60 && m<=120){
      histo->Fill(m);
      i++;
    }
  }while( i< N);
} 


enum MuTag { globalMu, trackerMu, standaloneMu, undefinedMu };

MuTag mu(double effTrk, double effSa, TRandom3 * eventGenerator) {
  double _isTraker     = eventGenerator->Rndm();
  double _isStandAlone = eventGenerator->Rndm();
  if((_isTraker< effTrk)   && (_isStandAlone< effSa)) return globalMu;
  if((_isStandAlone< effSa)&& (_isTraker >effTrk)) return standaloneMu;
  if((_isTraker < effTrk)  && (_isStandAlone > effSa)) return trackerMu;
  return undefinedMu;
}


bool efficiencyTag(double eff, TRandom3 * eventGenerator ){
 return (eventGenerator->Rndm()< eff);
}


class BkgShape {
public:
  BkgShape(double min, double max, double slope, double a0, double a1, double a2) :
    norm_(1), min_(min), max_(max), fmax_(0),
    slope_(slope), a0_(a0), a1_(a1), a2_(a2) { 
    normalize();
  }
  double operator()(double x) const {
    if(x < min_ || x > max_) return 0;
    return norm_* exp(-slope_*x)*(a0_ + (a1_ + a2_*x)*x);
  }
  double rndm(TRandom3 * eventGenerator) const {
    double x, f;
    do {
      x = eventGenerator->Uniform(min_, max_);
      f = operator()(x);
    } while(eventGenerator->Uniform(0, fmax_) > f);
    return x;
  }
  double integral() const { return norm_; }

private:
  void normalize() {
    static const unsigned int steps = 1000;
    double s = 0, x, f;
    double base = max_ - min_;
    double dx = base/steps;
    for(unsigned int n = 0; n < steps; ++n) {
      x = min_ * n * dx;
      s += (f = operator()(x));
      if(f > fmax_) fmax_ = f;
    }
    fmax_ *= 1.001;
    norm_ = s * base;
  }
  double norm_, min_, max_, fmax_;
  double slope_, a0_, a1_, a2_;
};

int main(int argc, char * argv[]){
  TRandom3 *eventGenerator = new TRandom3();
  cout<<"int main()"<<endl;
  int o;
  char* endPtr;
  char* pdf("analysis_Z_133pb_trackIso_3.root");
  double yield(50550), effTrk(.9982), effSa(.989626), effHlt(.915496), effIso(.978575);
  double slopeMuTk(0.015556), a0MuTk(0.00035202), a1MuTk(2.99663), a2MuTk(-0.0211138);
  double slopeMuMuNonIso(0.0246876), a0MuMuNonIso(0.884777), a1MuMuNonIso(6.67684), a2MuMuNonIso(-0.0523693);
  BkgShape zMuTkBkgPdf(60, 120, slopeMuTk, a0MuTk, a1MuTk, a2MuTk);
  BkgShape zMuMuNonIsoBkgPdf(60, 120, slopeMuMuNonIso, a0MuMuNonIso, a1MuMuNonIso, a2MuMuNonIso);
 
  int expt(1), seed(1);
  
  while ((o = getopt(argc, argv,"p:n:s:y:T:S:H:I:h"))!=EOF) {
    switch(o) {
    case 'p':
      pdf  = optarg;
      break;
    case 'n':
      expt  = strtoul(optarg,&endPtr,0);
      break;
    case 's':
      seed = strtoul(optarg,&endPtr,0);
      break;
    case 'y':
      yield = strtoul(optarg,&endPtr,0);
      break;
    case 'T':
      effTrk  = strtod(optarg,&endPtr);
      break;
    case 'S':
      effSa = strtod(optarg,&endPtr);
      break;
    case 'H':
      effHlt = strtod(optarg,&endPtr);
      break;
    case 'I':
      effIso  = strtod(optarg,&endPtr);
      break;
    case 'h':
      cout<< " -p : input root file for pdf"<<endl <<" -n : number of experiment (default 1)"<<endl <<" -s : seed for generator (default 1)"<<endl <<" -T : efficiency of track (default 0.99883)"<<endl <<" -S : efficiency of standAlone(default 0.9896)"<< endl <<" -I : efficiency of Isolation (default 0.9786)" << endl << " -H : efficiency of HLT (default 0.9155)" <<endl << " -y : yield (default 50550)"<<endl;
      break;
    default:
      break;
    }
  }
  MuTag mu1,mu2;
  eventGenerator->SetSeed(seed);
  int count = 0; 
  //PDF
  TFile *inputfile = new TFile(pdf);
  TH1F *pdfzmm = (TH1F*)inputfile->Get("goodZToMuMuPlots/zMass");//pdf signal Zmumu(1hlt,2hlt), ZMuMunotIso, ZmuTk
  TH1F *pdfzmsa = (TH1F*)inputfile->Get("zmumuSaMassHistogram/zMass");//pdf signal ZmuSa
  cout<<"take pdf"<<endl;
  for(int j = 1; j <=expt; ++j){//loop on number of experiments  
    int N0 = eventGenerator->Poisson(yield);
    int nMuTkBkg = eventGenerator->Poisson(zMuTkBkgPdf.integral());
    int nMuMuNonIsoBkg = eventGenerator->Poisson(zMuMuNonIsoBkgPdf.integral());
    cout<<"loop on experiment"<<endl;

    int Nmumu = 0;
    int N2HLT = 0;
    int N1HLT = 0;
    int NISO = 0;
    int NSa = 0;
    int NTk = 0;
    for(int i = 0; i < N0; ++i){//loop on Z Yield
      mu1=mu(effTrk,effSa, eventGenerator);
      mu2=mu(effTrk,effSa, eventGenerator);
      bool iso1 =  efficiencyTag(effIso,eventGenerator);
      bool iso2 =  efficiencyTag(effIso,eventGenerator);
      bool trig1 = efficiencyTag(effHlt,eventGenerator);
      bool trig2 = efficiencyTag(effHlt,eventGenerator);
   
      if(mu1 == globalMu && mu2 == globalMu){
	 if(iso1 && iso2){//two global mu isolated
	   if(trig1 && trig2) N2HLT++;//two trigger
	   else if((trig1 && !trig2)||(!trig1 && trig2)) N1HLT++;//one trigger
	 }
	 else if(!iso1 || !iso2){//at least one not iso
	   if( trig1 || trig2) NISO++;//at least one trigger
	 }
       }//end global
       else if(((mu1 == globalMu && trig1) &&  mu2 == standaloneMu ) 
	       ||((mu2 == globalMu && trig2) && mu1 == standaloneMu )){
	 //cout<<"find standAlone"<<endl;
	 if(iso1 && iso2) {
	   NSa++;
	   // cout<<"fill standAlone"<<endl;
	 }
       }//end mu sa
      else if(((mu1 == globalMu && trig1) && mu2 == trackerMu) 
	       || ((mu2 == globalMu && trig2) && mu1 == trackerMu)){
	 if(iso1 && iso2) NTk++;
       }//end mu tk
     

    }//end of generation given the yield
    cout<<"logic end"<<endl;
    
    Nmumu = N2HLT + N1HLT;
    
    //Define signal Histo
    TH1F *zMuMu = new TH1F("zMass_golden","zMass",200,0,200);
    TH1F *zMuMu2HLT = new TH1F("zMass_2hlt","zMass",200,0,200);
    TH1F *zMuMu1HLT = new TH1F("zMass_1hlt","zMass",200,0,200);
    TH1F *zMuMuNotIso= new TH1F("zMass_noIso","zMass",200,0,200);
    TH1F *zMuSa = new TH1F("zMass_sa","zMass",200,0,200);
    TH1F *zMuTk = new TH1F("zMass_tk","zMass",200,0,200);
    pdfzmsa->SetName("zMass_safromGolden");
  
    //Fill signal Histo
   
    cout<<"Yield ="<< yield <<endl;
    cout<<"N0 ="<< N0 <<endl;
    cout<<"Nmumu = "<< Nmumu <<endl;
    cout<<"N2HLT= "<< N2HLT <<endl;
    cout<<"N1HLT = "<< N1HLT <<endl;
    cout<<"NISO = "<< NISO <<endl;
    cout<<"NSa = "<< NSa <<endl;
    cout<<"NTk = "<< NTk <<endl;

    
    fillRandom(Nmumu,pdfzmm,zMuMu);
    fillRandom(N2HLT, pdfzmm,zMuMu2HLT);
    fillRandom(N1HLT, pdfzmm,zMuMu1HLT);
    fillRandom(NISO,pdfzmm,zMuMuNotIso);
    fillRandom(NSa,pdfzmsa,zMuSa);
    fillRandom(NTk, pdfzmm,zMuTk);
    //cout<<"Signal filled"<<endl;
    
    //output	
    char head[30];
    sprintf(head,"zmm_%d",j);
    string tail =".root";
    string title = head + tail;
    
    TFile *outputfile = new TFile(title.c_str(),"RECREATE");
    
    //Hierarchy directory  
    
    TDirectory * goodZToMuMu = outputfile->mkdir("goodZToMuMuPlots");
    TDirectory * goodZToMuMu2HLT = outputfile->mkdir("goodZToMuMu2HLTPlots");
    TDirectory * goodZToMuMu1HLT = outputfile->mkdir("goodZToMuMu1HLTPlots");
    TDirectory * nonIsolatedZToMuMu = outputfile->mkdir("nonIsolatedZToMuMuPlots");
    TDirectory * goodZToMuMuOneStandAloneMuon = outputfile->mkdir("goodZToMuMuOneStandAloneMuonPlots");
    TDirectory * zmumuSaMassHistogram = outputfile->mkdir("zmumuSaMassHistogram");
    TDirectory * goodZToMuMuOneTrack = outputfile->mkdir("goodZToMuMuOneTrackPlots");
  
    goodZToMuMu->cd();
    zMuMu->Write();
    
    goodZToMuMu2HLT->cd();
    zMuMu2HLT->Write();
    
    goodZToMuMu1HLT->cd();
    zMuMu1HLT->Write();
    
    nonIsolatedZToMuMu->cd();
    zMuMuNotIso->Write();
    
    goodZToMuMuOneStandAloneMuon->cd();
    zMuSa->Write();

    zmumuSaMassHistogram->cd();
    pdfzmsa->Write();


    goodZToMuMuOneTrack->cd();
    zMuTk->Write();
    
    
    outputfile->Write();
    outputfile->Close();
    
    delete zMuMu;
    delete zMuMu2HLT;
    delete zMuMu1HLT;
    delete zMuMuNotIso;
    delete zMuSa;
    delete zMuTk;
    
    delete inputfile;
    
    
    //Define Background Histo
    TH1F *zMuMuBkg = new TH1F("zMass_golden","zMass",200,0,200);
    TH1F *zMuMu2HLTBkg = new TH1F("zMass_2hlt","zMass",200,0,200);
    TH1F *zMuMu1HLTBkg = new TH1F("zMass_1hlt","zMass",200,0,200);
    TH1F *zMuSaBkg = new TH1F("zMass_sa","zMass",200,0,200);
    TH1F *zMuSafromGoldenBkg = new TH1F("zMass_safromGolden","zMass",200,0,200);
    TH1F *zMuMuNotIsoBkg= new TH1F("zMass_noIso","zMass",200,0,200);
    TH1F *zMuTkBkg = new TH1F("zMass_tk","zMass",200,0,200);
    
    
    
    //Fill >Bkg Histograms 
    for(int i = 0; i < nMuTkBkg; ++i) {
      zMuTkBkg->Fill(zMuTkBkgPdf.rndm(eventGenerator));
    }
    for(int i = 0; i < nMuMuNonIsoBkg; ++i) {
      zMuMuNotIsoBkg->Fill(zMuMuNonIsoBkgPdf.rndm(eventGenerator));
    }
    
    char head2[30];
    sprintf(head2,"bkg_%d",j);
    string title2 = head2 + tail;
    TFile *outputfile2 = new TFile(title2.c_str(),"RECREATE");
    
    //Hierarchy directory  
    TDirectory * goodZToMuMu2 = outputfile2->mkdir("goodZToMuMuPlots");
    TDirectory * goodZToMuMu2HLT2 = outputfile2->mkdir("goodZToMuMu2HLTPlots");
    TDirectory * goodZToMuMu1HLT2 = outputfile2->mkdir("goodZToMuMu1HLTPlots");
    TDirectory * nonIsolatedZToMuMu2 = outputfile2->mkdir("nonIsolatedZToMuMuPlots");
    TDirectory * goodZToMuMuOneStandAloneMuon2 = outputfile2->mkdir("goodZToMuMuOneStandAloneMuonPlots");
    TDirectory * zmumuSaMassHistogram2 = outputfile2->mkdir("zmumuSaMassHistogram");
    TDirectory * goodZToMuMuOneTrack2 = outputfile2->mkdir("goodZToMuMuOneTrackPlots");
    
    // cout<<" bkg Hierarchy "<<endl;
    goodZToMuMu2->cd();
    zMuMuBkg->Write();
    
    goodZToMuMu2HLT2->cd();
    zMuMu2HLTBkg->Write();
    
    goodZToMuMu1HLT2->cd();
    zMuMu1HLTBkg->Write();
    
    nonIsolatedZToMuMu2->cd();
    zMuMuNotIsoBkg->Write();
    
    goodZToMuMuOneStandAloneMuon2->cd();
    zMuSaBkg->Write();

    zmumuSaMassHistogram2->cd();
    zMuSafromGoldenBkg->Write();

    goodZToMuMuOneTrack2->cd();
    zMuTkBkg->Write();
    
    outputfile2->Write();
    outputfile2->Close();
    
    delete zMuMuBkg;
    delete zMuMu2HLTBkg;
    delete zMuMu1HLTBkg;
    delete zMuMuNotIsoBkg;
    delete zMuSafromGoldenBkg;
    delete zMuSaBkg;
    delete zMuTkBkg;
    
    
    cout<<count<<"\n";
    count++;
  }//end of experiments 
     
  return 0;
  
}
