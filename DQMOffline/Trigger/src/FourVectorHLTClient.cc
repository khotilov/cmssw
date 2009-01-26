/*

   FourVectorHLTClient.cc
   author:  Vladimir Rekovic, U Minn. 
   date of first version: Sept 2008

*/
//$Id: FourVectorHLTClient.cc,v 1.5 2009/01/20 11:18:04 rekovic Exp $

#include "DQMOffline/Trigger/interface/FourVectorHLTClient.h"

#include "DQMServices/Core/interface/QTest.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DQMServices/Core/interface/QReport.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "TRandom.h"
#include <TF1.h>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <memory>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include "TROOT.h"


using namespace edm;
using namespace std;

FourVectorHLTClient::FourVectorHLTClient(const edm::ParameterSet& ps)
{
  parameters_=ps;
  initialize();

}

FourVectorHLTClient::~FourVectorHLTClient(){
  LogDebug("FourVectorHLTClient")<< "FourVectorHLTClient: ending...." ;
}

//--------------------------------------------------------
void FourVectorHLTClient::initialize(){ 

  counterLS_=0; 
  counterEvt_=0; 
  
  // get back-end interface
  dbe_ = Service<DQMStore>().operator->();
  

  // base folder for the contents of this job
  sourceDir_ = TString(parameters_.getUntrackedParameter<string>("hltSourceDir",""));
  LogDebug("FourVectorHLTClient")<< "Monitor dir = " << sourceDir_ << endl;
    
  clientDir_ = TString(parameters_.getUntrackedParameter<string>("hltClientDir",""));
  LogDebug("FourVectorHLTClient")<< "Client dir = " << clientDir_ << endl;
    
  prescaleLS_ = parameters_.getUntrackedParameter<int>("prescaleLS", -1);
  LogDebug("FourVectorHLTClient")<< "DQM lumi section prescale = " << prescaleLS_ << " lumi section(s)"<< endl;
  
  prescaleEvt_ = parameters_.getUntrackedParameter<int>("prescaleEvt", -1);
  LogDebug("FourVectorHLTClient")<< "DQM event prescale = " << prescaleEvt_ << " events(s)"<< endl;
  
  std::vector<edm::ParameterSet> effpaths = parameters_.getParameter<std::vector<edm::ParameterSet> >("effpaths");

  std::pair<std::string, std::string> custompathnamepair;
  for(std::vector<edm::ParameterSet>::iterator pathconf = effpaths.begin() ; pathconf != effpaths.end(); 
      pathconf++) {
       custompathnamepair.first =pathconf->getParameter<std::string>("pathname"); 
       custompathnamepair.second = pathconf->getParameter<std::string>("denompathname");   
       custompathnamepairs_.push_back(custompathnamepair);
       //    customdenompathnames_.push_back(pathconf->getParameter<std::string>("denompathname"));  
       // custompathnames_.push_back(pathconf->getParameter<std::string>("pathname"));  
    }


      
}

//--------------------------------------------------------
void FourVectorHLTClient::beginJob(const EventSetup& context){


  LogDebug("FourVectorHLTClient")<<"[FourVectorHLTClient]: beginJob" << endl;
  // get backendinterface  
  dbe_ = Service<DQMStore>().operator->();



  

}

//--------------------------------------------------------
void FourVectorHLTClient::beginRun(const Run& r, const EventSetup& context) {

  LogDebug("FourVectorHLTClient")<<"[FourVectorHLTClient]: beginRun" << endl;
	/*
  TString summaryFolder = clientDir_ + TString("/reportSummaryContents/");
  TString summaryPath = summaryFolder + TString("reportSummary");

  dbe_->setCurrentFolder(summaryFolder.Data());
  

  reportSummary_ = dbe_->get(summaryPath.Data());
  if ( reportSummary_ ) {
      dbe_->removeElement(reportSummary_->getName()); 
  }

  reportSummary_ = dbe_->bookFloat("reportSummary");

  //initialize reportSummary to 1
  if (reportSummary_) reportSummary_->Fill(-999);
	*/

}

//--------------------------------------------------------
void FourVectorHLTClient::beginLuminosityBlock(const LuminosityBlock& lumiSeg, const EventSetup& context) {
   // optionally reset histograms here
}

void FourVectorHLTClient::endLuminosityBlock(const edm::LuminosityBlock& lumiSeg, const edm::EventSetup& c){   


} 



//--------------------------------------------------------
void FourVectorHLTClient::analyze(const Event& e, const EventSetup& context){
   
  counterEvt_++;
  if (prescaleEvt_<1) return;
  if (prescaleEvt_>0 && counterEvt_%prescaleEvt_ != 0) return;
  
  LogDebug("FourVectorHLTClient")<<"analyze..." << endl;

  


  //reportSummary = average of report summaries of each system
  
 
}

//--------------------------------------------------------
void FourVectorHLTClient::endRun(const Run& r, const EventSetup& context){
 LogDebug("FourVectorHLTClient")<<"FourVectorHLTClient:: endLuminosityBlock "  << endl;

 // QTests
 ////////////////////////////
 Comp2RefKolmogorov * kl_test_ = new Comp2RefKolmogorov("my_kolm");

 // HLT paths
 TObjArray *hltPathNameColl = new TObjArray(100);
 hltPathNameColl->SetOwner();

 if(! dbe_->containsAnyMonitorable(sourceDir_.Data())) {

   LogDebug("FourVectorHLTClient")<<"[FourVectorHLTClient] endLuminosityBlock: sourceDir_(" << sourceDir_.Data() << ") has no MEs " << endl;
    return;
 }

 LogDebug("FourVectorHLTClient")<<"[FourVectorHLTClient] endLuminosityBlock: Interesting dir has MEs " << endl;

  hltMEs = dbe_->getContents(sourceDir_.Data());
  
 // For easier work, build vector<TString> to hold HLT ME names
 for(unsigned int i=0;i<hltMEs.size();i++) {

   if(hltMEs[i]->getName().find("HLT") == 0) hltMEName.push_back(TString(hltMEs[i]->getName()));

 }

 // loop over all HLT DQM MEs to get HLT path names
 //////////////////////////////////////////////////
 LogDebug("FourVectorHLTClient")<<"[FourVectorHLTClient] HLT DQM MEs are: "  << endl;

 for(unsigned int i=0;i<hltMEName.size();i++) {

   LogDebug("FourVectorHLTClient")<< hltMEName[i]  << endl;

   // collect HLT path names from MEs names, removing " *" from the end
	 // e.g, of a ME name: "HLT_Mu3 HLT_Jet30_NOn"
   ////////////////////////////////////////////////////////////////
   TString tempName = hltMEName[i];
   tempName.Remove(tempName.Last(' '),tempName.Length());
   hltPathName.push_back(tempName);

   TObjString *tempHltName = new TObjString(tempName.Data());

   if( !hltPathNameColl->Contains(tempHltName)) {

     hltPathNameColl->Add(tempHltName);

   } 

 } //end for
  
 LogDebug("FourVectorHLTClient")<< "Number of derived HLT paths from DQM source directory: " << hltPathNameColl->GetEntriesFast()  << endl;
 
 ///////////////////////////
 // Loop over all path names
 ///////////////////////////
 for (int j=0;j<hltPathNameColl->GetEntriesFast();j++) {


   TObjString *curHltPath = (TObjString*)hltPathNameColl->At(j);
   TString pathPrefis = curHltPath->String() + TString("_");
   LogDebug("FourVectorHLTClient")<< " path " << curHltPath->String() << endl;

   TString pathFolder = clientDir_ + TString("/paths/") + curHltPath->String() + TString("/any");
   dbe_->setCurrentFolder(pathFolder.Data());

   ///////////////////////////////////////////
   // KolmogorovTest for each DQM histogram
   ///////////////////////////////////////////
   for(unsigned int i=0;i<hltMEName.size();i++) {

     TString tempHistName = hltMEName[i];

     // Only consider histos belonging to current HLT path
		 // Such a histo should have "pathName" at the very begining
		 // of its name.  E.g.  histogram "HLT_Mu3 HLT_Jet30_NOn" 
		 // belongs to the path HLT_Mu3, and not to the path HLT_Jet30.
     //if( !tempHistName.Contains(curHltPath->String()) ) continue;
		 TString mainPath = curHltPath->String()+" ";
     //if( tempHistName.Index(curHltPath->String().Data(), curHltPath->String().Length(),0) !=0 ) continue;
     if( tempHistName.Index(mainPath.Data(), mainPath.Length(),0,TString::kExact) !=0 ) continue;

     // Some string gymnastics:
     //
     // there are paths like HLT_Jet180 and HLT_Jet180_MET60
     // make sure we dont cross use histograms of such paths
     // histogram name can only have one character "_" in addition to the pathname
     TString radical = tempHistName;
     //radical.ReplaceAll(curHltPath->String()+"_","");
     radical.Remove(0,tempHistName.Last('_'));
      LogDebug("FourVectorHLTClient")<< " tempHistName is " << tempHistName << "     radical is: " << radical << endl;

     if(radical.CountChar('_')>1) continue;


		 TString diagName = curHltPath->String() +" "+ curHltPath->String();
		 if(tempHistName.Contains(diagName,TString::kExact)) {
   	   pathFolder = clientDir_ + TString("/paths/") + curHltPath->String() + TString("/distributions");
		 }
		 else {
   	  pathFolder = clientDir_ + TString("/paths/") + curHltPath->String() + TString("/cross_distributions");
		 }

     dbe_->setCurrentFolder(pathFolder.Data());
     
     TString histName = sourceDir_+TString("/")+tempHistName;
     TString klmgrvName = tempHistName+TString("_klmgrTest");

     //klmgrvTest_ = dbe_->bookFloat(klmgrvName.Data());
     //if (klmgrvTest_) klmgrvTest_->Fill(-999);

     MonitorElement *histME = dbe_->get(histName.Data());


     //std::cout << " Running test " << chi2_test_->getName() 
               //<< " (Algorithm: " << chi2_test_->getAlgoName() << ") "
               //<< std::endl;
     //float prob = chi2_test_->runTest(my_test);


      LogDebug("FourVectorHLTClient")<< " histME pathName:" << histME->getPathname() << ",   kind: " << histME->kind() << "  (TH1F=" << MonitorElement::DQM_KIND_TH1F << ")" << endl;

     if(histME->kind() == MonitorElement::DQM_KIND_TH1F) {

       TH1F *hist = new TH1F();
       hist = histME->getTH1F();
       dbe_->book1D(hist->GetName(),hist);

       // must be non empty hist to do KolmogorovTest
       if(hist->GetEntries() == 0) continue;
       if(hist->Integral() == 0) continue;

       LogDebug("FourVectorHLTClient")<< "endRun:   histName " << hist->GetName() << endl;
       LogDebug("FourVectorHLTClient")<< "endRun:   hist entries:" << hist->GetEntries() << endl;

       float kl_prob = kl_test_->runTest(histME);

       LogDebug("FourVectorHLTClient")<< "endRun:   KLMGRV test = " << kl_prob << endl;
     
       //if (klmgrvTest_) klmgrvTest_->Fill(kl_prob);

		 }

     if(histME->kind() == MonitorElement::DQM_KIND_TH2F) {

       TH2F *hist = new TH2F();
       hist = histME->getTH2F();
       dbe_->book2D(hist->GetName(),hist);

			 /*
       // must be non empty hist to do KolmogorovTest
       if(hist->GetEntries() == 0) continue;
       if(hist->Integral() == 0) continue;

       LogDebug("FourVectorHLTClient")<< "endRun:   histName " << hist->GetName() << endl;
       LogDebug("FourVectorHLTClient")<< "endRun:   hist entries:" << hist->GetEntries() << endl;

       float kl_prob = kl_test_->runTest(histME);

       LogDebug("FourVectorHLTClient")<< "endRun:   KLMGRV test = " << kl_prob << endl;
     
       if (klmgrvTest_) klmgrvTest_->Fill(kl_prob);
			 */

		 }



   } // end for

    ///////////////////////////
    // Efficiencies
    ///////////////////////////
    pathFolder = clientDir_ + TString("/paths/") + curHltPath->String() + TString("/efficiencies");

    dbe_->setCurrentFolder(pathFolder.Data());

    for (std::vector<std::pair<std::string, std::string> >::iterator custompathnamepair = custompathnamepairs_.begin(); custompathnamepair != custompathnamepairs_.end(); ++custompathnamepair)
    {

      TString numPathName=TString(custompathnamepair->first);
      if(!curHltPath->String().Contains(numPathName,TString::kExact)) continue;

			TString denPathName=TString(custompathnamepair->second);


    vector<TString> vObj;
    vObj.push_back(TString("offPhi"));
    vObj.push_back(TString("offEta"));
    vObj.push_back(TString("offEt"));

    for (unsigned int k=0; k<vObj.size();k++) {


      //TString oldHistPathNum = sourceDir_+TString("/")+curHltPath->String()+" "+curHltPath->String()+TString("_")+vObj[k]+TString("L1On"); 
      //TString oldHistPathDen = sourceDir_+TString("/")+curHltPath->String()+" "+curHltPath->String()+TString("_")+vObj[k]+TString("L1"); 
      TString oldHistPathNum = sourceDir_+TString("/")+numPathName+" "+denPathName+TString("_")+vObj[k]+TString("OnOff"); 
      TString oldHistPathDen = sourceDir_+TString("/")+numPathName+" "+denPathName+TString("_")+vObj[k]+TString("Off"); 
  
  
      TString newHistName = curHltPath->String()+TString("_"+vObj[k]+"_Eff_HLTtoL1");
      TString newHistTitle = numPathName+" given " + denPathName +TString("  "+vObj[k]+" Eff  HLTtoL1");
      TString newHistPath = pathFolder+newHistName;
  
      MonitorElement *numME = dbe_->get(oldHistPathNum.Data());
      MonitorElement *denME = dbe_->get(oldHistPathDen.Data());
  
      //check if HLTOffline histogram exist
      if ( numME &&  denME ) {
      
        //check if booked HLTClient histogram exist
        if ( dbe_->get(newHistPath.Data()) ) {
          dbe_->removeElement(newHistPath.Data());
        }
      
      
        TH1F* numHist = numME->getTH1F();
        TH1F* denHist = denME->getTH1F();
  
        TH1F* effHist = (TH1F*) numHist->Clone(newHistName.Data());
        effHist->SetTitle(newHistTitle.Data());
        effHist->Sumw2();
        effHist->Divide(denHist);
        
  
        //reportSummaryMap_ = dbe_->book1D(numHist->Divide(denHist));
        
        dbe_->book1D(newHistName.Data(), effHist);
  
      } // end if
      else {
  
        LogDebug("FourVectorHLTClient")<< "Cannot find NUM and DEN histograms to derive " << newHistName <<  "."<< endl;
  

      }
    } //end for obj k
   } // end for custompathpair

  } //end for 

  hltPathNameColl->Delete();
  delete kl_test_;


}

//--------------------------------------------------------
void FourVectorHLTClient::endJob(){
}



TH1F * FourVectorHLTClient::get1DHisto(string meName, DQMStore * dbi)
{

  MonitorElement * me_ = dbi->get(meName);

  if (!me_) { 
    LogDebug("FourVectorHLTClient")<< "ME NOT FOUND." << endl;
    return NULL;
  }

  return me_->getTH1F();
}

TH2F * FourVectorHLTClient::get2DHisto(string meName, DQMStore * dbi)
{


  MonitorElement * me_ = dbi->get(meName);

  if (!me_) { 
    LogDebug("FourVectorHLTClient")<< "ME NOT FOUND." << endl;
    return NULL;
  }

  return me_->getTH2F();
}



TProfile2D *  FourVectorHLTClient::get2DProfile(string meName, DQMStore * dbi)
{


  MonitorElement * me_ = dbi->get(meName);

  if (!me_) { 
    LogDebug("FourVectorHLTClient")<< "ME NOT FOUND." << endl;
   return NULL;
  }

  return me_->getTProfile2D();
}


TProfile *  FourVectorHLTClient::get1DProfile(string meName, DQMStore * dbi)
{


  MonitorElement * me_ = dbi->get(meName);

  if (!me_) { 
    LogDebug("FourVectorHLTClient")<< "ME NOT FOUND." << endl;
    return NULL;
  }

  return me_->getTProfile();
}








