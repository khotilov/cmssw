/*
 * \file DQMEventInfo.cc
 * \author M. Zanetti - CERN PH
 * Last Update:
 * $Date: 2009/12/11 13:11:16 $
 * $Revision: 1.4 $
 * $Author: ameyer $
 *
 */

#include "DQMProvInfo.h"
#include <TSystem.h>
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

// Framework


#include <stdio.h>
#include <sstream>
#include <math.h>

using namespace edm;
using namespace std;

DQMProvInfo::DQMProvInfo(const ParameterSet& ps){
  
  parameters_ = ps;

  dbe_ = edm::Service<DQMStore>().operator->();

  provinfofolder_ = parameters_.getUntrackedParameter<string>("provInfoFolder", "ProvInfo") ;
  subsystemname_ = parameters_.getUntrackedParameter<string>("subSystemFolder", "YourSubsystem") ;
  makedcsinfo_ = parameters_.getUntrackedParameter<bool>("makeDcsInfo",true);
  dcsword_=0xffff; // set true and switch off in case a single event in a given LS has all subsys off.
  physDecl_=true; // set true and switch off in case a single event in a given LS does not have it set.
  for (int i=0;i<24;i++) dcs24[i]=true;
}

DQMProvInfo::~DQMProvInfo(){
}

void 
DQMProvInfo::beginRun(const edm::Run& r, const edm::EventSetup &c ) {

  makeProvInfo();
  
  dbe_->setCurrentFolder(subsystemname_ +"/EventInfo/");

  if (makedcsinfo_)
  {
    reportSummaryMap_ = dbe_->book2D("reportSummaryMap",
                	"Info Report Summary Map", 15, 1., 16., 10, 1., 11.);
    reportSummaryMap_->setBinLabel(10, "Physics", 2);
    reportSummaryMap_->setBinLabel( 9, "CSC",     2);
    reportSummaryMap_->setBinLabel( 8, "DT",      2);
    reportSummaryMap_->setBinLabel( 7, "EB",      2);
    reportSummaryMap_->setBinLabel( 6, "EE",      2);
    reportSummaryMap_->setBinLabel( 5, "ES",      2);
    reportSummaryMap_->setBinLabel( 4, "Hcal",    2);
    reportSummaryMap_->setBinLabel( 3, "Pixel",   2);
    reportSummaryMap_->setBinLabel( 2, "RPC",     2);
    reportSummaryMap_->setBinLabel( 1, "Strip",   2);
    for (int i=1;i<16;i++)
      for (int j=1;j<11;j++)
	reportSummaryMap_->setBinContent(i,j,-1.);

    dbe_->cd();
    dbe_->setCurrentFolder( subsystemname_ + "/Conditions") ;
    dcsVsLumi_ = dbe_->book2D("dcsVsLumi",
                       "DCS vs Lumi", 200, 0., 200., 24, 0., 24.);
    dcsVsLumi_->setBinLabel(1,"CSCp",2);   
    dcsVsLumi_->setBinLabel(2,"CSCm",2);   
    dcsVsLumi_->setBinLabel(3,"DT0",2);    
    dcsVsLumi_->setBinLabel(4,"DTp",2);    
    dcsVsLumi_->setBinLabel(5,"DTm",2);    
    dcsVsLumi_->setBinLabel(6,"EBp",2);    
    dcsVsLumi_->setBinLabel(7,"EBm",2);    
    dcsVsLumi_->setBinLabel(8,"EEp",2);    
    dcsVsLumi_->setBinLabel(9,"EEm",2);    
    dcsVsLumi_->setBinLabel(10,"ESp",2);    
    dcsVsLumi_->setBinLabel(11,"ESm",2);   
    dcsVsLumi_->setBinLabel(12,"HBHEa",2); 
    dcsVsLumi_->setBinLabel(13,"HBHEb",2); 
    dcsVsLumi_->setBinLabel(14,"HBHEc",2); 
    dcsVsLumi_->setBinLabel(15,"HF",2);    
    dcsVsLumi_->setBinLabel(16,"HO",2);    
    dcsVsLumi_->setBinLabel(17,"BPIX",2);  
    dcsVsLumi_->setBinLabel(18,"FPIX",2);  
    dcsVsLumi_->setBinLabel(19,"RPC",2);   
    dcsVsLumi_->setBinLabel(20,"TIBTID",2);
    dcsVsLumi_->setBinLabel(21,"TOB",2);   
    dcsVsLumi_->setBinLabel(22,"TECp",2);  
    dcsVsLumi_->setBinLabel(23,"TECm",2);  
    dcsVsLumi_->setBinLabel(24,"CASTOR",2);
    
  }
  else
  { 
    reportSummaryMap_ = dbe_->book2D("reportSummaryMap",
                	"Info Report Summary Map", 1, 0., 1., 1, 0., 1.);
    reportSummaryMap_->setBinContent(1,1,1.);
  }
  
  reportSummary_=dbe_->bookFloat("reportSummary");
  reportSummary_->Fill(1.); // to be refined based on some algorithm

  dcsword_=0xffff;
  physDecl_=true;
  for (int i=0;i<24;i++) dcs24[i]=true;
} 

void DQMProvInfo::analyze(const Event& e, const EventSetup& c){
 
//  makeProvInfo();

  if (makedcsinfo_) 
  { 
     makeDcsInfo(e);
     makeGtInfo(e);
  }
  return;
}

void
DQMProvInfo::endLuminosityBlock(const edm::LuminosityBlock& l, const edm::EventSetup& c)
{

 
  if (!makedcsinfo_) return;

  int nlumi = l.id().luminosityBlock();

  // put dcsword into reportSummary 
  // FIXME consider making renderplugin "dynamic"
  // first move everything to the left by 1
  for (int i=2;i<16;i++) 
  {
    for (int j=1;j<11;j++) 
    {
      float cont = reportSummaryMap_->getBinContent(i,j);
      reportSummaryMap_->setBinContent(i-1,j,cont);
    }
  }

  // fill last bin 15 for detector HV
  for (int j=0;j<9;j++) 
  {
     float cont = 0.;
     if (dcsword_&(0x1<<j)) cont = 1.;
     cout << j << " " << (0x1<<j) << " " << cont << endl;
     reportSummaryMap_->setBinContent(15,9-j,cont);
  }
  // reset
  dcsword_=0xffff;

  // fill dcs vs lumi
  for (int i=0;i<24;i++)
  {
    if (dcs24[i])
      dcsVsLumi_->setBinContent(nlumi,i+1,1.);
    else
      dcsVsLumi_->setBinContent(nlumi,i+1,0.);
   
    dcs24[i]=true;
  }

  // fill physics decl. bit in y bin 10.
  if (physDecl_) 
    reportSummaryMap_->setBinContent(15,10,1.);
  else
    reportSummaryMap_->setBinContent(15,10,0.);
  // reset   
  physDecl_=true;  


  // set labels 
  char label[10];
  sprintf(label, "lumi %d", nlumi);
  reportSummaryMap_->setBinLabel( 15, label ,     1);
  if (nlumi>15) 
  {
    sprintf(label, "lumi %d", nlumi-14);
    reportSummaryMap_->setBinLabel(1, label , 1);
  }

  return;
  
}

// run showtag command line
std::string 
DQMProvInfo::getShowTags(void)
{
   TString out;
   FILE *pipe = gSystem->OpenPipe("showtags u -t", "r");

   TString line;
   while (line.Gets(pipe,true)) {
     if (line.Contains("Test Release")) continue;
     if (line.Contains("Base Release")) continue;
     if (line.Contains("Test release")) continue;
     if (line.Contains("--- Tag ---")) continue;
     if (line.Contains(" ")) line.Replace(line.First(" "),1,":");
     line.ReplaceAll(" ","");
     out = out + line + ";";
     if (line.Contains("-------------------")) break;
     if (out.Length()>2000) break;
   }
   out.ReplaceAll("--","");
   out.ReplaceAll(";-",";");
   out.ReplaceAll(";;",";");
   out.ReplaceAll("\n","");

   Int_t r = gSystem->ClosePipe(pipe);
   if (r) {
     gSystem->Error("ShowTags","problem running command showtags -u -t");
   }

   std::string str(out);
   if (str.length()>2000) str.resize(2000);

   std::string safestr =
     "/ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-;:";
   size_t found=str.find_first_not_of(safestr);
   if (found!=std::string::npos)
   {
     std::cout << "DQMFileSaver::ShowTags: Illegal character found: " 
               << str[found] 
               << " at position " 
               << int(found) << std::endl;
     return "notags";
   }   
   return str;
}


void 
DQMProvInfo::makeProvInfo()
{
    dbe_->cd() ;
    dbe_->setCurrentFolder( subsystemname_ + "/" +  provinfofolder_) ;

    // if (dbe_->get("ProvInfo/CMSSW")) return ;
    
    versCMSSW_     = dbe_->bookString("CMSSW",edm::getReleaseVersion().c_str() );
    hostName_      = dbe_->bookString("hostName",gSystem->HostName());
    workingDir_    = dbe_->bookString("workingDir",gSystem->pwd());
    processId_     = dbe_->bookInt("processID"); processId_->Fill(gSystem->GetPid());

    //versDataset_   = dbe_->bookString("Dataset",workflow_);
    versGlobaltag_ = dbe_->bookString("Globaltag","global tag"); // FIXME
    versTaglist_   = dbe_->bookString("Taglist",getShowTags()); 

    isComplete_ = dbe_->bookInt("runIsComplete"); 
    //isComplete_->Fill((runIsComplete_?1:0));
    fileVersion_ = dbe_->bookInt("fileVersion");
    //fileVersion_->Fill(version_);
    
    return ;
}
void 
DQMProvInfo::makeDcsInfo(const edm::Event& e)
{

  edm::Handle<DcsStatusCollection> dcsStatus;
  e.getByLabel("scalersRawToDigi", dcsStatus);
  for (DcsStatusCollection::const_iterator dcsStatusItr = dcsStatus->begin(); 
                                         dcsStatusItr != dcsStatus->end(); ++dcsStatusItr) 
  {
      // cout << "here " << dcsStatusItr->ready() << endl;
      if (!dcsStatusItr->ready(DcsStatus::CSCp) && 
          !dcsStatusItr->ready(DcsStatus::CSCm))   
     	                              dcsword_ = dcsword_ & 0xfffe;

      if (!dcsStatusItr->ready(DcsStatus::DT0) &&
          !dcsStatusItr->ready(DcsStatus::DTp) &&
	  !dcsStatusItr->ready(DcsStatus::DTm))   
	                              dcsword_ = dcsword_ & 0xfffd; 

      if (!dcsStatusItr->ready(DcsStatus::EBp) &&
          !dcsStatusItr->ready(DcsStatus::EBm))   
	                              dcsword_ = dcsword_ & 0xfffb;

      if (!dcsStatusItr->ready(DcsStatus::EEp) &&
          !dcsStatusItr->ready(DcsStatus::EEm))
	                              dcsword_ = dcsword_ & 0xfff7;

      if (!dcsStatusItr->ready(DcsStatus::ESp) &&
          !dcsStatusItr->ready(DcsStatus::ESm)) 
	                              dcsword_ = dcsword_ & 0xffef;

      if (!dcsStatusItr->ready(DcsStatus::HBHEa) &&
          !dcsStatusItr->ready(DcsStatus::HBHEb) &&
	  !dcsStatusItr->ready(DcsStatus::HBHEc) &&
	  !dcsStatusItr->ready(DcsStatus::HF) &&
	  !dcsStatusItr->ready(DcsStatus::HO)) 
	                              dcsword_ = dcsword_ & 0xffdf; 

      if (!dcsStatusItr->ready(DcsStatus::BPIX) &&
          !dcsStatusItr->ready(DcsStatus::FPIX))  
	                              dcsword_ = dcsword_ & 0xffbf;

      if (!dcsStatusItr->ready(DcsStatus::RPC))   
                                      dcsword_ = dcsword_ & 0xff7f;

      if (!dcsStatusItr->ready(DcsStatus::TIBTID) &&
          !dcsStatusItr->ready(DcsStatus::TOB) && 
	  !dcsStatusItr->ready(DcsStatus::TECp) &&
	  !dcsStatusItr->ready(DcsStatus::TECm))  
	                              dcsword_ = dcsword_ & 0xfeff;

      if (!dcsStatusItr->ready(DcsStatus::CSCp))   dcs24[0]=false;
      if (!dcsStatusItr->ready(DcsStatus::CSCm))   dcs24[1]=false;   
      if (!dcsStatusItr->ready(DcsStatus::DT0))    dcs24[2]=false;
      if (!dcsStatusItr->ready(DcsStatus::DTp))    dcs24[3]=false;
      if (!dcsStatusItr->ready(DcsStatus::DTm))    dcs24[4]=false;
      if (!dcsStatusItr->ready(DcsStatus::EBp))    dcs24[5]=false;
      if (!dcsStatusItr->ready(DcsStatus::EBm))    dcs24[6]=false;
      if (!dcsStatusItr->ready(DcsStatus::EEp))    dcs24[7]=false;
      if (!dcsStatusItr->ready(DcsStatus::EEm))    dcs24[8]=false;
      if (!dcsStatusItr->ready(DcsStatus::ESp))    dcs24[9]=false;
      if (!dcsStatusItr->ready(DcsStatus::ESm))    dcs24[10]=false; 
      if (!dcsStatusItr->ready(DcsStatus::HBHEa))  dcs24[11]=false;
      if (!dcsStatusItr->ready(DcsStatus::HBHEb))  dcs24[12]=false;
      if (!dcsStatusItr->ready(DcsStatus::HBHEc))  dcs24[13]=false; 
      if (!dcsStatusItr->ready(DcsStatus::HF))     dcs24[14]=false;
      if (!dcsStatusItr->ready(DcsStatus::HO))     dcs24[15]=false;
      if (!dcsStatusItr->ready(DcsStatus::BPIX))   dcs24[16]=false;
      if (!dcsStatusItr->ready(DcsStatus::FPIX))   dcs24[17]=false;
      if (!dcsStatusItr->ready(DcsStatus::RPC))    dcs24[18]=false;
      if (!dcsStatusItr->ready(DcsStatus::TIBTID)) dcs24[19]=false;
      if (!dcsStatusItr->ready(DcsStatus::TOB))    dcs24[20]=false;
      if (!dcsStatusItr->ready(DcsStatus::TECp))   dcs24[21]=false;
      if (!dcsStatusItr->ready(DcsStatus::TECm))   dcs24[22]=false;
      if (!dcsStatusItr->ready(DcsStatus::CASTOR)) dcs24[23]=false;

      // cout << hex << dcsword_ << endl;
  }
      
  return ;
}

void 
DQMProvInfo::makeGtInfo(const edm::Event& e)
{

  edm::Handle<L1GlobalTriggerReadoutRecord> gtrr_handle;
  e.getByLabel("gtDigis", gtrr_handle);
  L1GlobalTriggerReadoutRecord const* gtrr = gtrr_handle.product();
  L1GtFdlWord fdlWord ; 
  if (gtrr)
    fdlWord = gtrr->gtFdlWord();
  else
  {
    cout << "phys decl. bit not accessible !!!"  << endl;
    return;
  }
    
  // gtrr->print(std::cout);

  // cout << "phys decl. bit =" << static_cast<int>(fdlWord.physicsDeclared()) << endl;
  if (fdlWord.physicsDeclared() !=1) physDecl_=false;

  return;
}
