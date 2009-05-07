
/*
 * \file DTDataIntegrityTest.cc
 * 
 * $Date: 2009/03/27 13:22:21 $
 * $Revision: 1.27 $
 * \author S. Bolognesi - CERN
 *
 */

#include <DQM/DTMonitorClient/src/DTDataIntegrityTest.h>

//Framework
#include <DataFormats/FEDRawData/interface/FEDNumbering.h>
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <string>


using namespace std;
using namespace edm;


DTDataIntegrityTest::DTDataIntegrityTest(const ParameterSet& ps) : nevents(0) {
  
  LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") << "[DTDataIntegrityTest]: Constructor";

  //Number of bin in time histo
  nTimeBin =  ps.getUntrackedParameter<int>("nTimeBin", 10);
 //If you want info VS time histos
  doTimeHisto =  ps.getUntrackedParameter<bool>("doTimeHisto", false);
  // switch to write histos to file
  writeHisto = ps.getUntrackedParameter<bool>("writeHisto", false);
  // prefix of the name of the root file (lumi# and run# will be appended)
  outputFile = ps.getUntrackedParameter<string>("outputFile", "DTDataIntegrityTest");
  // prescale on the # of LS to update the test
  prescaleFactor = ps.getUntrackedParameter<int>("diagnosticPrescale", 1);


}


DTDataIntegrityTest::~DTDataIntegrityTest(){

  LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") << "DataIntegrityTest: analyzed " << nupdates << " updates";

}


void DTDataIntegrityTest::beginJob(const EventSetup& context){

  LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") << "[DTDataIntegrityTest]: BeginJob";

  //nSTAEvents = 0;
  nupdates = 0;
  run=0;

  dbe = Service<DQMStore>().operator->();
  
  // book the summary histogram
  dbe->setCurrentFolder("DT/00-DataIntegrity");
  summaryHisto = dbe->book2D("DataIntegritySummary","Summary Data Integrity",12,1,13,5,770,775);
  summaryHisto->setAxisTitle("ROS",1);
  summaryHisto->setBinLabel(1,"FED770",2);
  summaryHisto->setBinLabel(2,"FED771",2);
  summaryHisto->setBinLabel(3,"FED772",2);
  summaryHisto->setBinLabel(4,"FED773",2);
  summaryHisto->setBinLabel(5,"FED774",2);

  dbe->setCurrentFolder("DT/00-DataIntegrity");
  glbSummaryHisto = dbe->book2D("DataIntegrityGlbSummary","Summary Data Integrity",12,1,13,5,770,775);
  glbSummaryHisto->setAxisTitle("ROS",1);
  glbSummaryHisto->setBinLabel(1,"FED770",2);
  glbSummaryHisto->setBinLabel(2,"FED771",2);
  glbSummaryHisto->setBinLabel(3,"FED772",2);
  glbSummaryHisto->setBinLabel(4,"FED773",2);
  glbSummaryHisto->setBinLabel(5,"FED774",2);

}



void DTDataIntegrityTest::beginLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {

  LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]: Begin of LS transition";

  // Get the run number
  run = lumiSeg.run();

}



void DTDataIntegrityTest::analyze(const Event& e, const EventSetup& context){
  // count the analyzed events
  nevents++;
  if(nevents%1000 == 0)
    LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest")
      << "[DTDataIntegrityTest]: "<<nevents<<" events";
}



void DTDataIntegrityTest::endLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {

  // counts number of lumiSegs 
  nLumiSegs = lumiSeg.id().luminosityBlock();
  stringstream nLumiSegs_s; nLumiSegs_s << nLumiSegs;
  
  // prescale factor
  if (nLumiSegs%prescaleFactor != 0) return;
  
  LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest")
    <<"[DTDataIntegrityTest]: End of LS " << nLumiSegs << ", performing client operations";


  // counts number of updats 
  nupdates++;

  if(writeHisto && nupdates%nTimeBin == 0){
    LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]: saving all histos";
    stringstream runNumber; runNumber << run;
    stringstream lumiNumber; lumiNumber << nLumiSegs;
    string rootFile = outputFile + "_" + lumiNumber.str() + "_" + runNumber.str() + ".root";
    dbe->save(rootFile);
  }
  
  //Counter for x bin in the timing histos
  counter++;

  //Loop on FED id
  for (int dduId=FEDNumbering::MINDTFEDID; dduId<FEDNumbering::MAXDTFEDID; ++dduId){
    LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest")
      <<"[DTDataIntegrityTest]:FED Id: "<<dduId;
 
    //Each nTimeBin onUpdate remove timing histos and book a new bunch of them
    stringstream dduId_s; dduId_s << dduId;
    if(doTimeHisto && nupdates%nTimeBin == 1) {
      LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]: booking a new bunch of time histos";
      //if(nupdates>nTimeBin)
      //dbe->rmdir("DT/Tests/DTDataIntegrity/FED" + dduId_s.str() + "/TimeInfo"); //FIXME: it doesn't work anymore
      //    (dduVectorHistos.find("TTSVSTime")->second).find(dduId)->second.clear();
      bookTimeHistos("TTSVSTime",dduId, nLumiSegs);
      bookTimeHistos("ROSVSTime",dduId, nLumiSegs);
      bookTimeHistos("EvLenghtVSTime",dduId,nLumiSegs);
      bookTimeHistos("FIFOVSTime",dduId,nLumiSegs);
    }

    string histoType;

//     //1D histo: % of tts values 
//     MonitorElement * tts_histo = dbe->get(getMEName("TTSValues",dduId));
//     if (tts_histo) {
//       // FIXME: test?
//     }

    //Check if the list of ROS is compatible with the channels enabled
    MonitorElement * ros_histo = dbe->get(getMEName("ROSStatus",dduId));
    if (ros_histo) {
        LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:histo DDUChannelStatus found";

 	for(int i=1;i<13;i++){
	  if(ros_histo->getBinContent(1,i) != ros_histo->getBinContent(9,i))
	    LogError ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:WARNING: ROS"<<i<<" in "
					    <<ros_histo->getBinContent(9,i)<<" events"<<endl
					    <<"               but channel"<<i<<" enabled in "
					    <<ros_histo->getBinContent(1,i)<<" events";
	  //FIXME: how to notify this warning in a LogFile?
	}
    }
    //Monitor the number of ROS VS time
     MonitorElement * rosNumber_histo = dbe->get(getMEName("ROSList",dduId));
    if (rosNumber_histo && doTimeHisto) {
      LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:histo DDUROSList found";

      double rosNumber_mean = rosNumber_histo->getMean();
      //Fill timing histos and set x label with luminosity block number
      histoType = "ROSVSTime";
      if (dduHistos[histoType].find(dduId) == dduHistos[histoType].end()) {
	bookTimeHistos(histoType,dduId,nLumiSegs);
      }
      (dduHistos.find(histoType)->second).find(dduId)->second->setBinContent(counter,rosNumber_mean);
      (dduHistos.find(histoType)->second).find(dduId)->second->setBinLabel(counter, nLumiSegs_s.str(), 1);
    }
    
    //Monitor the event lenght VS time
     MonitorElement * evLenght_histo = dbe->get(getMEName("EventLenght",dduId));
     if (evLenght_histo && doTimeHisto) {
       LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:histo DDUEventLenght found";
       
       double evLenght_mean = evLenght_histo->getMean();
       //Fill timing histos and set x label with luminosity block number
       histoType = "EvLenghtVSTime";
       if (dduHistos[histoType].find(dduId) == dduHistos[histoType].end()) {
	 bookTimeHistos(histoType,dduId,nLumiSegs);
       }
       (dduHistos.find(histoType)->second).find(dduId)->second->setBinContent(counter,evLenght_mean);
       (dduHistos.find(histoType)->second).find(dduId)->second->setBinLabel(counter, nLumiSegs_s.str(), 1);
       
     }
     
     //Monitor the FIFO occupancy VS time 
     MonitorElement * fifo_histo = dbe->get(getMEName("FIFOStatus",dduId));
     if (fifo_histo && doTimeHisto) {
       LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:histo DDUFIFOStatus found";
       
       //Fill timing histos and set x label with luminosity block number
       histoType = "FIFOVSTime";
       if (dduVectorHistos[histoType].find(dduId) == dduVectorHistos[histoType].end()) {
	 bookTimeHistos(histoType,dduId,nLumiSegs);
       }
       for(int i=1;i<8;i++){
	 (dduVectorHistos.find("FIFOVSTime")->second).find(dduId)->second[i-1]->
	   setBinContent(counter,(fifo_histo->getBinContent(i,1) + 2*(fifo_histo->getBinContent(i,2)))/fifo_histo->getEntries());
	 (dduVectorHistos.find("FIFOVSTime")->second).find(dduId)->second[i-1]->
	   setBinLabel(counter, nLumiSegs_s.str(), 1);
       }
     }

     // Fill the summary histo   
     // Get the error summary histo
     string wheelSummaryName = "DT/00-DataIntegrity/FED" + dduId_s.str() + "_ROSSummary";
     MonitorElement * FED_ROSSummary = dbe->get(wheelSummaryName);
     // Get the histos for FED integrity
     string fedIntegrityFolder = "DT/FEDIntegrity_SM/";
     MonitorElement * hFEDEntry = dbe->get(fedIntegrityFolder+"FEDEntries");
     MonitorElement * hFEDFatal = dbe->get(fedIntegrityFolder+"FEDFatal");
     MonitorElement * hFEDNonFatal = dbe->get(fedIntegrityFolder+"FEDNonFatal");

     if(FED_ROSSummary) {
       TH2F * histo_FEDSummary = FED_ROSSummary->getTH2F();
       for(int rosNumber = 1; rosNumber <= 12; ++rosNumber) { // loop on the ROS
	 int result = -2;
	 float nErrors = histo_FEDSummary->Integral(1,11,rosNumber,rosNumber);
	 if(nErrors == 0) { // no errors
	   result = 0;
	 } else { // there are errors
	   result = 2;
	 }
	 summaryHisto->setBinContent(rosNumber,dduId-769,result);
	 // FIXME: different errors should have different weights
	 float sectPerc = max((float)0., ((float)nevents-nErrors)/(float)nevents);
	 glbSummaryHisto->setBinContent(rosNumber,dduId-769,sectPerc);
       }
       // Check that the FED is in the ReadOut using the FEDIntegrity histos
       if(hFEDEntry->getBinContent(dduId-769) == 0 &&
	  hFEDFatal->getBinContent(dduId-769) == 0 &&
	  hFEDNonFatal->getBinContent(dduId-769) == 0) {
	 // no data in this FED: it is off
	 for(int rosNumber = 1; rosNumber <= 12; ++rosNumber) {
	   summaryHisto->setBinContent(rosNumber,dduId-769,1);
	   glbSummaryHisto->setBinContent(rosNumber,dduId-769,0);
	 }
       }

     } else { // no data in this FED: it is off
       for(int rosNumber = 1; rosNumber <= 12; ++rosNumber) {
	 summaryHisto->setBinContent(rosNumber,dduId-769,1);
	 glbSummaryHisto->setBinContent(rosNumber,dduId-769,0);
       }
     }
  
  }


     


  
}



void DTDataIntegrityTest::endJob(){

  LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest] endjob called!";

//   dbe->rmdir("DT/DTDataIntegrity");
}



string DTDataIntegrityTest::getMEName(string histoType, int FEDId){
  //Use the DDU name to find the ME
  stringstream dduID_s; dduID_s << FEDId;

  string folderName = "DT/00-DataIntegrity/FED" + dduID_s.str(); 

  string histoName = folderName + "/FED" + dduID_s.str() + "_" + histoType;
  return histoName;
}



void DTDataIntegrityTest::bookHistos(string histoType, int dduId){
  stringstream dduId_s; dduId_s << dduId;
  dbe->setCurrentFolder("DT/00-DataIntegrity/FED" + dduId_s.str());
  string histoName;

}

void DTDataIntegrityTest::bookTimeHistos(string histoType, int dduId, int nLumiSegs){
  stringstream dduId_s; dduId_s << dduId;
  stringstream nLumiSegs_s; nLumiSegs_s << nLumiSegs;
  string histoName;
  LogTrace ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"Booking time histo "<<histoType<<" for ddu "<<dduId<<" from luminosity block "<<nLumiSegs;

  //Counter for x bin in the timing histos
  counter = 1;//assuming synchronized booking for all histo VS time

  if(histoType == "TTSVSTime"){
    dbe->setCurrentFolder("DT/00-DataIntegrity/FED" + dduId_s.str()+ "/TimeInfo/TTSVSTime");
    histoName = "FED" + dduId_s.str() + "_" + histoType + "_disconn_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + histoType + "_overflow_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + histoType + "_outSynch_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + histoType + "_busy_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + histoType + "_ready_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + histoType + "_error_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + histoType + "_disconnect_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
  }
  else if(histoType == "ROSVSTime"){
    dbe->setCurrentFolder("DT/00-DataIntegrity/FED" + dduId_s.str()+ "/TimeInfo/ROSVSTime");
    histoName = "FED" + dduId_s.str() + "_" + histoType + "_LumBlock" + nLumiSegs_s.str();
    (dduHistos[histoType])[dduId] = dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin);
  }
  else if(histoType == "EvLenghtVSTime"){
    dbe->setCurrentFolder("DT/00-DataIntegrity/FED" + dduId_s.str()+ "/TimeInfo/EvLenghtVSTime");
    histoName = "FED" + dduId_s.str() + "_" + histoType + "_LumBlock" +  nLumiSegs_s.str();
    (dduHistos[histoType])[dduId] = dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin);
  }
  else if(histoType == "FIFOVSTime"){
    dbe->setCurrentFolder("DT/00-DataIntegrity/FED" + dduId_s.str()+ "/TimeInfo/FIFOVSTime");
    histoName = "FED" + dduId_s.str() + "_" + histoType + "_Input1_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + "_" + histoType + "_Input2_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + "_" + histoType + "_Input3_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + "_" + histoType + "_L1A1_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + "_" + histoType + "_L1A2_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + "_" + histoType + "_L1A3_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
    histoName = "FED" + dduId_s.str() + "_" + histoType + "_Output_LumBlock" + nLumiSegs_s.str();
    ((dduVectorHistos[histoType])[dduId]).push_back(dbe->book1D(histoName,histoName,nTimeBin,nLumiSegs,nLumiSegs+nTimeBin));
  }
} 

