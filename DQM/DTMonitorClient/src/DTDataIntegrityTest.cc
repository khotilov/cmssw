
/*
 * \file DTDataIntegrityTest.cc
 * 
 * $Date: 2008/07/02 14:19:27 $
 * $Revision: 1.22 $
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
  
  LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") << "[DTDataIntegrityTest]: Constructor";

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

  LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") << "DataIntegrityTest: analyzed " << nupdates << " updates";

}


void DTDataIntegrityTest::beginJob(const EventSetup& context){

  LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") << "[DTDataIntegrityTest]: BeginJob";

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
}



void DTDataIntegrityTest::beginLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {

  LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]: Begin of LS transition";

  // Get the run number
  run = lumiSeg.run();

}



void DTDataIntegrityTest::analyze(const Event& e, const EventSetup& context){

  nevents++;
  LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") << "[DTDataIntegrityTest]: "<<nevents<<" events";

}



void DTDataIntegrityTest::endLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {

  //nSTAEvents++;
 // if running in standalone perform diagnostic only after a reasonalbe amount of events
  //if ( parameters.getUntrackedParameter<bool>("runningStandalone", false) && 
  //   nSTAEvents%parameters.getUntrackedParameter<int>("diagnosticPrescale", 1000) != 0 ) return;
 
  LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest")
    <<"[DTDataIntegrityTest]: End of LS transition, performing the DQM client operation";

  // counts number of lumiSegs 
  nLumiSegs = lumiSeg.id().luminosityBlock();
  stringstream nLumiSegs_s; nLumiSegs_s << nLumiSegs;

  // prescale factor
  if ( nLumiSegs%prescaleFactor != 0 ) return;


  // counts number of updats 
  nupdates++;
 
  LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]: "<<nupdates<<" updates";
  LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]: "<<nLumiSegs<<" luminosity block number";

  if(writeHisto && nupdates%nTimeBin == 0){
    LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]: saving all histos";
    stringstream runNumber; runNumber << run;
    stringstream lumiNumber; lumiNumber << nLumiSegs;
    string rootFile = outputFile + "_" + lumiNumber.str() + "_" + runNumber.str() + ".root";
    dbe->save(rootFile);
  }
  
  //Counter for x bin in the timing histos
  counter++;

  //Loop on FED id
  for (int dduId=FEDNumbering::getDTFEDIds().first; dduId<=FEDNumbering::getDTFEDIds().second; ++dduId){
    LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:FED Id: "<<dduId;
 
    //Each nTimeBin onUpdate remove timing histos and book a new bunch of them
    stringstream dduId_s; dduId_s << dduId;
    if(doTimeHisto && nupdates%nTimeBin == 1){
      LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]: booking a new bunch of time histos";
      //if(nupdates>nTimeBin)
      //dbe->rmdir("DT/Tests/DTDataIntegrity/FED" + dduId_s.str() + "/TimeInfo"); //FIXME: it doesn't work anymore
      //    (dduVectorHistos.find("TTSVSTime")->second).find(dduId)->second.clear();
      bookTimeHistos("TTSVSTime",dduId, nLumiSegs);
      bookTimeHistos("ROSVSTime",dduId, nLumiSegs);
      bookTimeHistos("EvLenghtVSTime",dduId,nLumiSegs);
      bookTimeHistos("FIFOVSTime",dduId,nLumiSegs);
    }

    string histoType;
    //1D histo: % of tts values 
    MonitorElement * tts_histo = dbe->get(getMEName("TTSValues",dduId));
    if (tts_histo) {
        LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:histo DDUTTSValues found";

	histoType = "TTSValues_Percent";   
	if (dduHistos[histoType].find(dduId) == dduHistos[histoType].end()) {
	  bookHistos(histoType,dduId);
	} 
	//Loop on possible tts values
	for(int i=1;i<8;i++){
	  (dduHistos.find(histoType)->second).find(dduId)->second->
	    setBinContent(i,tts_histo->getBinContent(i)/tts_histo->getEntries());

	  if(doTimeHisto){
	    //Fill timing histos and set x label with luminosity block number
	    if( dduVectorHistos["TTSVSTime"].find(dduId) == dduVectorHistos["TTSVSTime"].end() ){
	      bookTimeHistos("TTSVSTime",dduId,nLumiSegs); 
	    }
	    (dduVectorHistos.find("TTSVSTime")->second).find(dduId)->second[i-1]->
	      setBinContent(counter,tts_histo->getBinContent(i)/tts_histo->getEntries());
	    (dduVectorHistos.find("TTSVSTime")->second).find(dduId)->second[i-1]->
	      setBinLabel(counter, nLumiSegs_s.str(), 1);
	  }
	}

	//Check if there are too many events with wrong tts value
	double alert_tts1 = 0.5, alert_tts4 = 0.5, alert_tts20 = 0.5;
	if((dduHistos.find(histoType)->second).find(dduId)->second->getBinContent(2) > alert_tts1)
	  LogWarning ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:WARNING: "<<
	    (dduHistos.find(histoType)->second).find(dduId)->second->getBinContent(2)<<" % events with warning overflow";

   	if(((dduHistos.find(histoType)->second).find(dduId)->second->getBinContent(1) +
	    (dduHistos.find(histoType)->second).find(dduId)->second->getBinContent(3)) > alert_tts20 )
	  LogWarning ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:WARNING: "<<
	    ((dduHistos.find(histoType)->second).find(dduId)->second->getBinContent(1) +
	     (dduHistos.find(histoType)->second).find(dduId)->second->getBinContent(3))<<" % events with out of synch or disconnected";

	if((dduHistos.find(histoType)->second).find(dduId)->second->getBinContent(4) > alert_tts4)
	  LogWarning ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:WARNING: "<<
	    (dduHistos.find(histoType)->second).find(dduId)->second->getBinContent(4)<<" % events with busy";
	//FIXME: how to notify this warning in a LogFile?
         }

    //Check if the list of ROS is compatible with the channels enabled
    MonitorElement * ros_histo = dbe->get(getMEName("ROSStatus",dduId));
    if (ros_histo) {
        LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:histo DDUChannelStatus found";

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
      LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:histo DDUROSList found";

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
       LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:histo DDUEventLenght found";
       
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
       LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest]:histo DDUFIFOStatus found";
       
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
     string fedIntegrityFolder = "DT/FEDIntegrity/";
     MonitorElement * hFEDEntry = dbe->get(fedIntegrityFolder+"FEDEntries");
     MonitorElement * hFEDFatal = dbe->get(fedIntegrityFolder+"FEDFatal");
     MonitorElement * hFEDNonFatal = dbe->get(fedIntegrityFolder+"FEDNonFatal");

     if(FED_ROSSummary) {
       TH2F * histo_FEDSummary = FED_ROSSummary->getTH2F();
       for(int rosNumber = 1; rosNumber <= 12; ++rosNumber) { // loop on the ROS
	 int result = -2;
	 if(histo_FEDSummary->Integral(1,11,rosNumber,rosNumber) == 0) { // no errors
	   result = 0;
	 } else { // there are errors
	   result = 2;
	 }
	 summaryHisto->setBinContent(rosNumber,dduId-769,result);
       }
       // Check that the FED is in the ReadOut using the FEDIntegrity histos
       if(hFEDEntry->getBinContent(dduId-769) == 0 &&
	  hFEDFatal->getBinContent(dduId-769) == 0 &&
	  hFEDNonFatal->getBinContent(dduId-769) == 0) {
	 // no data in this FED: it is off
	 for(int rosNumber = 1; rosNumber <= 12; ++rosNumber) {
	   summaryHisto->setBinContent(rosNumber,dduId-769,1);
	 }
       }

     } else { // no data in this FED: it is off
       for(int rosNumber = 1; rosNumber <= 12; ++rosNumber) {
	 summaryHisto->setBinContent(rosNumber,dduId-769,1);
       }
     }
  
  }


     


  
}



void DTDataIntegrityTest::endJob(){

  LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"[DTDataIntegrityTest] endjob called!";

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

  if(histoType == "TTSValues_Percent"){
    histoName = "FED" + dduId_s.str() + histoType;
    (dduHistos[histoType])[dduId] = dbe->book1D(histoName,histoName,7,0,7);
    ((dduHistos[histoType])[dduId])->setBinLabel(1,"disconnected",1);	
    ((dduHistos[histoType])[dduId])->setBinLabel(2,"warning overflow",1);	
    ((dduHistos[histoType])[dduId])->setBinLabel(3,"out of synch",1);	
    ((dduHistos[histoType])[dduId])->setBinLabel(4,"busy",1);	
    ((dduHistos[histoType])[dduId])->setBinLabel(5,"ready",1);	
    ((dduHistos[histoType])[dduId])->setBinLabel(6,"error",1);	
    ((dduHistos[histoType])[dduId])->setBinLabel(7,"disconnected",1);	
  }
}

void DTDataIntegrityTest::bookTimeHistos(string histoType, int dduId, int nLumiSegs){
  stringstream dduId_s; dduId_s << dduId;
  stringstream nLumiSegs_s; nLumiSegs_s << nLumiSegs;
  string histoName;
  LogVerbatim ("DTDQM|DTRawToDigi|DTMonitorClient|DTDataIntegrityTest") <<"Booking time histo "<<histoType<<" for ddu "<<dduId<<" from luminosity block "<<nLumiSegs;

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

