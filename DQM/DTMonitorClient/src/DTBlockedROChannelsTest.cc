
/*
 * \file DTBlockedROChannelsTest.cc
 * 
 * $Date: 2009/10/19 16:18:47 $
 * $Revision: 1.6 $
 * \author G. Cerminara - University and INFN Torino
 *
 */

#include <DQM/DTMonitorClient/src/DTBlockedROChannelsTest.h>

//Framework
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "CondFormats/DTObjects/interface/DTReadOutMapping.h"
#include "CondFormats/DataRecord/interface/DTReadOutMappingRcd.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DQM/DTMonitorModule/interface/DTTimeEvolutionHisto.h"
#include <iostream>
#include <string>


using namespace std;
using namespace edm;


DTBlockedROChannelsTest::DTBlockedROChannelsTest(const ParameterSet& ps) : nevents(0),
									   neventsPrev(0),
									   prevNLumiSegs(0),
									   prevTotalPerc(0),
									   hSystFractionVsLS(0) {
  LogTrace("DTDQM|DTRawToDigi|DTMonitorClient|DTBlockedROChannelsTest")
    << "[DTBlockedROChannelsTest]: Constructor";

  // prescale on the # of LS to update the test
  prescaleFactor = ps.getUntrackedParameter<int>("diagnosticPrescale", 1);

  offlineMode = ps.getUntrackedParameter<bool>("offlineMode", true);
}



DTBlockedROChannelsTest::~DTBlockedROChannelsTest() {
  LogTrace("DTDQM|DTRawToDigi|DTMonitorClient|DTBlockedROChannelsTest")
    << "DataIntegrityTest: analyzed " << nupdates << " updates";
}



// book histos
void DTBlockedROChannelsTest::beginJob(const EventSetup& context) {
  LogTrace("DTDQM|DTRawToDigi|DTMonitorClient|DTBlockedROChannelsTest")
    << "[DTBlockedROChannelsTest]: BeginJob";

  //nSTAEvents = 0;
  nupdates = 0;
  run=0;

  dbe = Service<DQMStore>().operator->();
  
  // book the summary histogram
  dbe->setCurrentFolder("DT/00-ROChannels");
  summaryHisto = dbe->book2D("ROChannelSummary","Summary Blocked RO Channels",12,1,13,5,-2,3);
  summaryHisto->setAxisTitle("Sector",1);
  summaryHisto->setAxisTitle("Wheel",2);

  for(int wheel = -2; wheel != 3; ++wheel) {
    stringstream namestream;  namestream << "ROChannelSummary_W" << wheel;
    stringstream titlestream; titlestream << "Blocked RO Channels (Wh " << wheel << ")";
    wheelHitos[wheel] = dbe->book2D(namestream.str().c_str(),titlestream.str().c_str(),12,1,13,4,1,5);
    wheelHitos[wheel]->setAxisTitle("Sector",1);
    wheelHitos[wheel]->setBinLabel(1,"MB1",2);
    wheelHitos[wheel]->setBinLabel(2,"MB2",2);
    wheelHitos[wheel]->setBinLabel(3,"MB3",2);
    wheelHitos[wheel]->setBinLabel(4,"MB4",2);
  }

  if(!offlineMode) {
    hSystFractionVsLS = new DTTimeEvolutionHisto(dbe, "EnabledROChannelsVsLS", "% RO channels",
						 500, 5, true, 3);
  }

}




void DTBlockedROChannelsTest::beginRun(const Run& run, const EventSetup& context) {
  // get the RO mapping
  context.get<DTReadOutMappingRcd>().get(mapping);
  
  // fill the map of the robs per chamber
  for(int dduId = FEDNumbering::MINDTFEDID; dduId<=FEDNumbering::MAXDTFEDID; ++dduId) { //loop over DDUs
    for(int ros = 1; ros != 13; ++ros) { // loop over ROSs
      for(int rob = 1; rob != 26; ++rob) { // loop over ROBs	
	int wheel = 0;
	int station = 0;
	int sector = 0;
	int dummy = 0;
	if(!mapping->readOutToGeometry(dduId,ros,rob-1,0,2,wheel,station,sector,dummy,dummy,dummy) ||
	   !mapping->readOutToGeometry(dduId,ros,rob-1,0,16,wheel,station,sector,dummy,dummy,dummy)) {
	  DTChamberId chId(wheel, station, sector);
	  if(chamberMap.find(chId) == chamberMap.end()) {
	    chamberMap[chId] = DTRobBinsMap(dduId, ros, dbe);
	    chamberMap[chId].addRobBin(rob);
	  } 
	  chamberMap[chId].addRobBin(rob);
	} else {
	   LogTrace("DTDQM|DTRawToDigi|DTMonitorClient|DTBlockedROChannelsTest")
	     << "[DTRobBinsMap] FED: " << dduId << " ROS " << ros << " ROB: " << rob-1
	     << " not in the mapping!" << endl;
	}
      }
    }
  }
}



void DTBlockedROChannelsTest::beginLuminosityBlock(LuminosityBlock const& lumiSeg,
						   EventSetup const& context) {
  LogTrace("DTDQM|DTRawToDigi|DTMonitorClient|DTBlockedROChannelsTest")
    <<"[DTBlockedROChannelsTest]: Begin of LS transition";

  // Get the run number
  run = lumiSeg.run();
  nevents = 0;

//   // loop over all chambers and read the values at the beginning of the LS
//   for(map<DTChamberId, DTRobBinsMap>::iterator chAndRobs = chamberMap.begin();
//       chAndRobs != chamberMap.end(); ++chAndRobs) {
//     (*chAndRobs).second.readNewValues(); 
//   }

}



void DTBlockedROChannelsTest::analyze(const Event& e, const EventSetup& context){
  // count the analyzed events
  nevents++;
  if(nevents%1000 == 0)
    LogTrace("DTDQM|DTRawToDigi|DTMonitorClient|DTBlockedROChannelsTest")
      << "[DTBlockedROChannelsTest]: "<<nevents<<" events";
}



void DTBlockedROChannelsTest::endLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {

  // counts number of lumiSegs 
  nLumiSegs = lumiSeg.id().luminosityBlock();
  
  // prescale factor
  if (nLumiSegs%prescaleFactor != 0) return;
  
  LogTrace("DTDQM|DTRawToDigi|DTMonitorClient|DTBlockedROChannelsTest")
    <<"[DTBlockedROChannelsTest]: End of LS " << nLumiSegs << ", performing client operations";


  // counts number of updats 
  nupdates++;

  // skip empty LSs
  

  if(nevents == 0) { // hack to work also in offline DQM
    MonitorElement *procEvt =  dbe->get("DT/EventInfo/processedEvents");
    if(procEvt != 0) {
      int procEvents = procEvt->getIntValue();
      nevents = procEvents - neventsPrev;
      neventsPrev = procEvents;
    }
  }

  double totalPerc = prevTotalPerc;
  // check again!
  if(nevents != 0) { // skip the computration if no events in the last LS

    // reset the histos
    summaryHisto->Reset();
    for(int wheel = -2; wheel != 3; ++wheel) {
      wheelHitos[wheel]->Reset();
    }
  
    totalPerc = 0.;

    // loop over all chambers and fill the wheel plots
    for(map<DTChamberId, DTRobBinsMap>::iterator chAndRobs = chamberMap.begin();
	chAndRobs != chamberMap.end(); ++chAndRobs) {
      DTChamberId chId = (*chAndRobs).first;
      double scale = 1.;
      int sectorForPlot = chId.sector();
      if(sectorForPlot == 13 || (sectorForPlot == 4 && chId.station() ==4)) {
	sectorForPlot = 4;
	scale = 0.5;
      } else if(sectorForPlot == 14 || (sectorForPlot == 10 && chId.station() ==4)) {
	sectorForPlot = 10;
	scale = 0.5;
      }

      // NOTE: can be called only ONCE per event per each chamber
      double chPercent = (*chAndRobs).second.getChamberPercentage(); 
      wheelHitos[chId.wheel()]->Fill(sectorForPlot, chId.station(),
				     scale*chPercent);
      totalPerc += chPercent*scale*1./250.;
//       if(chPercent != 1.) {
// 	cout << "Ch: " << (*chAndRobs).first << endl;
// 	cout << "      perc: " << chPercent << endl;
//       }
      // Fill the summary
      summaryHisto->Fill(sectorForPlot, chId.wheel(), 0.25*scale*chPercent);
    }
  }

  
  // this part is executed even if no events were processed in order to include the last LS 
  if(offlineMode) { // save the results in a map and draw them in the end-run
    if(resultsPerLumi.size() == 0) { // the first 2 LS are analyzed together
//       cout << "LS: " << nLumiSegs << " total %: " << totalPerc << endl;
      resultsPerLumi[nLumiSegs] = totalPerc;
    } else {
//       cout << "LS: " << nLumiSegs << " total %: " << prevTotalPerc << endl;
      resultsPerLumi[nLumiSegs] = prevTotalPerc;
    }
    prevTotalPerc = totalPerc;
    prevNLumiSegs = nLumiSegs;

  } else { // directly fill the histo
    hSystFractionVsLS->accumulateValueTimeSlot(totalPerc);
    hSystFractionVsLS->updateTimeSlot(nLumiSegs, nevents);
    prevTotalPerc = totalPerc;
  }
}



void DTBlockedROChannelsTest::endJob(){
  LogTrace("DTDQM|DTRawToDigi|DTMonitorClient|DTBlockedROChannelsTest")
    <<"[DTBlockedROChannelsTest] endjob called!";
}




void DTBlockedROChannelsTest::endRun(edm::Run const& run, edm::EventSetup const& eSetup) {

//     Commeneted to fix almost infinite loop on MC
//   if(offlineMode) {
//     // fill a trend plot based on the results stored in the map
//     float lBin = resultsPerLumi.rbegin()->first;
//     dbe->setCurrentFolder("DT/00-ROChannels");
  
//     //   MonitorElement* hSystFractionVsLS =  dbe->book1D("EnabledROChannelsVsLS", "% RO channels vs LS", nBins,fBin,lBin);

//     hSystFractionVsLS = new DTTimeEvolutionHisto(dbe, "EnabledROChannelsVsLS", "% RO channels",
// 						 (int)lBin, 1, false, 2);
								 
//     for(map<int, double>::const_iterator bin = resultsPerLumi.begin();
// 	bin != resultsPerLumi.end(); ++bin) {
//       hSystFractionVsLS->setTimeSlotValue((*bin).second, (*bin).first);
//     }
//   }
}



int DTBlockedROChannelsTest::readOutToGeometry(int dduId, int ros, int& wheel, int& sector){

  int dummy;
  return mapping->readOutToGeometry(dduId,ros,2,2,2,wheel,dummy,sector,dummy,dummy,dummy);

}





DTBlockedROChannelsTest::DTRobBinsMap::DTRobBinsMap(const int fed, const int ros, const DQMStore* dbe) {
  // get the pointer to the corresondig histo
  stringstream mename; mename << "DT/00-DataIntegrity/FED" << fed << "/ROS" << ros
			      << "/FED" << fed << "_ROS" << ros << "_ROSError";
  hName = mename.str();
  meROS = dbe->get(mename.str());
  theDbe = dbe;
}




DTBlockedROChannelsTest::DTRobBinsMap::DTRobBinsMap() {
  meROS = 0;
}



DTBlockedROChannelsTest::DTRobBinsMap::~DTRobBinsMap() {}



// add a rob to the set of robs
void DTBlockedROChannelsTest::DTRobBinsMap::addRobBin(int robBin) {
  robsAndValues[robBin] = getValueRobBin(robBin);
}



    
int DTBlockedROChannelsTest::DTRobBinsMap::getValueRobBin(int robBin) const {
  int value = 0;
  if(meROS) {
    value += (int)meROS->getBinContent(9,robBin);
    value += (int)meROS->getBinContent(11,robBin);
  }
  return value;
}




bool DTBlockedROChannelsTest::DTRobBinsMap::robChanged(int robBin) {
  // check that this is a valid ROB for this map (= it has been added!)
  if(robsAndValues.find(robBin) == robsAndValues.end()) {
     LogWarning("DTDQM|DTRawToDigi|DTMonitorClient|DTBlockedROChannelsTest")
       << "[DTRobBinsMap]***Error: ROB: " << robBin << " is not valid" << endl;
    return false;
  }

  int newValue = getValueRobBin(robBin);
  if(newValue > robsAndValues[robBin]) {
    robsAndValues[robBin] = newValue;
    return true;
  }
  return false;
}




double DTBlockedROChannelsTest::DTRobBinsMap::getChamberPercentage() {
  meROS = theDbe->get(hName);
  int nChangedROBs = 0;
  for(map<int, int>::const_iterator robAndValue = robsAndValues.begin();
      robAndValue != robsAndValues.end(); ++robAndValue) {
    if(robChanged((*robAndValue).first)) nChangedROBs++;
  }
  return 1.-((double)nChangedROBs/(double)robsAndValues.size());
}


void DTBlockedROChannelsTest::DTRobBinsMap::readNewValues() {
  meROS = theDbe->get(hName);
  for(map<int, int>::const_iterator robAndValue = robsAndValues.begin();
      robAndValue != robsAndValues.end(); ++robAndValue) {
    robChanged((*robAndValue).first);
  }
  
}




// FIXME: move to SealModule
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_ANOTHER_FWK_MODULE(DTBlockedROChannelsTest);
