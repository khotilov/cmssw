#include "DQMOffline/CalibMuon/interface/DTtTrigDBValidation.h"

// Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTTopology.h"

// DataFormats
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"

// tTrig
#include "CondFormats/DTObjects/interface/DTTtrig.h"
#include "CondFormats/DataRecord/interface/DTTtrigRcd.h"

#include <stdio.h>
#include <sstream>
#include <math.h>
#include "TFile.h"

using namespace edm;
using namespace std;

DTtTrigDBValidation::DTtTrigDBValidation(const ParameterSet& pset) {

  cout << "[DTtTrigDBValidation] Constructor called!" << endl;

  // Get the DQM needed services
  dbe = edm::Service<DQMStore>().operator->();
  dbe->setCurrentFolder("DT/DTDBValidation");

  // Get dataBase label
  labelDBRef = pset.getUntrackedParameter<string>("labelDBRef");
  labelDB = pset.getUntrackedParameter<string>("labelDB");

  parameters = pset;
}


DTtTrigDBValidation::~DTtTrigDBValidation(){}


void DTtTrigDBValidation::beginJob(const EventSetup& setup) {


  metname = "tTrigdbValidation";
  LogTrace(metname)<<"[DTtTrigDBValidation] Parameters initialization";
 
  outputFileName = parameters.getUntrackedParameter<std::string>("OutputFileName");

  ESHandle<DTTtrig> tTrig_Ref;
  setup.get<DTTtrigRcd>().get(labelDBRef, tTrig_Ref);
  const DTTtrig* DTTtrigRefMap = &*tTrig_Ref;
  LogTrace(metname)<<"[DTtTrigDBValidation] reference Ttrig version: " << tTrig_Ref->version();

  ESHandle<DTTtrig> tTrig;
  setup.get<DTTtrigRcd>().get(labelDB, tTrig);
  const DTTtrig* DTTtrigMap = &*tTrig;
  LogTrace(metname)<<"[DTtTrigDBValidation] Ttrig to validate version: " << tTrig->version();

  // Get the geometry
  setup.get<MuonGeometryRecord>().get(dtGeom);

  // Loop over Ref DB entries
  for(DTTtrig::const_iterator it = DTTtrigRefMap->begin();
      it != DTTtrigRefMap->end(); ++it) {
    DTSuperLayerId slId((*it).first.wheelId,
		        (*it).first.stationId,
		        (*it).first.sectorId,
		        (*it).first.slId);
    float tTrigmean = (*it).second.tTrig;
    float tTrigrms = (*it).second.tTrms;
    LogTrace(metname)<< "Ref Superlayer: " <<  slId <<endl
		     << " Ttrig mean (TDC counts): " << tTrigmean
		     << " Ttrig rms (TDC counts): " << tTrigrms;

    //t0RefMap[wireId].push_back(t0mean);
    //t0RefMap[wireId].push_back(t0rms);
    tTrigRefMap[slId] = std::make_pair<float,float>(tTrigmean,tTrigrms);
  }

  // Loop over Ref DB entries
  for(DTTtrig::const_iterator it = DTTtrigMap->begin();
      it != DTTtrigMap->end(); ++it) {
    DTSuperLayerId slId((*it).first.wheelId,
		        (*it).first.stationId,
		        (*it).first.sectorId,
		        (*it).first.slId);
    float tTrigmean = (*it).second.tTrig;
    float tTrigrms = (*it).second.tTrms;
    LogTrace(metname)<< "SuperLayer: " <<  slId <<endl
		     << " Ttrig mean (TDC counts): " << tTrigmean
		     << " Ttrig rms (TDC counts): " << tTrigrms;

    //t0Map[wireId].push_back(t0mean);
    //t0Map[wireId].push_back(t0rms);
    tTrigMap[slId] = std::make_pair<float,float>(tTrigmean,tTrigrms);
  }

  double difference=0;
  for(map<DTSuperLayerId, pair<float,float> >::const_iterator it = tTrigRefMap.begin();
      it != tTrigRefMap.end();
      ++it) {  
      if(tTrigMap.find((*it).first) != tTrigMap.end()) {

      // compute the difference
      difference = tTrigMap[(*it).first].first - (*it).second.first;

      //book histo
      /*DTLayerId layerId = (*theMap).first.layerId();
      if(t0DiffHistos.find(layerId) == t0DiffHistos.end()) {
	const DTTopology& dtTopo = dtGeom->layer(layerId)->specificTopology();
	const int firstWire = dtTopo.firstChannel();
	const int lastWire = dtTopo.lastChannel();
	bookHistos(layerId, firstWire, lastWire);
      }*/

      int wheel = (*it).first.chamberId().wheel();
      int sector = (*it).first.chamberId().sector();	
      if(tTrigDiffHistos.find(make_pair(wheel,sector)) == tTrigDiffHistos.end()) bookHistos(wheel,sector);
			
      cout<< "Filling the histo for super-layer: "<<(*it).first
	  <<"  difference: "<<difference<<endl;
      //t0DiffHistos[layerId]->Fill((*theMap).first.wire(),difference);

      // Fill the test histos
      int entry=-1;
      int station = (*it).first.chamberId().station();	
      if(station == 1) entry=0;
      if(station == 2) entry=3;
      if(station == 3) entry=6;
      if(station == 4) entry=9;

      int BinNumber = entry + (*it).first.superLayer();
      if(BinNumber == 12) BinNumber=11;	
	
      tTrigDiffHistos[make_pair(wheel,sector)]->setBinContent(BinNumber, difference);	

    }
  } // Loop over the tTrig map reference
  
  
}


void DTtTrigDBValidation::endJob() {

  //check the histos
  string testCriterionName = parameters.getUntrackedParameter<string>("tTrigTestName","tTrigDifferenceInRange"); 
  for(map<pair<int,int>, MonitorElement*>::const_iterator hDiff = tTrigDiffHistos.begin();
      hDiff != tTrigDiffHistos.end();
      hDiff++) {
      const QReport * theDiffQReport = (*hDiff).second->getQReport(testCriterionName);
      if(theDiffQReport) {
        vector<dqm::me_util::Channel> badChannels = theDiffQReport->getBadChannels();
        for (vector<dqm::me_util::Channel>::iterator channel = badChannels.begin(); 
	    channel != badChannels.end(); channel++) {
	  //cout << "layer:"<<(*hDiff).first<<" Bad mean channels: "<<(*channel).getBin()<<"  Contents : "<<(*channel).getContents()<<endl;
	  cout << "Bad mean channel: wh: " << (*hDiff).first.first
				           << " st: " << stationFromBin((*channel).getBin())
				           << " sect: " << (*hDiff).first.second
				           << " sl: " << slFromBin((*channel).getBin())
				           << " mean (cm): " << (*channel).getContents() << endl;
        }
        cout << "-------- Wheel, Sector: "<< (*hDiff).first.first << ", " << (*hDiff).first.second << "  " << theDiffQReport->getMessage() << " ------- " << theDiffQReport->getStatus() << endl; 
	
    }
  }

  // write the histos on a file
  dbe->save(outputFileName);

}

// Book a set of histograms for a given Layer
/*void DTtTrigDBValidation::bookHistos(DTLayerId lId, int firstWire, int lastWire) {
  
  LogTrace(metname)<< "   Booking histos for L: " << lId;

  // Compose the chamber name
  stringstream wheel; wheel << lId.superlayerId().chamberId().wheel();	
  stringstream station; station << lId.superlayerId().chamberId().station();	
  stringstream sector; sector << lId.superlayerId().chamberId().sector();	
  stringstream superLayer; superLayer << lId.superlayerId().superlayer();	
  stringstream layer; layer << lId.layer();

  string lHistoName =
    "_W" + wheel.str() +
    "_St" + station.str() +
    "_Sec" + sector.str() +
    "_SL" + superLayer.str()+
    "_L" + layer.str();
  
  dbe->setCurrentFolder("DT/t0Validation/Wheel" + wheel.str() +
			   "/Station" + station.str() +
			   "/Sector" + sector.str() +
			   "/SuperLayer" +superLayer.str());
  // Create the monitor elements
  MonitorElement * hDifference;
  hDifference = dbe->book1D("hDifference"+lHistoName, "difference between the two t0 values",lastWire-firstWire+1, firstWire-0.5, lastWire+0.5);
  
  t0DiffHistos[lId] = hDifference;
}*/

void DTtTrigDBValidation::bookHistos(int wheel, int sector) {

  LogTrace(metname)<< "   Booking histos for Wheel, Sector: " << wheel << ", " << sector;

  // Compose the chamber name
  stringstream str_wheel; str_wheel << wheel;
  stringstream str_sector; str_sector << sector;

  string lHistoName = "_W" + str_wheel.str() + "_Sec" + str_sector.str();

  dbe->setCurrentFolder("DT/tTrigValidation/Wheel" + str_wheel.str());

  // Create the monitor elements
  MonitorElement * hDifference;
  hDifference = dbe->book1D("htTrigDifference"+lHistoName, "difference between the two tTrig values",11,0,11);

  pair<int,int> mypair(wheel,sector);
  tTrigDiffHistos[mypair] = hDifference;

  (tTrigDiffHistos[mypair])->setBinLabel(1,"MB1_SL1",1);
  (tTrigDiffHistos[mypair])->setBinLabel(2,"MB1_SL2",1);
  (tTrigDiffHistos[mypair])->setBinLabel(3,"MB1_SL3",1);
  (tTrigDiffHistos[mypair])->setBinLabel(4,"MB2_SL1",1);
  (tTrigDiffHistos[mypair])->setBinLabel(5,"MB2_SL2",1);
  (tTrigDiffHistos[mypair])->setBinLabel(6,"MB2_SL3",1);
  (tTrigDiffHistos[mypair])->setBinLabel(7,"MB3_SL1",1);
  (tTrigDiffHistos[mypair])->setBinLabel(8,"MB3_SL2",1);
  (tTrigDiffHistos[mypair])->setBinLabel(9,"MB3_SL3",1);
  (tTrigDiffHistos[mypair])->setBinLabel(10,"MB4_SL1",1);
  (tTrigDiffHistos[mypair])->setBinLabel(11,"MB4_SL3",1);
}

int DTtTrigDBValidation::stationFromBin(int bin) const {
  return (int) (bin /3.1)+1;
}
 
int DTtTrigDBValidation::slFromBin(int bin) const {
  int ret = bin%3;
  if(ret == 0 || bin == 11) ret = 3;
  
  return ret;
}
