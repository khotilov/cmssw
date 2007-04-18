
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2007/04/12 07:40:04 $
 *  $Revision: 1.3 $
 *  \author G. Mila - INFN Torino
 */


#include <DQM/DTMonitorClient/src/DTDeadChannelTest.h>

// Framework
#include <FWCore/Framework/interface/Event.h>
#include "DataFormats/Common/interface/Handle.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <DQMServices/Core/interface/MonitorElementBaseT.h>

// Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTTopology.h"

#include "CondFormats/DataRecord/interface/DTStatusFlagRcd.h"
#include "CondFormats/DTObjects/interface/DTStatusFlag.h"

#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>


using namespace edm;
using namespace std;


DTDeadChannelTest::DTDeadChannelTest(const edm::ParameterSet& ps){

  debug = ps.getUntrackedParameter<bool>("debug", "false");
  if(debug)
    cout<<"[DTDeadChannelTest]: Constructor"<<endl;

  parameters = ps;

  dbe = edm::Service<DaqMonitorBEInterface>().operator->();
  dbe->setVerbose(1);

}

DTDeadChannelTest::~DTDeadChannelTest(){

  if(debug)
    cout << "DTDeadChannelTest: analyzed " << nevents << " events" << endl;

}

void DTDeadChannelTest::endJob(){

  if(debug)
    cout<<"[DTDeadChannelTest] endjob called!"<<endl;

  dbe->rmdir("DT/Tests/DTDeadChannel");

}

void DTDeadChannelTest::beginJob(const edm::EventSetup& context){

  if(debug)
    cout<<"[DTDeadChannelTest]: BeginJob"<<endl;

  nevents = 0;

  // Get the geometry
  context.get<MuonGeometryRecord>().get(muonGeom);

}

void DTDeadChannelTest::bookHistos(const DTLayerId & lId, int firstWire, int lastWire) {

  stringstream wheel; wheel << lId.superlayerId().wheel();
  stringstream station; station << lId.superlayerId().station();	
  stringstream sector; sector << lId.superlayerId().sector();
  stringstream superLayer; superLayer << lId.superlayerId().superlayer();
  stringstream layer; layer << lId.layer();

  string HistoName = "W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str() +  "_SL" + superLayer.str() +  "_L" + layer.str();
  string OccupancyDiffHistoName =  "OccupancyDiff_" + HistoName; 

  dbe->setCurrentFolder("DT/Tests/DTDeadChannel/Wheel" + wheel.str() +
			   "/Station" + station.str() +
			   "/Sector" + sector.str());

  OccupancyDiffHistos[HistoName] = dbe->book1D(OccupancyDiffHistoName.c_str(),OccupancyDiffHistoName.c_str(),lastWire-firstWire+1, firstWire-0.5, lastWire+0.5);

}

void DTDeadChannelTest::analyze(const edm::Event& e, const edm::EventSetup& context){
  
  nevents++;
  if (nevents%1 == 0 && debug) 
    cout<<"[DTDeadChannelTest]: "<<nevents<<" updates"<<endl;

  vector<DTChamber*>::const_iterator ch_it = muonGeom->chambers().begin();
  vector<DTChamber*>::const_iterator ch_end = muonGeom->chambers().end();

  cout<<endl;
  cout<<"[DTDeadChannelTest]: Occupancy tests results"<<endl;

  // Loop over the chambers
  for (; ch_it != ch_end; ++ch_it) {
    DTChamberId chID = (*ch_it)->id();
    vector<const DTSuperLayer*>::const_iterator sl_it = (*ch_it)->superLayers().begin(); 
    vector<const DTSuperLayer*>::const_iterator sl_end = (*ch_it)->superLayers().end();

    stringstream wheel; wheel << chID.wheel();
    stringstream station; station << chID.station();
    stringstream sector; sector << chID.sector();
    
    string HistoName = "W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str(); 

    // Get the ME produced by DigiTask Source
    MonitorElement * noise_histo = dbe->get(getMEName("OccupancyNoise_perCh", chID));	
    MonitorElement * hitInTime_histo = dbe->get(getMEName("OccupancyInTimeHits_perCh", chID));

    // ME -> TH2F
    if(noise_histo && hitInTime_histo) {	  
      MonitorElementT<TNamed>* occNoise = dynamic_cast<MonitorElementT<TNamed>*>(noise_histo);
      MonitorElementT<TNamed>* occInTime = dynamic_cast<MonitorElementT<TNamed>*>(hitInTime_histo);
      
      if (occNoise && occInTime) {
	TH2F * noise_histo_root = dynamic_cast<TH2F*> (occNoise->operator->());
	TH2F * hitInTime_histo_root = dynamic_cast<TH2F*> (occInTime->operator->());

	if (noise_histo_root && hitInTime_histo_root) {	      
	  
	  // Loop over the SuperLayers
	  for(; sl_it != sl_end; ++sl_it) {
	    DTSuperLayerId slID = (*sl_it)->id();
	    vector<const DTLayer*>::const_iterator l_it = (*sl_it)->layers().begin();
	    vector<const DTLayer*>::const_iterator l_end = (*sl_it)->layers().end();
      
	    // Loop over the layers
	    for(; l_it != l_end; ++l_it) {
	      DTLayerId lID = (*l_it)->id();

	      //Parameters to fill histos
	      stringstream superLayer; superLayer << slID.superlayer();
	      stringstream layer; layer << lID.layer();
	      string HistoNameTest = "W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str() +  "_SL" + superLayer.str() +  "_L" + layer.str();

	      const int firstWire = muonGeom->layer(lID)->specificTopology().firstChannel();
	      const int lastWire = muonGeom->layer(lID)->specificTopology().lastChannel();

	      int entry=-1;
	      if(slID.superlayer() == 1) entry=0;
	      if(slID.superlayer() == 2) entry=4;
	      if(slID.superlayer() == 3) entry=8;
	      int YBinNumber = entry+lID.layer();
	      

	      // Loop over the TH2F bin and fill the ME to be used for the Quality Test
	      for(int bin=firstWire; bin <= lastWire; bin++) {
		if (OccupancyDiffHistos.find(HistoNameTest) == OccupancyDiffHistos.end()) bookHistos(lID, firstWire, lastWire);
		float difference = hitInTime_histo_root->GetBinContent(bin, YBinNumber) - noise_histo_root->GetBinContent(bin, YBinNumber);
		OccupancyDiffHistos.find(HistoNameTest)->second->setBinContent(bin, difference);
	      }
	    } // loop on layers
	  } // loop on superlayers
	}
      }
    }
  } // loop on chambers

  // Occupancy Difference test 
  string OccupancyDiffCriterionName = parameters.getUntrackedParameter<string>("OccupancyDiffTestName","OccupancyDiffInRange"); 
  for(map<string, MonitorElement*>::const_iterator hOccDiff = OccupancyDiffHistos.begin();
      hOccDiff != OccupancyDiffHistos.end();
      hOccDiff++) {
    const QReport * theOccupancyDiffQReport = (*hOccDiff).second->getQReport(OccupancyDiffCriterionName);
    if(theOccupancyDiffQReport) {
      vector<dqm::me_util::Channel> badChannels = theOccupancyDiffQReport->getBadChannels();
      for (vector<dqm::me_util::Channel>::iterator channel = badChannels.begin(); 
	   channel != badChannels.end(); channel++) {
	cout<<"Bad occupancy difference channels: "<<(*channel).getBin()<<" "<<(*channel).getContents()<<endl;
      }
      cout<<"-------- "<<theOccupancyDiffQReport->getMessage()<<" ------- "<<theOccupancyDiffQReport->getStatus()<<endl;
    }
  }

  if (nevents%parameters.getUntrackedParameter<int>("resultsSavingRate",10) == 0){
    if ( parameters.getUntrackedParameter<bool>("writeHisto", true) ) 
      dbe->save(parameters.getUntrackedParameter<string>("outputFile", "DTDeadChannelTest.root"));
  }
}


string DTDeadChannelTest::getMEName(string histoTag, const DTChamberId & chId) {

  stringstream wheel; wheel << chId.wheel();
  stringstream station; station << chId.station();
  stringstream sector; sector << chId.sector();

  string folderName = 
    "Collector/FU0/DT/DTDigiTask/Wheel" +  wheel.str() +
    "/Station" + station.str() +
    "/Sector" + sector.str() + 
    "/Occupancies" + "/";
  
  string histoname = folderName + histoTag  
    + "_W" + wheel.str() 
    + "_St" + station.str() 
    + "_Sec" + sector.str();
  
  return histoname;
  
}
