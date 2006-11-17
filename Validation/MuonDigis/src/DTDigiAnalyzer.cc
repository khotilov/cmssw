/** \class DTDigiAnalyzer
 *  Analyse the the muon-drift-tubes digitizer. 
 *  
 *  $Date: 2006/08/07 15:47:36 $
 *  $Revision: 1.5 $
 *  \authors: R. Bellan
 * $Revision: 1.6 $
 *   B. de la Cruz 
 *  Date 2006/11
 *  Addition of some plots
 */



#include <FWCore/Framework/interface/Event.h>
#include <DataFormats/DTDigi/interface/DTDigiCollection.h>
#include "Validation/MuonDigis/interface/DTDigiAnalyzer.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/Vector/interface/LocalPoint.h"

#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"

#include <iostream>
#include <string>

#include "TFile.h"

#include "SimMuon/DTDigitizer/test/Histograms.h"

hDigis hDigis_global("Global");
hDigis hDigis_W0("Wheel0");
hDigis hDigis_W1("Wheel1");
hDigis hDigis_W2("Wheel2");

hHits hAllHits("AllHits");


using namespace edm;
using namespace std;

DTDigiAnalyzer:: DTDigiAnalyzer(const ParameterSet& pset){
//   MCStatistics = new DTMCStatistics();
//   MuonDigiStatistics = new DTMuonDigiStatistics();
//   HitsAnalysis = new DTHitsAnalysis();
  label = pset.getUntrackedParameter<string>("label");  
  file = new TFile("DTDigiPlots.root","RECREATE");
  file->cd();
  DigiTimeBox = new TH1F("DigiTimeBox","Digi Time Box",2048,0,1600);
  DigiTimeBox_wheel2m = new TH1F("DigiTimeBox_wheel2m","Digi Time Box wheel -2",512,0,1600);
  
  DigiTimeBox_wheel1m = new TH1F("DigiTimeBox_wheel1m","Digi Time Box wheel -1",512,0,1600);
  DigiTimeBox_wheel0 = new TH1F("DigiTimeBox_wheel0","Digi Time Box wheel 0",512,0,1600);
  DigiTimeBox_wheel1p = new TH1F("DigiTimeBox_wheel1p","Digi Time Box wheel 1",512,0,1600);
  DigiTimeBox_wheel2p = new TH1F("DigiTimeBox_wheel2p","Digi Time Box wheel 2",512,0,1600);
  DigiEfficiencyMu = new TH1F("DigiEfficiencyMu", "Ratio (#Digis Mu)/(#SimHits Mu)",100, 0., 5.);
  DigiEfficiency = new TH1F("DigiEfficiency", "Ratio (#Digis)/(#SimHits)",100, 0., 5.); 
  DoubleDigi = new TH1F("Number_Digi_per_layer","Number_Digi_per_layer",10,0.,10.);  
  SimvsDigi = new TH2F("Number_simhit_vs_digi","Number_simhit_vs_digi",70, 0., 70., 70, 0., 70.);
  Wire_DoubleDigi = new TH1F("Wire_Number_with_double_Digi","Wire_Number_with_double_Digi",100,0.,100.);

  MB1_sim_occup = new TH1F("Simhit_occupancy_MB1", "Simhit_occupancy_MB1", 53, 0., 53. );
  MB1_digi_occup = new TH1F("Digi_occupancy_MB1", "Digi_occupancy_MB1", 53, 0., 53. );

  MB2_sim_occup = new TH1F("Simhit_occupancy_MB2", "Simhit_occupancy_MB2", 63, 0., 63. );
  MB2_digi_occup = new TH1F("Digi_occupancy_MB2", "Digi_occupancy_MB2", 63, 0., 63. );

  MB3_sim_occup = new TH1F("Simhit_occupancy_MB3", "Simhit_occupancy_MB3", 75, 0., 75. );
  MB3_digi_occup = new TH1F("Digi_occupancy_MB3", "Digi_occupancy_MB3", 75, 0., 75. );

  MB4_sim_occup = new TH1F("Simhit_occupancy_MB4", "Simhit_occupancy_MB4", 99, 0., 99. );
  MB4_digi_occup = new TH1F("Digi_occupancy_MB4", "Digi_occupancy_MB4", 99, 0., 99. );

  char stringcham[40];

  for ( int slnum = 0; slnum < 62; ++slnum )
  {
   sprintf(stringcham, "DigiTimeBox_slid_%d", slnum) ;
   DigiHisto =  new TH1F(stringcham, stringcham, 100,0,1600);
   DigiTimeBox_SL.push_back(DigiHisto);
  }


  if(file->IsOpen()) cout<<"file open!"<<endl;
  else cout<<"*** Error in opening file ***"<<endl;
}

DTDigiAnalyzer::~DTDigiAnalyzer(){
  //  cout<<"Number of analyzed event: "<<nevts<<endl;
  //HitsAnalysis->Report();
  file->cd();
  DigiTimeBox->Write();
  DigiTimeBox_wheel2m->Write();
  DigiTimeBox_wheel1m->Write();
  DigiTimeBox_wheel0->Write();
  DigiTimeBox_wheel1p->Write();
  DigiTimeBox_wheel2p->Write();
  DigiEfficiency->Write();
  DigiEfficiencyMu->Write();
  DoubleDigi->Write();
  SimvsDigi->Write();
  Wire_DoubleDigi->Write();

  MB1_sim_occup->Write();
  MB1_digi_occup->Write();
  MB2_sim_occup->Write();
  MB2_digi_occup->Write();
  MB3_sim_occup->Write();
  MB3_digi_occup->Write();
  MB4_sim_occup->Write();
  MB4_digi_occup->Write();

   for ( int slnum = 0; slnum < 61; ++slnum )
   {
   DigiTimeBox_SL[slnum]->Write();
   }
  
  //cout << "--- OK DigiTimeBox- " << endl;
  // hDigis_global.Write();
  // hDigis_W0.Write();
  // hDigis_W1.Write();
  // hDigis_W2.Write();
  // hAllHits.Write();
  file->Close();
  //    delete file;
  // delete DigiTimeBox;
}

void  DTDigiAnalyzer::analyze(const Event & event, const EventSetup& eventSetup){
  cout << "--- Run: " << event.id().run()
       << " Event: " << event.id().event() << endl;
  
  Handle<DTDigiCollection> dtDigis;
  event.getByLabel(label, dtDigis);
   
 // cout << " OK for the moment " << endl;

 
  Handle<PSimHitContainer> simHits; 
//   event.getByLabel("g4SimHits","MuonDTHits",simHits);    
//   For the PSimHit part generated by G4, for the 'simevent.root' type files
//   event.getByLabel("SimG4Object","muonDTHits",simHits);  
     event.getByLabel("SimG4Object","MuonDTHits",simHits);
// event.getByType(simHits);

  
  ESHandle<DTGeometry> muonGeom;
  eventSetup.get<MuonGeometryRecord>().get(muonGeom);

  int num_simhits;
  int num_mudigis;
  int num_musimhits;
  int num_digis;
  int num_digis_layer;
  int cham_num ;
  int wire_touched; 
  num_digis = 0;
  num_mudigis = 0;
  num_musimhits = 0;
  DTWireIdMap wireMap;     
 
   num_simhits = simHits->size();
 //  cout << "num simhits " << num_simhits << endl;

  for(vector<PSimHit>::const_iterator hit = simHits->begin();
      hit != simHits->end(); hit++){    
    // Create the id of the wire, the simHits in the DT known also the wireId
     DTWireId wireId(hit->detUnitId());
   //  cout << " PSimHits wire id " << wireId << " part type " << hit->particleType() << endl;
    // Fill the map
    wireMap[wireId].push_back(&(*hit));

    LocalPoint entryP = hit->entryPoint();
    LocalPoint exitP = hit->exitPoint();
    int partType = hit->particleType();
    if ( abs(partType) == 13 ) num_musimhits++;    
  
    if ( wireId.station() == 1 && abs(partType) == 13 ) MB1_sim_occup->Fill(wireId.wire());
    if ( wireId.station() == 2 && abs(partType) == 13 ) MB2_sim_occup->Fill(wireId.wire());
    if ( wireId.station() == 3 && abs(partType) == 13 ) MB3_sim_occup->Fill(wireId.wire());
    if ( wireId.station() == 4 && abs(partType) == 13 ) MB4_sim_occup->Fill(wireId.wire());


    float path = (exitP-entryP).mag();
    float path_x = fabs((exitP-entryP).x());
    
    hAllHits.Fill(entryP.x(),exitP.x(),
		   entryP.y(),exitP.y(),
		   entryP.z(),exitP.z(),
		   path , path_x, 
		   partType, hit->processType(),
		  hit->pabs());
  }

  // cout << "num muon simhits " << num_musimhits << endl;
 
  DTDigiCollection::DigiRangeIterator detUnitIt;
  for (detUnitIt=dtDigis->begin();
       detUnitIt!=dtDigis->end();
       ++detUnitIt){
    
    const DTLayerId& id = (*detUnitIt).first;
    const DTDigiCollection::Range& range = (*detUnitIt).second;
    
    // DTLayerId print-out
 //   cout<<"--------------"<<endl;
 //   cout<<"id: "<<id;
  
   // num_mudigis = 0;  
    num_digis_layer = 0 ;  
    cham_num = 0 ;
    wire_touched = 0;

   // Loop over the digis of this DetUnit
    for (DTDigiCollection::const_iterator digiIt = range.first;
	 digiIt!=range.second;
	 ++digiIt){
  //    cout<<" Wire: "<<(*digiIt).wire()<<endl
  //	  <<" digi time (ns): "<<(*digiIt).time()<<endl;
      
      num_digis++;
      num_digis_layer++;
      if (num_digis_layer > 1 )
      {
       if ( (*digiIt).wire() == wire_touched )
        {
          Wire_DoubleDigi->Fill((*digiIt).wire());
    //      cout << "antiguo wire y nuevo " << wire_touched << " " << (*digiIt).wire() << endl; 
        }
      }
      wire_touched = (*digiIt).wire();   
 
      DigiTimeBox->Fill((*digiIt).time());
      if (id.wheel() == -2 ) DigiTimeBox_wheel2m->Fill((*digiIt).time());      
      if (id.wheel() == -1 ) DigiTimeBox_wheel1m->Fill((*digiIt).time());
      if (id.wheel() == 0 ) DigiTimeBox_wheel0->Fill((*digiIt).time());
      if (id.wheel() == 1 ) DigiTimeBox_wheel1p->Fill((*digiIt).time());
      if (id.wheel() == 2 ) DigiTimeBox_wheel2p->Fill((*digiIt).time());
  
   //   Superlayer number and fill histo with digi timebox
 
      cham_num = (id.wheel() +2)*12 + (id.station() -1)*3 + id.superlayer();
   //   cout << " Histo number " << cham_num << endl;

      DigiTimeBox_SL[cham_num]->Fill((*digiIt).time());

    //  cout << " size de digis " << (*digiIt).size() << endl;
      
      DTWireId wireId(id,(*digiIt).wire());
      if (wireId.station() == 1 ) MB1_digi_occup->Fill((*digiIt).wire());
      if (wireId.station() == 2 ) MB2_digi_occup->Fill((*digiIt).wire());
      if (wireId.station() == 3 ) MB3_digi_occup->Fill((*digiIt).wire());
      if (wireId.station() == 4 ) MB4_digi_occup->Fill((*digiIt).wire());

      int mu=0;
      float theta = 0;
      
      for(vector<const PSimHit*>::iterator hit = wireMap[wireId].begin();
	  hit != wireMap[wireId].end(); hit++)
	if( abs((*hit)->particleType()) == 13){
	  theta = atan( (*hit)->momentumAtEntry().x()/ (-(*hit)->momentumAtEntry().z()) )*180/M_PI;
    //	  cout<<"momentum x: "<<(*hit)->momentumAtEntry().x()<<endl
    //	      <<"momentum z: "<<(*hit)->momentumAtEntry().z()<<endl
    //	      <<"atan: "<<theta<<endl;
	  mu++;
	}
     
  //     cout << " mu number " << mu << "  num_digis_layer " << num_digis_layer << endl;
     if(mu ) num_mudigis++;  
     if(mu && theta){
	hDigis_global.Fill((*digiIt).time(),theta,id.superlayer());
	//filling digi histos for wheel and for RZ and RPhi
	WheelHistos(id.wheel())->Fill((*digiIt).time(),theta,id.superlayer());
      }
	  
    }// for digis in layer
    DoubleDigi->Fill( (float)num_digis_layer );
  }// for layers
 // cout << " num digis " << num_digis << " num_mudigis  " << num_mudigis << endl;
  DigiEfficiencyMu->Fill( (float)num_mudigis/(float)num_musimhits );
  DigiEfficiency->Fill( (float)num_digis/(float)num_musimhits );
  SimvsDigi->Fill( (float)num_musimhits, (float)num_digis ) ;
//  cout<<"--------------"<<endl;
}

hDigis* DTDigiAnalyzer::WheelHistos(int wheel){
  switch(abs(wheel)){

  case 0: return  &hDigis_W0;
  
  case 1: return  &hDigis_W1;
    
  case 2: return  &hDigis_W2;
     
  default: return NULL;
  }
}


#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DTDigiAnalyzer)
