/** \file
 * 
 * 
 * Authors N Ippolito and S Durkin
 *
*/
#include <FWCore/Utilities/interface/Exception.h>
#include <FWCore/Framework/interface/EventProcessor.h>
#include <FWCore/Utilities/interface/ProblemTracker.h>
#include <CondFormats/CSCObjects/interface/CSCReadoutMappingFromFile.h>
#include <CondFormats/CSCObjects/interface/CSCReadoutMappingForSliceTest.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>


#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include <DataFormats/FEDRawData/interface/FEDRawDataCollection.h>
#include <DataFormats/FEDRawData/interface/FEDHeader.h>
#include <DataFormats/FEDRawData/interface/FEDTrailer.h>
#include <DataFormats/FEDRawData/interface/FEDNumbering.h>

#include <EventFilter/CSCRawToDigi/interface/CSCDCCHeader.h>
#include <EventFilter/CSCRawToDigi/interface/CSCDCCEventData.h>
#include <EventFilter/CSCRawToDigi/interface/CSCDDUHeader.h>
#include <EventFilter/CSCRawToDigi/interface/CSCDDUEventData.h>
#include <DataFormats/CSCDigi/interface/CSCStripDigi.h>
#include <DataFormats/CSCDigi/interface/CSCWireDigi.h>
#include "EventFilter/CSCRawToDigi/interface/CSCEventData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCAnodeData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCALCTHeader.h"
#include "EventFilter/CSCRawToDigi/interface/CSCCLCTData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDDUTrailer.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDMBHeader.h"

#include "RecoLocalMuon/CSCStandAlone/interface/Strip_Fit_Constants.h"
#include "RecoLocalMuon/CSCStandAlone/interface/CatTrkFnd.h"
#include "RecoLocalMuon/CSCStandAlone/interface/AnoTrkFnd.h"
#include "RecoLocalMuon/CSCStandAlone/interface/EvtDisp.h"
#include "RecoLocalMuon/CSCStandAlone/interface/TrkFit3D.h"
#include "RecoLocalMuon/CSCStandAlone/interface/CatTrkDisp.h"
#include "RecoLocalMuon/CSCStandAlone/interface/TrkLink3D.h"
#include "RecoLocalMuon/CSCStandAlone/interface/EvtDumpTrack.h"
#include "RecoLocalMuon/CSCStandAlone/interface/RootWriter.h"

#include </afs/cern.ch/cms/external/lcg/external/root/5.10.00a/slc3_ia32_gcc323/root/include/TFile.h>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

using namespace edm;
using namespace std;

namespace test{
  
 
  
  class CSCTrackLink: public EDAnalyzer{

   
  public:
    CSCTrackLink(const ParameterSet& pset)
    {
      writer.setup();
      init=0;
    }
    
    ~CSCTrackLink()
    {
      
      writer.done();
    }
    
    void analyze(const Event & e, const EventSetup& c){ 
      
      
      // start durkin stuff
      if(init==0){
        init=1;
        printf(" Initialize \n");
        theMapping=CSCReadoutMappingFromFile("/home/ippolito/CMSSW_0_6_0/src/RecoLocalMuon/CSCStandAlone/test/csc_slice_test_map.txt");
	chamb_const= Strip_Fit_Constants();  
        // register all chambers
        int dmb_val[9]={1,2,3,4,5,7,8,9,10};
        for(int vmecrate=0;vmecrate<4;vmecrate++){
	  for(int i=0;i<9;i++){
	    int dmb=dmb_val[i];  
	    chamb_const.Register(theMapping,vmecrate,dmb);
	    
	  }
        }
      }

      
      // end durkin stuff 
      
      Handle<FEDRawDataCollection> rawdata;
      e.getByLabel("DaqSource", rawdata);
      for (int id=FEDNumbering::getCSCFEDIds().first;
	   id<=FEDNumbering::getCSCFEDIds().second; ++id){ //for each of our DCCs
	const FEDRawData& data = rawdata->FEDData(id);
	if(size_t size=data.size()) {
	  cout << "FED# " << id << " " << size << endl;
          if(id==750){
            unsigned short * buf = (unsigned short *)data.data();
            CSCDCCEventData dccEvent(buf);
            std::vector<CSCDDUEventData> & ddudata = dccEvent.dduData();
            
	    if(ddudata.size()==0)continue;
            // start durkin stuff 
            TrkLink3D trklink;

            printf(" DDU data size %d \n",ddudata.size());
            for(unsigned iddu=0;iddu<ddudata.size();iddu++){
              const std::vector<CSCEventData> & cscData = ddudata[iddu].cscData();
              int lvl1num=ddudata[iddu].header().lvl1num();
	      printf(" Level 1 Number %d \n",lvl1num);
              printf(" ddu number %d cscData.size() %d \n",iddu,cscData.size());                                                        
	      for (unsigned k=0; k<cscData.size(); ++k) {
	        AnoTrkFnd atrk(cscData[k]);              // find Anode tracks
	        CatTrkFnd ctrk(theMapping,cscData[k],chamb_const);   // find Cathode tracks (1/2 strip)
	        ctrk.CatTrkTime(theMapping,cscData[k],chamb_const);             // Fit time Cathode track hits
	        TrkFit3D triD(theMapping,cscData[k],ctrk,atrk,chamb_const);  // full blown gatti+
	        // least squares
	        //evtdump.dump3(theMapping,cscData[k],ctrk,atrk,triD);
	        evtdump.print_head(lvl1num,theMapping,cscData[k],ctrk,atrk,triD);
	        evtdump.print_gatti(theMapping,cscData[k],ctrk,atrk,triD);
	        evtdump.print_3Dtrack_fits(theMapping,cscData[k],ctrk,atrk,triD);
	        trklink.TrkLink3D_fill(theMapping,cscData[k],ctrk,atrk,triD);
		
	      }
            }
	    evtdump.TrkLink3D_dump(trklink);
	    // end durkin stuff
	    int *a = evtdump.getID1();
	    float* b= evtdump.getCSLP1();
	    float* c=evtdump.getCINT1();
	    float* d=evtdump.getASLP1();
	    float* e=evtdump.getAINT1();
	    int *f = evtdump.getID2();
	    float* g= evtdump.getCSLP2();
	    float* h=evtdump.getCINT2();
	    float *i=evtdump.getASLP2();
	    float *j=evtdump.getAINT2();
	    writer.save1(a,b,c,d,e,f,g,h,i,j);
          }
	}
	
      }
      
    }
   
  private:
    
    int init;
    Strip_Fit_Constants chamb_const;
    CSCReadoutMappingFromFile theMapping; 
    EvtDump evtdump;
    RootWriter writer; 
   
  };
  DEFINE_FWK_MODULE(CSCTrackLink)
    }
