/*
 * \file L1TCSCTF.cc
 *
 * $Date: 2008/04/24 13:10:17 $
 * $Revision: 1.20 $
 * \author J. Berryhill
 *
 */

#include "DQM/L1TMonitor/interface/L1TCSCTF.h"
#include "DQMServices/Core/interface/DQMStore.h"

// KK_start: includes to fetch all reguired data products from the edm::Event
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCStatusDigiCollection.h"
// Also remember geometry classes, which are needed to SR LUTs to function
#include "L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
// KK_end

// JAG start: grab reco info
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/MuonDetId/interface/CSCDetId.h"
//#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"
//JAG end

using namespace std;
using namespace edm;

L1TCSCTF::L1TCSCTF(const ParameterSet& ps)
  // KK_start: if some piece of data is absent - configure corresponding source with 'null:'
//  : csctfSource_( ps.getParameter< InputTag >("csctfSource") )
    : gmtProducer( ps.getParameter< InputTag >("gmtProducer") ),
      lctProducer( ps.getParameter< InputTag >("lctProducer") ),
      trackProducer( ps.getParameter< InputTag >("trackProducer") ),
      statusProducer( ps.getParameter< InputTag >("statusProducer") )
  // KK_end
{

  	// verbosity switch
  	verbose_ = ps.getUntrackedParameter<bool>("verbose", false);

  	if(verbose_) cout << "L1TCSCTF: constructor...." << endl;


  	dbe = NULL;
  	if ( ps.getUntrackedParameter<bool>("DQMStore", false) )
  	{
    	dbe = Service<DQMStore>().operator->();
    	dbe->setVerbose(0);
  	}

  	outputFile_ = ps.getUntrackedParameter<string>("outputFile", "");
  	if ( outputFile_.size() != 0 ) 
	{
    	cout << "L1T Monitoring histograms will be saved to " << outputFile_.c_str() << endl;
  	}

  	bool disable = ps.getUntrackedParameter<bool>("disableROOToutput", false);
  	if(disable){
    	outputFile_="";
  	}


  	if ( dbe !=NULL ) 
	{
    	dbe->setCurrentFolder("L1T/L1TCSCTF");
  	}

  	// KK_start: instantiate standard on-fly SR LUTs from CSC TF emulator package
  	bzero(srLUTs_,sizeof(srLUTs_));
  	int endcap=1, sector=1; // assume SR LUTs are all same for every sector in either of endcaps
  	bool TMB07=true; // specific TMB firmware
  	// Create a dumy pset for SR LUTs
  	edm::ParameterSet srLUTset;
  	srLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
  	srLUTset.addUntrackedParameter<bool>("Binary",   false);
  	srLUTset.addUntrackedParameter<std::string>("LUTPath", "./");
  	for(int station=1,fpga=0; station<=4 && fpga<5; station++)
	{
    	if(station==1)
        	for(int subSector=0; subSector<2 && fpga<5; subSector++)
           		srLUTs_[fpga++] = new CSCSectorReceiverLUT(endcap, sector, subSector+1, station, srLUTset, TMB07);
     	else
        	srLUTs_[fpga++] = new CSCSectorReceiverLUT(endcap, sector, 0, station, srLUTset, TMB07);
  	}
  	// KK_end
}

L1TCSCTF::~L1TCSCTF()
{
}

void L1TCSCTF::beginJob(const EventSetup& c)
{

  	nev_ = 0;

  	// get hold of back-end interface
  	DQMStore* dbe = 0;
  	dbe = Service<DQMStore>().operator->();

  	if ( dbe ) {
    	dbe->setCurrentFolder("L1T/L1TCSCTF");
    	dbe->rmdir("L1T/L1TCSCTF");
  	}


	if ( dbe )
  	{
    	dbe->setCurrentFolder("L1T/L1TCSCTF");

    	//csctfetavalue[1] = dbe->book1D("CSCTF_eta_value","CSCTF eta value", 100, -2.5, 2.5 ) ;
    	//csctfetavalue[2] = dbe->book1D("CSCTF_eta_value_+1","CSCTF eta value bx +1", 100, -2.5, 2.5 ) ;
    	//csctfetavalue[0] = dbe->book1D("CSCTF_eta_value_-1","CSCTF eta value bx -1", 100, -2.5, 2.5 ) ;
    	//csctfphivalue[1] = dbe->book1D("CSCTF_phi_value","CSCTF phi value", 100, 0.0, 6.2832 ) ;
    	//csctfphivalue[2] = dbe->book1D("CSCTF_phi_value_+1","CSCTF phi value bx +1", 100, 0.0, 6.2832 ) ;
    	//csctfphivalue[0] = dbe->book1D("CSCTF_phi_value_-1","CSCTF phi value bx -1", 100, 0.0, 6.2832 ) ;
    	//csctfptvalue[1] = dbe->book1D("CSCTF_pt_value","CSCTF pt value", 160, -0.5, 159.5 ) ;
    	//csctfptvalue[2] = dbe->book1D("CSCTF_pt_value_+1","CSCTF pt value bx +1", 160, -0.5, 159.5 ) ;
    	//csctfptvalue[0] = dbe->book1D("CSCTF_pt_value_-1","CSCTF pt value bx -1", 160, -0.5, 159.5 ) ;
    	//csctfchargevalue[1] = dbe->book1D("CSCTF_charge_value","CSCTF charge value", 3, -1.5, 1.5 ) ;
    	//csctfchargevalue[2] = dbe->book1D("CSCTF_charge_value_+1","CSCTF charge value bx +1", 3, -1.5, 1.5 ) ;
    	//csctfchargevalue[0] = dbe->book1D("CSCTF_charge_value_-1","CSCTF charge value bx -1", 3, -1.5, 1.5 ) ;
    	//csctfquality[1] = dbe->book1D("CSCTF_quality","CSCTF quality", 20, -0.5, 19.5 ) ;
    	//csctfquality[2] = dbe->book1D("CSCTF_quality_+1","CSCTF quality bx +1", 20, -0.5, 19.5 ) ;
    	//csctfquality[0] = dbe->book1D("CSCTF_quality_-1","CSCTF quality bx -1", 20, -0.5, 19.5 ) ;
    	//csctfntrack = dbe->book1D("CSCTF_ntrack","CSCTF ntrack", 20, -0.5, 19.5 ) ;
    	
		
  		// KK_start: declaration of two monitoring histograms
  		//  Error counting histogram:
  		//  1) checks TF data integrity (error rate - first bin),
  		//  2) monitors sychronization on input links (4 errors types: SE/SM/BX/AF; ORed for all time bins, links, and SPs),
  		//  3) reports FMM status (if in any SP FMM status != "Ready" - fill the last bin)
  		csctferrors = dbe->book1D("CSCTF_errors","CSCTF Errors",6,0,6);
  		csctferrors->setAxisTitle("Error type",1);
  		csctferrors->setAxisTitle("Number of Errors",2);
  		csctferrors->setBinLabel(1,"Corruptions",1);
  		csctferrors->setBinLabel(2,"Synch. Err.",1);
  		csctferrors->setBinLabel(3,"Synch. Mod.",1);
  		csctferrors->setBinLabel(4,"BX mismatch",1);
  		csctferrors->setBinLabel(5,"Time misalign.",1);
  		csctferrors->setBinLabel(6,"FMM != Ready",1);
  		//  Occupancy histogram Eta x Y, where Y:
  		//  1) Phi_packed of input LCTs from 1st, 2nd, 3rd, and 4th stations
  		//  2) Phi_packed of output tracks
  		//  (all 12 SPs - 360 degree coveradge)
  		csctfoccupancies = dbe->book2D("CSCTF_occupancies","CSCTF Occupancies",100,0.8,2.5,1229,0,1.2);
  		csctfoccupancies->setAxisTitle("#eta",1);
  		csctfoccupancies->setAxisTitle("#phi of LCTs x station x endcap & #phi of tracks",2);
  		csctfoccupancies->setBinLabel(64,  "ME-1",2);
  		csctfoccupancies->setBinLabel(192, "ME+1",2);
  		csctfoccupancies->setBinLabel(320, "ME-2",2);
  		csctfoccupancies->setBinLabel(448, "ME+2",2);
  		csctfoccupancies->setBinLabel(576, "ME-3",2);
  		csctfoccupancies->setBinLabel(704, "ME+3",2);
  		csctfoccupancies->setBinLabel(832, "ME-4",2);
  		csctfoccupancies->setBinLabel(960, "ME+4",2);
  		csctfoccupancies->setBinLabel(1088,"Tracks",2);
  		// KK_end
  		
		//JAG
		haloDelEta23 = dbe->book1D("CSCTF_Halo_Eta23","Delta station 2 to station 3 Eta for Halo Muons", 40, -0.20,0.30);
		
		
		csctfTrackQ = dbe->book1D("CSCTF_Track_Q","CSC Track Quality", 16, -0.5, 15.5);
		csctfTrackQ->setAxisTitle("Track Type", 1);
		csctfTrackQ->setBinLabel(3,"ME1-2-3",1);
		csctfTrackQ->setBinLabel(7,"ME1-2",1);
		csctfTrackQ->setBinLabel(8,"ME1-3",1);
		csctfTrackQ->setBinLabel(9,"ME2-3",1);
		csctfTrackQ->setBinLabel(16,"Halo Trigger",1);
		
  		csctfChamberOccupancies = dbe->book2D("CSCTF_Chamber_Occupancies","CSCTF Chamber Occupancies", 54, -0.05, 5.35, 10, -5.5, 4.5);
		csctfChamberOccupancies->setAxisTitle("Sector (Endcap), (chambers 1-9 not labeled)",1);
  		csctfChamberOccupancies->setBinLabel(1,"ME-4",2);
		csctfChamberOccupancies->setBinLabel(2,"ME-3",2);
		csctfChamberOccupancies->setBinLabel(3,"ME-2",2);
		csctfChamberOccupancies->setBinLabel(4,"ME-1b",2);
		csctfChamberOccupancies->setBinLabel(5,"ME-1a",2);
		csctfChamberOccupancies->setBinLabel(6,"ME+1a",2);
		csctfChamberOccupancies->setBinLabel(7,"ME+1b",2);
		csctfChamberOccupancies->setBinLabel(8,"ME+2",2);
		csctfChamberOccupancies->setBinLabel(9,"ME+3",2);
		csctfChamberOccupancies->setBinLabel(10,"ME+4",2);
		csctfChamberOccupancies->setBinLabel(1, "1(+), 7(-)",1);
		csctfChamberOccupancies->setBinLabel(10,"2(+), 8(-)",1);
		csctfChamberOccupancies->setBinLabel(19,"3(+), 9(-)",1);
		csctfChamberOccupancies->setBinLabel(28,"4(+), 10(-)",1);
		csctfChamberOccupancies->setBinLabel(37,"5(+), 11(-)",1);
		csctfChamberOccupancies->setBinLabel(46,"6(+), 12(-)",1);
		
		csctfTrackPhi = dbe->book1D("CSCTF_Track_Phi", "CSCTF Track Phi", 144, 0, 360);
		csctfTrackPhi->setAxisTitle("Track #phi", 1);
		csctfTrackEta = dbe->book1D("CSCTF_Track_Eta", "CSCTF Track Eta", 32, 0.9, 2.5);
		csctfTrackEta->setAxisTitle("Track #eta", 1);
		
		csctfbx = dbe->bookProfile("CSCTF_bx","CSCTF bx", 36, 0.5, 36.5, 15, -5.5,9.5 ) ;
		csctfbx->setAxisTitle("Sector, Endcap, MPC Link", 1);
		csctfbx->setBinLabel(1,"1, +",1);
		csctfbx->setBinLabel(4,"2, +",1);
		csctfbx->setBinLabel(7,"3, +",1);
		csctfbx->setBinLabel(10,"4, +",1);
		csctfbx->setBinLabel(13,"5, +",1);
		csctfbx->setBinLabel(16,"6, +",1);
		csctfbx->setBinLabel(19,"1, -",1);
		csctfbx->setBinLabel(22,"2, -",1);
		csctfbx->setBinLabel(25,"3, -",1);
		csctfbx->setBinLabel(28,"4, -",1);
		csctfbx->setBinLabel(31,"5, -",1);
		csctfbx->setBinLabel(34,"6, -",1);
		
		cscTrackStubNumbers = dbe->book1D("CSCTF_TrackStubs", "Number of Stubs in CSC Tracks", 5, 0.5, 5.5);
		
		csctfntrack = dbe->book1D("CSCTF_ntrack","Number of CSCTracks found per event", 5, 0.5, 5.5 ) ;
  		//JAG
	}
}


void L1TCSCTF::endJob(void)
{

	if(verbose_) cout << "L1TCSCTF: end job...." << endl;
  	LogInfo("EndJob") << "analyzed " << nev_ << " events";

 	if ( outputFile_.size() != 0  && dbe ) dbe->save(outputFile_);

	return;
}

void L1TCSCTF::analyze(const Event& e, const EventSetup& c)
{
	int NumCSCTfTracksRep = 0;
	nev_++;
  	if(verbose_) cout << "L1TCSCTF: analyze...." << endl;

  	edm::Handle<L1MuGMTReadoutCollection> pCollection;

  	// KK_start ///////////////////////////////////
	
  	if( gmtProducer.label() != "null" )
	{ // GMT block
    	e.getByLabel(gmtProducer,pCollection);
  		// KK_end   ///////////////////////////////////

		if (!pCollection.isValid()) 
  		{
    		edm::LogInfo("DataNotFound") << "can't find L1MuGMTReadoutCollection with label ";  // << csctfSource_.label() ;
    		return;
  		}

  		L1MuGMTReadoutCollection const* gmtrc = pCollection.product();
  		vector<L1MuGMTReadoutRecord> gmt_records = gmtrc->getRecords();
  		vector<L1MuGMTReadoutRecord>::const_iterator RRItr;

  		int ncsctftrack = 0;
  		
      	//csctfntrack->Fill(ncsctftrack);
		
      	if (verbose_)
		{
     	cout << "\tCSCTFCand ntrack " << ncsctftrack << endl;
		}
  
  	// KK_start  ///////////////////////////////////
  	} // end of GMT block
  	// KK_end    ///////////////////////////////////


  	// KK_start ///////////////////////////////////
  	if( statusProducer.label() != "null" )
	{
    	edm::Handle<L1CSCStatusDigiCollection> status;
     	e.getByLabel(statusProducer.label(),statusProducer.instance(),status);
     	bool integrity=status->first, se=false, sm=false, bx=false, af=false, fmm=false;
     	for(std::vector<L1CSCSPStatusDigi>::const_iterator stat=status->second.begin(); stat!=status->second.end(); stat++)
		{
        	se |= stat->SEs()&0xFFF;
        	sm |= stat->SMs()&0xFFF;
        	bx |= stat->BXs()&0xFFF;
        	af |= stat->AFs()&0xFFF;
        	fmm|= stat->FMM()!=8;
     	}
     
	 	if(integrity) csctferrors->Fill(0.5);
     	if(se)        csctferrors->Fill(1.5);
     	if(sm)        csctferrors->Fill(2.5);
     	if(bx)        csctferrors->Fill(3.5);
     	if(af)        csctferrors->Fill(4.5);
     	if(fmm)       csctferrors->Fill(5.5);
  	}


  	if( lctProducer.label() != "null" )
	{
     	edm::ESHandle<CSCGeometry> pDD;
     	c.get<MuonGeometryRecord>().get( pDD );
     	CSCTriggerGeometry::setGeometry(pDD);

     	edm::Handle<CSCCorrelatedLCTDigiCollection> corrlcts;
     	e.getByLabel(lctProducer.label(),lctProducer.instance(),corrlcts);

     	for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc=corrlcts.product()->begin(); csc!=corrlcts.product()->end(); csc++)
		{
        	CSCCorrelatedLCTDigiCollection::Range range1 = corrlcts.product()->get((*csc).first);
        	for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range1.first; lct!=range1.second; lct++)
			{
           		int endcap  = (*csc).first.endcap()-1;
           		int station = (*csc).first.station()-1;
           		int sector  = (*csc).first.triggerSector()-1;
           		int subSector = CSCTriggerNumbering::triggerSubSectorFromLabels((*csc).first);
           		int cscId   = (*csc).first.triggerCscId()-1;
           		int fpga    = ( subSector ? subSector-1 : station+1 );
				
				//JAG
				int endcapAssignment = 1;
				int shift = 1;
				float sectorArg = sector;
				//float sectorArg = j;
				
				if( endcap == 1 ){
					endcapAssignment = -1;
					shift = 2;
					//sectorArg = sector - 6;
				}
				
				int signedStation = (station + shift)* endcapAssignment;
				if( (station == 0) && (endcap == 0)) signedStation = subSector - 1;
				if( (station == 0) && (endcap == 1)) signedStation = (-1)*subSector;
				
				float chamberArg1 = cscId * 0.1 + sectorArg;
				//float chamberArg1 = i*0.1 + sectorArg;
				//std::cout << "First" << i << " " << sectorArg << " " << chamberArg1 << std::endl;
				
				float chamberArg11 = chamberArg1;
				if(sectorArg == 1) chamberArg1 = chamberArg11 - 0.1;
				if(sectorArg == 2) chamberArg1 = chamberArg11 - 0.2;
				if(sectorArg == 3) chamberArg1 = chamberArg11 - 0.3;
				if(sectorArg == 4) chamberArg1 = chamberArg11 - 0.4;
				if(sectorArg == 5) chamberArg1 = chamberArg11 - 0.5;

				//std::cout << "cscId, station, sector, endcap, sectorArg, chamber Arg: " << cscId << ", " << station << ", " <<sector << ", " << endcap << ", " << chamberArg1 << ", " << signedStation << std::endl;			

				csctfChamberOccupancies->Fill(chamberArg1, signedStation); 
				int bunchX = ( (lct->getBX()) - 6 );
				
				int timingSectorArg = 3*(sector) + (lct->getMPCLink());
				if( endcap == 1) timingSectorArg = 3*(sector + 6) + (lct->getMPCLink());
				//std::cout << "Sector, MPCLink, TSA, endcap: " << sector << ", " << lct->getMPCLink() << ", " << timingSectorArg << ", " << endcap << std::endl;
				
				csctfbx->Fill(timingSectorArg, bunchX );
				
				//std::cout << "LCT'S, encap: " << endcap << ", station: " << station << ", sector: " << sector << ", subSector: " << subSector << ", cscId: " << cscId << std:: endl;

				//End JAG
				
				// Check if Det Id is within pysical range:
           		if( endcap<0||endcap>1 || sector<0||sector>6 || station<0||station>3 || cscId<0||cscId>8 || fpga<0||fpga>4)
				{
              		edm::LogError("L1CSCTF: CSC TP are out of range: ")<<"  endcap: "<<(endcap+1)<<"  station: "<<(station+1) <<"  sector: "<<(sector+1)<<"  subSector: "<<subSector<<"  fpga: "<<fpga<<"  cscId: "<<(cscId+1);
              		continue;
           		}
				
           		lclphidat lclPhi;
           		
				try {
             		lclPhi = srLUTs_[fpga]->localPhi(lct->getStrip(), lct->getPattern(), lct->getQuality(), lct->getBend());
           		} catch(...) { 
					bzero(&lclPhi,sizeof(lclPhi)); 
				}
				
           		gblphidat gblPhi;
				
           		try {
             		gblPhi = srLUTs_[fpga]->globalPhiME(lclPhi.phi_local, lct->getKeyWG(), cscId+1);
           		} catch(...) { 
					bzero(&gblPhi,sizeof(gblPhi)); 
				}
				
           		gbletadat gblEta;
				
           		try {
             		gblEta = srLUTs_[fpga]->globalEtaME(lclPhi.phi_bend_local, lclPhi.phi_local, lct->getKeyWG(), cscId+1);
           		} catch(...) { 
					bzero(&gblEta,sizeof(gblEta)); 
				}
           
		   		// SR LUT gives packed eta and phi values -> normilize them to 1 by scale them to 'max' and shift by 'min'
           		csctfoccupancies->Fill( gblEta.global_eta/127. * 1.5 + 0.9, (gblPhi.global_phi + ( sector + (endcap?0:6) )*4096 + station*4096*12) * 1./(4*4096*12) );
        	}//lct != range1.scond
     	}//csc!=corrlcts.product()->end()
  	}// lctProducer.label() != "null"

  	if( trackProducer.label() != "null" )
	{
		
     	edm::Handle<L1CSCTrackCollection> tracks;
     	e.getByLabel(trackProducer.label(),trackProducer.instance(),tracks);
     	for(L1CSCTrackCollection::const_iterator trk=tracks->begin(); trk<tracks->end(); trk++){
        	NumCSCTfTracksRep++;
			
			csctfoccupancies->Fill( trk->first.eta_packed()/32. * 1.5 + 0.9, trk->first.phi_packed()*0.2/32. + 1.);
			//JAG_START
			edm::Handle<CSCCorrelatedLCTDigiCollection> corrlcts;
     		e.getByLabel(lctProducer.label(),lctProducer.instance(),corrlcts);
			
			long LUTAdd = trk->first.ptLUTAddress();
			int trigMode = ( (LUTAdd)&0xf0000 ) >> 16;
			
			if( trigMode == 15 ){
				
				csctfTrackQ->Fill( trigMode );
				double haloVals[4][4];
				for( int i = 0; i < 4; i++)
				{
					haloVals[i][0] = 0;
				}
				
				edm::Handle<CSCCorrelatedLCTDigiCollection> corrlcts;
     			e.getByLabel(lctProducer.label(),lctProducer.instance(),corrlcts);
     			for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc=corrlcts.product()->begin(); csc!=corrlcts.product()->end(); csc++)
				{
        			CSCCorrelatedLCTDigiCollection::Range range1 = corrlcts.product()->get((*csc).first);
        			for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range1.first; lct!=range1.second; lct++)
					{
						int endcap  = (*csc).first.endcap()-1;
           				int station = (*csc).first.station()-1;
           				int sector  = (*csc).first.triggerSector()-1;
		           		int cscId   = (*csc).first.triggerCscId()-1;
           				int subSector = CSCTriggerNumbering::triggerSubSectorFromLabels((*csc).first);
						int fpga    = ( subSector ? subSector-1 : station+1 );
						
						if( (station == 1) || (station == 2) )
						{
							int modEnd;
							if( endcap == 0 ) modEnd = -1;
							if( endcap == 1 ) modEnd = 1;
							int indexHalo = modEnd + station;
							if(haloVals[indexHalo][0] == 1.0) haloVals[indexHalo][3] = 1.0;
							if(haloVals[indexHalo][0] == 0) haloVals[indexHalo][0] = 1.0;
							haloVals[indexHalo][1] = sector*1.0;
						
							lclphidat lclPhi;
           		   			lclPhi = srLUTs_[fpga]->localPhi(lct->getStrip(), lct->getPattern(), lct->getQuality(), lct->getBend());
           					
           					gblphidat gblPhi;
				   			gblPhi = srLUTs_[fpga]->globalPhiME(lclPhi.phi_local, lct->getKeyWG(), cscId+1);
           		
           					gbletadat gblEta;
				   			gblEta = srLUTs_[fpga]->globalEtaME(lclPhi.phi_bend_local, lclPhi.phi_local, lct->getKeyWG(), cscId+1);
							
							haloVals[indexHalo][2] = gblEta.global_eta/127. * 1.5 + 0.9;
						
						} //station1 or 2
					} //lct first to second
				} //corrlcts
				
				if( (haloVals[0][0] == 1.) && (haloVals[1][0] == 1.) && (haloVals[0][3] != 1.) && (haloVals[1][3] != 1.)  ){
					if( haloVals[0][1] == haloVals[1][1] ){
						double delEta23 = haloVals[1][2] - haloVals[0][2];
						haloDelEta23->Fill( delEta23 );
					}
				}
				
				if( (haloVals[2][0] == 1.) && (haloVals[3][0] == 1.) && (haloVals[2][3] != 1.) && (haloVals[3][3] != 1.)  ){
					if( haloVals[2][1] == haloVals[3][1] ){
						double delEta23 = haloVals[3][2] - haloVals[2][2];
						haloDelEta23->Fill( delEta23 );
					}
				}
				
			} //halo trigger
			
			csctfTrackPhi->Fill( (trk->first.phi_packed()) *2.5 + (trk->first.sector() - 1 )*61 );
			csctfTrackEta->Fill( ( (trk->first.eta_packed() )*0.0125) + 0.9 );
			
			int cscTrackStub = 0;

			CSCCorrelatedLCTDigiCollection lctsOfTracks=trk->second;

     		for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator trackStub=lctsOfTracks.begin(); trackStub!=lctsOfTracks.end(); trackStub++)
			{
        		CSCCorrelatedLCTDigiCollection::Range range2 = lctsOfTracks.get((*trackStub).first);
        		for(CSCCorrelatedLCTDigiCollection::const_iterator lct=range2.first; lct!=range2.second; lct++)
				{
					cscTrackStub++;
           			
				}
			}
			
			//
     		//edm::Handle<L1CSCTrack> tfInf;
			//e.getByLabel(trackProducer.label(), tfInf);
			
			cscTrackStubNumbers->Fill(cscTrackStub);
				
			//JAG_END
		}
  	}
	csctfntrack->Fill(NumCSCTfTracksRep);
  	// KK_end    ///////////////////////////////////
}
