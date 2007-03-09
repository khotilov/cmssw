#include <L1Trigger/CSCTrackFinder/src/CSCTFDTReceiver.h>
#include <L1Trigger/DTTrackFinder/src/L1MuDTTFConfig.h>
#include <L1Trigger/DTTraco/interface/DTConfigTraco.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/MuonDetId/interface/CSCTriggerNumbering.h>
#include <DataFormats/L1CSCTrackFinder/interface/CSCBitWidths.h>
#include <L1Trigger/CSCCommonTrigger/interface/CSCConstants.h>
#include <DataFormats/MuonDetId/interface/DTChamberId.h>

CSCTriggerContainer<csctf::TrackStub> CSCTFDTReceiver::process(const L1MuDTChambPhContainer* dttrig)
{
  dtstubs.clear();

  const int dt_minBX = L1MuDTTFConfig::getBxMin();
  const int dt_maxBX = L1MuDTTFConfig::getBxMax();
  const int dt_toffs = 0;// changed since DT tpg now centers around zero //abs(dt_maxBX - dt_minBX);

  // consider all BX
  for(int bx = dt_minBX + dt_toffs; bx <= dt_maxBX + dt_toffs; ++bx)
    for(int e = CSCDetId::minEndcapId(); e <= CSCDetId::maxEndcapId(); ++e)
      for(int s = CSCTriggerNumbering::minTriggerSectorId(); s <= CSCTriggerNumbering::maxTriggerSectorId(); ++s)
	{
	  int wheel = (e == 1) ? 2 : -2;
	  int sector = 2*s - 1;
	  int csc_bx = bx + CSCConstants::TIME_OFFSET;
	  
	  // combine two 30 degree DT sectors into a 60 degree CSC
	  // sector.
	  for(int is = sector; is <= sector+1; ++is)
	    {
	      int iss = (is == 12) ? 0 : is;
	      L1MuDTChambPhDigi* dtts[2];
	      
	      for(int stub = 0; stub < 2; ++stub)
		{
		  dtts[stub] = (stub == 0) ? dttrig->chPhiSegm1(wheel,1,iss,bx) :
		                             dttrig->chPhiSegm2(wheel,1,iss,bx);
		  if(dtts[stub])
		    {
		      // Convert stubs to CSC format (signed -> unsigned)
		      // phi was 12 bits (signed) for pi radians = 57.3 deg
		      // relative to center of 30 degree DT sector
		      double tmp = static_cast<const double> (dtts[stub]->phi()) /
			DTConfigTraco::RESOLPSIR * 180./M_PI + 15.;

		      int phi = static_cast<int> (tmp/60. * (1<<(CSCBitWidths::kGlobalPhiDataBitWidth)));
		      if (is>sector) phi = phi + (1<<(CSCBitWidths::kGlobalPhiDataBitWidth - 1));

		      // DT chambers may lie outside CSC sector boundary
		      // Eventually we need to extend CSC phi definition
		      phi = (phi>0) ? phi : 0;
		      phi = (phi<(1<<(CSCBitWidths::kGlobalPhiDataBitWidth))) ? phi : 
			(1<<(CSCBitWidths::kGlobalPhiDataBitWidth))-1;

		      // account for slope in DT/CSC comparison
		      phi = static_cast<int>(phi*(1.-40./4096.)) + 25;
		      phi = (phi>0) ? phi : 0;
		      phi = (phi<(1<<(CSCBitWidths::kGlobalPhiDataBitWidth))) ? phi : 
			(1<<(CSCBitWidths::kGlobalPhiDataBitWidth))-1;

		      // change phib from 10 bits to 6
		      int phib = (dtts[stub]->phiB() + DTConfigTraco::RESOLPSI) / 16;
		      int qual = dtts[stub]->code();
		      // barrel allows quality=0!
		      /// shift all by one and take mod 8, since DT quality of 7 is a null stub
		      qual = (qual + 1)%8;

		      CSCCorrelatedLCTDigi dtinfo(stub+1,1, qual, 0, 0, 0, phib, csc_bx, (stub+1) + 2*((is+1)%2));
		      DTChamberId dtid(wheel,1,is);
		      csctf::TrackStub tsCSC(dtinfo,dtid, phi, 0);

		      dtstubs.push_back(tsCSC);
		    }
		}
	    }
	}

  return dtstubs;
}

