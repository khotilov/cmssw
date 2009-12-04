#include "CondFormats/DataRecord/interface/L1MuGMTScalesRcd.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "L1Trigger/DTTrackFinder/interface/L1MuDTTrack.h"
#include "L1Trigger/GlobalMuonTrigger/src/L1MuGMTConfig.h"

#include "SLHCUpgradeSimulations/L1Trigger/interface/DTL1SimOperation.h"


void L1DTSimOperation::getDTSimTrigger(edm::Event& event, 
				       const edm::EventSetup& eventSetup)
{
  outAscii << "\nTrigger block\n-------------" << endl;

  // L1 Local Trigger Block ----------------------------------------------------
  bool BTI_done   = false;
  // bool TRACO_done = false;

  // BTI
  vector<DTBtiTrigData> btitrigs = theTrigger->BtiTrigs();
  vector<DTBtiTrigData>::const_iterator pbti;
  int ibti = 0;
  if(debug_bti) 
    outAscii << "\n[DTTrigger]: " << btitrigs.size() << " BTI triggers found" 
	     << endl;
  DTBtiTrigger* aDTBti;
  for ( pbti = btitrigs.begin(); pbti != btitrigs.end(); pbti++ ) 
    {
      Global3DPoint pos = theTrigger->CMSPosition(&(*pbti));
      Global3DVector dir = theTrigger->CMSDirection(&(*pbti));
      aDTBti = new DTBtiTrigger(*pbti, pos, dir);
      BtiTrigs->push_back(*aDTBti);
      ++ibti;
      if(debug_bti && ( ibti < 100 ) ) {
	outAscii << ibti << ")" << endl; 
	outAscii << aDTBti->sprint();
      }
      delete aDTBti;
    }
  BTI_done = true;  
  
  
  // TSPhi
  vector<DTChambPhSegm> theTSPhTrigs = theTrigger->TSPhTrigs();
  vector<DTChambPhSegm>::const_iterator tsphi;
  //DTTSphiTrigger* aTSphiTrig;
  if(debug_tsphi)
    outAscii << "\n[DTTrigger]: " << theTSPhTrigs.size() << " TSPhi triggers found" 
	     << endl;
  for (tsphi = theTSPhTrigs.begin(); tsphi != theTSPhTrigs.end(); tsphi++) {
    /*
    Global3DPoint pos = theTrigger->CMSPosition(&(*tsphi)); 
    Global3DVector dir = theTrigger->CMSDirection(&(*tsphi));
    aTSphiTrig = new DTTSphiTrigger(*tsphi, pos, dir);
    TSPhiTrigs->push_back(*aTSphiTrig);
    ++itsphi;
    if(debug_tsphi && (itsphi < 40 ) ) {
      outAscii << itsphi << ")" << endl;
      outAscii << aTSphiTrig->sprint();
    }
    */
    // *********************************
    // *** match TSphi-BTItheta (SV) *** 
    // *********************************
    if( BTI_done &&  (tsphi->station()==1 || tsphi->station()==2) ) {
      for(BtiTrigsCollection::iterator bti = BtiTrigs->begin();
	  bti != BtiTrigs->end(); bti++) 
	//if( match(*bti, *aTSphiTrig) )
	if( match(*bti, *tsphi) )
	  // Add DTMatch to DTStubCollection
	  //DTStubMatches->addDT(*bti, *aTSphiTrig, debug_dttrackmatch_extrapolation);
	  DTStubMatches->addDT(*bti, *tsphi, debug_dttrackmatch_extrapolation);
    }
    // *** end match TSphi-BTItheta ***
    //delete aTSphiTrig; 
  }
  if(debug_dtmatch) {
    outAscii 
      << "\n=========================================================="
      << "\nNumber of DTMatch: " <<  DTStubMatches->numDt()
      << "; in station 1: " << DTStubMatches->numDt(1)
      << "; in station 2: " << DTStubMatches->numDt(2)
      << endl;
    for (int idtmatch = 0; idtmatch < DTStubMatches->numDt(); idtmatch++)
      outAscii 
	<< " " << idtmatch << ")"
	<< " bx " << DTStubMatches->dtmatch(idtmatch)->bx()
	<< " st " << DTStubMatches->dtmatch(idtmatch)->station()
	<< " wh " << DTStubMatches->dtmatch(idtmatch)->wheel()
	<< " se " << DTStubMatches->dtmatch(idtmatch)->sector()
	<< " code " << DTStubMatches->dtmatch(idtmatch)->code()
	<< "\n    phi_glo " <<  DTStubMatches->dtmatch(idtmatch)->phi_glo()
	<< " phi_ts " <<  DTStubMatches->dtmatch(idtmatch)->phi_ts()
	<< " phib " <<  DTStubMatches->dtmatch(idtmatch)->phib_ts()
	<< " theta " << DTStubMatches->dtmatch(idtmatch)->theta_ts()
	<< endl;
    outAscii 
      << "==========================================================" 
      << endl;
  }

  /*
  // TSTheta
  vector<DTChambThSegm> theTSThTrigs = theTrigger->TSThTrigs();
  vector<DTChambThSegm>::const_iterator tsTheta;
  int itsTheta = 0; 
  if(debug_tstheta)
    outAscii << "\n[DTTrigger]: " << theTSThTrigs.size() << " TSTheta triggers found" 
	     << endl;
  DTTSthetaTrigger* aTSthetaTrig;
  for (tsTheta = theTSThTrigs.begin(); tsTheta != theTSThTrigs.end(); tsTheta++) {
    aTSthetaTrig = new DTTSthetaTrigger(*tsTheta);
    TSThetaTrigs->push_back(*aTSthetaTrig);
    ++itsTheta;
    if(debug_tstheta && (itsTheta < 40 ) ) {
      outAscii << itsTheta << ")" << endl;
      outAscii << aTSthetaTrig->sprint();
    }
    delete aTSthetaTrig; 
  }
  */
  return;



}
