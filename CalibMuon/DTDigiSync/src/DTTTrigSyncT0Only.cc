/*
 *  See header file for a description of this class.
 *
 *  $Date: 2007/02/19 11:45:21 $
 *  $Revision: 1.1 $
 *  \author G. Cerminara - INFN Torino
 */

#include "DTTTrigSyncT0Only.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "CondFormats/DTObjects/interface/DTT0.h"
#include "CondFormats/DataRecord/interface/DTT0Rcd.h"

#include <iostream>

using namespace std;
using namespace edm;


DTTTrigSyncT0Only::DTTTrigSyncT0Only(const ParameterSet& config){
  debug = config.getUntrackedParameter<bool>("debug");
}



DTTTrigSyncT0Only::~DTTTrigSyncT0Only(){}



void DTTTrigSyncT0Only::setES(const EventSetup& setup) {
  ESHandle<DTT0> t0;
  setup.get<DTT0Rcd>().get(t0);
  tZeroMap = &*t0;
  
  if(debug) {
    cout << "[DTTTrigSyncT0Only] T0 version: " << t0->version() << endl;
  }
}




double DTTTrigSyncT0Only::offset(const DTLayer* layer,
				  const DTWireId& wireId,
				  const GlobalPoint& globPos,
				  double& tTrig,
				  double& wirePropCorr,
				  double& tofCorr) {
  tTrig = offset(wireId);
  wirePropCorr = 0;
  tofCorr = 0;

  if(debug) {
    cout << "[DTTTrigSyncT0Only] Offset (ns): " << tTrig + wirePropCorr - tofCorr << endl
	 << "      various contributions are: " << endl
	 //<< "      tZero (ns):   " << t0 << endl
	 << "      Propagation along wire delay (ns): " <<  wirePropCorr << endl
	 << "      TOF correction (ns): " << tofCorr << endl
	 << endl;
  }
  //The global offset is the sum of various contributions
  return tTrig + wirePropCorr - tofCorr;
}

double DTTTrigSyncT0Only::offset(const DTWireId& wireId) {
  float t0 = 0;
  float t0rms = 0;
  tZeroMap->get(wireId,
		t0,
		t0rms,
		DTTimeUnits::ns);
  return t0;
}


// Set the verbosity level
bool DTTTrigSyncT0Only::debug;




