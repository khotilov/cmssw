#ifndef DQM_HCALMONITORTASKS_HCALBASEMONITOR_H
#define DQM_HCALMONITORTASKS_HCALBASEMONITOR_H

// Define number of eta, phi bins for histogram objects
#define ETABINS 87
#define PHIBINS 72

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CondFormats/HcalObjects/interface/HcalQIECoder.h"
#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "TH1F.h"
#include "TH2F.h"
#include <map>

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "FWCore/Utilities/interface/CPUTimer.h"

#include <iostream>

// Temporary fix:  Add this into base class until I figure why multiple inclusions are a problem -- Jeff, 23 May 2008
#include "CalibFormats/HcalObjects/interface/HcalCalibrationWidths.h"

using namespace std;
/** \class HcalBaseMonitor
  *  
  * $Date: 2008/10/15 20:08:41 $
  * $Revision: 1.13 $
  * \author W. Fisher - FNAL
  */
class HcalBaseMonitor {
public:
  HcalBaseMonitor(); 
  virtual ~HcalBaseMonitor(); 

  virtual void setup(const edm::ParameterSet& ps, DQMStore* dbe);
  virtual void done();
  virtual void clearME();

  void setVerbosity(int verb) { fVerbosity = verb; }
  int getVerbosity() const { return fVerbosity; }
  
  void setDiagnostics(bool myval) { makeDiagnostics=myval;}
  bool getDiagnostics() const { return makeDiagnostics;}

  bool vetoCell(HcalDetId id);
  bool validDetId(HcalSubdetector subdet, int tower_ieta, int tower_iphi, int depth); // determine whether ID is valid (disable at some point)


protected:
  
  int fVerbosity;
  bool showTiming; // controls whether to show timing diagnostic info
  int checkNevents_; // controls when histograms should be updated

  double etaMax_, etaMin_;
  double phiMax_, phiMin_;
  int etaBins_, phiBins_;
  double minErrorFlag_;
  bool checkHB_, checkHE_, checkHO_, checkHF_;

  edm::CPUTimer cpu_timer; // 
    
  bool makeDiagnostics; // controls whether to make diagnostic plots

  DQMStore* m_dbe;
  vector<string> hotCells_;
  string rootFolder_;
  string baseFolder_;

};

#endif
