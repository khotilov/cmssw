#ifndef EETriggerTowerTask_H
#define EETriggerTowerTask_H

/*
 * \file EETriggerTowerTask.h
 *
 * $Date: 2007/07/19 11:48:49 $
 * $Revision: 1.9 $
 * \author C. Bernet
 *
*/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// #include "boost/multi_array.hpp"

// #include <iostream>
// #include <fstream>
#include <vector>


class MonitorElement;
class DaqMonitorBEInterface;

class EETriggerTowerTask : public edm::EDAnalyzer {

 public:

  /// Constructor
  EETriggerTowerTask(const edm::ParameterSet& ps);

  /// Destructor
  virtual ~EETriggerTowerTask();

  /// number of trigger towers in eta
  static const int nTTEta; 

  /// number of trigger towers in phi
  static const int nTTPhi; 

  /// number of supermodules
  static const int nSM; 

 protected:

  /// Analyze
  void analyze(const edm::Event& e, 
	       const edm::EventSetup& c);

  /// BeginJob
  void beginJob(const edm::EventSetup& c);

  /// EndJob
  void endJob(void);

  /// Setup
  void setup(void);

  /// Cleanup
  void cleanup(void);

 private:
  
  /// 1D array
  typedef std::vector<MonitorElement*> array1;

  /// reserve an array to hold one histogram per supermodule
  void reserveArray( array1& array );

  /// process a collection of digis, either real or emulated
  void processDigis( const edm::Handle<EcalTrigPrimDigiCollection>& digis, 
		     array1& meEtMap,
		     array1& meVeto,
		     array1& meFlags,
		     const edm::Handle<EcalTrigPrimDigiCollection>& digis
		     = edm::Handle<EcalTrigPrimDigiCollection>());


  /// book monitor elements for real, or emulated digis
  void setup( DaqMonitorBEInterface* dbe,
	      const char* nameext,
	      const char* folder, 
	      bool emulated);
  

  /// local event counter
  int ievt_;

  /// Et vs ix vs iy, for each SM 
  array1 meEtMapReal_;

  /// fine grain veto vs ix vs iy, for each SM 
  array1 meVetoReal_;

  /// flag vs ix vs iy, for each SM   
  array1 meFlagsReal_;
  
  /// Emulated Et vs ix vs iy, for each SM 
  array1 meEtMapEmul_;

  /// Emulated fine grain veto vs ix vs iy, for each SM 
  array1 meVetoEmul_;

  /// Emulated flag vs ix vs iy, for each SM   
  array1 meFlagsEmul_;
  
  /// error flag vs ix vs iy, for each SM
  /// the error flag is set to true in case of a discrepancy between 
  /// the emulator and the real data
  array1 meEmulError_;
  array1 meVetoEmulError_;
  array1 meFlagEmulError_;

  /// not sure this is necessary
  bool init_;

  /// to find the input collection of real digis 
  edm::InputTag  realCollection_;

  /// to find the input collection of emulated digis
  edm::InputTag  emulCollection_;
  
  /// debug output root file. if empty, no output file created.
  std::string   outputFile_;
};

#endif
