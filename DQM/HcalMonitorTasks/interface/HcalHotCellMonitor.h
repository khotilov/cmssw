#ifndef DQM_HCALMONITORTASKS_HCALHOTCELLMONITOR_H
#define DQM_HCALMONITORTASKS_HCALHOTCELLMONITOR_H

#include "DQM/HcalMonitorTasks/interface/HcalBaseMonitor.h"
//#include "CalibFormats/HcalObjects/interface/HcalCalibrationWidths.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include <cmath>
#include <iostream>
#include <fstream>

/** \class HcalHotCellMonitor
  *
  * $Date: 2008/10/28 20:05:11 $
  * $Revision: 1.20 $
  * \author J. Temple - Univ. of Maryland
  */

struct hotNeighborParams{
  int DeltaIphi;
  int DeltaIeta;
  int DeltaDepth;
  double minCellEnergy; // cells below this threshold can never be considered "hot" by this algorithm
  double minNeighborEnergy; //neighbors must have some amount of energy to be counted
  double maxEnergy; //  a cell above this energy will always be considered hot
  double HotEnergyFrac; // a cell will be considered hot if neighbor energy/ cell energy is less than this value
};

class HcalHotCellMonitor: public HcalBaseMonitor {

 public:
  HcalHotCellMonitor();

  ~HcalHotCellMonitor();

  void setup(const edm::ParameterSet& ps, DQMStore* dbe);
  void setupNeighborParams(const edm::ParameterSet& ps, hotNeighborParams& N, char* type);
  void done(); // overrides base class function
  void clearME(); // overrides base class function
  void reset();

  void createMaps(const HcalDbService& cond);
  
  void processEvent(const HBHERecHitCollection& hbHits,
                    const HORecHitCollection& hoHits,
                    const HFRecHitCollection& hfHits,
		    //const ZDCRecHitCollection& zdcHits,
		    const HBHEDigiCollection& hbhedigi,
                    const HODigiCollection& hodigi,
                    const HFDigiCollection& hfdigi,
		    //const ZDCDigiCollection& zdcdigi, 
		    const HcalDbService& cond
		    );

  void processEvent_pedestal(const HBHEDigiCollection& hbhedigi,
			     const HODigiCollection& hodigi,
			     const HFDigiCollection& hfdigi,
			     //const ZDCDigiCollection& zdcdigi, 
			     const HcalDbService& cond
			     );

  void processEvent_rechitenergy( const HBHERecHitCollection& hbheHits,
				  const HORecHitCollection& hoHits,
				  const HFRecHitCollection& hfHits);

  void processEvent_rechitneighbors( const HBHERecHitCollection& hbheHits,
				     const HORecHitCollection& hoHits,
				     const HFRecHitCollection& hfHits);
  void fillHotHistosAtEndRun();

 private:
  void fillNevents_pedestal();
  void fillNevents_neighbor();
  void fillNevents_energy();
  void fillNevents_persistentenergy();

  void fillNevents_problemCells();

  bool doFCpeds_; //specify whether pedestals are in fC (if not, assume ADC)
  bool hotmon_makeDiagnostics_;

  // Booleans to control which of the three hot cell checking routines are used
  bool hotmon_test_pedestal_;
  bool hotmon_test_neighbor_;
  bool hotmon_test_energy_;
  bool hotmon_test_persistent_;

  int hotmon_checkNevents_;  // specify how often to check is cell is hot
  // Let each test have its own checkNevents value
  int hotmon_checkNevents_pedestal_;
  int hotmon_checkNevents_neighbor_;
  int hotmon_checkNevents_energy_;
  int hotmon_checkNevents_persistent_;

  double energyThreshold_, HBenergyThreshold_, HEenergyThreshold_, HOenergyThreshold_, HFenergyThreshold_, ZDCenergyThreshold_;
  double persistentThreshold_, HBpersistentThreshold_, HEpersistentThreshold_, HOpersistentThreshold_, HFpersistentThreshold_, ZDCpersistentThreshold_;

  MonitorElement* meEVT_;
  int ievt_;

  double hotmon_minErrorFlag_; // minimum error rate needed to dump out bad bin info 

  // Problem Histograms
  MonitorElement* ProblemHotCells;
  std::vector<MonitorElement*> ProblemHotCellsByDepth;
  
  double nsigma_;
  double HBnsigma_, HEnsigma_, HOnsigma_, HFnsigma_, ZDCnsigma_;
  std::vector<MonitorElement*>AboveNeighborsHotCellsByDepth;
  std::vector<MonitorElement*>AboveEnergyThresholdCellsByDepth;
  std::vector<MonitorElement*>AbovePersistentThresholdCellsByDepth; 
  std::vector<MonitorElement*>AbovePedestalHotCellsByDepth;
 
  // map of pedestals from database (in ADC)
  std::map<HcalDetId, float> pedestals_;
  std::map<HcalDetId, float> widths_;
  std::map<HcalDetId, float> pedestal_thresholds_;
  std::map<HcalDetId, double> rechitEnergies_;
  

  unsigned int abovepedestal[ETABINS][PHIBINS][4]; // filled when digi is above pedestal+nsigma
  unsigned int aboveneighbors[ETABINS][PHIBINS][4];
  unsigned int aboveenergy[ETABINS][PHIBINS][4]; // when rechit is above threshold energy
  unsigned int abovepersistent[ETABINS][PHIBINS][4]; // when rechit is consistently above some threshold
  unsigned int rechit_occupancy_sum[ETABINS][PHIBINS][4];
  float rechit_energy_sum[ETABINS][PHIBINS][4];

  // Diagnostic plots
  MonitorElement* d_HBnormped;
  MonitorElement* d_HEnormped;
  MonitorElement* d_HOnormped;
  MonitorElement* d_HFnormped;
  MonitorElement* d_ZDCnormped;

  MonitorElement* d_HBrechitenergy;
  MonitorElement* d_HErechitenergy;
  MonitorElement* d_HOrechitenergy;
  MonitorElement* d_HFrechitenergy;
  MonitorElement* d_ZDCrechitenergy;
 
  MonitorElement* d_HBenergyVsNeighbor;
  MonitorElement* d_HEenergyVsNeighbor;
  MonitorElement* d_HOenergyVsNeighbor;
  MonitorElement* d_HFenergyVsNeighbor;
  MonitorElement* d_ZDCenergyVsNeighbor;

  std::vector<MonitorElement*> d_avgrechitenergymap;
  
  hotNeighborParams defaultNeighborParams_, HBNeighborParams_, HENeighborParams_, HONeighborParams_, HFNeighborParams_, ZDCNeighborParams_;
};

#endif
